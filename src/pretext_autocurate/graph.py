"""hifiasm assembly graphs and YaHS AGPs: which contig ends belong together.

YaHS builds each scaffold from hifiasm's haplotype contigs (`--yahs-agp`, one AGP per
haplotype) and puts a gap between neighbours. hifiasm's contig graph (`--gfa`, the
`*.hapN.p_ctg.gfa` files; with or without sequence) links contig ends that overlap in the
reads (L lines; the overlap, 7-21 kb on Troglobisium, is the CIGAR length). Together they
say, for a gap in a scaffold:

  link     the graph links the two contig ends YaHS put next to each other: the join is
           right, and the gap can be closed by dropping the overlap from one side
  bridge   the only link out of the left contig and the only link into the right one
           both go to the same third contig, which YaHS left as a scaffold of its own:
           that contig fills the gap (typically a piece too small for Hi-C to place)

The contig sequence is in the assembly FASTA (through the AGP), so GFAs without sequence
are enough here; GFA sequence is read when present but not needed. Sequence changes (closing
a gap, inserting a bridge) cannot be expressed in a PretextView savestate; they are made on
the assembly files after pretext-to-asm (`close_gaps`, called by `asm`).
"""
import collections
import gzip
import os
import re
from typing import NamedTuple

FLIP = {"+": "-", "-": "+"}
_COMP = str.maketrans("ACGTNacgtnRYKMSWBDHVrykmswbdhv", "TGCANtgcanYRMKSWVHDByrmkswvhdb")


def revcomp(s):
    return s.translate(_COMP)[::-1]


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


# --------------------------------------------------------------------------- FASTA access

class Faidx:
    """Random access to a plain FASTA through its .fai (built if missing); for a gzipped
    FASTA the needed records are read whole into memory on first use."""

    def __init__(self, path):
        self.path = str(path)
        self.gz = self.path.endswith(".gz")
        self.idx, self._cache = {}, {}
        if self.gz:
            return
        fai = self.path + ".fai"
        if not os.path.exists(fai) or os.path.getmtime(fai) < os.path.getmtime(self.path):
            self._build(fai)
        with open(fai) as f:
            for line in f:
                n, ln, off, lb, lw = line.split("\t")[:5]
                self.idx[n] = (int(ln), int(off), int(lb), int(lw))

    def _build(self, fai):
        """samtools-style .fai (name, length, offset, bases per line, bytes per line)."""
        out, name = [], None
        with open(self.path, "rb") as f:
            off = 0
            ln = lb = lw = 0
            seq_off = 0
            for raw in f:
                if raw.startswith(b">"):
                    if name is not None:
                        out.append(f"{name}\t{ln}\t{seq_off}\t{lb}\t{lw}\n")
                    name = raw[1:].split()[0].decode()
                    ln = lb = lw = 0
                    seq_off = off + len(raw)
                else:
                    b = len(raw.rstrip(b"\r\n"))
                    if lb == 0:
                        lb, lw = b, len(raw)
                    ln += b
                off += len(raw)
            if name is not None:
                out.append(f"{name}\t{ln}\t{seq_off}\t{lb}\t{lw}\n")
        with open(fai, "w") as f:
            f.writelines(out)

    def length(self, name):
        if self.gz:
            return len(self._whole(name))
        return self.idx[name][0]

    def _whole(self, name):
        if name not in self._cache:
            cur, buf = None, []
            with _open(self.path) as f:
                for line in f:
                    if line.startswith(">"):
                        if cur == name:
                            break
                        cur = line[1:].split()[0]
                    elif cur == name:
                        buf.append(line.strip())
            self._cache[name] = "".join(buf)
        return self._cache[name]

    def fetch(self, name, start, end, upper=True):
        """Bases [start, end) (0-based), upper case unless `upper` is False."""
        if self.gz:
            s = self._whole(name)[start:end]
            return s.upper() if upper else s
        ln, off, lb, lw = self.idx[name]
        start, end = max(0, start), min(ln, end)
        if end <= start:
            return ""
        a = off + (start // lb) * lw + start % lb
        b = off + ((end - 1) // lb) * lw + (end - 1) % lb + 1
        with open(self.path, "rb") as f:
            f.seek(a)
            raw = f.read(b - a)
        s = raw.replace(b"\n", b"").replace(b"\r", b"").decode()
        return s.upper() if upper else s


# --------------------------------------------------------------------------- YaHS AGPs

class Placement(NamedTuple):
    """A contig in an assembly scaffold: [start, end) 0-based, orientation in the scaffold,
    and its index among the scaffold's contigs."""
    scaffold: str
    start: int
    end: int
    orient: str
    k: int


def parse_spec(spec):
    """'H1=path' -> ('H1', path); a bare path -> (prefix from 'hapN' in the file name, path)."""
    if "=" in spec and not os.path.exists(spec):
        p, path = spec.split("=", 1)
        return p, path
    m = re.search(r"hap(\d+)", os.path.basename(spec))
    return (f"H{m.group(1)}" if m else ""), spec


class Layout:
    """The YaHS scaffolds: contig placements and, per scaffold, its contigs in order."""

    def __init__(self, specs, names=None):
        """`specs`: AGP paths ('H1=path' or a path whose name contains hapN). YaHS scaffold
        names (scaffold_N) get the haplotype prefix (H1_scaffold_N) unless `names` (the
        assembly's scaffold names) already contains them unprefixed."""
        self.place, self.scaffolds = {}, collections.defaultdict(list)
        for spec in specs:
            prefix, path = parse_spec(spec)
            with _open(path) as f:
                for line in f:
                    if line.startswith("#") or not line.strip():
                        continue
                    c = line.rstrip("\n").split("\t")
                    if len(c) < 9 or c[4] not in "WD":   # W: a contig (component)
                        continue
                    sc = c[0]
                    if prefix and not (names is not None and sc in names) and not sc.startswith(prefix + "_"):
                        sc = f"{prefix}_{sc}"
                    lst = self.scaffolds[sc]
                    p = Placement(sc, int(c[1]) - 1, int(c[2]), c[8], len(lst))
                    lst.append(c[5])
                    self.place[c[5]] = p

    def neighbours(self):
        """(scaffold, left contig, its orientation, right contig, its orientation, gap
        start, gap end) for every two contigs YaHS put next to each other (bp, 0-based)."""
        for sc, lst in self.scaffolds.items():
            for a, b in zip(lst, lst[1:]):
                pa, pb = self.place[a], self.place[b]
                yield sc, a, pa.orient, b, pb.orient, pa.end, pb.start

    def check(self, lengths):
        """Scaffolds of the AGPs missing from `lengths` ({name: bp}) or of another length."""
        bad = []
        for sc, lst in self.scaffolds.items():
            end = self.place[lst[-1]].end
            if sc not in lengths:
                bad.append(f"{sc} not in the assembly")
            elif abs(lengths[sc] - end) > 1:
                bad.append(f"{sc}: {end} bp in the AGP, {lengths[sc]} in the assembly")
        return bad


# --------------------------------------------------------------------------- contig graph

class Graph:
    """Links between contig ends from hifiasm GFAs: out[(c, o)] = {(d, od, overlap)}: leaving
    contig c traversed in orientation o, the path enters d in orientation od, the two
    sharing `overlap` bp."""

    def __init__(self, paths, end_reads=8):
        """Also keeps, from the `A` lines (reads placed in each contig), the `end_reads`
        outermost reads at each contig end (`head`, `tail`: [(read, strand)]) and the contig
        of every read (`read_ctg`), to place contigs on the unitig graph."""
        self.out = collections.defaultdict(set)
        self.length, self.has_seq = {}, False
        self.read_ctg, places = {}, collections.defaultdict(list)
        for path in paths:
            with _open(path) as f:
                for line in f:
                    t = line[:2]
                    if t == "S\t":
                        c = line.rstrip("\n").split("\t", 3)
                        self.has_seq |= c[2] != "*"
                        tags = c[3] if len(c) > 3 else ""
                        m = re.search(r"LN:i:(\d+)", tags)
                        self.length[c[1]] = int(m.group(1)) if m else len(c[2])
                    elif t == "A\t":
                        c = line.split("\t", 7)
                        off = int(c[2])
                        places[c[1]].append((off, off + int(c[6]) - int(c[5]), c[3], c[4]))
                        self.read_ctg[c[4]] = c[1]
                    elif t == "L\t":
                        c = line.split("\t")
                        m = re.match(r"(\d+)M$", c[5])
                        ov = int(m.group(1)) if m else 0
                        self.out[(c[1], c[2])].add((c[3], c[4], ov))
                        self.out[(c[3], FLIP[c[4]])].add((c[1], FLIP[c[2]], ov))
        self.head, self.tail = {}, {}
        for ctg, v in places.items():
            self.head[ctg] = [(r, s) for _, _, s, r in sorted(v)[:end_reads]]
            self.tail[ctg] = [(r, s) for _, _, s, r in sorted(v, key=lambda x: -x[1])[:end_reads]]

    def link(self, a, oa, b, ob):
        """Overlap of the link a(oa) -> b(ob), or None."""
        for d, od, ov in self.out.get((a, oa), ()):
            if (d, od) == (b, ob):
                return ov
        return None

    def into(self, b, ob):
        """Links entering b(ob): {(a, oa, overlap)}."""
        return {(d, FLIP[od], ov) for d, od, ov in self.out.get((b, FLIP[ob]), ())}


# --------------------------------------------------------------------------- unitig graph

class UnitigGraph:
    """hifiasm's unitig graph (`*.p_utg.gfa` or `*.r_utg.gfa`, with or without sequence):
    oriented links with their overlaps, the unitig of every read (`A` lines), and, if the file
    has sequence, where each unitig's sequence is in the file (read on demand)."""

    def __init__(self, path):
        self.path = str(path)
        self.length, self.offset, self.depth = {}, {}, {}
        self.adj = collections.defaultdict(list)       # (u, o) -> [((v, ov_), overlap)]
        self.read = {}                                 # read -> (unitig, strand)
        self.reads = collections.defaultdict(list)     # unitig -> [read]
        if self.path.endswith(".gz"):
            raise ValueError("the unitig GFA must be uncompressed (its sequence is read by position)")
        pos = 0
        with open(self.path, "rb") as f:
            for raw in f:
                t = raw[:2]
                if t == b"S\t":
                    c = raw[2:].split(b"\t", 2)
                    name = c[0].decode()
                    tags = c[2] if len(c) > 2 else b""
                    m = re.search(rb"LN:i:(\d+)", tags)
                    seq_len = len(c[1].rstrip(b"\r\n"))
                    self.length[name] = int(m.group(1)) if m else seq_len
                    m = re.search(rb"rd:i:(\d+)", tags)
                    if m:
                        self.depth[name] = int(m.group(1))
                    if c[1].rstrip(b"\r\n") != b"*":
                        self.offset[name] = pos
                elif t == b"L\t":
                    c = raw.decode().split("\t")
                    m = re.match(r"(\d+)M$", c[5])
                    ov = int(m.group(1)) if m else 0
                    self.adj[(c[1], c[2])].append(((c[3], c[4]), ov))
                elif t == b"A\t":
                    c = raw.decode().split("\t", 6)
                    self.read[c[4]] = (c[1], c[3])
                    self.reads[c[1]].append(c[4])
                pos += len(raw)
        for k in self.adj:
            self.adj[k] = sorted(set(self.adj[k]))
        self.has_seq = bool(self.offset)
        # haploid read depth: the length-weighted mode of hifiasm's rd over long unitigs (the
        # lower of the two peaks: heterozygous sequence; collapsed homozygous unitigs sit at
        # about twice that). Troglobisium: 26 and 56.
        hist = collections.Counter()
        for u, d in self.depth.items():
            if self.length[u] >= 50000:
                hist[d] += self.length[u]
        smooth = {d: sum(hist.get(d + k, 0) for k in (-2, -1, 0, 1, 2)) for d in hist}
        peaks = [d for d in smooth if smooth[d] >= max(smooth.get(d - 3, 0), smooth.get(d + 3, 0))
                 and smooth[d] > 0.05 * max(smooth.values())] if smooth else []
        self.haploid = float(min(peaks)) if peaks else None

    def depth_class(self, units):
        """'haploid', 'collapsed' (about twice the haploid depth), 'repeat' (more) or
        'between', from the length-weighted depth of `units`; with that depth."""
        units = [u for u in units if u in self.depth]
        if not units or not self.haploid:
            return "unknown", 0.0
        L = sum(self.length[u] for u in units)
        d = sum(self.depth[u] * self.length[u] for u in units) / L
        r = d / self.haploid
        cls = "haploid" if r <= 1.4 else "collapsed" if 1.6 <= r <= 2.6 else "repeat" if r > 2.6 else "between"
        return cls, d

    def seq(self, u):
        """Sequence of unitig u (upper case), from the file."""
        with open(self.path, "rb") as f:
            f.seek(self.offset[u])
            return f.readline()[2:].split(b"\t", 2)[1].decode().upper()

    def overlap(self, a, b):
        return next(ov for n, ov in self.adj[a] if n == b)

    def end_node(self, graph, ctg, side):
        """The unitig (oriented along the contig) where contig `ctg` starts ('h') or ends ('t'),
        from its outermost reads that are in a unitig; None if none is."""
        for rd, st in (graph.head if side == "h" else graph.tail).get(ctg, []):
            if rd in self.read:
                u, ust = self.read[rd]
                return (u, "+" if st == ust else "-")
        return None

    def ends(self, graph, a, oa, b, ob):
        """(exit node of a as laid out, entry node of b) or (None, None)."""
        ex = self.end_node(graph, a, "t" if oa == "+" else "h")
        en = self.end_node(graph, b, "h" if ob == "+" else "t")
        if ex is None or en is None:
            return None, None
        return (ex if oa == "+" else (ex[0], FLIP[ex[1]])), (en if ob == "+" else (en[0], FLIP[en[1]]))

    def shortest(self, ex, en, limit):
        """Shortest oriented path ex .. en by bp of the unitigs between them, up to `limit`:
        (bp, [nodes]) or (None, None)."""
        import heapq
        if ex == en:
            return 0, [ex]
        seen, prev, pq = {ex: 0}, {}, [(0, ex)]
        while pq:
            d, n = heapq.heappop(pq)
            if n == en:
                path = [n]
                while path[-1] != ex:
                    path.append(prev[path[-1]])
                return d, path[::-1]
            if d > seen.get(n, 1e18) or d > limit:
                continue
            for m, _ in self.adj.get(n, ()):
                nd = d + (0 if m == en else self.length[m[0]])
                if nd < seen.get(m, 1e18):
                    seen[m], prev[m] = nd, n
                    heapq.heappush(pq, (nd, m))
        return None, None

    def paths(self, ex, en, best, slack=150000, cap=20, budget=300000):
        """Simple paths ex .. en no longer than the shortest + `slack` bp, up to `cap`, shortest
        first: [(bp, (nodes))]. Several: bubbles (alternative alleles) or repeats."""
        out, steps, stack = [], 0, [(ex, 0, (ex,))]
        while stack and len(out) < cap and steps < budget:
            steps += 1
            x, d, p = stack.pop()
            if x == en:
                out.append((d, p))
                continue
            for m, _ in self.adj.get(x, ()):
                nd = d + (0 if m == en else self.length[m[0]])
                if nd <= best + slack and m not in p:
                    stack.append((m, nd, p + (m,)))
        return sorted(out) or [(best, ())]

    def path_seq(self, path):
        """Sequence spelled by an oriented path (overlaps removed)."""
        out = []
        for i, (u, o) in enumerate(path):
            s = self.seq(u)
            s = s if o == "+" else revcomp(s)
            out.append(s if i == 0 else s[self.overlap(path[i - 1], (u, o)):])
        return "".join(out)


class UnitigPath(NamedTuple):
    """A unitig-graph path between the two contig ends at a gap."""
    scaffold: str
    gap_start: int
    gap_end: int
    left: str
    left_o: str
    right: str
    right_o: str
    bp: int             # bp of unitigs between the two end unitigs
    path: tuple         # oriented unitigs, first and last holding the contig ends
    n_paths: int        # 1: unique; more: bubbles (up to the count's cap)
    source: str         # where the reads of the inner unitigs went (contigs / none)
    alts: tuple = ()    # every path found, shortest first: ((bp, nodes), ...)


def unitig_paths(layout, graph, utg, limit=1_000_000, pairs=None):
    """A unitig path for every YaHS gap (or for `pairs`: [(scaffold, gap start, gap end, a,
    oa, b, ob)]) whose two contig ends are connected within `limit` bp."""
    out = []
    for sc, a, oa, b, ob, g0, g1 in (pairs if pairs is not None else layout.neighbours()):
        ex, en = utg.ends(graph, a, oa, b, ob)
        if ex is None:
            continue
        bp, path = utg.shortest(ex, en, limit)
        if path is None:
            continue
        home = collections.Counter()
        for u, _ in path[1:-1]:
            for rd in utg.reads.get(u, []):
                c = graph.read_ctg.get(rd)
                home["no contig" if c is None else c] += 1
        tot = sum(home.values())
        src = ", ".join(f"{k} {100 * v / tot:.0f}%" for k, v in home.most_common(3)) if tot else ""
        alts = utg.paths(ex, en, bp)
        if not alts[0][1]:
            alts = [(bp, tuple(path))]
        out.append(UnitigPath(sc, g0, g1, a, oa, b, ob, bp, tuple(path), len(alts), src, tuple(alts)))
    return out


def _hap_of(layout, ctg):
    """Haplotype prefix (H1, H2 ...) of the scaffold a contig is in, or ''."""
    p = layout.place.get(ctg)
    return p.scaffold.split("_")[0] if p else ""


def _along(utg, graph, c, nodes):
    """Orientation of contig c along a unitig path if the whole contig lies on it (its first
    and last unitig both on the path, in order): '+' / '-', else None."""
    h, t = utg.end_node(graph, c, "h"), utg.end_node(graph, c, "t")
    if h is None or t is None:
        return None
    pos = {n: i for i, n in enumerate(nodes)}
    if h in pos and t in pos and pos[h] <= pos[t]:
        return "+"
    hr, tr = (h[0], FLIP[h[1]]), (t[0], FLIP[t[1]])
    if hr in pos and tr in pos and pos[tr] <= pos[hr]:
        return "-"
    return None


def choose_fill(fa, layout, graph, utg, up, own_max=0.1):
    """What to put in a gap a unitig path crosses, without duplicating sequence the assembly
    already holds for that haplotype:
      - a path with no unitig between the two end unitigs: a direct join (overlap trimmed,
        or the few bases the end unitigs add, classed like a fill);
      - else the shortest path whose inner unitigs are at about twice the haploid depth
        (collapsed homozygous sequence: the genome has two copies, the assembly one) and
        whose reads are not in this haplotype's contigs (at most `own_max` of them): fill;
      - else nothing: a haploid path is one haplotype's sequence, in the other haplotype's
        contigs (filling would duplicate it), and a repeat is not resolved.
    Returns (action, value, note): ('overlap', bp), ('fill', sequence) or (None, why)."""
    hap = _hap_of(layout, up.left)              # the haplotype of the contigs around the gap
    flank = {up.left, up.right}                 # their reads reach into the unitigs past their ends

    def own_share(units):
        n = own = 0
        for u in units:
            for rd in utg.reads.get(u, []):
                c = graph.read_ctg.get(rd)
                if c is not None and c not in flank:
                    n += 1
                    own += _hap_of(layout, c) == hap
        return own / n if n else 0.0

    tried = []
    fills = [(bp, nodes, path_fill(fa, layout, utg, up._replace(path=nodes, bp=bp))) for bp, nodes in up.alts]
    # 1. contigs that overlap (or abut) along a path: nothing to add, close
    for bp, nodes, (kind, val) in fills:
        if kind == "overlap" or (kind == "fill" and not val):
            return "overlap", (val if kind == "overlap" else 0), (
                f"the contigs overlap along a unitig path ({len(nodes)} unitigs)" if kind == "overlap" and val
                else f"the contigs abut along a unitig path ({len(nodes)} unitigs)")
    # next to a circular contig (hifiasm's 'c': a tandem repeat unit), a fill is one more unit
    # of the repeat, with no way to know how many there are (Troglobisium: the 124 kb fills at
    # H1_scaffold_9 and H2_scaffold_2 and the one at H1_scaffold_46 are also in this haplotype)
    circ = [c for c in (up.left, up.right) if c.endswith("c")]
    if circ:
        return None, None, (f"next to the circular contig {circ[0]} (a tandem repeat unit): the number of units is "
                            "not known, not filled")
    # 2. a path through this haplotype's own contigs: its sequence for the gap is already in
    # the assembly; the contig holding it belongs here, not a copy
    for bp, nodes in up.alts:
        inner = [u for u, _ in nodes[1:-1]]
        if own_share(inner) > own_max:
            where = collections.Counter(graph.read_ctg[rd] for u in inner for rd in utg.reads.get(u, [])
                                        if rd in graph.read_ctg and graph.read_ctg[rd] not in flank
                                        and _hap_of(layout, graph.read_ctg[rd]) == hap)
            ctgs = ", ".join(f"{c} ({layout.place[c].scaffold})" if c in layout.place else c
                             for c, _ in where.most_common(2))
            note = (f"this haplotype's own sequence for the gap is in {ctgs} (a path through it, "
                    f"{bp / 1e3:.0f} kb): place that contig here rather than copy")
            # placeable: the leading contig lies wholly on the path and is a scaffold of its own
            c = where.most_common(1)[0][0]
            oc = _along(utg, graph, c, nodes)
            pc = layout.place.get(c)
            tangle = len(up.alts) >= 20                 # the path search's cap: a repeat tangle
            if oc and pc is not None and len(layout.scaffolds[pc.scaffold]) == 1 and not tangle:
                return "place", (c, oc), note
            why = ("it is not wholly on the path" if not oc else
                   "a repeat tangle (20 or more paths)" if tangle else
                   f"{pc.scaffold if pc else c} holds other contigs too")
            return None, None, note + f" (not placed automatically: {why})"
    # 3. collapsed sequence (two copies in the genome, one in the assembly): fill
    for bp, nodes, (kind, val) in fills:
        inner = [u for u, _ in nodes[1:-1]]
        if not inner:
            if kind == "fill":
                # the bases between the two contig ends come from the end unitigs themselves,
                # whose reads include the contigs' own: judge them by depth alone
                cls, d = utg.depth_class(sorted({nodes[0][0], nodes[-1][0]}))
                if cls == "collapsed":
                    return "fill", val, (f"direct join along the unitig graph; {len(val):,} bp from the end unitigs, "
                                         f"collapsed (depth {d:.0f})")
                tried.append(f"direct join with {len(val):,} bp from end unitigs of {cls} depth ({d:.0f})")
                continue
            tried.append(f"direct join: {val}")
            continue
        cls, d = utg.depth_class(inner)
        share = own_share(inner)
        if cls != "collapsed":
            tried.append(f"{bp / 1e3:.0f} kb {cls} (depth {d:.0f})")
            continue
        if share > own_max:
            tried.append(f"{bp / 1e3:.0f} kb collapsed but {100 * share:.0f}% of its reads in this haplotype's contigs")
            continue
        if kind == "fill":
            return "fill", val, (f"collapsed homozygous sequence (depth {d:.0f}, haploid {utg.haploid:.0f}), "
                                 f"{100 * share:.0f}% of its reads in this haplotype's contigs; "
                                 f"{len(up.alts)} path(s) in all")
        tried.append(f"{bp / 1e3:.0f} kb collapsed: {val}")
    return None, None, "; ".join(dict.fromkeys(tried)) or "no usable path"


def path_fill(fa, layout, utg, up, anchor=1000):
    """What lies between the two contig ends along a unitig path: ('fill', sequence),
    ('overlap', bp the right contig repeats) or (None, why). The contig ends (`anchor` bp)
    must each occur exactly once in the path's sequence."""
    if not utg.has_seq:
        return None, "the unitig GFA has no sequence"
    s = utg.path_seq(up.path)
    ta = contig_end(fa, layout.place[up.left], up.left_o, "tail", anchor)
    hb = contig_end(fa, layout.place[up.right], up.right_o, "head", anchor)
    pa, pb = s.find(ta), s.find(hb)
    if pa < 0 or pb < 0:
        return None, "a contig end is not in the path's sequence"
    if pb < pa:
        # the right contig starts before the left one ends: an overlap longer than the anchor,
        # kept only if the two contigs really share those bases
        L = pa + len(ta) - pb
        if contig_end(fa, layout.place[up.left], up.left_o, "tail", L) == \
                contig_end(fa, layout.place[up.right], up.right_o, "head", L):
            return "overlap", L
        return None, "the contig ends are in the wrong order along the path"
    if s.find(ta, pa + 1) >= 0 or s.find(hb, pb + 1) >= 0:
        return None, "a contig end occurs more than once in the path's sequence"
    gap = pb - (pa + len(ta))
    return ("overlap", -gap) if gap < 0 else ("fill", s[pa + len(ta):pb])


# --------------------------------------------------------------------------- evidence

def contig_end(fa, place, orient, side, n):
    """The n bp at the 'head' or 'tail' of a contig traversed in `orient` (relative to
    hifiasm's contig), read from the assembly through its placement."""
    p = place
    # head / tail of the contig as hifiasm wrote it (+), in assembly coordinates
    if (side == "tail") == (orient == "+"):        # tail of +, or head of - (= revcomp of + tail)
        s = fa.fetch(p.scaffold, p.end - n, p.end) if p.orient == "+" else revcomp(fa.fetch(p.scaffold, p.start, p.start + n))
        return s if orient == "+" else revcomp(s)
    s = fa.fetch(p.scaffold, p.start, p.start + n) if p.orient == "+" else revcomp(fa.fetch(p.scaffold, p.end - n, p.end))
    return s if orient == "+" else revcomp(s)


def verify_overlap(fa, layout, a, oa, b, ob, ov, slack=2000):
    """Does the end of a(oa) equal the start of b(ob) over the link's overlap? Returns the
    exact overlap length (searched within +-slack bp of the graph's, as a suffix-prefix
    match) or None. hifiasm's overlap comes from reads and may differ from the contigs'
    consensus by a little."""
    pa, pb = layout.place.get(a), layout.place.get(b)
    if pa is None or pb is None:
        return None
    n = ov + slack
    ta = contig_end(fa, pa, oa, "tail", n)
    hb = contig_end(fa, pb, ob, "head", n)
    if ta[-ov:] == hb[:ov]:
        return ov
    for d in range(1, slack + 1):
        for L in (ov - d, ov + d):
            if 0 < L <= min(len(ta), len(hb)) and ta[-L:] == hb[:L]:
                return L
    return None


class Junction(NamedTuple):
    """Graph evidence at a YaHS gap."""
    kind: str           # 'link' or 'bridge'
    scaffold: str
    gap_start: int      # bp, 0-based, in the assembly scaffold
    gap_end: int
    left: str           # contig names and orientations as laid out
    left_o: str
    right: str
    right_o: str
    overlap: int        # link: a-b; bridge: a-c
    via: str = ""       # bridge: the contig that fills the gap
    via_o: str = ""
    overlap2: int = 0   # bridge: c-b
    alternatives: int = 0   # other links at the ends involved
    via_scaffold: str = ""  # bridge: the assembly scaffold that is the bridging contig


def junctions(layout, graph):
    """Graph evidence at every YaHS gap: a direct link, or a bridge through a contig that is a
    scaffold of its own (see module docstring)."""
    out = []
    for sc, a, oa, b, ob, g0, g1 in layout.neighbours():
        ov = graph.link(a, oa, b, ob)
        if ov is not None:
            alt = len(graph.out.get((a, oa), ())) + len(graph.into(b, ob)) - 2
            out.append(Junction("link", sc, g0, g1, a, oa, b, ob, ov, alternatives=alt))
            continue
        nxt, prv = graph.out.get((a, oa), set()), graph.into(b, ob)
        if len(nxt) != 1 or len(prv) != 1:
            continue
        (c, oc, ov1), = nxt
        (c2, oc2, ov2), = prv
        if (c, oc) != (c2, oc2) or c in (a, b):
            continue
        pc = layout.place.get(c)
        if pc is None or len(layout.scaffolds[pc.scaffold]) != 1:
            continue                    # the bridge must be free to move: a scaffold of its own
        alt = len(graph.out.get((c, oc), ())) + len(graph.into(c, oc)) - 2
        out.append(Junction("bridge", sc, g0, g1, a, oa, b, ob, ov1, c, oc, ov2, alt, pc.scaffold))
    return out


def load(gfa, yahs_agp, fasta=None, names=None):
    """Junctions from the given files, each with its overlap(s) as found in the sequence
    (with `fasta`; None where not checked, 0 where the sequence disagrees):
    [(Junction, (overlap, overlap2))]. Raises ValueError if the AGPs do not fit the assembly
    (`names`: its scaffold names; with `fasta`, lengths are checked too)."""
    lay = Layout(yahs_agp, names=set(names) if names is not None else None)
    fa = Faidx(fasta) if fasta else None
    if fa is not None and not fa.gz:
        bad = lay.check({n: fa.length(n) for n in fa.idx})
    elif names is not None:
        bad = [f"{s} not in the assembly" for s in lay.scaffolds if s not in set(names)]
    else:
        bad = []
    if bad:
        raise ValueError("the YaHS AGPs do not match the assembly: " + "; ".join(bad[:5]))
    gr = Graph(gfa)
    out = []
    for j in junctions(lay, gr):
        if fa is None:
            out.append((j, None))
        elif j.kind == "link":
            out.append((j, (verify_overlap(fa, lay, j.left, j.left_o, j.right, j.right_o, j.overlap) or 0, 0)))
        else:
            out.append((j, (verify_overlap(fa, lay, j.left, j.left_o, j.via, j.via_o, j.overlap) or 0,
                            verify_overlap(fa, lay, j.via, j.via_o, j.right, j.right_o, j.overlap2) or 0)))
    return out, lay, gr


# --------------------------------------------------------------------------- joins at scaffold ends

def outward_ends(layout):
    """{(contig, orientation)}: leaving the contig that way leaves its YaHS scaffold (the
    scaffold's two outer ends)."""
    E = set()
    for lst in layout.scaffolds.values():
        first, last = layout.place[lst[0]], layout.place[lst[-1]]
        E.add((lst[-1], last.orient))
        E.add((lst[0], FLIP[first.orient]))
    return E


def end_joins(layout, graph, utg=None, limit=1_000_000, min_big=1_000_000, margin=3.0):
    """Joins only the graph sees: two outer scaffold ends of one haplotype that the contig
    graph links (the only link at each end), or, with `utg`, that are each other's nearest
    outer end along the unitig graph with no other within `margin` times the distance. At
    least one of the two scaffolds must be `min_big` bp or more (the join extends a
    chromosome). Returns [dict(kind, a, oa, b, ob, dist, overlap, sa, sb, size_a, size_b)]:
    leaving contig a as oa joins contig b left as ob (both written as exits)."""
    E = outward_ends(layout)
    size = {sc: layout.place[lst[-1]].end for sc, lst in layout.scaffolds.items()}
    hap = lambda c: _hap_of(layout, c)
    out, seen = [], set()

    def add(kind, a, oa, b, ob, dist, ov):
        k = frozenset([(a, oa), (b, ob)])
        sa, sb = layout.place[a].scaffold, layout.place[b].scaffold
        if k in seen or sa == sb or max(size[sa], size[sb]) < min_big or hap(a) != hap(b):
            return
        seen.add(k)
        out.append(dict(kind=kind, a=a, oa=oa, b=b, ob=ob, dist=dist, overlap=ov, sa=sa, sb=sb,
                        size_a=size[sa], size_b=size[sb]))
    for (c, o) in sorted(E):
        nxt = graph.out.get((c, o), set())
        if len(nxt) != 1:
            continue
        (d, od, ov), = nxt
        if (d, FLIP[od]) in E and len(graph.into(d, od)) == 1:
            add("link", c, o, d, FLIP[od], 0, ov)
    if utg is None:
        return out
    import heapq
    entry = collections.defaultdict(list)
    for (c, o) in E:
        n = utg.end_node(graph, c, "t" if o == "+" else "h")
        if n is not None:
            ex = n if o == "+" else (n[0], FLIP[n[1]])
            entry[(ex[0], FLIP[ex[1]])].append((c, FLIP[o]))
    near = {}
    for (c, o) in sorted(E):
        n = utg.end_node(graph, c, "t" if o == "+" else "h")
        if n is None:
            continue
        ex = n if o == "+" else (n[0], FLIP[n[1]])
        seen_, pq, found = {ex: 0}, [(0, ex)], []
        while pq and len(found) < 2:
            dd, x = heapq.heappop(pq)
            if dd > seen_.get(x, 1e18) or dd > limit:
                continue
            for d, od in entry.get(x, ()):
                if d != c and hap(d) == hap(c) and x != ex and all(f[1] != d for f in found):
                    found.append((dd, d, od))
            for m, _ in utg.adj.get(x, ()):
                nd = dd + utg.length[m[0]]
                if nd < seen_.get(m, 1e18):
                    seen_[m] = nd
                    heapq.heappush(pq, (nd, m))
        if found:
            dd, d, od = found[0]
            second = found[1][0] if len(found) > 1 else None
            if second is None or second > margin * max(dd, 50000):
                near[(c, o)] = ((d, FLIP[od]), dd)
    for (c, o), ((d, od), dd) in near.items():
        if near.get((d, od), (None,))[0] == (c, o):          # mutual
            add("path", c, o, d, od, dd, 0)
    return out


def layout_end_joins(pt, state, layout):
    """The joins between contig ends in the layout a savestate gives (as PretextView would
    export it): {frozenset({(contig, exit orientation), (contig, exit orientation)})}."""
    import tempfile
    from .layout import write_agp
    fd, tmp = tempfile.mkstemp(suffix=".agp")
    os.close(fd)
    try:
        write_agp(pt, state, tmp)
        _, objs = read_agp(tmp)
    finally:
        os.remove(tmp)
    tol = 1.5 * pt.bp_per_pix
    out = set()
    for rows in objs.values():
        w = [r for r in rows if r.kind == "W"]
        for x_, y_ in zip(w, w[1:]):
            x, y = _snap_end(layout, x_, "right", tol), _snap_end(layout, y_, "left", tol)
            if x and y:
                out.add(frozenset([(x[0], x[1]), (y[0], FLIP[y[1]])]))
    return out


def join_rows(pt, state, layout, joins):
    """Review rows (`graph_join_candidate`) for the graph joins (`end_joins`) the layout of
    `state` does not make, with the Hi-C end contact of the two scaffold ends against the
    best other end of a large scaffold for the smaller one, and where the layout put the
    smaller scaffold."""
    from . import allhaps
    from .mapstate import MapState
    made = layout_end_joins(pt, state, layout)
    sg = allhaps._signal(pt)
    g = sg.decay_long(230)
    idx = {n: i for i, n in enumerate(pt.names)}
    npx = lambda sc: int(pt.scaffold_npix[idx[sc]]) if sc in idx else 0
    rng = lambda sc: (int(pt.scaffold_start[idx[sc]]), int(pt.scaffold_start[idx[sc]] + npx(sc) - 1))

    def contact(c, o, d, od):
        A, B = layout.place[c].scaffold, layout.place[d].scaffold
        if npx(A) < 2 or npx(B) < 2:
            return None
        return allhaps._adjacency(sg, g, rng(A), o != layout.place[c].orient, rng(B), od != layout.place[d].orient)
    big_ends = [(c, o) for (c, o) in outward_ends(layout) if npx(layout.place[c].scaffold) * pt.bp_per_pix >= 1e6]
    ms = MapState(pt).replay(state.edits)
    paint = dict(state.scaffs)
    painted = {}
    for i, c in enumerate(ms.contigs):
        if paint.get(i):
            painted.setdefault(pt.names[c.scaffold], paint[i])
    rows = []
    for j in joins:
        k = frozenset([(j["a"], j["oa"]), (j["b"], j["ob"])])
        if k in made:
            continue
        big, small = (("a", "b") if j["size_a"] >= j["size_b"] else ("b", "a"))
        c, o, d, od = j[big], j["o" + big], j[small], j["o" + small]
        v = contact(c, o, d, FLIP[od])
        alts = [x for e in big_ends if e[0] not in (c, d)
                for x in [contact(e[0], e[1], d, FLIP[od])] if x is not None]
        best = max(alts) if alts else None
        if v is None:
            hic = "too small for Hi-C (under 2 px)"
        else:
            hic = (f"Hi-C end contact {v:.2f}" + (f", the best of any large scaffold end (next {best:.2f})"
                                                  if best is None or v > best else
                                                  f"; another large scaffold end touches it more ({best:.2f})"))
        what = (f"hifiasm's contig graph links {j['a']}{j['oa']} to {j['b']}{FLIP[j['ob']]} (overlap {j['overlap']:,} bp; "
                "the only link at both ends)" if j["kind"] == "link" else
                f"a unitig path of {j['dist'] / 1e3:.0f} kb joins {j['a']}{j['oa']} to {j['b']}{FLIP[j['ob']]}, each the "
                "other's nearest scaffold end with no close alternative")
        where = ("; the tool put " + j["s" + small] + " in a chromosome elsewhere" if painted.get(j["s" + small])
                 and painted.get(j["s" + small]) != painted.get(j["s" + big]) else
                 "; " + j["s" + small] + " is not placed" if not painted.get(j["s" + small]) else "")
        rows.append(dict(action="graph_join_candidate", scaffold=j["s" + big], length_mb=round(j["size_" + big] / 1e6, 2),
                         partner=j["s" + small], strand="", tag="", score=round(v, 2) if v is not None else "",
                         next_best=round(best, 2) if best is not None else "",
                         evidence=f"{what}: the end of {j['s' + big]} ({j['size_' + big] / 1e6:.1f} Mb) continues into "
                                  f"{j['s' + small]} ({j['size_' + small] / 1e6:.2f} Mb), a join no layout makes; "
                                  f"{hic}{where}; not made (review)"))
    return rows


# --------------------------------------------------------------------------- AGP editing

class AgpRow(NamedTuple):
    obj: str
    kind: str           # 'W' (component) or a gap type ('U', 'N')
    name: str           # component name; gap length (as text) for gaps
    beg: int            # component range, 1-based inclusive (0 for gaps)
    end: int
    orient: str
    extra: tuple        # gap columns 7-9, or component columns 10+ (tags)


def read_agp(path):
    """(header lines, {object: [AgpRow]} in file order)."""
    head, objs = [], {}
    with _open(path) as f:
        for line in f:
            if line.startswith("#"):
                head.append(line)
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 9:
                continue
            if c[4] in ("U", "N"):
                r = AgpRow(c[0], c[4], c[5], 0, 0, "", tuple(c[6:]))
            else:
                r = AgpRow(c[0], "W", c[5], int(c[6]), int(c[7]), c[8], tuple(c[9:]))
            objs.setdefault(c[0], []).append(r)
    return head, objs


def write_agp_rows(path, head, objs):
    """Write objects with object coordinates and part numbers recomputed."""
    with open(path, "w") as f:
        f.writelines(head)
        for obj, rows in objs.items():
            pos = 1
            for k, r in enumerate(rows, 1):
                ln = int(r.name) if r.kind != "W" else r.end - r.beg + 1
                if r.kind == "W":
                    cols = [obj, pos, pos + ln - 1, k, "W", r.name, r.beg, r.end, r.orient, *r.extra]
                else:
                    cols = [obj, pos, pos + ln - 1, k, r.kind, r.name, *r.extra]
                f.write("\t".join(map(str, cols)) + "\n")
                pos += ln


def _scaffold_orient(layout, contig, o):
    """Orientation of a contig's assembly scaffold that traverses the contig as `o`."""
    return o if layout.place[contig].orient == "+" else FLIP[o]


def insert_bridges(agp_in, agp_out, juncs, layout):
    """Put each bridge contig into its gap in a PretextView-style AGP (before pretext-to-asm):
    the line holding the gap is split there, and the bridge's own line (a scaffold of its
    own) moves in between, with the host's tags. Returns [(junction, done, note)]."""
    head, objs = read_agp(agp_in)
    report = []
    for j in [j for j in juncs if j.kind == "bridge"]:
        pc = layout.place[j.via]
        # the bridge must stand alone: the only component of its object
        home = [o for o, rows in objs.items() if any(r.kind == "W" and r.name == pc.scaffold for r in rows)]
        if len(home) != 1 or len([r for r in objs[home[0]] if r.kind == "W"]) != 1:
            report.append((j, False, f"{pc.scaffold} is not a scaffold of its own in the layout (placed by the curator?)"))
            continue
        host = None
        for o, rows in objs.items():
            for i, r in enumerate(rows):
                if r.kind == "W" and r.name == j.scaffold and r.beg <= j.gap_start and j.gap_end + 1 <= r.end:
                    host = (o, i)
        if host is None:
            report.append((j, False, f"the gap at {j.scaffold}:{j.gap_start + 1} is not inside one piece of the layout"))
            continue
        o, i = host
        r = objs[o][i]
        left = r._replace(end=j.gap_start)                 # sequence before the gap (1-based end)
        right = r._replace(beg=j.gap_end + 1)              # after it
        so = _scaffold_orient(layout, j.via, j.via_o)      # bridge scaffold orientation, host forward
        mid = AgpRow(o, "W", pc.scaffold, 1, pc.end, so if r.orient == "+" else FLIP[so], r.extra)
        gap = AgpRow(o, "U", "100", 0, 0, "", ("scaffold", "yes", "proximity_ligation"))
        new = [left, gap, mid, gap, right] if r.orient == "+" else [right, gap, mid, gap, left]
        objs[o][i:i + 1] = new
        del objs[home[0]]
        report.append((j, True, f"{pc.scaffold} ({pc.end} bp) put into the gap of {j.scaffold} at {j.gap_start + 1}"))
    write_agp_rows(agp_out, head, objs)
    return report


def _snap_end(layout, row, side, tol):
    """As `_end_contig`, for a PretextView-style line whose coordinates come from pixels:
    the contig end nearest the piece's end, within `tol` bp; (contig, orientation) or None."""
    lst = layout.scaffolds.get(row.name)
    if not lst:
        return None
    at_end = (side == "right") == (row.orient == "+")
    pos = row.end if at_end else row.beg - 1
    best = None
    for c in lst:
        p = layout.place[c]
        d = abs((p.end if at_end else p.start) - pos)
        if d <= tol and (best is None or d < best[0]):
            best = (d, c, p)
    if best is None:
        return None
    _, c, p = best
    return c, (p.orient if row.orient == "+" else FLIP[p.orient])


def place_between(agp_path, layout, graph, utg, fa):
    """Contigs of the gap's own haplotype that lie wholly on a unitig path across a junction
    between two pieces of a PretextView-style AGP (joins the tool or the curator made, not
    YaHS gaps): the contig's own line (a scaffold of its own) moves between them, with the
    host's tags. Edits the file in place; returns [(object, left, right, contig, note)]."""
    head, objs = read_agp(agp_path)
    m = next((re.search(r"RESOLUTION: ([\d.]+)", h) for h in head if "RESOLUTION" in h), None)
    tol = 1.5 * float(m.group(1)) if m else 0
    lone = {rows[0].name: o for o, rows in objs.items() if len(rows) == 1 and rows[0].kind == "W"}
    todo = []
    for o, rows in objs.items():
        for i in range(len(rows) - 2):
            a, g, b = rows[i], rows[i + 1], rows[i + 2]
            if a.kind != "W" or g.kind == "W" or b.kind != "W":
                continue
            x, y = _snap_end(layout, a, "right", tol), _snap_end(layout, b, "left", tol)
            if not (x and y) or _hap_of(layout, x[0]) != _hap_of(layout, y[0]):
                continue
            ups = unitig_paths(layout, graph, utg, pairs=[(o, x[0], x[1], y[0], y[1], 0, 0)])
            if not ups:
                continue
            kind, val, note = choose_fill(fa, layout, graph, utg, ups[0])
            if kind == "place" and layout.place[val[0]].scaffold in lone:
                todo.append((o, i, x, y, val, note))
    named = collections.Counter(v[0] for *_, v, _ in todo)
    report, shift = [], collections.Counter()
    for o, i, x, y, (c, oc), note in todo:
        if named[c] != 1:
            continue
        pc = layout.place[c]
        home = lone[pc.scaffold]
        host = objs[o][i + shift[o]]
        so = _scaffold_orient(layout, c, oc)
        mid = AgpRow(o, "W", pc.scaffold, 1, pc.end, so, host.extra)
        gap = AgpRow(o, "U", "100", 0, 0, "", ("scaffold", "yes", "proximity_ligation"))
        k = i + shift[o] + 1                        # after the left piece: gap, contig, (old gap), right
        objs[o][k:k] = [gap, mid]
        shift[o] += 2
        del objs[home]
        report.append((o, f"{x[0]}{x[1]}", f"{y[0]}{y[1]}", c, note))
    if report:
        write_agp_rows(agp_path, head, objs)
    return report


def _end_contig(layout, row, side):
    """The contig at the 'left' or 'right' end of a component as laid out, and how it is
    traversed there (relative to hifiasm's contig): (contig, orientation), or None if the
    component does not end at a contig end."""
    lst = layout.scaffolds.get(row.name)
    if not lst:
        return None
    # the scaffold coordinate at that end: the component's end for its right side when
    # forward (or its left side when reversed), else its start
    at_end = (side == "right") == (row.orient == "+")
    for c in lst:
        p = layout.place[c]
        if (at_end and p.end == row.end) or (not at_end and p.start + 1 == row.beg):
            return c, (p.orient if row.orient == "+" else FLIP[p.orient])
    return None


def close_gaps(agp_path, fa_path, fa_in, layout, graph, report_rows, utg=None, fill=False, fills=None,
               limit=1_000_000):
    """Close gaps of a pretext-to-asm output (AGP + FASTA) whose two sides are contig ends
    that belong together in hifiasm's graphs:
      - linked in the contig graph: the right-hand component loses the overlap (checked in
        the sequence) and the gap line goes;
      - with `utg` (unitig graph with sequence), joined by a single unitig path along which
        the two contigs overlap: the same, with the overlap the path shows;
      - with `fill` as well, joined by a single unitig path with sequence between them: the
        gap line is replaced by that sequence, a component named gapfill_N whose sequence is
        added to `fills` ({name: sequence}, to be written next to the AGP).
    Edited objects are rebuilt from the input assembly (and `fills`); the others are copied.
    Appends to `report_rows`; returns the number of gaps closed or filled."""
    head, objs = read_agp(agp_path)
    fills = {} if fills is None else fills
    changed, closed = {}, 0
    base = os.path.basename(agp_path)
    for obj, rows in objs.items():
        i, edited = 1, False
        while i < len(rows) - 1:
            a, g, b = rows[i - 1], rows[i], rows[i + 1]
            if not (a.kind == "W" and g.kind != "W" and b.kind == "W"):
                i += 1
                continue
            x, y = _end_contig(layout, a, "right"), _end_contig(layout, b, "left")
            if not (x and y):
                i += 1
                continue
            rep = dict(file=base, object=obj, left=f"{x[0]}{x[1]}", right=f"{y[0]}{y[1]}")
            ov = graph.link(x[0], x[1], y[0], y[1])
            L, new = None, None
            if ov is not None:
                L = verify_overlap(fa_in, layout, x[0], x[1], y[0], y[1], ov)
                report_rows.append(dict(rep, graph_overlap=ov, overlap=L or "",
                                        action="gap_closed" if L else "kept: the overlap is not in the sequence"))
            elif utg is not None:
                ups = unitig_paths(layout, graph, utg, limit, pairs=[(obj, x[0], x[1], y[0], y[1], 0, 0)])
                if ups:
                    up = ups[0]
                    kind, val, note = choose_fill(fa_in, layout, graph, utg, up)
                    if kind == "overlap" and val:
                        L = val
                        report_rows.append(dict(rep, graph_overlap="", overlap=L, action=f"gap_closed: {note}"))
                    elif kind == "overlap" or (kind == "fill" and fill):
                        if val:
                            name = f"gapfill_{len(fills) + 1}"
                            fills[name] = val
                            new = [a, AgpRow(obj, "W", name, 1, len(val), "+", ()), b]
                            what = f"gap_filled: {name}, {len(val):,} bp; {note}"
                        else:
                            new = [a, b]
                            what = f"gap_closed: {note} (the contigs abut)"
                        report_rows.append(dict(rep, graph_overlap="", overlap=f"fill {len(val) if val else 0}",
                                                action=what))
                    elif kind == "fill":
                        report_rows.append(dict(rep, graph_overlap="", overlap=f"fill {len(val)}",
                                                action=f"kept: {len(val):,} bp could fill it (--fill-gaps); {note}"))
                    else:
                        report_rows.append(dict(rep, graph_overlap="", overlap="",
                                                action=f"kept, not filled: {note}"))
            if L:
                b2 = b._replace(beg=b.beg + L) if b.orient == "+" else b._replace(end=b.end - L)
                new = [a, b2]
            if new is not None:
                rows[i - 1:i + 2] = new
                edited, closed = True, closed + 1
                i += len(new) - 1
                continue
            i += 1
        if edited:
            changed[obj] = rows
    if not changed:
        return 0
    # rebuild the edited objects' sequence from the input assembly; copy the rest
    width = 60
    with _open(fa_path) as f:
        for line in f:
            if not line.startswith(">"):
                width = len(line.rstrip("\n")) or 60
                break
    tmp_fa, tmp_agp = fa_path + ".tmp", agp_path + ".tmp"
    opener = (lambda p: gzip.open(p, "wt")) if fa_path.endswith(".gz") else (lambda p: open(p, "w"))
    with _open(fa_path) as fin, opener(tmp_fa) as fout:
        skip = False
        for line in fin:
            if line.startswith(">"):
                name = line[1:].split()[0]
                skip = name in changed
                fout.write(line)
                if skip:
                    seq = []
                    for r in changed[name]:
                        if r.kind == "W":
                            s = (fills[r.name][r.beg - 1:r.end] if r.name in fills
                                 else fa_in.fetch(r.name, r.beg - 1, r.end, upper=False))
                            seq.append(s if r.orient == "+" else revcomp(s))
                        else:
                            seq.append("N" * int(r.name))
                    s = "".join(seq)
                    fout.writelines(s[k:k + width] + "\n" for k in range(0, len(s), width))
            elif not skip:
                fout.write(line)
    write_agp_rows(tmp_agp, head, objs)
    os.replace(tmp_fa, fa_path)
    os.replace(tmp_agp, agp_path)
    return closed
