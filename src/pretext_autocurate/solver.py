"""Adjacency solver: an alternative to the rule-based allhaps pass (`allhaps.curate_allhaps`),
built alongside it (`run --mode allhaps --engine adjacency`) to be compared on the same genomes.

The whole curation is one decision: which contig ends sit next to each other.

1. Nodes are contigs: the stretches of each scaffold between gap pixels (the gap track). Cuts
   happen only at gaps, as the curator asks.
2. Every pair of contig ends that might be adjacent is an edge, scored as a sum of named parts:
     yahs      the assembly put them together (a prior; `--level` sets how strong)
     hic       end-to-end contact against what adjacent sequence shows (log2, from one
               binned contact matrix), plus how specific it is against the next best end
     homology  the two ends land next to each other on the other haplotype, facing the right way
               (+), or the other haplotype says they are far apart (-)
     graph     hifiasm's contig graph links them (a YaHS join it backs is fixed; a join only the
               graph sees gets a bonus)
     telomere  a telomere at an end counts against joining it (supporting evidence only)
3. Paths: edges are taken best first, each end used once, no cycles, a YaHS join only above 0
   and a new join only above a level-dependent bar. A YaHS join not taken is a cut; an edge
   that is not a YaHS join is a join.
4. The paths of the two haplotypes are paired by sequence alignment and laid out side by side,
   the second turned to run the same way; small leftovers go to the end of the chromosome
   their contact (or alignment) points to, tagged Unloc.

Every cut and join is a row with its score parts, so each edit traces to the evidence behind
it. Returns what `curate_allhaps` returns: (rows, savestate, info)."""
import bisect
import collections
import math
import re

import numpy as np

from .allhaps.common import _signal
from .mapstate import MapState, plan_edits
from .savestate import BREAK_DELTA, Edit, SaveState

# Weights at levels 0, 0.5, 1 (interpolated): (prior for a YaHS join, bar for a new join)
LEVELS = {"prior": (2.5, 1.5, 0.75), "join_min": (3.0, 2.0, 1.2), "specific_min": (2.0, 1.5, 1.25)}
# R_REF: contact ratio halfway (in log) between background (~0.06 on the toy and real maps) and
# adjacent sequence (~1), so a misjoin and a true join score about -3 and +3 with W_HIC
W_HIC, W_SPEC, W_HOM, W_GRAPH, W_TEL, R_REF = 1.5, 1.5, 2.0, 3.0, 1.0, 0.25
# a telomere facing outwards at an end (strand tracks): that end is a chromosome end
W_TEL_OUT = 2.5
# homology through two contigs of the other haplotype that YaHS joined is only as good as that
# join (the other haplotype was scaffolded with the same Hi-C and can share a misjoin)
HOM_ACROSS = 0.25


def _at(name, level):
    v0, v5, v1 = LEVELS[name]
    return v0 + (v5 - v0) * level / 0.5 if level <= 0.5 else v5 + (v1 - v5) * (level - 0.5) / 0.5


class Node:
    __slots__ = ("i", "s", "a", "b", "hap", "hic")

    def __init__(self, i, s, a, b, hap):
        self.i, self.s, self.a, self.b, self.hap, self.hic = i, s, a, b, hap, False

    @property
    def n(self):
        return self.b - self.a + 1


def contigs(pt, hap_regex, f=8, min_bins=2, min_px=2, breaks=()):
    """Nodes (contigs between gap pixels, in map order) and the YaHS joins between them:
    [(left node, right node, pixel q: the cut would go after q)]. `breaks`: further pixels to
    split after, inside contigs (sharp breaks: contact at background with no gap)."""
    brk = {int(q) for q in breaks}
    gap = pt.graphs.get("gap")
    gap = np.zeros(pt.npix) if gap is None else gap
    nodes, yahs = [], []
    for s in np.argsort(pt.scaffold_start):
        s, n = int(s), int(pt.scaffold_npix[s])
        if n == 0:
            continue
        a0 = int(pt.scaffold_start[s])
        m = re.match(hap_regex, pt.names[s])
        hap = m.group(1) if m else ""
        g = gap[a0:a0 + n] > 0
        ends, k = [], 0
        while k < n:                      # a run of gap pixels: the piece ends at its last pixel
            if g[k]:
                j = k
                while j + 1 < n and g[j + 1]:
                    j += 1
                if 0 < k and j < n - 1:
                    ends.append(a0 + j)
                k = j + 1
            else:
                k += 1
        ends = sorted(set(ends) | {q for q in brk if a0 <= q < a0 + n - 1})
        bounds = [a0 - 1] + ends + [a0 + n - 1]
        pieces = [[bounds[i] + 1, bounds[i + 1]] for i in range(len(bounds) - 1)]
        # a piece too short to judge by contact (under min_bins whole bins of f px) joins the piece
        # before it (or after, if first): its gaps cannot be cut on evidence, and left alone it
        # would only be squeezed out by a join between its neighbours (icLarUrsu H1_scaffold_5
        # 47.36-48.13 Mb)
        ok = lambda p: (p[1] + 1) // f - 1 - (-(-p[0] // f)) + 1 >= min_bins and p[1] - p[0] + 1 >= min_px
        merged = []
        for p in pieces:
            if merged and (not ok(p) or not ok(merged[-1])):
                merged[-1][1] = p[1]
            else:
                merged.append(p)
        first = len(nodes)
        for a, b in merged:
            nodes.append(Node(len(nodes), s, a, b, hap))
        for k in range(first, len(nodes) - 1):
            yahs.append((k, k + 1, nodes[k].b))
    return nodes, yahs


def sharp_breaks(pt, level, end_px=10):
    """Pixels inside contigs where contact drops to background on both scales (the rules pass's
    sharp breaks, `misjoin.detect_misjoins` at the level's `sharp_max`); dense maps only."""
    from .allhaps.common import Thresholds
    from .misjoin import Params, detect_misjoins
    sg = _signal(pt)
    if sg.decay[min(20, len(sg.decay) - 1)] < 0.25 * sg.decay[0]:
        return []
    smax = Thresholds(level).sharp_max
    out = []
    for c in detect_misjoins(pt, Params(cut_max=smax), sg):
        q = int(c.after_px)
        s = int(pt.orig_scaffold[q])
        a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        if c.apply and not c.at_gap and c.far < smax and q - a >= end_px and a + n - 1 - q >= end_px:
            out.append(q)
    return out


def _decay_bins(M, sc, D):
    g = np.zeros(D)
    n = M.shape[0]
    for d in range(D):
        v = np.diagonal(M, d)
        same = sc[:n - d] == sc[d:]
        g[d] = float(v[same].mean()) if same.any() else 0.0
    return g


def end_contact(pt, nodes, f=8, wmax=8):
    """Contact between every pair of contig ends, as observed / expected for adjacent ends
    (E x E, E = 2 x nodes; end 2i = head of node i, 2i+1 = tail). Ends of nodes under two bins
    have no usable window (nan); `node.hic` says whether a node has one."""
    sg = _signal(pt)
    while pt.tex_res % f:
        f //= 2
    M = sg.binned(f)
    nb = M.shape[0]
    sc = pt.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, pt.npix - 1)]
    g = _decay_bins(M, sc, 2 * wmax + 2)
    E = 2 * len(nodes)
    W = np.zeros((nb, E), dtype=np.float32)
    win = np.zeros(E, dtype=int)
    for v in nodes:
        fb, lb = -(-v.a // f), (v.b + 1) // f - 1            # bins wholly inside the node
        nbins = lb - fb + 1
        if nbins < 2:
            continue
        v.hic = True
        w = min(wmax, nbins // 2)
        W[fb:fb + w, 2 * v.i] = 1.0
        W[lb - w + 1:lb + 1, 2 * v.i + 1] = 1.0
        win[2 * v.i] = win[2 * v.i + 1] = w
    O = W.T @ (M @ W)
    i = np.arange(2 * wmax)
    cache = {}

    def expect(wa, wb):
        if (wa, wb) not in cache:
            cache[(wa, wb)] = float(g[i[:wa, None] + i[None, :wb] + 1].sum())
        return cache[(wa, wb)]
    R = np.full((E, E), np.nan)
    ok = np.nonzero(win)[0]
    for x in ok:
        for y in ok:
            if x // 2 != y // 2:
                R[x, y] = O[x, y] / max(expect(win[x], win[y]), 1e-9)
    return R


def _projections(pt, nodes, blocks, min_mapq=30, end_bp=2_000_000):
    """Per end: (other haplotype's scaffold, position bp, direction): where the end lands on the
    other haplotype and which way it moves along it as one walks out of the node through that
    end. Also, per node, aligned bp to each node of the other haplotype."""
    idx = {n: i for i, n in enumerate(pt.names)}
    bpp = pt.bp_per_pix
    by_s = collections.defaultdict(list)                   # scaffold -> nodes (in order)
    for v in nodes:
        by_s[v.s].append(v)
    rec = collections.defaultdict(list)                    # scaffold -> (s0, e0, other, os, oe, same)
    for bk in blocks or []:
        if bk.mapq < min_mapq or bk.q not in idx or bk.t not in idx:
            continue
        rec[idx[bk.q]].append((bk.qs, bk.qe, idx[bk.t], bk.ts, bk.te, bk.same_strand))
        rec[idx[bk.t]].append((bk.ts, bk.te, idx[bk.q], bk.qs, bk.qe, bk.same_strand))

    def node_at(s, bp):
        lst = by_s.get(s)
        if not lst:
            return None
        k = bisect.bisect_right([v.a for v in lst], int(pt.scaffold_start[s] + bp // bpp)) - 1
        return lst[max(k, 0)]
    olen_bp = lambda s: int(pt.scaffold_npix[s]) * bpp

    def landing(O, p, d):
        """(scaffold, position, direction out, contig index, bp of that contig beyond p)."""
        w = node_at(O, min(max(p, 0), olen_bp(O) - 1))
        if w is None:
            return (O, p, d, -1, 0.0)
        a0 = int(pt.scaffold_start[O])
        lo, hi = (w.a - a0) * bpp, (w.b + 1 - a0) * bpp
        return (O, p, d, w.i, max(0.0, (hi - p) if d > 0 else (p - lo)))
    proj, share = {}, collections.defaultdict(lambda: collections.defaultdict(float))
    for v in nodes:
        lo = (v.a - int(pt.scaffold_start[v.s])) * bpp
        hi = (v.b + 1 - int(pt.scaffold_start[v.s])) * bpp
        mine = []
        for s0, e0, o, os_, oe, same in rec.get(v.s, ()):
            a, b = max(s0, lo), min(e0, hi)
            if b <= a:
                continue
            f0, f1 = (a - s0) / (e0 - s0), (b - s0) / (e0 - s0)
            o0, o1 = (os_ + f0 * (oe - os_), os_ + f1 * (oe - os_)) if same else \
                (oe - f1 * (oe - os_), oe - f0 * (oe - os_))
            mine.append((a, b, o, o0, o1, same))
            w = node_at(o, (o0 + o1) / 2)
            if w is not None:
                share[v.i][w.i] += b - a
        if not mine:
            continue
        tot = collections.Counter()
        for a, b, o, *_ in mine:
            tot[o] += b - a
        O = tot.most_common(1)[0][0]
        mo = [m for m in mine if m[2] == O]
        span = min(end_bp, (hi - lo) / 2)
        head = [m for m in mo if m[0] < lo + span]
        tail = [m for m in mo if m[1] > hi - span]
        # the largest block near each end, not the outermost: near a chromosome end the outermost
        # is often a short subtelomeric repeat that also matches another chromosome's end
        # (bAegFun H1_scaffold_1 134.6 Mb: chromosome B's start landed on the end of A's homolog);
        # its projection is carried on to the node's end
        if head:
            a, b, _, o0, o1, same = max(head, key=lambda m: m[1] - m[0])
            d = -1 if same else 1          # walking out through the head goes down O if same strand
            p = (o0 if same else o1) + d * (a - lo)
            proj[(v.i, "h")] = landing(O, p, d)
        if tail:
            a, b, _, o0, o1, same = max(tail, key=lambda m: m[1] - m[0])
            d = 1 if same else -1
            p = (o1 if same else o0) + d * (hi - b)
            proj[(v.i, "t")] = landing(O, p, d)
    return proj, share


def _homology(p, q, olen=None, near=1_000_000, far=5_000_000, slack=500_000):
    """What the other haplotype says about joining ends p and q, from where they land on it
    (scaffold, position, direction out, contig, contig bp beyond it): +1 .. 0 if next to each
    other and facing; -1 if far apart, facing the wrong way, or on different scaffolds while one
    lands well inside its scaffold (its neighbour there is known, and it is not the other end);
    None if it cannot say. Evidence that rests on the other haplotype's contigs counts fully;
    evidence that rests on its YaHS joins (two different contigs, or room only beyond a contig's
    end) counts HOM_ACROSS: that haplotype was scaffolded with the same Hi-C and can share the
    misjoin (bAegFun H1_scaffold_1 / H2_scaffold_1 at 134.6 Mb). `olen`: bp per scaffold."""
    if p is None or q is None:
        return None
    if p[0] != q[0]:
        if olen is None:
            return None

        def room(e):                      # sequence of the other haplotype beyond where e lands
            return olen[e[0]] - e[1] if e[2] > 0 else e[1]
        if max(p[4], q[4]) > far:
            return -1.0                   # its own contig continues there
        return -HOM_ACROSS if max(room(p), room(q)) > far else None
    w = 1.0 if p[3] == q[3] and p[3] >= 0 else HOM_ACROSS
    d = (q[1] - p[1]) * p[2]
    if q[2] != -p[2] or d < -slack:
        return -w
    if d <= near:
        return w
    if d <= far:
        return w * (1.0 - (d - near) / (far - near))
    return -w


def telomere_ends(pt, nodes, n=3):
    """Per end: 'out' when a telomere at that end faces outwards (telomere strand tracks: a
    chromosome's left end reads the motif reversed, `telomere_5p`; its right end forward,
    `telomere_3p`): that end is a chromosome end. Without strand tracks, 'any' for telomere
    signal at the end (which way it faces unknown). A pixel spans tens of kb, so each end's window
    includes the pixel across the gap (as `sites.strand_at`)."""
    t5, t3 = pt.graphs.get("telomere_5p"), pt.graphs.get("telomere_3p")
    out = {}
    if t5 is not None and t3 is not None:
        from .tracks import telomere_threshold
        thr = telomere_threshold(t5 + t3)
        for v in nodes:
            if float(t5[max(0, v.a - 1):v.a + n].max()) >= thr:
                out[(v.i, "h")] = "out"
            if float(t3[max(0, v.b - n + 1):v.b + 2].max()) >= thr:
                out[(v.i, "t")] = "out"
        return out
    tel = pt.graphs.get("telomere")
    if tel is None:
        return out
    thr = max(50.0, 5 * float(np.percentile(tel, 99)))
    for v in nodes:
        if float(tel[v.a:v.a + n].max()) >= thr:
            out[(v.i, "h")] = "any"
        if float(tel[max(0, v.b - n + 1):v.b + 1].max()) >= thr:
            out[(v.i, "t")] = "any"
    return out


def solve(pt, blocks=None, level=0.5, hap_regex=r"^(H\d+)_", graph_junctions=None, end_joins=None,
          unloc_ratio=2.0, chrom_min_px=20, repeat_cov=1.5):
    """The adjacency solver (module docstring). Returns (rows, savestate, info)."""
    bpp = pt.bp_per_pix
    f = 8
    while pt.tex_res % f:
        f //= 2
    sharp = set(sharp_breaks(pt, level))
    nodes, yahs = contigs(pt, hap_regex, f=f, breaks=sharp)
    R = end_contact(pt, nodes, f=f)
    proj, share = _projections(pt, nodes, blocks)
    E = lambda v, e: 2 * v + (e == "t")
    mb = lambda v: nodes[v].n * bpp / 1e6
    name = lambda v: f"{pt.names[nodes[v].s]}:{(nodes[v].a - int(pt.scaffold_start[nodes[v].s])) * bpp / 1e6:.2f}" \
                     f"-{(nodes[v].b + 1 - int(pt.scaffold_start[nodes[v].s])) * bpp / 1e6:.2f}Mb"
    prior, join_min = _at("prior", level), _at("join_min", level)
    # graph evidence: YaHS joins it backs (fixed), joins between scaffold ends only it sees
    fixed = set()
    q_of = {q: (u, v) for u, v, q in yahs}
    idx = {n: i for i, n in enumerate(pt.names)}
    for j, _ in graph_junctions or []:
        s = idx.get(j.scaffold)
        if s is None:
            continue
        q0 = int(pt.scaffold_start[s] + j.gap_start // bpp)
        q1 = int(pt.scaffold_start[s] + max(j.gap_end - 1, j.gap_start) // bpp)
        for q in range(q0 - 2, q1 + 3):
            if q in q_of:
                fixed.add(q_of[q])
    first, last = {}, {}
    for v in nodes:
        first.setdefault(v.s, v.i)
        last[v.s] = v.i
    graph_edges = set()
    for ej in end_joins or []:
        sa, sb = idx.get(ej["sa"]), idx.get(ej["sb"])
        if sa is None or sb is None or sa not in first or sb not in first:
            continue
        x = (first[sa], "h") if ej["ea"] == "h" else (last[sa], "t")
        y = (first[sb], "h") if ej["eb"] == "h" else (last[sb], "t")
        graph_edges.add(frozenset((x, y)))
    # candidate edges: YaHS joins, graph joins, the best Hi-C partners of each end, homology
    cand = {}
    for u, v, q in yahs:
        cand[frozenset(((u, "t"), (v, "h")))] = q
    for g_ in graph_edges:
        cand.setdefault(g_, None)
    ends = [(v.i, e) for v in nodes for e in "ht"]
    same_hap = lambda x, y: nodes[x[0]].hap == nodes[y[0]].hap and x[0] != y[0]
    for x in ends:
        if not nodes[x[0]].hic:
            continue
        row = R[E(*x)]
        order = [k for k in np.argsort(-np.nan_to_num(row, nan=-1))[:12] if np.isfinite(row[k])]
        for k in [k for k in order if same_hap(x, (k // 2, "ht"[k % 2]))][:4]:
            cand.setdefault(frozenset((x, (k // 2, "ht"[k % 2]))), None)
    by_o = collections.defaultdict(list)
    for x, p in proj.items():
        by_o[p[0]].append((p[1], x))
    for lst in by_o.values():
        lst.sort()
        for k, (pos, x) in enumerate(lst):
            for pos2, y in lst[k + 1:k + 6]:
                if pos2 - pos > 5_000_000:
                    break
                if same_hap(x, y) and _homology(proj[x], proj[y]) == 1.0:
                    cand.setdefault(frozenset((x, y)), None)

    end_hap = np.array([nodes[k // 2].hap for k in range(2 * len(nodes))])

    def best_other(x, y):
        """Best contact ratio of end x with any same-haplotype end other than y."""
        row = np.where(end_hap == nodes[x[0]].hap, R[E(*x)], np.nan)
        row[E(*y)] = np.nan
        if np.isfinite(row).any():
            return float(np.nanmax(row))
        return float("nan")
    tel = telomere_ends(pt, nodes)
    spec_min = math.log2(_at("specific_min", level))
    covr = _signal(pt).coverage_ratio()
    rep = {v.i: (float(np.median(covr[v.a:v.b + 1])) if covr is not None else 1.0) for v in nodes}
    olen = pt.scaffold_npix * bpp
    edges = []
    for key, q in cand.items():
        x, y = sorted(key)
        parts = {}
        is_yahs = q is not None
        if is_yahs:
            parts["yahs"] = prior
        r = R[E(*x), E(*y)]
        if np.isfinite(r):
            parts["hic"] = W_HIC * float(np.clip(math.log2(max(r, 1e-6) / R_REF), -3, 3))
            nxt = max(best_other(x, y), best_other(y, x))
            if np.isfinite(nxt) and nxt > 0 and r > 0:
                parts["specific"] = W_SPEC * float(np.clip(math.log2(r / nxt), -2, 1))
        h = _homology(proj.get(x), proj.get(y), olen)
        if h is not None:
            parts["homology"] = W_HOM * h
        if key in graph_edges:
            parts["graph"] = W_GRAPH
        t = -sum(W_TEL_OUT if tel.get(z) == "out" else W_TEL if tel.get(z) == "any" else 0.0 for z in (x, y))
        if t:
            parts["telomere"] = t
        total = sum(parts.values())
        if is_yahs and (x[0], y[0]) in fixed or (is_yahs and (y[0], x[0]) in fixed):
            total, parts["graph_link"] = math.inf, math.inf
        # a new join on Hi-C alone (no homology, no graph) must be specific, and not touch a
        # collapsed repeat (MQ0 multimapping reads contact everything)
        veto = ""
        if not is_yahs and "graph" not in parts and parts.get("homology", 0) <= 0:
            if parts.get("specific", -W_SPEC) / W_SPEC < spec_min:
                veto = f"Hi-C alone, not specific (needs {2 ** spec_min:.2f}x the next end)"
            elif max(rep[x[0]], rep[y[0]]) >= repeat_cov:
                veto = f"Hi-C alone, on a collapsed repeat (coverage {max(rep[x[0]], rep[y[0]]):.1f}x)"
        edges.append(dict(x=x, y=y, q=q, parts=parts, total=total, r=r, veto=veto))
    edges.sort(key=lambda e: -e["total"])
    # paths: best edge first, each end once, no cycles
    used, parent = {}, list(range(len(nodes)))

    def root(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a
    taken = []
    for e in edges:
        bar = 0.0 if e["q"] is not None else join_min
        if e["veto"] or e["total"] <= bar or e["x"] in used or e["y"] in used or root(e["x"][0]) == root(e["y"][0]):
            continue
        used[e["x"]], used[e["y"]] = e["y"], e["x"]
        parent[root(e["x"][0])] = root(e["y"][0])
        taken.append(e)
    # walk the paths: [(node, reversed)], starting from an end with no edge
    paths, seen = [], set()
    for v in nodes:
        if v.i in seen:
            continue
        start = None
        for e0 in "ht":
            if (v.i, e0) not in used:
                start = (v.i, e0)
                break
        if start is None:
            continue
        path, cur = [], start
        while True:
            u, e_in = cur
            seen.add(u)
            path.append((u, e_in == "t"))           # entered through its tail: laid out reversed
            out = (u, "t" if e_in == "h" else "h")
            if out not in used:
                break
            cur = used[out]
        # keep the assembly's direction where it does not matter: first node forward if possible
        if path[0][1] and len(path) == 1:
            path = [(path[0][0], False)]
        paths.append(path)
    plen = lambda p: sum(nodes[u].n for u, _ in p)
    big = [p for p in paths if plen(p) >= chrom_min_px and any(nodes[u].hic for u, _ in p)]
    small = [p for p in paths if p not in big]
    # pair the big paths of the two haplotypes by aligned bp
    pid = {}
    for k, p in enumerate(big):
        for u, _ in p:
            pid[u] = k
    S = collections.defaultdict(float)
    for u, d in share.items():
        for w, bp in d.items():
            if u in pid and w in pid and pid[u] != pid[w] and nodes[u].hap != nodes[w].hap:
                S[(pid[u], pid[w])] += bp
    best = {}
    for (a, b), bp in S.items():
        if bp > best.get(a, (None, 0))[1]:
            best[a] = (b, bp)
    pairs = []
    for a, (b, bp) in best.items():
        if a < b and best.get(b, (None,))[0] == a:
            second = max([v for (x, y), v in S.items() if x == a and y != b] + [0])
            if bp >= 2 * second:
                pairs.append((a, b, bp))
    # turn the second copy to run like the first: compare node positions along both paths
    pos_in = {}
    for k, p in enumerate(big):
        off = 0
        for u, r in p:
            pos_in[u] = (off, r)
            off += nodes[u].n

    def same_dir(a, b):
        xs, ys, ws = [], [], []
        for u, _ in big[a]:
            for w, bp in share.get(u, {}).items():
                if pid.get(w) == b:
                    xs.append(pos_in[u][0])
                    ys.append(pos_in[w][0])
                    ws.append(bp)
        if len(xs) < 2:
            return True
        c = np.cov(xs, ys, aweights=ws)[0, 1] if len(set(xs)) > 1 and len(set(ys)) > 1 else 1.0
        return c >= 0
    for a, b, _ in pairs:
        if not same_dir(a, b):
            big[b] = [(u, not r) for u, r in reversed(big[b])]
    # small leftovers: to the big path of their haplotype their contact (else alignment) points to
    f = 8
    while pt.tex_res % f:
        f //= 2
    M = _signal(pt).binned(f)
    nb = M.shape[0]
    pbins = {}
    for k, p in enumerate(big):
        bins = np.zeros(nb, dtype=bool)
        for u, _ in p:
            bins[nodes[u].a // f:nodes[u].b // f + 1] = True
        pbins[k] = bins
    unloc = collections.defaultdict(list)
    small_rows = []
    for p in small:
        hp = nodes[p[0][0]].hap
        bins = np.zeros(nb, dtype=bool)
        for u, _ in p:
            bins[nodes[u].a // f:nodes[u].b // f + 1] = True
        rowsum = M[bins].sum(axis=0)
        dens = {k: float(rowsum[pb].sum() / pb.sum()) for k, pb in pbins.items()
                if nodes[big[k][0][0]].hap == hp and pb.sum()}
        ranked = sorted(dens.items(), key=lambda kv: -kv[1])
        al = collections.Counter()
        for u, _ in p:
            for w, bp in share.get(u, {}).items():
                if w in pid:
                    al[pid[w]] += bp
        # alignment to a chromosome of the other haplotype points to its paired copy
        mate = {a: b for a, b, _ in pairs} | {b: a for a, b, _ in pairs}
        al_own = collections.Counter()
        for k, bp in al.items():
            if k in mate and nodes[big[mate[k]][0][0]].hap == hp:
                al_own[mate[k]] += bp
        target, why = None, ""
        if ranked and (len(ranked) == 1 or ranked[0][1] >= unloc_ratio * max(ranked[1][1], 1e-9)):
            target = ranked[0][0]
            why = (f"contact {ranked[0][1]:.2f} per bin with that chromosome, "
                   f"{ranked[0][1] / max(ranked[1][1], 1e-9) if len(ranked) > 1 else float('inf'):.1f}x the next")
            if al_own and al_own.most_common(1)[0][0] != target:
                target, why = None, why + "; but its sequence aligns to another chromosome: left unplaced"
        elif al_own:
            k, bp = al_own.most_common(1)[0]
            if bp >= 0.5 * plen(p) * bpp:
                target, why = k, f"aligns ({bp / 1e6:.2f} Mb) to the other haplotype's copy of that chromosome"
        if target is not None:
            unloc[target].append(p)
        small_rows.append((p, target, why))
    # layout: pairs (longer copy first, HAP1), then unpaired big paths, then the rest
    in_pair = {a for a, b, _ in pairs} | {b for a, b, _ in pairs}
    copies = []
    for a, b, bp in sorted(pairs, key=lambda t: -max(plen(big[t[0]]), plen(big[t[1]]))):
        h1, h2 = (a, b) if plen(big[a]) >= plen(big[b]) else (b, a)
        copies += [(h1, "HAP1", (a, b, bp)), (h2, "HAP2", (a, b, bp))]
    for k in sorted(set(range(len(big))) - in_pair, key=lambda k: -plen(big[k])):
        copies.append((k, None, None))
    order, paint, tag = [], {}, {}
    for n_, (k, t, _) in enumerate(copies, 1):
        for u, r in big[k]:
            order.append((u, r))
            paint[u], tag[u] = n_, t
        for p in unloc.get(k, []):
            for u, r in p:
                order.append((u, r))
                paint[u], tag[u] = n_, t
    placed = {u for u, _ in order}
    rest = [(v.i, False) for v in nodes if v.i not in placed]
    perm = np.concatenate([np.arange(nodes[u].a, nodes[u].b + 1)[::-1] if r else np.arange(nodes[u].a, nodes[u].b + 1)
                           for u, r in order + rest])
    covered = np.zeros(pt.npix, dtype=bool)
    covered[perm] = True
    if len(perm) != pt.npix or not covered.all():
        raise AssertionError("solver layout does not cover every pixel once")
    where = np.empty(pt.npix, dtype=np.int64)
    where[perm] = np.arange(pt.npix)
    kept = {frozenset((e["x"], e["y"])) for e in taken}
    cuts = [q for u, v, q in yahs if frozenset(((u, "t"), (v, "h"))) not in kept]
    st = SaveState.blank()
    st.edits = [Edit(q, q, BREAK_DELTA) for q in cuts if where[q + 1] - where[q] == 1] + plan_edits(perm)
    ms = MapState(pt).replay(st.edits)
    if not np.array_equal(ms.pix, perm):
        raise AssertionError("planned edits do not reproduce the solver layout")
    node_of = np.empty(pt.npix, dtype=np.int64)
    for v in nodes:
        node_of[v.a:v.b + 1] = v.i
    unloc_nodes = {u for ps in unloc.values() for p in ps for u, _ in p}
    for i, c in enumerate(ms.contigs):
        u = int(node_of[ms.pix[c.start]])
        if u in paint:
            st.scaffs.append((i, paint[u]))
        mask = 0
        if tag.get(u):
            mask |= 1 << st.tag_id(tag[u])
        if u in unloc_nodes:
            mask |= 1 << st.tag_id("Unloc")
        if mask:
            st.meta_flags.append((i, mask))
    # decision rows
    fmt = lambda parts: ", ".join(f"{k} {v:+.2f}" if math.isfinite(v) else f"{k} fixed" for k, v in parts.items())
    rows = []
    eby = {frozenset((e["x"], e["y"])): e for e in edges}
    for u, v, q in yahs:
        e = eby[frozenset(((u, "t"), (v, "h")))]
        pos = (q + 1 - int(pt.scaffold_start[nodes[u].s])) * bpp / 1e6
        if frozenset((e["x"], e["y"])) in kept:
            rows.append(dict(action="kept", scaffold=pt.names[nodes[u].s], length_mb=round(pos, 2), partner="",
                             strand="", tag="", score=round(e["total"], 2) if math.isfinite(e["total"]) else "fixed",
                             next_best="", evidence=f"YaHS join kept: score {e['total']:.2f} ({fmt(e['parts'])})"))
            continue
        rival = [(x, used.get(x)) for x in (e["x"], e["y"]) if used.get(x) is not None]
        why = (("no gap here (inside a contig): contact drops to background; " if q in sharp else "gap inside the scaffold; ")
               + f"join not kept: score {e['total']:.2f} ({fmt(e['parts'])})"
               + ("".join(f"; {name(x[0])} {x[1]} joined to {name(y[0])} {y[1]} instead" for x, y in rival)
                  if rival else "; below 0" if e["total"] <= 0 else
                  "; it would close a loop: the two pieces are already in one path through other joins"))
        rows.append(dict(action="cut", scaffold=pt.names[nodes[u].s], length_mb=round(pos, 2), partner="", strand="",
                         tag="", score=round(e["total"], 2), next_best="", evidence=why))
    for e in taken:
        if e["q"] is not None:
            continue
        (u, eu), (v, ev) = e["x"], e["y"]
        rows.append(dict(action="join", scaffold=f"{name(u)} {eu}", length_mb=round(mb(u), 2), partner=f"{name(v)} {ev}",
                         strand="", tag="", score=round(e["total"], 2), next_best="",
                         evidence=f"new adjacency, score {e['total']:.2f} ({fmt(e['parts'])}); bar {join_min:.2f} at "
                                  f"level {level:g}"))
    for a, b, bp in pairs:
        rows.append(dict(action="pair", scaffold=" + ".join(name(u) for u, _ in big[a][:3]) + (" ..." if len(big[a]) > 3 else ""),
                         length_mb=round(plen(big[a]) * bpp / 1e6, 2),
                         partner=" + ".join(name(u) for u, _ in big[b][:3]) + (" ..." if len(big[b]) > 3 else ""),
                         strand="", tag="", score=round(bp / 1e6, 2), next_best="",
                         evidence=f"{bp / 1e6:.1f} Mb aligned between the two paths"))
    for p, target, why in small_rows:
        act = "unloc" if target is not None else "unplaced"
        rows.append(dict(action=act, scaffold=" + ".join(name(u) for u, _ in p), length_mb=round(plen(p) * bpp / 1e6, 2),
                         partner=(name(big[target][0][0]) if target is not None else ""), strand="",
                         tag="Unloc" if target is not None else "", score="", next_best="",
                         evidence=why or "no chromosome stands out by contact or alignment"))
    for k in sorted(set(range(len(big))) - in_pair):
        rows.append(dict(action="unpaired", scaffold=name(big[k][0][0]), length_mb=round(plen(big[k]) * bpp / 1e6, 2),
                         partner="", strand="", tag="", score="", next_best="",
                         evidence="no path of the other haplotype aligns to it at 2x the next"))
    info = dict(nodes=len(nodes), yahs=len(yahs), cuts=len(cuts), joins=sum(e["q"] is None for e in taken),
                pairs=len(pairs), paths=len(big), unlocs=len(unloc_nodes))
    return rows, st, info
