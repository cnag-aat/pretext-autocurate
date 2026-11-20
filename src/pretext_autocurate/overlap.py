"""Contig-end overlaps: the end of one scaffold is another copy (the other haplotype) of the
end of another scaffold, and the copy runs to that scaffold's end.

hifiasm sometimes leaves both haplotypes of a stretch at two contig ends (icLarUrsu.pri:
scaffold_16 0.26-6.32 Mb is a 0.1-3% divergent copy of scaffold_22 16.82-23.07 Mb, to its end;
all of scaffold_32 is a copy of scaffold_14 0-6.1 Mb). The two ends then touch strongly and a
join of them looks right, but the copy has to go: the curator trimmed it and joined the other
scaffold's end to where the copy stops.

On the map an overlap of D px pairs depth k from one end with depth D-1-k from the other: a
line of contact as bright as the diagonal. Mere adjacency of the two ends puts those pairs D px
apart, so it gives that line only the decay at D. The overlap length is the D with the largest
excess. Inside a contig there is no gap to cut at, and a copy line alone does not tell a
haplotypic copy from a paralog, so a call is applied only when the sequence confirms it (the two
ends align in order, at haplotype-level divergence, `align_overlap`); the alignment also gives
the cut to the base.
"""
import os
import subprocess
import tempfile
from dataclasses import dataclass

import numpy as np

from .signal import Signal


@dataclass
class Params:
    w: int = 10               # px: end windows for picking partner ends
    top: int = 3              # partner ends examined per end
    min_unit: int = 20        # px: shorter pieces are not examined
    min_px: int = 20          # px: shortest overlap
    max_px: int = 400         # px: longest overlap examined
    line_min: float = 0.8     # copy line / diagonal
    excess_min: float = 0.35  # copy line minus what adjacency at that length gives
    cov_max: float = 1.3      # a copy with more coverage than this is a repeat (subtelomeric arrays)


@dataclass
class EndOverlap:
    a: int                    # scaffold (piece) index and end ('h' / 't') of each copy
    a_end: str
    b: int
    b_end: str
    overlap_px: int
    line: float
    expected: float
    cov_a: float
    cov_b: float

    @property
    def excess(self):
        return self.line - self.expected


def _depth_px(pt, s, e, k):
    a, L = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
    return a + k if e == "h" else a + L - 1 - k


def detect_end_overlaps(pt, params=None, signal=None, group=None, skip=None):
    """End overlaps in `pt` (Pretext or SegmentedMap), strongest first, each end in one call at
    most. `group(name)`: only ends of the same group are compared (haplotype prefix in allhaps);
    `skip(u, v)`: pieces never compared with each other."""
    p = params or Params()
    sg = signal or Signal(pt)
    cov = sg.coverage_ratio()
    g = np.array(sg.decay_long(max(p.max_px + 2, 2 * p.w + 2)), dtype=np.float64)
    g = g[:max(1, int(np.max(np.nonzero(g > 0)[0])) + 1)] if (g > 0).any() else g
    g0 = float(g[0]) or 1.0
    units = [s for s in range(len(pt.names)) if pt.scaffold_npix[s] >= p.min_unit]
    ends = [(s, e) for s in units for e in "ht"]
    k = np.arange(p.w)
    strip = {x: sg.strip(int(min(_depth_px(pt, *x, 0), _depth_px(pt, *x, p.w - 1))),
                         int(max(_depth_px(pt, *x, 0), _depth_px(pt, *x, p.w - 1)))).astype(np.float64) for x in ends}
    grp = (lambda s: group(pt.names[s])) if group else (lambda s: 0)
    seen, found = set(), []
    for x in ends:
        sc = sorted(((float(strip[x][:, _depth_px(pt, *y, k)].mean()), y) for y in ends
                     if y[0] != x[0] and grp(y[0]) == grp(x[0]) and not (skip and skip(x[0], y[0]))),
                    reverse=True)
        for _, y in sc[:p.top]:
            if (y, x) in seen or (x, y) in seen:
                continue
            seen.add((x, y))
            D = min(p.max_px, int(pt.scaffold_npix[x[0]]), int(pt.scaffold_npix[y[0]]))
            if D < p.min_px:
                continue
            kx, ky = _depth_px(pt, *x, np.arange(D)), _depth_px(pt, *y, np.arange(D))
            B = sg.block(int(kx.min()), int(kx.max()), int(ky.min()), int(ky.max())).astype(np.float64)
            F = B[kx - kx.min()][:, ky - ky.min()][:, ::-1]      # overlap d: diagonal D - d of F
            best = None
            for d in range(p.min_px, D + 1):
                line = float(np.diagonal(F, D - d).mean()) / g0
                exp = float(g[min(d, len(g) - 1)]) / g0
                if best is None or line - exp > best[0]:
                    best = (line - exp, d, line, exp)
            _, d, line, exp = best
            ca = float(np.median(cov[kx[:d]])) if cov is not None else float("nan")
            cb = float(np.median(cov[ky[:d]])) if cov is not None else float("nan")
            if line >= p.line_min and line - exp >= p.excess_min and not (ca > p.cov_max or cb > p.cov_max):
                found.append(EndOverlap(x[0], x[1], y[0], y[1], d, line, exp, ca, cb))
    found.sort(key=lambda o: -o.excess)
    used, out = set(), []
    for o in found:
        if (o.a, o.a_end) in used or (o.b, o.b_end) in used:
            continue
        used |= {(o.a, o.a_end), (o.b, o.b_end)}
        out.append(o)
    return out


@dataclass
class Side:
    """One copy's end in an input scaffold: `edge` is the 0-based coordinate of the end (a head
    starts there, a tail ends just before it); the copy runs inwards from it."""
    name: str
    edge: int
    head: bool

    def depth(self, x):
        return x - self.edge if self.head else self.edge - x

    def region(self, window, length):
        return (self.edge, min(length, self.edge + window)) if self.head else (max(0, self.edge - window), self.edge)


@dataclass
class OverlapAlignment:
    aligned_bp: int           # bases of the trimmed copy in alignment blocks (union)
    divergence: float         # alignment-length-weighted gap-compressed divergence (minimap2 de)
    copy_bp: int              # depth (bp from the trimmed end) where the copy stops: the cut
    keep_reach_bp: int        # how close (bp) the copy comes to the kept end
    blocks: int
    # Gapless switch: a stretch of >= `match_bp` identical bases nearest the kept end. The output
    # can run along the kept copy to keep_switch (0-based, included) and go on in the trimmed
    # scaffold after trim_switch (the base matching it), with no gap. None when there is none.
    keep_switch: int = None
    trim_switch: int = None
    match_bp: int = 0


def _copy_chain(items, max_gap, near_bp=100_000, keep_share=0.8):
    """The copy as a collinear chain of alignment blocks in depth coordinates: items (a0, a1, b0,
    b1, length, de), a = depth from the trimmed end, b = depth from the kept end. Along the copy
    a grows while b shrinks (depth k pairs with D-1-k); jumps under `max_gap` bp in each, so
    indels between the haplotypes do not break it (icLarUrsu scaffold_14 carries 1.5 Mb more
    than scaffold_32 over their 6 Mb copy) while repeats aligning elsewhere fall out.

    The copy ends where it first meets the kept end: of the blocks within `near_bp` of the
    smallest b, whose best chain (from the trimmed end) weighs at least `keep_share` of the
    heaviest, the one at the smallest depth a. (scaffold_14 holds scaffold_32's last Mb twice,
    at 6.73-7.62 and 8.21-9.12 Mb.) Returns the chain, trimmed-end first."""
    bl = sorted(items, key=lambda b: b[0])
    best, prev = [b[4] for b in bl], [-1] * len(bl)
    slack = 5000
    for j, bj in enumerate(bl):
        for i in range(j):
            bi = bl[i]
            da, db = bj[0] - bi[1], bi[2] - bj[3]          # a forwards, b backwards
            if -slack <= da <= max_gap and -slack <= db <= max_gap and best[i] + bj[4] > best[j]:
                best[j], prev[j] = best[i] + bj[4], i
    near = min(b[2] for b in bl) + near_bp
    cand = [j for j, b in enumerate(bl) if b[2] <= near]
    top = max(best[j] for j in cand)
    j = min((j for j in cand if best[j] >= keep_share * top), key=lambda j: bl[j][1])
    out = []
    while j >= 0:
        out.append(bl[j])
        j = prev[j]
    return out[::-1]


def _match_runs(q0, q1, t0, strand, cs, min_len):
    """Runs of >= min_len identical bases in a minimap2 --cs (short) string: [(q, t, n)], base
    q + k pairs with t + k ('+'), or q - k with t + k ('-', q the run's highest query base)."""
    out, qi, ti = [], (q0 if strand == "+" else q1 - 1), t0
    step = 1 if strand == "+" else -1
    i, n = 0, len(cs)
    while i < n:
        op = cs[i]
        j = i + 1
        while j < n and cs[j] not in ":*+-~":
            j += 1
        arg = cs[i + 1:j]
        if op == ":":
            m = int(arg)
            if m >= min_len:
                out.append((qi, ti, m))
            qi, ti = qi + step * m, ti + m
        elif op == "*":
            qi, ti = qi + step, ti + 1
        elif op == "+":                     # bases in the query only
            qi += step * len(arg)
        elif op == "-":                     # bases in the target only
            ti += len(arg)
        i = j
    return out


def align_overlap(fa, trim, keep, window_bp, minimap2="minimap2", threads=4, min_block=10000,
                  max_gap=2_000_000, match_bp=1000):
    """Align `window_bp` inwards from the trimmed copy's end (a Side) against the same from the
    kept copy's end (minimap2 asm20, base level, --cs). Copies pair depth k in one with depth
    D-1-k in the other, so the expected strand is '+' for a head and a tail, '-' for two heads
    or two tails. `fa` is a graph.Faidx. Returns an OverlapAlignment, or None when nothing
    aligns in that strand."""
    ta, tb = trim.region(window_bp, fa.length(trim.name))
    ka, kb = keep.region(window_bp, fa.length(keep.name))
    want = "+" if trim.head != keep.head else "-"
    with tempfile.TemporaryDirectory() as d:
        q, t = os.path.join(d, "q.fa"), os.path.join(d, "t.fa")
        with open(q, "w") as f:
            f.write(f">q\n{fa.fetch(trim.name, ta, tb)}\n")
        with open(t, "w") as f:
            f.write(f">t\n{fa.fetch(keep.name, ka, kb)}\n")
        r = subprocess.run([minimap2, "-x", "asm20", "-c", "--cs", "--secondary=no", "-t", str(threads), t, q],
                           capture_output=True, text=True)
    if r.returncode != 0:
        raise RuntimeError(f"minimap2 failed: {r.stderr.strip()[-300:]}")
    raw = []
    for line in r.stdout.splitlines():
        c = line.split("\t")
        if len(c) < 12 or c[4] != want or int(c[10]) < min_block:
            continue
        de = next((float(x[5:]) for x in c[12:] if x.startswith("de:f:")), float("nan"))
        cs = next((x[5:] for x in c[12:] if x.startswith("cs:Z:")), "")
        q0, q1, t0, t1 = int(c[2]) + ta, int(c[3]) + ta, int(c[7]) + ka, int(c[8]) + ka
        a = sorted((trim.depth(q0), trim.depth(q1)))
        b = sorted((keep.depth(t0), keep.depth(t1)))
        raw.append((a[0], a[1], b[0], b[1], int(c[10]), de, (q0, q1, t0, cs)))
    if not raw:
        return None
    blocks = _copy_chain(raw, max_gap)
    copy, reach = blocks[-1][1], blocks[-1][2]      # the cut (depth from the trimmed end); distance to the kept end
    # gapless switch: identical run nearest the kept end, from the end of the chain inwards
    sw = None
    for x in reversed(blocks):
        q0, q1, t0, cs = x[6]
        runs = _match_runs(q0, q1, t0, want, cs, match_bp)
        if runs:
            qs, ts, m = min(runs, key=lambda r: keep.depth(r[1] + r[2] // 2))
            k = m // 2
            sw = (ts + k, qs + k if want == "+" else qs - k, m)
            break
    # union of the trimmed copy's aligned bases
    iv = sorted((x[0], x[1]) for x in blocks)
    tot, cur = 0, None
    for s, e in iv:
        if cur is None or s > cur[1]:
            if cur:
                tot += cur[1] - cur[0]
            cur = [s, e]
        else:
            cur[1] = max(cur[1], e)
    tot += cur[1] - cur[0]
    w = sum(b[4] for b in blocks)
    div = sum(b[4] * b[5] for b in blocks) / w if w else float("nan")
    return OverlapAlignment(tot, div, copy, reach, len(blocks), *(sw if sw else (None, None, 0)))


def resolve_overlaps(raw, vt, fasta=None, minimap2=None, threads=4, group=None, params=None):
    """Detect end overlaps among the pieces of `vt` (a SegmentedMap of `raw`; `group(name)`
    limits them to ends of one group, the haplotype in allhaps) and decide each with the
    sequence. The copy on the longer piece is the one cut off (the curator, icLarUrsu scaffold_16
    and 14). Returns dicts: ov, trim / keep (piece, end), tside / kside (Side), al
    (OverlapAlignment or None), apply, q (cut after this raw pixel), cut_bp, why (when not
    applied), exact (breaks rows: the cut, and the gapless switch when there is one)."""
    fa, mm = None, None
    if fasta:
        from .graph import Faidx
        from .homology import find_tool
        fa, mm = Faidx(fasta), find_tool("minimap2", minimap2)
    lens_bp = raw.scaffold_len_bp()
    N, G = float(raw.npix), float(raw.total_len)

    def side(u, e):
        """The piece end as a Side of its input scaffold (bp from the pixels at a cut)."""
        p = int(vt.scaffold_start[u]) if e == "h" else int(vt.scaffold_start[u] + vt.scaffold_npix[u] - 1)
        s = int(raw.orig_scaffold[p])
        if e == "h":
            edge = 0 if raw.rel_coord[p] == 0 else int(round(raw.rel_coord[p] * raw.bp_per_pix))
        else:
            last = p + 1 >= raw.npix or raw.orig_scaffold[p + 1] != s
            edge = int(lens_bp[s]) if last else int(round((raw.rel_coord[p] + 1) * raw.bp_per_pix))
        return Side(raw.names[s], edge, e == "h"), s
    out = []
    # Two pieces of one input scaffold are not compared: their ends meet at a cut the misjoin
    # rules made, and rejoining them past a copy undid that cut (xgTylRafi scaffold_2 and 3,
    # 3-3.6% divergence: junctions 45/68 -> 43/68). A copy around one break is a gap's business.
    same = lambda u, v: raw.orig_scaffold[int(vt.scaffold_start[u])] == raw.orig_scaffold[int(vt.scaffold_start[v])]
    for o in detect_end_overlaps(vt, params or Params(), group=group, skip=same):
        A, B = (o.a, o.a_end), (o.b, o.b_end)
        tr, kp = (A, B) if vt.scaffold_npix[o.a] >= vt.scaffold_npix[o.b] else (B, A)
        (ts, tsi), (ks, _) = side(*tr), side(*kp)
        D_bp = o.overlap_px * raw.bp_per_pix
        rec = dict(ov=o, trim=tr, keep=kp, tside=ts, kside=ks, al=None, apply=False, q=None, cut_bp=None, exact=[])
        why = ""
        if fa is None or mm is None:
            why = "no --fasta (or no minimap2): proposed, the sequence decides where the copy ends"
        else:
            al = align_overlap(fa, ts, ks, int(1.5 * D_bp + 20 * raw.bp_per_pix), minimap2=mm, threads=threads)
            rec["al"] = al
            pa, pb = int(vt.scaffold_start[tr[0]]), int(vt.scaffold_start[tr[0]] + vt.scaffold_npix[tr[0]] - 1)
            if al is None:
                why = "the two ends do not align in the expected orientation: proposed"
            elif al.aligned_bp < 0.5 * al.copy_bp or al.divergence > 0.05:
                why = (f"the sequence aligns over {al.aligned_bp / max(al.copy_bp, 1):.0%} of the copy at "
                       f"{al.divergence:.1%} divergence (needs 50%, <= 5%): proposed")
            elif al.keep_reach_bp > max(200_000, 2 * raw.bp_per_pix):
                why = f"the copy stops {al.keep_reach_bp / 1e3:.0f} kb short of the other end: proposed"
            elif abs(al.copy_bp - D_bp) > 0.3 * D_bp:
                why = (f"copy {al.copy_bp / 1e6:.2f} Mb in the sequence but {D_bp / 1e6:.2f} Mb on the map "
                       "(repeats or duplications in the copy): proposed")
            else:
                cut_bp = ts.edge + al.copy_bp if ts.head else ts.edge - al.copy_bp
                q = int(raw.scaffold_start[tsi]) + int(round(cut_bp / raw.bp_per_pix)) - 1
                if not pa <= q < pb:
                    why = "the copy covers the whole piece: proposed"
                else:
                    rec.update(apply=True, q=q, cut_bp=cut_bp)
                    rec["exact"].append(dict(kind="cut", scaffold=ts.name, px_bp=int((raw.rel_coord[q] + 1) * G / N),
                                             bp=cut_bp, partner="", partner_bp=0, method="alignment",
                                             evidence=f"copy of {ks.name} ends here ({al.divergence:.2%} divergence)"))
                    if al.keep_switch is not None:
                        rec["exact"].append(dict(
                            kind="switch", scaffold=ks.name, px_bp=int(ks.edge), bp=al.keep_switch + 1, partner=ts.name,
                            partner_bp=al.trim_switch + 1, method="alignment",
                            evidence=f"{al.match_bp:,} identical bases around the switch; {ks.name} "
                                     f"{'start' if ks.head else 'end'} to {ts.name} past its copy"))
        rec["why"] = why
        out.append(rec)
    return out


def overlap_evidence(rec, bpp):
    """The evidence text of an overlap decision row."""
    o, ts, ks, al = rec["ov"], rec["tside"], rec["kside"], rec["al"]
    wh = lambda sd: f"{sd.name} {'start' if sd.head else 'end'}"
    ev = (f"{wh(ts)} is another copy of {wh(ks)}: copy line {o.line:.2f} of the diagonal over "
          f"{o.overlap_px * bpp / 1e6:.2f} Mb (adjacency alone gives {o.expected:.2f}); coverage "
          f"{o.cov_a:.2f} / {o.cov_b:.2f}")
    if al is not None:
        ev += (f"; sequence: {al.aligned_bp / 1e6:.2f} Mb aligned in order at {al.divergence:.2%} divergence, the "
               f"copy ends {al.copy_bp / 1e6:.2f} Mb from the {'start' if ts.head else 'end'}, "
               f"{al.keep_reach_bp / 1e3:.0f} kb from {wh(ks)}")
        if al.keep_switch is not None:
            ev += (f"; gapless join possible: {ks.name} to {al.keep_switch + 1:,} then {ts.name} from "
                   f"{al.trim_switch + 2 if ts.head else al.trim_switch:,} ({al.match_bp:,} identical bases there)")
    return ev


def gap_overlaps(raw, fasta, minimap2=None, threads=4, window=1_000_000, min_copy=10_000):
    """Gaps (from the gap track, placed to the base in the sequence) whose two sides are copies:
    the right contig's start aligns in order with the left contig's end and runs into it to
    within max(20 kb, 10% of the copy), >= `min_copy`, half of it aligned at <= 5% divergence,
    with an identical stretch to switch at. Returns [dict(scaffold, gap, al, switch row or None)]
    in gap order. The gaps are aligned `threads` at a time (one minimap2 thread each: the runs
    are small, and on a cluster node more cores simply means more gaps at once)."""
    import re
    from concurrent.futures import ThreadPoolExecutor
    from .graph import Faidx
    from .homology import find_tool
    mm = find_tool("minimap2", minimap2)
    if mm is None:
        return []
    fa, gap, bpp = Faidx(fasta), raw.graphs["gap"], raw.bp_per_pix
    lens = raw.scaffold_len_bp()
    gaps, seen = [], set()
    for q in np.nonzero(gap > 0)[0]:                 # the gaps, placed to the base in the sequence
        s = int(raw.orig_scaffold[q])
        name = raw.names[s]
        a = max(0, int((raw.rel_coord[q] - 1) * bpp))
        b = min(int(lens[s]), int((raw.rel_coord[q] + 2) * bpp))
        for m in re.finditer(r"[Nn]{10,}", fa.fetch(name, a, b, upper=False)):
            g0, g1 = a + m.start(), a + m.end()
            if (name, g0) in seen or g0 == 0 or g1 >= lens[s]:
                continue
            seen.add((name, g0))
            gaps.append((name, g0, g1))

    def one(g):
        name, g0, g1 = g
        left, right = Side(name, g0, False), Side(name, g1, True)
        al = align_overlap(fa, right, left, window, minimap2=mm, threads=1, min_block=5000)
        if al is not None and al.copy_bp > 0.9 * window:            # longer than the window
            al = align_overlap(fa, right, left, 3 * window, minimap2=mm, threads=1, min_block=5000)
        ok = (al is not None and al.copy_bp >= min_copy and al.aligned_bp >= 0.5 * al.copy_bp
              and al.divergence <= 0.05 and al.keep_reach_bp <= max(20_000, 0.1 * al.copy_bp)
              and al.keep_switch is not None)
        sw = None
        if ok:
            sw = dict(kind="switch", scaffold=name, px_bp=g0, bp=al.keep_switch + 1, partner=name,
                      partner_bp=al.trim_switch + 1, method="alignment", src="gap",
                      evidence=f"gap at {g0 + 1:,}: the {al.copy_bp / 1e3:.0f} kb after it is another copy of the "
                               f"sequence before it ({al.divergence:.2%} divergence); {al.match_bp:,} identical "
                               f"bases around the switch")
        return dict(scaffold=name, gap=(g0, g1), al=al, ok=ok, switch=sw)
    with ThreadPoolExecutor(max_workers=max(1, int(threads))) as ex:     # minimap2 runs outside the GIL
        return list(ex.map(one, gaps))


@dataclass
class CopyCall:
    """What a candidate haplotig is, from its alignment to the partner region its Hi-C copy line
    points at (`classify_candidates`)."""
    piece: str
    partner: str
    kind: str                 # gap_fill | copy | none
    strand: str               # '+' / '-': the piece's direction along the partner
    aligned: float            # share of the piece in the chain of alignment blocks
    divergence: float
    gap: tuple = None         # gap_fill: the partner's gap (0-based [g0, g1))
    left: tuple = None        # gap_fill: the piece's copy at its start ([a, b), piece bp)
    right: tuple = None       # gap_fill: its copy at its end
    fill: tuple = None        # gap_fill: the piece's own sequence that goes into the gap
    # gap_fill: where each copy can switch to the flank without a gap, 0-based (piece base,
    # identical flank base), in an identical stretch nearest the fill; None when there is none
    left_switch: tuple = None
    right_switch: tuple = None


def _chain(blocks, strand, max_gap):
    """Heaviest collinear chain of (q0, q1, t0, t1, len, de) in one strand, sorted by q."""
    tt = (lambda b: (b[2], b[3])) if strand == "+" else (lambda b: (-b[3], -b[2]))
    bl = sorted(blocks, key=lambda b: b[0])
    if not bl:
        return [], 0
    best, prev = [b[4] for b in bl], [-1] * len(bl)
    for j, bj in enumerate(bl):
        tj0 = tt(bj)[0]
        for i in range(j):
            ti1 = tt(bl[i])[1]
            dq, dt = bj[0] - bl[i][1], tj0 - ti1
            if -5000 <= dq <= max_gap and -5000 <= dt <= max_gap and best[i] + bj[4] > best[j]:
                best[j], prev[j] = best[i] + bj[4], i
    j = int(np.argmax(best))
    w, out = best[j], []
    while j >= 0:
        out.append(bl[j])
        j = prev[j]
    return out[::-1], w


def classify_candidates(fasta, cands, minimap2=None, threads=4, flank=1_000_000, max_gap=2_000_000, min_copy=10_000):
    """Each candidate (piece, partner, partner_bp 0-based) aligned against the partner around
    partner_bp (the piece's length plus `flank` each way), minimap2 asm20, the heaviest collinear
    chain in either strand. Returns {piece: CopyCall}:

      gap_fill  the chain starts at one end of the piece and ends at the other, and jumps across
                a gap of the partner between two blocks: the piece's ends are copies of the two
                flanks (each >= min_copy, reaching the gap within max(20 kb, 10%) of it, <= 5%
                divergence) and the stretch between them fills the gap (xgTylRafi scaffold_29: its
                first 236 kb copy scaffold_1 just after the 18.23 Mb gap, its last 77 kb just
                before it, both on the - strand, and 392 kb between them fill it);
      copy      >= 70% of the piece in the chain at <= 5% divergence: a haplotypic copy;
      none      otherwise (the copy line was something else, or the copy is too diverged).

    The pieces are aligned `threads` at a time, one minimap2 thread each."""
    import re
    from concurrent.futures import ThreadPoolExecutor
    from .graph import Faidx
    from .homology import find_tool
    mm = find_tool("minimap2", minimap2)
    if mm is None:
        return {}
    fa = Faidx(fasta)
    tol = lambda n: max(20_000, 0.1 * n)

    def one(c):
        piece, partner, pbp = c
        L, Lp = fa.length(piece), fa.length(partner)
        r0, r1 = max(0, pbp - L - flank), min(Lp, pbp + L + flank)
        tseq = fa.fetch(partner, r0, r1)
        with tempfile.TemporaryDirectory() as d:
            q, t = os.path.join(d, "q.fa"), os.path.join(d, "t.fa")
            with open(q, "w") as f:
                f.write(f">q\n{fa.fetch(piece, 0, L)}\n")
            with open(t, "w") as f:
                f.write(f">t\n{tseq}\n")
            r = subprocess.run([mm, "-x", "asm20", "-c", "--cs", "--secondary=no", "-t", "1", t, q], capture_output=True,
                               text=True)
        bl = {"+": [], "-": []}
        for line in r.stdout.splitlines():
            x = line.split("\t")
            if len(x) < 12 or int(x[10]) < 5000:
                continue
            de = next((float(y[5:]) for y in x[12:] if y.startswith("de:f:")), float("nan"))
            cs = next((y[5:] for y in x[12:] if y.startswith("cs:Z:")), "")
            bl[x[4]].append((int(x[2]), int(x[3]), int(x[7]) + r0, int(x[8]) + r0, int(x[10]), de, cs))
        (ch, w), strand = max(((_chain(bl[s], s, max_gap), s) for s in "+-"), key=lambda z: z[0][1])
        if not ch:
            return piece, CopyCall(piece, partner, "none", "", 0.0, float("nan"))
        tot, cur = 0, None
        for a, b in sorted((b[0], b[1]) for b in ch):
            if cur is None or a > cur[1]:
                tot += (cur[1] - cur[0]) if cur else 0
                cur = [a, b]
            else:
                cur[1] = max(cur[1], b)
        tot += cur[1] - cur[0]
        div = sum(b[4] * b[5] for b in ch) / max(sum(b[4] for b in ch), 1)
        call = CopyCall(piece, partner, "none", strand, tot / L, div)
        gaps = [(r0 + m.start(), r0 + m.end()) for m in re.finditer("N{10,}", tseq)]
        if div <= 0.05 and ch[0][0] <= tol(L) and L - ch[-1][1] <= tol(L):
            for i in range(len(ch) - 1):
                a, b = ch[i], ch[i + 1]
                for g0, g1 in gaps:
                    if strand == "+":           # a before the gap, b after it
                        ok = a[3] <= g0 + 5000 and b[2] >= g1 - 5000 and g0 - a[3] <= tol(a[1]) and b[2] - g1 <= tol(L - b[0])
                    else:                       # the piece runs backwards: a after the gap, b before it
                        ok = a[2] >= g1 - 5000 and b[3] <= g0 + 5000 and a[2] - g1 <= tol(a[1]) and g0 - b[3] <= tol(L - b[0])
                    if ok and a[1] >= min_copy and L - b[0] >= min_copy and b[0] >= a[1]:
                        call.kind, call.gap = "gap_fill", (g0, g1)
                        call.left, call.fill, call.right = (0, a[1]), (a[1], b[0]), (b[0], L)

                        def switch(blocks, inner):
                            """Identical run nearest the fill in these copy blocks: its middle base pair."""
                            runs = [r_ for x in blocks for r_ in _match_runs(x[0], x[1], x[2], strand, x[6], 1000)]
                            if not runs:
                                return None
                            mid = lambda r_: r_[0] + (r_[2] // 2 if strand == "+" else -(r_[2] // 2))
                            qs, ts_, m = min(runs, key=lambda r_: abs(mid(r_) - inner))
                            return (mid((qs, ts_, m)), ts_ + m // 2)
                        call.left_switch = switch(ch[:i + 1], a[1])
                        call.right_switch = switch(ch[i + 1:], b[0])
                        return piece, call
        if call.aligned >= 0.7 and div <= 0.05:
            call.kind = "copy"
        return piece, call
    with ThreadPoolExecutor(max_workers=max(1, int(threads))) as ex:
        return dict(ex.map(one, cands))
