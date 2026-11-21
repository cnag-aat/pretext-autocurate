"""Displaced stretches: a stretch of a scaffold with little contact with its neighbours on the
map and a clear off-diagonal line to one place elsewhere, where it belongs (curator,
2026-11-21: "no contact signal near the diagonal; instead a horizontal or vertical line way
off diagonal, specifically and clearly indicating its proper position somewhere else").

xgTylRafi: scaffold_9 14.82-15.67 Mb has contact 0.04 / 0.07 with its two sides and 0.65 with
scaffold_10 at 56.43 Mb, where its last 340 kb are a haplotypic copy of scaffold_10 just after
the gap at 56.72 Mb; scaffold_10's own 56.48-56.72 Mb (0.06 / 0.04 with its sides) lines up
with scaffold_11 at 40.85 Mb (0.82).

`row_peaks` makes one pass over the map (every texture once): for each row, its strongest
contact in each block of columns away from the row's own neighbourhood, and the columns'
totals (bins that contact everything are masked). `scan` groups rows whose peak keeps to one
place into stretches and tests each one as a whole:

  home    contact of the stretch with `W` px on each side of it, over what adjacent sequence
          shows at those distances (the side with more contact counts);
  line    contact with `W` px on each side of the place its line points at, over the same;
  second  the same at the next best place, away from the first.

A stretch is displaced when line >= line_min, line >= move_ratio * home and line >= spec *
second. Its place: the junction (gap, cut, or weakest contact) nearest where its line
starts, and its strand from which end touches which side there.
"""
from dataclasses import dataclass, field

import numpy as np


@dataclass
class Params:
    min_px: int = 3               # shortest stretch
    max_px: int = 120             # longest stretch (a larger misplacement is a cut and join)
    excl_px: int = 150            # columns this near the row are its own neighbourhood
    bin: int = 4                  # px per bin for the row scan
    sub: int = 256                # column block (px) for the per-row maxima
    row_min: float = 0.12         # a row's peak (over the map's own contact at distance 0) to be a candidate
    drift_px: int = 40            # how far the peak may move from one row to the next within a stretch
    hole_px: int = 1              # rows below row_min tolerated inside a stretch
    W: int = 10                   # side windows for home / line contact
    line_min: float = 0.25        # line contact over what adjacent sequence shows
    move_ratio: float = 3.0       # line over home
    spec: float = 1.5             # line over the next best place
    promiscuous: float = 4.0      # column total over the median: a bin that contacts everything
    weak: float = 0.3             # contact across (5 px) below this at a local minimum: a junction


@dataclass
class Displaced:
    s: int                        # scaffold (index in the map)
    a: int                        # first pixel of the stretch
    b: int                        # last pixel
    t: int                        # the place: after this pixel
    strand: str                   # '+' or '-' as it goes in there
    home: float
    line: float
    second: float
    orient: float                 # orientation ratio (>= 1; 1 = no idea)
    left: float = 0.0             # home contact on each side
    right: float = 0.0
    reason: str = ""
    extra: dict = field(default_factory=dict)


_PEAKS = {}


def cached_peaks(pt, p=None, groups=None):
    """row_peaks of `pt`, kept for the run: a --sweep runs every level on one map, and the
    pass reads every texture (about 100 s on a 65,536 px map)."""
    p = p or Params()
    key = (getattr(pt, "path", id(pt)), pt.npix, p.bin, p.sub, p.excl_px, groups is not None)
    if key not in _PEAKS:
        _PEAKS.clear()                     # one map at a time: the arrays are large
        _PEAKS[key] = row_peaks(pt, p, groups)
    return _PEAKS[key]


def row_peaks(pt, p=None, groups=None):
    """One pass over every texture, in bins of `f` x `f` px (single pixels are too noisy: on a
    saturated MQ0 map nearly every row has some pixel at 0.9 of the diagonal somewhere).
    Returns (val, col, colsum, f): val[r, k] the highest 3 x 3-bin average of bin-row r in
    column block k, away from the row's own neighbourhood (`excl_px`, same scaffold), in raw
    map units; col[r, k] where (a bin); colsum: each bin column's total over rows outside its
    neighbourhood. `groups` (one integer per scaffold, e.g. the haplotype): a row looks only at
    columns of its own group (allhaps: every stretch has a line to its homolog in the other
    haplotype; the place it belongs is in its own)."""
    p = p or Params()
    R, n, nt = pt.tex_res, pt.npix, pt.n_tex
    f = p.bin if R % p.bin == 0 else 1
    rb, nb = R // f, n // f
    sub = max(1, min(p.sub // f, rb))
    kb = nb // sub
    val = np.zeros((nb, kb), dtype=np.float32)
    col = np.zeros((nb, kb), dtype=np.int32)
    colsum = np.zeros(nb, dtype=np.float64)
    sc = pt.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, n - 1)]
    grp = np.asarray(groups)[sc] if groups is not None else None
    ex = p.excl_px // f

    def shift_sum(m, axis):
        """m plus its neighbours along `axis` (edges: the missing neighbour counts as m)."""
        out = m.copy()
        lo = [slice(None)] * 2
        hi = [slice(None)] * 2
        lo[axis], hi[axis] = slice(1, None), slice(None, -1)
        out[tuple(lo)] += m[tuple(hi)]
        out[tuple(hi)] += m[tuple(lo)]
        e0 = [slice(None)] * 2
        e1 = [slice(None)] * 2
        e0[axis], e1[axis] = 0, -1
        out[tuple(e0)] += m[tuple(e0)]
        out[tuple(e1)] += m[tuple(e1)]
        return out / 3

    def take(rows0, cols0, m):
        """Bin-rows rows0.. of binned block m (bin-columns cols0..): update val / col / colsum."""
        r = np.arange(rows0, rows0 + m.shape[0])
        c = np.arange(cols0, cols0 + m.shape[1])
        if abs(rows0 - cols0) <= ex + m.shape[0]:
            near = (np.abs(r[:, None] - c[None, :]) <= ex) & (sc[r][:, None] == sc[c][None, :])
            if near.any():
                m = np.where(near, 0, m)
        if grp is not None:
            m = np.where(grp[r][:, None] == grp[c][None, :], m, 0)
        colsum[c] += m.sum(axis=0)
        sm = shift_sum(shift_sum(m, 0), 1)          # a line, not a bright pixel
        if abs(rows0 - cols0) <= ex + m.shape[0] + 1:
            sm = np.where(np.abs(r[:, None] - c[None, :]) <= ex + 1, 0, sm)
        blk = sm.reshape(m.shape[0], m.shape[1] // sub, sub)
        a = blk.argmax(axis=2)
        k0 = cols0 // sub
        val[r, k0:k0 + blk.shape[1]] = np.take_along_axis(blk, a[:, :, None], 2)[:, :, 0]
        col[r, k0:k0 + blk.shape[1]] = cols0 + np.arange(blk.shape[1])[None, :] * sub + a

    for x in range(nt):
        for y in range(x, nt):
            m = pt.read_tex(x, y).astype(np.float32).reshape(rb, f, rb, f).mean(axis=(1, 3))
            take(x * rb, y * rb, m)
            if x != y:
                take(y * rb, x * rb, m.T)
    return val, col, colsum, f


def _runs(best, bcol, sc_bin, p, f):
    """Bin-rows whose peak reaches row_min, grouped into runs: consecutive (holes of at most
    `hole_px` bins), within one scaffold, the peak in one target scaffold and moving at most
    `drift_px` from one row to the next. Returns [(first bin, last bin)]."""
    cand = best >= p.row_min
    tsc = sc_bin[bcol]
    out, cur, holes = [], None, 0
    drift = max(1, p.drift_px // f)
    for i in range(len(best)):
        if cand[i]:
            if cur is not None and sc_bin[i] == sc_bin[cur[0]] and tsc[i] == tsc[cur[1]] \
                    and abs(int(bcol[i]) - int(bcol[cur[1]])) <= drift * (1 + holes):
                cur[1] = i
            else:
                if cur is not None:
                    out.append(tuple(cur))
                cur = [i, i]
            holes = 0
        elif cur is not None:
            holes += 1
            if holes > p.hole_px or sc_bin[i] != sc_bin[cur[0]]:
                out.append(tuple(cur))
                cur, holes = None, 0
    if cur is not None:
        out.append(tuple(cur))
    return out


def scan(pt, sg, p=None, peaks=None, junctions=None, gap=None, groups=None):
    """Displaced stretches of `pt` (see the module). `sg`: its Signal; `peaks`: row_peaks
    output if already computed; `junctions`: extra pixels after which a stretch may go in
    (cuts already made); `gap`: the gap track (default: the map's). Returns [Displaced]."""
    p = p or Params()
    if peaks is None:
        peaks = row_peaks(pt, p, groups)
    val, col, colsum, f = peaks
    n = pt.npix
    nb = n // f
    g0 = float(sg.decay[0])
    sc = pt.orig_scaffold
    sc_bin = sc[np.minimum(np.arange(nb) * f + f // 2, n - 1)]
    med = float(np.median(colsum[colsum > 0])) if (colsum > 0).any() else 0.0
    prom = colsum > p.promiscuous * med
    v = np.where(prom[col], 0.0, val) / g0
    k = v.argmax(axis=1)
    best = v[np.arange(nb), k]
    bcol = col[np.arange(nb), k]
    gap = pt.graphs.get("gap") if gap is None else gap
    ins = sg.insulation(5, 0)
    junc = set(junctions or [])
    W = p.W
    g = sg.decay_long(p.max_px + 2 * W + 4)
    start, npx = pt.scaffold_start, pt.scaffold_npix

    def expect(L, w):
        """Mean contact a block of L rows x w columns shows next to them, from the decay."""
        i, j = np.arange(L)[:, None], np.arange(w)[None, :]
        return float(g[np.minimum(i + j + 1, len(g) - 1)].mean())

    def is_junction(q):
        """Pixel q ends a scaffold, sits at a gap, has an applied cut after it, or is a local
        minimum of contact across below `weak` (a misjoin the cuts have not taken)."""
        s = sc[q]
        if q == start[s] + npx[s] - 1 or q in junc or (gap is not None and (gap[q] > 0 or (q + 1 < n and gap[q + 1] > 0))):
            return True
        x = ins[max(0, q - 2):q + 3]
        return bool(np.isfinite(ins[q]) and ins[q] < p.weak and ins[q] <= np.nanmin(x))

    def snap_end(q, lo, hi):
        """A stretch boundary (after pixel q) moved to a gap within lo..hi, else to the point of
        least contact across there."""
        cands = [x for x in range(lo, hi + 1) if 0 <= x < n - 1 and sc[x] == sc[q] and sc[x + 1] == sc[q]]
        if not cands:
            return q
        gp = [x for x in cands if gap is not None and gap[x] > 0]
        if gp:
            return min(gp, key=lambda x: abs(x - q))
        w = np.array([ins[x] if np.isfinite(ins[x]) else np.inf for x in cands])
        return cands[int(w.argmin())] if np.isfinite(w).any() else q

    def place(a, b, w0, w1, gone):
        """The junction in columns w0..w1 where rows a..b touch most (both sides, W px each,
        skipping the pixels of `gone` stretches): (after pixel, line contact over expectation,
        (strand, orientation ratio)); (None, 0, None) without a junction there."""
        L = b - a + 1
        lo, hi = max(0, w0 - 3 * W - L), min(n - 1, w1 + 3 * W + L)
        keep = np.ones(hi - lo + 1, bool)
        for ga, gb in gone:
            keep[max(ga, lo) - lo:min(gb, hi) - lo + 1] = False
        idx = lo + np.nonzero(keep)[0]                 # columns left once `gone` is out
        B = sg.block(a, b, lo, hi).astype(np.float64)[:, keep]
        cs = np.concatenate([np.zeros((B.shape[0], 1)), B.cumsum(axis=1)], axis=1)
        ex = expect(L, W)
        h = max(1, L // 3)
        best_t, best_line, best_o = None, 0.0, None
        for c in np.nonzero((idx >= w0) & (idx <= w1))[0]:
            t = int(idx[c])
            if not is_junction(t) or c + 1 - W < 0 or c + 1 + W > B.shape[1]:
                continue
            lw = (cs[:, c + 1] - cs[:, c + 1 - W])     # per row: the W columns up to t
            rw = (cs[:, c + 1 + W] - cs[:, c + 1])     # the W columns after it
            ln = (lw.mean() + rw.mean()) / (2 * W) / ex
            if ln > best_line:
                hl, tl = lw[:h].mean(), lw[-h:].mean()
                hr, tr = rw[:h].mean(), rw[-h:].mean()
                plus, minus = hl + tr, hr + tl
                best_t, best_line = t, ln
                best_o = ("+" if plus >= minus else "-", max(plus, minus) / max(min(plus, minus), 1e-9))
        return best_t, best_line, best_o

    out = []
    for r0, r1 in _runs(best, bcol, sc_bin, p, f):
        s = int(sc_bin[r0])
        s0, s1 = int(start[s]), int(start[s] + npx[s] - 1)
        a = snap_end(max(s0, r0 * f) - 1, r0 * f - f - 1, r0 * f + f - 1) + 1
        b = snap_end(min(s1, r1 * f + f - 1), r1 * f, r1 * f + 2 * f - 1)
        L = b - a + 1
        if L < p.min_px or L > p.max_px or a - s0 < W or s1 - b < W:
            continue                       # too short, too long, or at a scaffold end (joins handle those)
        ex = expect(L, W)
        left = float(sg.block(a, b, a - W, a - 1).mean()) / ex
        right = float(sg.block(a, b, b + 1, b + W).mean()) / ex
        home = max(left, right)
        # where the line points: the peak columns of the run, with room either side
        cols = bcol[r0:r1 + 1] * f
        ts = int(sc[int(np.median(cols))])
        t0, t1 = int(start[ts]), int(start[ts] + npx[ts] - 1)
        w0, w1 = max(t0, int(cols.min()) - L - W), min(t1, int(cols.max()) + f + L + W)
        if ts == s and not (w1 < a - W or w0 > b + W):
            w0, w1 = (w0, min(w1, a - W - 1)) if cols.mean() < a else (max(w0, b + W + 1), w1)
        w0, w1 = max(w0, W), min(w1, n - 1 - W)     # the side windows stay on the map
        if w1 - w0 < 2 * W:
            continue
        best_t, best_line, best_o = place(a, b, w0, w1, ())
        if best_t is None:
            continue
        # the next best place: the rows' peaks outside the window, against inside it
        rows = v[r0:r1 + 1]
        cb = col[r0:r1 + 1] * f
        inside = (cb >= w0 - W) & (cb <= w1 + W)
        if ts == s:
            inside |= (cb >= a - p.excl_px) & (cb <= b + p.excl_px)
        rin = np.where(inside, rows, 0).max(axis=1).mean()
        rout = np.where(inside, 0, rows).max(axis=1).mean()
        second = best_line * rout / max(rin, 1e-9)
        d = Displaced(s=s, a=a, b=b, t=best_t, strand=best_o[0], home=home, line=best_line, second=second,
                      orient=best_o[1], left=left, right=right)
        d.extra = dict(row_in=float(rin), row_out=float(rout), window=(w0, w1))
        if best_line >= p.line_min and best_line >= p.move_ratio * home and rin >= p.spec * rout:
            out.append(d)
    # Places again, with every displaced stretch out of the way: a stretch that belongs where
    # another one sits now meets that one's neighbours (xgTylRafi: scaffold_9 14.82-15.70 goes
    # where scaffold_10 56.48-56.74 is, and that one belongs in scaffold_11)
    gone = [(d.a, d.b) for d in out]
    for d in out:
        w0, w1 = d.extra["window"]
        t, ln, o = place(d.a, d.b, w0, w1, [x for x in gone if x != (d.a, d.b)])
        if t is not None:
            d.t, d.line, d.strand, d.orient = t, ln, o[0], o[1]
    return out


@dataclass
class StretchAlignment:
    """A displaced stretch aligned to where its line points (`align_stretch`). Positions are
    0-based bp of the input scaffolds; copies are [start, end) on the stretch's scaffold with
    their [start, end) on the target scaffold."""
    strand: str               # '+' / '-': the stretch along the target ('' if nothing aligned)
    weight: int               # bp in the chain
    aligned: float            # share of the stretch in the chain
    divergence: float
    head: tuple = None        # (q0, q1, t0, t1): a copy at the stretch's start, or None
    tail: tuple = None        # the same at its end
    head_switch: tuple = None  # (stretch base, identical target base) nearest the copy's inner end
    tail_switch: tuple = None


def align_stretch(fasta, name, a_bp, b_bp, target, t0_bp, t1_bp, minimap2=None, max_gap=2_000_000,
                  min_copy=10_000, end_tol=20_000, max_div=0.05):
    """The stretch name:[a_bp, b_bp) aligned (minimap2 asm20) to target:[t0_bp, t1_bp), the
    heaviest collinear chain in either strand (overlap._chain). A copy at an end: chain blocks
    from within `end_tol` of that end, at most `max_div` diverged, `min_copy` bp or more.
    Returns a StretchAlignment, or None without minimap2."""
    import os
    import subprocess
    import tempfile
    from .graph import Faidx
    from .homology import find_tool
    from .overlap import _chain, _match_runs
    mm = find_tool("minimap2", minimap2)
    if mm is None:
        return None
    fa = Faidx(fasta)
    t0_bp, t1_bp = max(0, int(t0_bp)), min(fa.length(target), int(t1_bp))
    a_bp, b_bp = max(0, int(a_bp)), min(fa.length(name), int(b_bp))
    L = b_bp - a_bp
    with tempfile.TemporaryDirectory() as d:
        q, t = os.path.join(d, "q.fa"), os.path.join(d, "t.fa")
        with open(q, "w") as f:
            f.write(f">q\n{fa.fetch(name, a_bp, b_bp)}\n")
        with open(t, "w") as f:
            f.write(f">t\n{fa.fetch(target, t0_bp, t1_bp)}\n")
        r = subprocess.run([mm, "-x", "asm20", "-c", "--cs", "--secondary=no", "-t", "1", t, q],
                           capture_output=True, text=True)
    bl = {"+": [], "-": []}
    for line in r.stdout.splitlines():
        x = line.split("\t")
        if len(x) < 12 or int(x[10]) < 5000:
            continue
        de = next((float(y[5:]) for y in x[12:] if y.startswith("de:f:")), float("nan"))
        cs = next((y[5:] for y in x[12:] if y.startswith("cs:Z:")), "")
        bl[x[4]].append((int(x[2]), int(x[3]), int(x[7]) + t0_bp, int(x[8]) + t0_bp, int(x[10]), de, cs))
    (ch, w), strand = max(((_chain(bl[s], s, max_gap), s) for s in "+-"), key=lambda z: z[0][1])
    if not ch:
        return StretchAlignment("", 0, 0.0, float("nan"))
    cov = 0
    cur = None
    for x0, x1 in sorted((c[0], c[1]) for c in ch):
        if cur is None or x0 > cur[1]:
            cov += (cur[1] - cur[0]) if cur else 0
            cur = [x0, x1]
        else:
            cur[1] = max(cur[1], x1)
    cov += cur[1] - cur[0]
    div = sum(c[4] * c[5] for c in ch) / max(sum(c[4] for c in ch), 1)
    out = StretchAlignment(strand, int(w), cov / max(L, 1), div)

    def end_copy(blocks):
        """Blocks running on from the first one with gaps under 20 kb: (q0, q1, t0, t1) or None."""
        if not blocks:
            return None, []
        run = [blocks[0]]
        for c in blocks[1:]:
            prev = run[-1]
            # the next block on along both sequences (blocks come in either order)
            dq = max(c[0] - prev[1], prev[0] - c[1])
            dt = max(c[2] - prev[3], prev[2] - c[3])
            if dq > 20_000 or dt > 20_000:
                break
            run.append(c)
        q0, q1 = min(c[0] for c in run), max(c[1] for c in run)
        t0, t1 = min(c[2] for c in run), max(c[3] for c in run)
        dv = sum(c[4] * c[5] for c in run) / max(sum(c[4] for c in run), 1)
        if q1 - q0 < min_copy or dv > max_div:
            return None, []
        return (q0 + a_bp, q1 + a_bp, t0, t1), run

    def switch(blocks, inner):
        """Identical run nearest the stretch base `inner` (local): (stretch base, target base)."""
        runs = [r_ for x in blocks for r_ in _match_runs(x[0], x[1], x[2], strand, x[6], 1000)]
        if not runs:
            return None
        mid = lambda r_: r_[0] + (r_[2] // 2 if strand == "+" else -(r_[2] // 2))
        qs, ts_, m = min(runs, key=lambda r_: abs(mid(r_) - inner))
        return (mid((qs, ts_, m)) + a_bp, ts_ + m // 2)

    if ch[0][0] <= end_tol:
        out.head, run = end_copy(ch)
        if out.head:
            out.head_switch = switch(run, out.head[1] - a_bp)
    if L - ch[-1][1] <= end_tol:
        out.tail, run = end_copy(ch[::-1])
        if out.tail:
            out.tail_switch = switch(run, out.tail[0] - a_bp)
    return out
