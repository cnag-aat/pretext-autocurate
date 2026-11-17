"""Order and orientation of the pieces of a chromosome copy: end contact
(`arrange_pieces`, `_adjacency`) and the other haplotype's alignment (`_project`,
`_aligned_dir`, `_contig_conflicts`)."""
import numpy as np


def _adjacency(sg, g, P, rP, Q, rQ, window=100):
    """Contact between the right end of piece P and the left end of piece Q as laid out
    (pieces are inclusive pixel ranges; r = reversed), relative to what adjacent sequence
    shows at the same distances. The outermost pixels are skipped."""
    def win(piece, tail):
        """End window of a piece: (first px, last px, skipped px K, at the tail). K skips the
        outermost pixels (gap and edge effects); W is up to `window` px, half the piece."""
        L = piece[1] - piece[0] + 1
        K = min(10, L // 10)
        W = max(1, min(window, L // 2 - K))
        return (piece[1] - K - W + 1, piece[1] - K, K, True) if tail else (piece[0] + K, piece[0] + K + W - 1, K, False)
    a0, a1, k1, t1 = win(P, not rP)
    b0, b1, k2, t2 = win(Q, rQ)
    B = sg.block(a0, a1, b0, b1).astype(np.float64)
    # i, j: each pixel's distance from the joining end of its piece, so i + j + 1 + K's is
    # its distance across the junction if the pieces were adjacent
    i, j = np.arange(a1 - a0 + 1), np.arange(b1 - b0 + 1)
    i, j = (i[::-1] if t1 else i), (j[::-1] if t2 else j)
    expect = g[np.minimum(i[:, None] + j[None, :] + 1 + k1 + k2, len(g) - 1)]
    return float(B.sum() / max(expect.sum(), 1e-9))


def _homolog_dir(sg, a, b, cols, max_rows=120, k=6):
    """Direction of the homology line of rows a..b within the column ranges `cols` (the other
    copy of the chromosome): +1 if the homolog runs the same way, -1 if reversed, 0 if no
    line stands out. Sums the map along diagonals and anti-diagonals and takes the sharpest
    one (a line is a few pixels wide; cis contact is broad and cancels in the contrast)."""
    L = b - a + 1
    if L > max_rows:
        a += (L - max_rows) // 2
        b = a + max_rows - 1
    W = b - a + 1
    St = sg.strip(a, b).astype(np.float32)
    best = {1: 0.0, -1: 0.0}
    noise = False
    for lo, hi in cols:
        M = hi - lo + 1 - W + 1                 # offsets at which W rows fit along the column range
        if M <= 4 * k:
            continue
        # F[o]: mean along the diagonal starting at column lo+o; R[o]: the anti-diagonal
        F, R = np.zeros(M, np.float32), np.zeros(M, np.float32)
        for i in range(W):
            F += St[i, lo + i:lo + i + M]
            R += St[i, lo + W - 1 - i:lo + W - 1 - i + M]
        for d, D in ((1, F / W), (-1, R / W)):
            c = D[k:-k] - 0.5 * (D[:-2 * k] + D[2 * k:])     # contrast with diagonals k px away
            best[d] = max(best[d], float(c.max()))
            noise = True
    if not noise:
        return 0
    d = 1 if best[1] >= best[-1] else -1
    if best[d] < 2.0 * best[-d]:                # the other direction shows what "no line" gives here
        return 0
    return d


def _block_index(blocks):
    """{scaffold name: [(start, end, other name, other start, other end, same strand)]} with
    every alignment entered from both of its sides."""
    idx = {}
    for bk in blocks:
        if bk.mapq < 30:
            continue
        idx.setdefault(bk.q, []).append((bk.qs, bk.qe, bk.t, bk.ts, bk.te, bk.same_strand))
        idx.setdefault(bk.t, []).append((bk.ts, bk.te, bk.q, bk.qs, bk.qe, bk.same_strand))
    return idx


def _aligned_dir(raw, idx, piece, cols, min_bp=100000, min_frac=0.1, clear=0.7):
    """Direction of a piece relative to the other copy from sequence alignments: +1 if at
    least 70% of its bases aligned to the pixel ranges `cols` are on the same strand, -1 if
    at most 30%, 0 if too little aligns or the strands are mixed."""
    bpp = raw.bp_per_pix
    name = raw.names[raw.orig_scaffold[piece[0]]]
    lo, hi = raw.rel_coord[piece[0]] * bpp, (raw.rel_coord[piece[1]] + 1) * bpp
    targets = {}
    for a, b in cols:
        targets.setdefault(raw.names[raw.orig_scaffold[a]], []).append((raw.rel_coord[a] * bpp, (raw.rel_coord[b] + 1) * bpp))
    tot = same = 0.0
    for s0, e0, other, os_, oe_, st in idx.get(name, []):
        ov = min(e0, hi) - max(s0, lo)
        if ov <= 0 or other not in targets:
            continue
        # the part of the other side that corresponds to the overlap (alignments are long)
        f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
        if not st:
            f0, f1 = 1 - f1, 1 - f0
        o0, o1 = os_ + f0 * (oe_ - os_), os_ + f1 * (oe_ - os_)
        if any(min(o1, t1) - max(o0, t0) > 0 for t0, t1 in targets[other]):
            tot += ov
            same += ov if st else 0.0
    if tot < max(min_bp, min_frac * (hi - lo)):
        return 0
    return 1 if same >= clear * tot else -1 if same <= (1 - clear) * tot else 0


def _aligned_bp(raw, idx, a, b):
    """Bases of pixel range `a` (inclusive) that align into pixel range `b`, from a block
    index (`_block_index`). Each alignment counts for its overlap with `a` when the
    corresponding stretch of the other side overlaps `b` (alignments are long chains)."""
    bpp = raw.bp_per_pix
    na, nb = raw.names[raw.orig_scaffold[a[0]]], raw.names[raw.orig_scaffold[b[0]]]
    lo, hi = raw.rel_coord[a[0]] * bpp, (raw.rel_coord[a[1]] + 1) * bpp
    t0, t1 = raw.rel_coord[b[0]] * bpp, (raw.rel_coord[b[1]] + 1) * bpp
    tot = 0.0
    for s0, e0, other, os_, oe_, st in idx.get(na, []):
        ov = min(e0, hi) - max(s0, lo)
        if other != nb or ov <= 0:
            continue
        f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
        if not st:
            f0, f1 = 1 - f1, 1 - f0
        if min(os_ + f1 * (oe_ - os_), t1) > max(os_ + f0 * (oe_ - os_), t0):
            tot += ov
    return tot


def _project(raw, idx, piece, layout, min_bp=100000, min_frac=0.1, clear=0.7):
    """Where a piece lies along the other copy as it is laid out, from sequence alignments.
    `layout`: the other copy's pieces in order, (first pixel, last pixel, reversed).
    Returns (position in layout pixels, direction +1 / -1, the aligned stretches of the layout)
    or None when too little of the piece aligns there or the strands are mixed."""
    bpp = raw.bp_per_pix
    # the piece in bp of its original scaffold (alignments are in scaffold coordinates)
    name = raw.names[raw.orig_scaffold[piece[0]]]
    lo, hi = raw.rel_coord[piece[0]] * bpp, (raw.rel_coord[piece[1]] + 1) * bpp
    # ref: per scaffold of the other copy, its pieces in the layout as (start bp, end bp,
    # offset of the piece in layout pixels, reversed)
    ref, off = {}, 0
    for a, b, rev in layout:
        ref.setdefault(raw.names[raw.orig_scaffold[a]], []).append((raw.rel_coord[a] * bpp, (raw.rel_coord[b] + 1) * bpp, off, rev))
        off += b - a + 1
    tot = same = pos = 0.0
    segs = []
    for s0, e0, other, os_, oe_, st in idx.get(name, []):
        ov = min(e0, hi) - max(s0, lo)
        if ov <= 0 or other not in ref:
            continue
        # the stretch of the other scaffold that corresponds to the overlap, by linear
        # interpolation along the block (mirrored on the reverse strand)
        f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
        if not st:
            f0, f1 = 1 - f1, 1 - f0
        o0, o1 = os_ + f0 * (oe_ - os_), os_ + f1 * (oe_ - os_)
        for t0, t1, off_k, rev_k in ref[other]:
            w = min(o1, t1) - max(o0, t0)
            if w <= 0:
                continue
            # bp of the other scaffold -> layout pixels, counting from the reversed end if the
            # layout shows that piece reversed
            to_px = lambda x: off_k + ((t1 - x) if rev_k else (x - t0)) / bpp
            mid = (max(o0, t0) + min(o1, t1)) / 2
            pos += w * to_px(mid)                  # position: mean weighted by aligned bp
            tot += w
            # runs the same way as the layout: same strand on a forward piece, or opposite
            # strand on a reversed one
            same += w if st != rev_k else 0.0
            segs.append(tuple(sorted((to_px(max(o0, t0)), to_px(min(o1, t1))))))
    # too little aligned (under the larger of min_bp and min_frac of the piece), or strands
    # mixed: no answer
    if tot < max(min_bp, min_frac * (hi - lo)):
        return None
    d = 1 if same >= clear * tot else -1 if same <= (1 - clear) * tot else 0
    return (pos / tot, d, segs) if d else None


def _contig_conflicts(raw, idx, pieces, layout, min_bp=50000):
    """How much contiguous sequence of one copy contradicts the other copy's layout.

    Curator's rule (2026-11-14): contigs outrank scaffold joins. The pieces of a copy meet at
    gaps (scaffold joins or cut sites), so a layout's piece order is only as good as YaHS or
    this tool; a contig is sequence the assembler put together. Each contig (gap-free stretch,
    from the map's `gap` track) of `pieces` (one copy, inclusive pixel ranges) is followed
    through its alignments into `layout` (the other copy, [(first px, last px, reversed)]).
    Where one contig runs across two pieces of the layout, it must meet them in layout order
    and orientation; if it does not, the layout is contradicted there.

    Returns (contradicted bp, [(contig as 'scaffold:start-end Mb', layout pieces involved)]).
    Only alignment stretches of at least `min_bp` count (repeats and small indels are noise);
    rearrangements inside one layout piece are not counted (they say nothing about the
    order of pieces)."""
    gap = raw.graphs.get("gap")
    bpp = raw.bp_per_pix
    # the layout per scaffold of the other copy, as in `_project`, with each piece's rank k
    ref, off = {}, 0
    for k, (a, b, rev) in enumerate(layout):
        ref.setdefault(raw.names[raw.orig_scaffold[a]], []).append(
            (raw.rel_coord[a] * bpp, (raw.rel_coord[b] + 1) * bpp, off, rev, k))
        off += b - a + 1
    bad, where = 0.0, []
    for a, b in pieces:
        name = raw.names[raw.orig_scaffold[a]]
        # contigs: runs of non-gap pixels inside the piece
        stops = [a - 1] + ([q for q in range(a, b + 1) if gap[q] > 0] if gap is not None else []) + [b + 1]
        for c0, c1 in zip(stops, stops[1:]):
            if c1 - c0 < 2:
                continue
            lo, hi = raw.rel_coord[c0 + 1] * bpp, (raw.rel_coord[c1 - 1] + 1) * bpp
            pts = []                    # (position along the contig bp, layout px, piece k, same way, bp)
            for s0, e0, other, os_, oe_, st in idx.get(name, []):
                if min(e0, hi) - max(s0, lo) < min_bp or other not in ref:
                    continue
                f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
                if not st:
                    f0, f1 = 1 - f1, 1 - f0
                o0, o1 = os_ + f0 * (oe_ - os_), os_ + f1 * (oe_ - os_)
                for t0, t1, off_k, rev_k, k in ref[other]:
                    w = min(o1, t1) - max(o0, t0)
                    if w < min_bp:
                        continue
                    mid = (max(o0, t0) + min(o1, t1)) / 2
                    g = (mid - os_) / max(oe_ - os_, 1)         # back along the block to this contig
                    d_pos = s0 + (g if st else 1 - g) * (e0 - s0)
                    l_px = off_k + ((t1 - mid) if rev_k else (mid - t0)) / bpp
                    pts.append((d_pos, l_px, k, st != rev_k, w))
            if len({x[2] for x in pts}) < 2:
                continue                # the contig meets one layout piece: nothing to say about order
            pts.sort()
            # runs: consecutive stretches in the same layout piece
            runs = []
            for d_pos, l_px, k, fwd, w in pts:
                if runs and runs[-1]["k"] == k:
                    r = runs[-1]
                    r["w"] += w
                    r["fw"] += w if fwd else 0.0
                    r["last"] = l_px
                else:
                    runs.append(dict(k=k, w=w, fw=w if fwd else 0.0, first=l_px, last=l_px))
            # the contig's direction along the layout: by most aligned bp
            tot = sum(r["w"] for r in runs)
            sign = 1 if sum(r["fw"] for r in runs) >= tot / 2 else -1
            hit = 0.0
            for r in runs:              # a piece turned against the contig's direction
                if (1 if r["fw"] >= r["w"] / 2 else -1) != sign:
                    hit += r["w"]
            for r1, r2 in zip(runs, runs[1:]):  # two pieces met out of layout order
                if r1["k"] != r2["k"] and (r2["first"] - r1["last"]) * sign < 0:
                    hit += min(r1["w"], r2["w"])
            if hit > 0:
                bad += hit
                where.append((f"{name}:{lo / 1e6:.1f}-{hi / 1e6:.1f} Mb",
                              sorted({raw.names[raw.orig_scaffold[layout[r['k']][0]]] for r in runs})))
    return bad, where


def arrange_pieces(sg, pieces, init, max_n=9, orient_gain=0.12, guide=None):
    """Order and orientation of the pieces of one chromosome copy that give the most contact
    between neighbouring ends (exact search over all orders and orientations).
    `pieces`: inclusive pixel ranges; `init`: the arrangement to compare with, a list of
    (index, reversed). Returns (best arrangement, its score, score of init, link scores,
    pieces whose orientation was too close to call and was left as it was, pieces oriented
    by `guide`).

    Scores are sums of `_adjacency` over neighbouring ends (contact relative to adjacent
    sequence). `max_n`: above this many pieces init is returned unsearched (the search is
    O(2^n n^2)); `orient_gain`: least loss of score from flipping a piece for its
    orientation to count as clear; `guide(i)`: +1 / -1 / 0, direction of piece i's homolog."""
    n = len(pieces)
    g = sg.decay_long(230)
    # link score for every ordered pair of piece ends, in both orientations
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                for ri in (False, True):
                    for rj in (False, True):
                        adj[(i, ri, j, rj)] = _adjacency(sg, g, pieces[i], ri, pieces[j], rj)
    total = lambda arr: sum(adj[(a, ra, b, rb)] for (a, ra), (b, rb) in zip(arr, arr[1:]))
    if n > max_n or n < 2:
        return list(init), total(init), total(init), adj, [], []
    # Held-Karp style DP over subsets: the best path through the pieces in `mask` ending with
    # piece `last` in a given orientation. Masks are visited in increasing order, so every
    # subset is complete before it is extended.
    best = {}                                   # (mask, last, reversed) -> (score, previous state)
    for i in range(n):
        for r in (False, True):
            best[(1 << i, i, r)] = (0.0, None)
    for mask in range(1, 1 << n):
        for i in range(n):
            if not mask >> i & 1:
                continue
            for r in (False, True):
                cur = best.get((mask, i, r))
                if cur is None:
                    continue
                for j in range(n):
                    if mask >> j & 1:
                        continue
                    for rj in (False, True):
                        key = (mask | 1 << j, j, rj)
                        v = cur[0] + adj[(i, r, j, rj)]
                        if key not in best or v > best[key][0]:
                            best[key] = (v, (mask, i, r))
    full = (1 << n) - 1
    end = max(((best[(full, i, r)][0], i, r) for i in range(n) for r in (False, True)))
    arr, state = [], (full, end[1], end[2])     # walk back from the best full path
    while state is not None:
        arr.append((state[1], state[2]))
        state = best[state][1]
    arr.reverse()
    # a path and its reverse score the same (contact is symmetric); pick the one that leaves
    # the largest piece as it was, so a copy is not turned round for nothing
    big = max(range(n), key=lambda i: pieces[i][1] - pieces[i][0])
    if dict(arr)[big] != dict(init)[big]:       # keep the largest piece the way it was
        arr = [(i, not r) for i, r in reversed(arr)]
    # A piece whose two orientations score about the same: end contact cannot orient it.
    # Then the other haplotype is the guide (`guide(i)`: +1 / -1 / 0, direction of piece i's
    # homology line in the other copy): the piece is turned so that its homolog runs the same
    # way as that of the nearest piece whose orientation is clear. Without a guide it keeps
    # the orientation it had.
    r_of = dict(arr)
    flip_cost, guided = {}, []
    # flip_cost: score lost by turning only piece i round in place
    for k_, (i, r) in enumerate(arr):
        alt = arr[:k_] + [(i, not r)] + arr[k_ + 1:]
        flip_cost[i] = total(arr) - total(alt)
    unclear = [i for i, _ in arr if flip_cost[i] < orient_gain]
    if unclear and len(unclear) == n:           # nothing is clear: the largest piece is the reference
        unclear.remove(big)
    order_ = [i for i, _ in arr]
    for i in list(unclear):
        ref = min((j for j in order_ if j not in unclear), key=lambda j: abs(order_.index(j) - order_.index(i)), default=None)
        di, dj = (guide(i), guide(ref)) if guide is not None and ref is not None else (0, 0)
        if di and dj:
            # homologs running opposite ways: the piece is turned relative to the reference
            r_of[i] = r_of[ref] ^ (di != dj)
            unclear.remove(i)                   # settled by the homolog
            guided.append(i)
        else:
            r_of[i] = dict(init)[i]
    arr = [(i, r_of[i]) for i in order_]
    if dict(arr)[big] != dict(init)[big]:       # the reference may have been the one turned
        arr = [(i, not r) for i, r in reversed(arr)]
    return arr, total(arr), total(init), adj, unclear, guided
