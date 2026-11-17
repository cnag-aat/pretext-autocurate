"""Pairs that are parts of one chromosome by chromosome-level contact alone (curator,
2026-11-17): two arms whose unique sequence is pushed apart by centromeric repeats touch only
weakly, over the whole of both, not at their ends; common in beetles, hard to see in
PretextView, and missed by YaHS. Small chromosomes that cluster in the nucleus, and repeats
shared between chromosomes, give contact of their own, so the test is specificity:

  - repeat bins masked: coverage at least `cov_max` times a single copy (collapsed repeats),
    telomere signal, and the `promisc` % of bins with the most contact to other scaffolds of
    their haplotype (repeats shared by many) -- the MQ0 map is all there is (no MQ10 maps
    for allhaps);
  - each chromosome copy (with its unlocs) against every other copy of its haplotype: mean
    contact over the whole block, over the geometric mean of each one's mean with all others;
  - mutual best, `margin` times either one's second best;
  - spread: at least `spread_min` of the block's rows and of its columns carry contact (a
    hot spot -- a shared repeat, interstitial telomeric sequence -- does not count);
  - not between pieces the pass cut apart (that cut was decided on contact across), and each
    copy with its unlocs at least `min_mb` (small chromosomes cluster in the nucleus);
  - in the other haplotype the two copies' homologs are each other's best too, by
    `other_margin` (one haplotype is often more fragmented, its signal partly in pieces not
    yet paired); either haplotype may be the one that passes in full.

Troglobisium (measured 2026-11-17 on whole scaffolds): every mutual-best pair with spread
>= 37% was a curator's join (or the curator's cross-haplotype pair), every one <= 10% a hot
spot. In the pass at level 0.5 it merges H2_scaffold_18 / H1_scaffold_23 with
H2_scaffold_66 / H1_scaffold_62 (curator's chromosome; the other joins are found earlier by
end contact, whole-block contact or the graph). On the owl, before the cut and size guards it
rejoined a cut the curator also made and two clustering microchromosomes; with them, nothing."""
import numpy as np

from .common import _level


def candidates(pt, pairs, hap, extra, margin=2.0, spread_min=0.25, cov_max=1.5, promisc=2.0, min_px=60, debug=None,
               other_margin=1.5, min_mb=8.0):
    """Returns [dict(P, Q, x, y, ratio, second, spread, ratio_other, spread_other)]: pairs P
    and Q (kind 'pair') whose copies x (of P) and y (of Q), one haplotype, and their mates in
    the other, pass every test. `extra`: {copy: [units counted with it]} (its unlocs, the
    units already merged into it). Each pair appears at most once (mutual best)."""
    lv = _level(pt)
    base = getattr(pt, "pt", pt)
    if ("_mip", lv) not in base.__dict__:
        base.__dict__[("_mip", lv)] = base.matrix(lv).astype(np.float32)
    M = base.__dict__[("_mip", lv)]
    nb = M.shape[0]
    f = pt.npix // nb
    start, npx = pt.scaffold_start, pt.scaffold_npix
    unit_of = pt.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, pt.npix - 1)]
    bhap = np.array([hap[int(u)] for u in unit_of])
    # repeat bins
    mask = np.zeros(nb, dtype=bool)
    tel, cov = pt.graphs.get("telomere"), pt.graphs.get("coverage")
    if tel is not None:
        thr = max(50.0, 5 * float(np.percentile(tel, 99)))
        mask |= tel[:nb * f].reshape(nb, f).max(1) >= thr
    if cov is not None:
        c_ = cov[:nb * f].reshape(nb, f).astype(float).mean(1)
        large = npx[unit_of] >= 0.2 * npx.max()
        mask |= c_ / max(float(np.median(c_[large])), 1e-9) >= cov_max
    if promisc > 0:
        trans = np.zeros(nb)
        for h in set(bhap.tolist()):
            if not h:
                continue
            sel = np.nonzero(bhap == h)[0]
            sub = M[np.ix_(sel, sel)]
            same = unit_of[sel][:, None] == unit_of[sel][None, :]
            trans[sel] = np.where(same, 0.0, sub).sum(1)
        mask |= trans >= np.percentile(trans, 100 - promisc)
    keep = ~mask

    def bins(u):
        """Unmasked bins wholly inside unit u (rounded inwards)."""
        a = -(-int(start[u]) // f)
        b = max(int(start[u] + npx[u]) // f, a + 1)
        r = np.arange(a, min(b, nb))
        return r[keep[r]]
    main = [p for p in pairs if p.kind == "pair" and min(npx[p.a], npx[p.b]) >= min_px and hap[p.a] != hap[p.b]]
    mate = {}
    for p in main:
        mate[p.a], mate[p.b] = p.b, p.a
    pair_of = {u: p for p in main for u in (p.a, p.b)}
    group = {c: np.concatenate([bins(c)] + [bins(u) for u in extra.get(c, [])]) for c in mate}
    copies = {}
    for c in mate:
        copies.setdefault(hap[c], []).append(c)
    R, S, best = {}, {}, {}
    for h, cs in copies.items():
        mean = {}
        for i, x in enumerate(cs):
            for y in cs[i + 1:]:
                gx, gy = group[x], group[y]
                v = float(M[np.ix_(gx, gy)].mean()) if len(gx) and len(gy) else 0.0
                mean[(x, y)] = mean[(y, x)] = v
        bg = {x: max(float(np.mean([mean[(x, y)] for y in cs if y != x] or [0.0])), 1e-9) for x in cs}
        for (x, y), v in mean.items():
            R[(x, y)] = v / np.sqrt(bg[x] * bg[y])
        for x in cs:
            c = sorted(((R[(x, y)], y) for y in cs if y != x), reverse=True)
            if len(c) >= 2:
                best[x] = (c[0][1], c[0][0], c[1][0])

    def spread(x, y):
        gx, gy = group[x], group[y]
        if not len(gx) or not len(gy):
            return 0.0
        b = M[np.ix_(gx, gy)] > 0
        return float(min(b.any(1).mean(), b.any(0).mean()))

    def mutual(x, y, m):
        return (best.get(x, (None,))[0] == y and best.get(y, (None,))[0] == x
                and R[(x, y)] >= m * max(best[x][2], best[y][2]))

    def passes(x, y):
        """All tests in this haplotype."""
        if not mutual(x, y, margin):
            return False
        S[(x, y)] = spread(x, y)
        return S[(x, y)] >= spread_min

    def confirms(x, y):
        """The other haplotype: mutual best with a smaller margin (one haplotype is often more
        fragmented, its signal partly in unpaired pieces: Troglobisium H1_scaffold_62)."""
        S.setdefault((x, y), spread(x, y))
        return mutual(x, y, min(margin, other_margin))
    if debug is not None:                                  # each copy: (best, ratio, second best ratio)
        debug.update(best=best, R=R, spread=spread, mate=mate)
    raw = getattr(pt, "pt", pt)
    members = lambda c: [c] + list(extra.get(c, []))
    origin = lambda c: {int(raw.orig_scaffold[int(start[u])]) for u in members(c)}
    size_mb = lambda c: sum(int(npx[u]) for u in members(c)) * pt.bp_per_pix / 1e6
    out, seen = [], set()
    for x in mate:
        y = best.get(x, (None,))[0]
        if y is None or pair_of[x] is pair_of[y]:
            continue
        # never against a cut the pass made (it was made on contact across: two chromosomes;
        # owl H1/H2_scaffold_1 at 134.5 Mb), and not between small copies, which cluster in the
        # nucleus (owl microchromosomes H1_scaffold_32 + H1_scaffold_24:2, 4.6 and 2.6 Mb)
        if origin(x) & origin(y) or origin(mate[x]) & origin(mate[y]):
            continue
        if min(size_mb(x), size_mb(y), size_mb(mate[x]), size_mb(mate[y])) < min_mb:
            continue
        key = frozenset((id(pair_of[x]), id(pair_of[y])))       # each pair of pairs once, from either haplotype
        if key in seen:
            continue
        seen.add(key)
        mx, my = mate[x], mate[y]
        if not ((passes(x, y) and confirms(mx, my)) or (passes(mx, my) and confirms(x, y))):
            continue
        out.append(dict(P=pair_of[x], Q=pair_of[y], x=x, y=y, ratio=R[(x, y)], second=max(best[x][2], best[y][2]),
                        spread=S[(x, y)], ratio_other=R[(mate[x], mate[y])], spread_other=S[(mate[x], mate[y])]))
    return out
