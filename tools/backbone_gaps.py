"""Does pooling the Hi-C of both copies of a chromosome ("shared backbone", hapscotch)
separate the gaps a curator cut from the gaps they kept better than one copy alone?

Usage: python3 tools/backbone_gaps.py MAP.pretext H2_vs_H1.paf CURATOR.savestate [OUT.tsv]

For every gap inside a scaffold (both sides >= 2*K+10 px) the contact across it is measured
as observed / expected at the same distances (expected: the map's decay within scaffolds),
in windows of up to W px on each side, skipping K px next to the gap:

  own      the gap's own copy
  other    the homologous windows in the other copy (through the alignment), when both
           windows map colinearly into one scaffold of the other haplotype
  min      the weaker of own / other (what "seen in both haplotypes" asks)
  pooled   (observed own + other) / (expected own + other): the shared backbone
  pooled4  pooled plus the two cross blocks (own left x other right and the reverse), each
           with its own expected level (at MQ0 reads of identical sequence fall in either copy)

The curator's call: the gap is kept if its two sides are still neighbours, the same way
round, in the curator's layout; cut otherwise. Reports, per measure, how well it ranks cut
gaps below kept ones (AUC: 1 = perfect, 0.5 = chance), on all gaps with a mapped homolog and
on those where the homologous point in the other copy is also a gap (both copies joined
there by scaffolding, so the alignment cannot confirm the join). Exploratory, not part of
the package.
"""
import os
import sys

import numpy as np

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate import Pretext, SaveState, MapState          # noqa: E402
from pretext_autocurate.homology import read_paf                     # noqa: E402
from pretext_autocurate.signal import Signal                         # noqa: E402

K, W, DMAX = 2, 60, 240

pt = Pretext(sys.argv[1])
blocks = read_paf(sys.argv[2], min_mapq=30)
st = SaveState.load(sys.argv[3])
out_path = sys.argv[4] if len(sys.argv) > 4 else None
bpp, N = pt.bp_per_pix, pt.npix
sg = Signal(pt)
g = sg.decay_long(DMAX)
idx_of = {n: i for i, n in enumerate(pt.names)}
start, npx = pt.scaffold_start, pt.scaffold_npix
hap = lambda s: pt.names[s].split("_")[0]

# --- pixel map between the copies: for every pixel, (other pixel, same strand) from the
# unique alignment blocks (both directions); -1 where nothing aligns
other_px = np.full(N, -1, dtype=np.int64)
other_fw = np.zeros(N, dtype=bool)
best_len = np.zeros(N)
for b in blocks:
    for (s, s0, s1), (t, t0, t1) in (((b.q, b.qs, b.qe), (b.t, b.ts, b.te)), ((b.t, b.ts, b.te), (b.q, b.qs, b.qe))):
        si, ti = idx_of[s], idx_of[t]
        p0, p1 = int(start[si] + s0 // bpp), int(start[si] + (s1 - 1) // bpp)
        L = s1 - s0
        for p in range(p0, p1 + 1):
            if L <= best_len[p]:
                continue
            f = ((p - start[si]) * bpp + bpp / 2 - s0) / L         # fraction along the block
            f = min(max(f, 0.0), 1.0)
            tb = t0 + f * (t1 - t0) if b.same_strand else t1 - f * (t1 - t0)
            other_px[p] = int(start[ti] + min(tb // bpp, npx[ti] - 1))
            other_fw[p] = b.same_strand
            best_len[p] = L

# --- the curator's layout: position of every raw pixel and its direction
ms = MapState(pt).replay(st.edits)
pos = np.empty(N, dtype=np.int64)
pos[ms.pix] = np.arange(N)

gap = pt.graphs.get("gap")
gp = np.nonzero(gap > 0)[0]
runs = np.split(gp, np.nonzero(np.diff(gp) > 1)[0] + 1) if len(gp) else []


def obs_exp(L, Rr, cross=None):
    """Observed / expected contact between pixel lists L (left side, nearest the gap last)
    and Rr (right side, nearest first). Expected from the decay at their distance across
    the gap (i + j + 1 + 2K + gap width) unless `cross` gives a curve."""
    a0, a1, b0, b1 = min(L), max(L), min(Rr), max(Rr)
    if a0 > b0:                       # keep the block in the upper triangle
        B = sg.block(b0, b1, a0, a1).astype(np.float64).T
    else:
        B = sg.block(a0, a1, b0, b1).astype(np.float64)
    o = B[np.ix_(np.asarray(L) - a0, np.asarray(Rr) - b0)]
    i = np.arange(len(L))[::-1][:, None]
    j = np.arange(len(Rr))[None, :]
    d = np.minimum(i + j + 1 + 2 * K, DMAX - 1)
    e = (cross if cross is not None else g)[d]
    return float(o.sum()), float(e.sum())


cross_cache = {}


def cross_curve(s, t):
    """Expected contact between a pixel of copy s and the pixel d away (in s) as seen in
    copy t, through the pixel map: the level of the homolog block at each distance."""
    if (s, t) in cross_cache:
        return cross_cache[(s, t)]
    a, n = int(start[s]), int(npx[s])
    ok = np.nonzero((other_px[a:a + n] >= 0) & (pt.orig_scaffold[np.maximum(other_px[a:a + n], 0)] == t))[0] + a
    tot, cnt = np.zeros(DMAX), np.zeros(DMAX)
    rng = np.random.default_rng(0)
    for p in rng.choice(ok, size=min(300, len(ok)), replace=False) if len(ok) else []:
        qs = np.arange(p, min(p + DMAX, a + n))
        qs = qs[other_px[qs] >= 0]
        if len(qs) < 2:
            continue
        tp = other_px[qs]
        lo, hi = int(tp.min()), int(tp.max())
        row = sg.block(p, p, lo, hi)[0].astype(np.float64) if p <= lo else sg.block(lo, hi, p, p)[:, 0].astype(np.float64)
        dd = qs - p
        tot[dd] += row[tp - lo]
        cnt[dd] += 1
    c = tot / np.maximum(cnt, 1)
    cross_cache[(s, t)] = c
    return c


rows = []
for r in runs:
    g0, g1 = int(r[0]), int(r[-1])
    s = int(pt.orig_scaffold[g0])
    a, b = int(start[s]), int(start[s] + npx[s] - 1)
    left_n, right_n = g0 - a, b - g1
    if min(left_n, right_n) < 2 * K + 10:
        continue
    w = min(W, (left_n - K) // 2, (right_n - K) // 2)
    L = list(range(g0 - K - w, g0 - K))
    Rr = list(range(g1 + 1 + K, g1 + 1 + K + w))
    ow, ew = obs_exp(L, Rr)
    # curator: sides still neighbours, same way round?
    kept = abs(int(pos[g1 + 1 + K]) - int(pos[g0 - K - 1]) - (g1 - g0 + 2 + 2 * K)) <= 2 or \
        abs(int(pos[g0 - K - 1]) - int(pos[g1 + 1 + K]) - (g1 - g0 + 2 + 2 * K)) <= 2
    rec = dict(scaffold=pt.names[s], mb=round((g0 - a) * bpp / 1e6, 2), kept=kept, w=w, own=ow / max(ew, 1e-9),
               other=np.nan, min=np.nan, pooled=np.nan, pooled4=np.nan, partner="", partner_gap=False)
    # homologous windows in the other copy
    lo, ro = other_px[L], other_px[Rr]
    if (lo >= 0).mean() >= 0.7 and (ro >= 0).mean() >= 0.7:
        t = np.bincount(pt.orig_scaffold[np.concatenate([lo[lo >= 0], ro[ro >= 0]])]).argmax()
        Lm = [x for x in lo if x >= 0 and pt.orig_scaffold[x] == t]
        Rm = [x for x in ro if x >= 0 and pt.orig_scaffold[x] == t]
        colinear = len(Lm) >= 0.7 * w and len(Rm) >= 0.7 * w and hap(t) != hap(s)
        if colinear:
            fw = np.median(Lm[-5:]) < np.median(Rm[:5])
            sideL = (max(Lm) < min(Rm)) if fw else (min(Lm) > max(Rm))
            colinear = sideL
        if colinear:
            # pairs (i, j) kept where both pixels map; distances measured in the own copy
            Li = [k for k, x in enumerate(lo) if x >= 0 and pt.orig_scaffold[x] == t]
            Rj = [k for k, x in enumerate(ro) if x >= 0 and pt.orig_scaffold[x] == t]
            o_sum = e_sum = 0.0
            for i_ in Li:                  # one row at a time: the map is not monotone pixel by pixel
                p = int(lo[i_])
                cols = ro[Rj]
                c0, c1 = int(cols.min()), int(cols.max())
                row = sg.block(p, p, c0, c1)[0] if p <= c0 else sg.block(c0, c1, p, p)[:, 0]
                o_sum += float(row[cols - c0].sum())
                dd = np.minimum((w - 1 - i_) + np.array(Rj) + 1 + 2 * K, DMAX - 1)
                e_sum += float(g[dd].sum())
            rec.update(partner=pt.names[t], other=o_sum / max(e_sum, 1e-9))
            rec["min"] = min(rec["own"], rec["other"])
            rec["pooled"] = (ow + o_sum) / max(ew + e_sum, 1e-9)
            # cross blocks: own left x other right, other left x own right
            cc = cross_curve(s, t)
            ocx = ecx = 0.0
            for (A, Bl, idxA, idxB) in ((L, ro, range(len(L)), Rj), (lo, Rr, Li, range(len(Rr)))):
                for i_ in idxA:
                    p = int(A[i_])
                    cols = np.array([Bl[k] for k in idxB])
                    c0, c1 = int(cols.min()), int(cols.max())
                    row = sg.block(p, p, c0, c1)[0] if p <= c0 else sg.block(c0, c1, p, p)[:, 0]
                    ocx += float(row[cols - c0].sum())
                    dd = np.minimum((w - 1 - i_) + np.array(list(idxB)) + 1 + 2 * K, DMAX - 1)
                    ecx += float(cc[dd].sum())
            # the cross blocks' expected is the homolog block's own level at each distance
            rec["pooled4"] = (ow + o_sum + ocx) / max(ew + e_sum + ecx, 1e-9)
            # is the homologous point a gap in the other copy too (scaffold join in both)?
            mid = sorted([int(Lm[-1]), int(Rm[0])])
            rec["partner_gap"] = bool((gap[max(mid[0] - 3, 0):mid[1] + 4] > 0).any())
    rows.append(rec)


def auc(vals, kept):
    """P(a kept gap scores above a cut gap); ties count half."""
    v, k = np.asarray(vals, float), np.asarray(kept, bool)
    ok = ~np.isnan(v)
    v, k = v[ok], k[ok]
    pos_, neg = v[k], v[~k]
    if not len(pos_) or not len(neg):
        return np.nan, len(pos_), len(neg)
    gt = (pos_[:, None] > neg[None, :]).sum() + 0.5 * (pos_[:, None] == neg[None, :]).sum()
    return gt / (len(pos_) * len(neg)), len(pos_), len(neg)


print(f"gaps measured: {len(rows)}; cut by the curator: {sum(not r['kept'] for r in rows)}")
for label, sel in (("all gaps", lambda r: True),
                   ("homolog mapped", lambda r: r["partner"] != ""),
                   ("homolog mapped, other copy also a gap", lambda r: r["partner_gap"]),
                   ("homolog mapped, other copy a contig", lambda r: r["partner"] != "" and not r["partner_gap"]),
                   ("no homolog mapped", lambda r: r["partner"] == "")):
    sub = [r for r in rows if sel(r)]
    line = f"{label:42} n={len(sub):4}"
    for m in ("own", "other", "min", "pooled", "pooled4"):
        a_, nk, nc = auc([r[m] for r in sub], [r["kept"] for r in sub])
        line += f"  {m} {a_:.3f}" if not np.isnan(a_) else f"  {m}   -  "
    line += f"  (kept {sum(r['kept'] for r in sub)}, cut {sum(not r['kept'] for r in sub)})"
    print(line)
if out_path:
    cols = ["scaffold", "mb", "kept", "w", "own", "other", "min", "pooled", "pooled4", "partner", "partner_gap"]
    with open(out_path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(round(r[c], 3)) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
