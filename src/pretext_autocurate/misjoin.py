"""Misjoin detection from contact across a boundary ("insulation"), cuts preferred at gaps.

A misjoin is a sharp dip: contact between the sequence left and right of a boundary falls
to a small fraction of what an intact scaffold shows at the same distances, while each side
is normally connected within itself. That last condition is what separates a misjoin from
compartment (checkerboard) boundaries, which weaken but do not remove contact, and from
regions with no signal at all, where everything is low over a wide stretch.
"""
from dataclasses import dataclass, asdict

import numpy as np

from .signal import Signal


@dataclass
class Params:
    w_near: int = 10          # px: window touching the boundary
    w_far: int = 30           # px: window set back by `skip` (ignores the saturated diagonal)
    skip: int = 10
    cut_max: float = 0.2      # both windows below this fraction of expectation
    strong_max: float = 0.1   # needed to cut where there is no gap nearby
    side_min: float = 0.5     # each side must be internally connected at least this well
    side_reach: int = 30      # px searched on each side for that
    snap: int = 10            # px: a cut is moved to a gap within this distance
    min_piece: int = 5        # px: minima closer than this are merged


@dataclass
class Cut:
    scaffold: str
    after_px: int             # cut after this original pixel
    bp: int
    at_gap: bool
    near: float
    far: float
    apply: bool
    reason: str


def detect_misjoins(pt, params=None, signal=None):
    p = params or Params()
    sg = signal or Signal(pt)
    near = sg.insulation(p.w_near, 0)
    far = sg.insulation(p.w_far, p.skip)
    gap = pt.graphs.get("gap")
    gpos = np.nonzero(gap > 0)[0] if gap is not None else np.array([], dtype=int)
    n = pt.npix
    # close to scaffold ends only the near window exists
    score = np.where(np.isfinite(far), np.maximum(near, far), near)
    low = np.isfinite(score) & (score < p.cut_max)
    cuts, q = [], 0
    idx = np.nonzero(low)[0]
    groups = np.split(idx, np.nonzero(np.diff(idx) > p.min_piece)[0] + 1) if len(idx) else []
    for grp in groups:
        q = int(grp[np.argmin(score[grp])])
        s = int(pt.orig_scaffold[q])
        a0, a1 = int(pt.scaffold_start[s]), int(pt.scaffold_start[s] + pt.scaffold_npix[s] - 1)
        lo, hi = max(a0, int(grp[0]) - p.side_reach), min(a1, int(grp[-1]) + p.side_reach)
        left = near[lo:max(lo, int(grp[0]) - p.w_near + 1)]
        right = near[min(hi, int(grp[-1]) + p.w_near):hi + 1]
        sides_ok = (np.isfinite(left).any() and np.nanmax(left) >= p.side_min and
                    np.isfinite(right).any() and np.nanmax(right) >= p.side_min)
        cand = gpos[(gpos >= max(a0, int(grp[0]) - p.snap)) & (gpos <= min(a1 - 1, int(grp[-1]) + p.snap))]
        at_gap = len(cand) > 0
        if at_gap:
            q = int(cand[np.argmin(np.abs(cand - q))])
        nv, fv = float(near[grp].min()), float(np.nanmin(far[grp])) if np.isfinite(far[grp]).any() else float("nan")
        if not sides_ok:
            ok, why = False, "low contact over a wide region, not a sharp break: flagged only"
        elif at_gap:
            ok, why = True, "contact across the gap is at background level"
        elif np.nanmin([nv, fv]) < p.strong_max:
            ok, why = True, "no gap nearby, but contact across this point is below %.0f%% of expectation" % (100 * p.strong_max)
        else:
            ok, why = True, "no gap nearby; contact across this point is below %.0f%% of expectation on both scales, with normal contact on each side" % (100 * p.cut_max)
        cuts.append(Cut(pt.names[s], q, int(round((pt.rel_coord[q] + 1) * pt.bp_per_pix)), at_gap, nv, fv, ok, why))
    # Small foreign pieces: two sharp dips of the near window close together, while the far
    # window (which reaches over the piece) is normal. Both ends are cut, at gaps only.
    used = np.zeros(n, dtype=bool)
    for grp in groups:
        used[max(0, grp[0] - p.min_piece):grp[-1] + p.min_piece + 1] = True
    idx = np.nonzero(np.isfinite(near) & (near < p.cut_max) & ~used)[0]
    dips = []
    for grp in (np.split(idx, np.nonzero(np.diff(idx) > p.min_piece)[0] + 1) if len(idx) else []):
        q = int(grp[np.argmin(near[grp])])
        cand = gpos[np.abs(gpos - q) <= p.snap]
        cand = cand[pt.orig_scaffold[cand] == pt.orig_scaffold[q]] if len(cand) else cand
        dips.append((int(cand[np.argmin(np.abs(cand - q))]) if len(cand) else q, len(cand) > 0, float(near[grp].min()), q))
    reach = 2 * (p.w_far + p.skip)
    for (q1, g1, v1, r1), (q2, g2, v2, r2) in zip(dips, dips[1:]):
        if pt.orig_scaffold[q1] != pt.orig_scaffold[q2] or not p.min_piece <= q2 - q1 <= reach:
            continue
        for q, g, v in ((q1, g1, v1), (q2, g2, v2)):
            if any(c.after_px == q for c in cuts):
                continue
            s = int(pt.orig_scaffold[q])
            ok = True
            why = (f"piece of {q2 - q1} px between two breaks: no contact with either side" +
                   ("" if g else "; this end is not at a gap"))
            cuts.append(Cut(pt.names[s], q, int(round((pt.rel_coord[q] + 1) * pt.bp_per_pix)), g, v,
                            float(far[q]) if np.isfinite(far[q]) else float("nan"), ok, why))
    cuts.sort(key=lambda c: c.after_px)
    return cuts


def weak_junctions(pt, signal, taken=(), near_max=0.5, far_max=0.6, big_px=400, merge=5):
    """Gaps where contact across is clearly reduced but not at background.

    These are not cut outright. They are candidate breakpoints for rearrangements: the
    caller cuts them tentatively and keeps a cut only if the pieces fit better elsewhere.
    Returns [(pixel, near, far)].
    """
    gap = pt.graphs.get("gap")
    if gap is None:
        return []
    near, far = signal.insulation(10, 0), signal.insulation(30, 10)
    taken = np.array(sorted(taken), dtype=int)
    out = []
    for q in np.nonzero(gap > 0)[0]:
        q = int(q)
        if pt.scaffold_npix[pt.orig_scaffold[q]] <= big_px or (len(taken) and np.abs(taken - q).min() <= merge):
            continue
        w, f = near[max(0, q - 2):q + 3], far[max(0, q - 2):q + 3]
        if not np.isfinite(w).any():
            continue
        v = float(np.nanmin(w))
        fv = float(np.nanmin(f)) if np.isfinite(f).any() else float("nan")
        if v < near_max and not fv >= far_max:
            if out and q - out[-1][0] <= merge:
                if v < out[-1][1]:
                    out[-1] = (q, v, fv)
            else:
                out.append((q, v, fv))
    return out


def write_cuts(rows, path):
    cols = list(asdict(rows[0])) if rows else ["scaffold"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in asdict(r).values()) + "\n")
