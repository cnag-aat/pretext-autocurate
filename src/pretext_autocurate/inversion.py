"""Inversions inside a scaffold, from the bow-tie they leave on the contact map.

A segment [i, j) assembled in the wrong orientation keeps its own diagonal intact, so contact
across its two ends need not drop much on a saturated map (the misjoin rules miss it). What
gives it away is which ends touch. With L the window just left of i, R the window just right of
j, and B0 / B1 the segment's first and last windows:

  as assembled   L.B0 and B1.R are adjacent sequence;
  reversed       L.B1 and B0.R are adjacent sequence (the two corners of the bow-tie).

Each corner is a w x w block divided by what two abutting windows give (the adjacent level,
from the map's decay curve). A real inversion needs both: the reversed corners at the adjacent
level (gain) AND the assembled ones down to what sequence the segment's length apart shows
(drop). Reversed corners high while the assembled ones stay high is what a heterozygous
inversion (both orientations in the reads) or a repeat gives, and is not called. Measured on
icLarUrsu.pri.yahs (curated): scaffold_13 39.7-46.1 Mb, gain 1.24 / 1.24, drop 0.12 / 0.29;
scaffold_4 68.8-78.2 Mb (left as assembled by the curator), gain 1.32, drop 0.95.

A segment that runs to a scaffold end has one junction, so half a bow-tie (curator, 2026-11-19):
it is judged on that junction alone, with stricter thresholds.
"""
from dataclasses import dataclass

import numpy as np

from .signal import Signal


@dataclass
class Params:
    w: int = 20               # px: corner window (~1 Mb at 54 kb/px)
    min_px: int = 40          # px: shortest segment
    max_expect: float = 0.75  # contact expected at the segment's length must be below this
    gain_min: float = 0.8     # reversed corners: (corner - expected) / (1 - expected), each
    drop_max: float = 0.4     # assembled corners: the same scale, each
    end_gain_min: float = 1.0     # half a bow-tie (segment at a scaffold end)
    end_drop_max: float = 0.25
    max_corner: float = 1.4   # a reversed corner above this (x adjacent level): repeat, not a junction
    max_bins: int = 3000      # scaffolds longer than this many px are scanned on 2, 3, .. px bins


@dataclass
class Inversion:
    scaffold: str
    start_px: int             # first pixel of the segment (original map)
    end_px: int               # last pixel of the segment
    at_start: bool            # segment begins at the scaffold start (one junction)
    at_end: bool              # segment ends at the scaffold end (one junction)
    gain: float               # weakest reversed corner, 1 = adjacent, 0 = as expected apart
    drop: float               # strongest assembled corner, same scale
    expected: float           # contact expected at the segment's length (adjacent = 1)

    @property
    def score(self):
        return self.gain - self.drop


def _scan(M, w, min_n, g, adj, p):
    """Every segment [i, j) of one scaffold's contact block M (n x n, already binned); returns
    [(score, i, j, gain, drop, expected)] passing the thresholds."""
    n = M.shape[0]
    I = np.zeros((n + 1, n + 1))
    I[1:, 1:] = M.cumsum(0).cumsum(1)

    def S(r0, c0):          # w x w blocks at rows r0.., cols c0.. (arrays), / adjacent level
        return (I[r0 + w, c0 + w] - I[r0, c0 + w] - I[r0 + w, c0] + I[r0, c0]) / adj
    k = np.arange(w)
    # expected corner (scale: adjacent = 1) for windows whose starts are d apart (d >= w)
    D = np.arange(n + 1)
    expect = np.array([g[np.minimum(k[:, None] - k[None, :] + d, len(g) - 1)].sum() for d in range(w, n + 1)]) / adj \
        if n >= w else np.zeros(0)
    ex = lambda d: expect[np.minimum(d, n) - w]
    out = []
    for j in range(min_n + w, n + 1):
        i = np.arange(0, j - min_n + 1)
        jj = np.full_like(i, j)
        e = ex(j - i)                   # L=[i-w, i) and B1=[j-w, j) start j-i apart
        sc = 1.0 - e
        inner_l = i >= w
        inner_r = j + w <= n
        lb0 = np.where(inner_l, S(np.maximum(i - w, 0), i), np.nan)
        lb1 = np.where(inner_l, S(np.maximum(i - w, 0), jj - w), np.nan)
        if inner_r:
            b1r, b0r = S(jj - w, jj), S(i, jj)
        else:
            b1r = b0r = np.full(len(i), np.nan)
        gl, dl = (lb1 - e) / sc, (lb0 - e) / sc
        gr, dr = (b0r - e) / sc, (b1r - e) / sc
        both = inner_l & inner_r
        gain = np.where(both, np.fmin(gl, gr), np.where(inner_l, gl, gr))
        drop = np.where(both, np.fmax(dl, dr), np.where(inner_l, dl, dr))
        one = ~both
        ok = (e <= p.max_expect) & np.isfinite(gain) & np.isfinite(drop) & \
            np.where(one, (gain >= p.end_gain_min) & (drop <= p.end_drop_max), (gain >= p.gain_min) & (drop <= p.drop_max))
        # one junction only at a scaffold end itself: a segment stopping a few px short of the end
        # is the end folding back (bAegFun: both haplotypes at the same places), not an inversion
        ok &= (inner_l | (i == 0)) & (inner_r | (j == n))
        ok &= inner_l | inner_r          # a whole scaffold reversed is no inversion
        # a reversed corner far above adjacent contact is a repeat or a copy, not a junction
        ok &= ~(np.fmax(lb1, b0r) > p.max_corner)
        for x in np.nonzero(ok)[0]:
            out.append((float(gain[x] - drop[x]), int(i[x]), j, float(gain[x]), float(drop[x]), float(e[x])))
    return out


def detect_inversions(pt, params=None, signal=None):
    """Inversion calls in `pt` (a Pretext or SegmentedMap), best first per scaffold, segments of
    one scaffold not overlapping. Breakpoints are at the corner windows' edges (± a few px)."""
    p = params or Params()
    sg = signal or Signal(pt)
    # decay out to one texture (beyond it decay_long has no samples); further out its last
    # value stands in, which overstates what is expected there and so only makes calls harder
    g = np.array(sg.decay_long(max(4 * p.w + 2, int(pt.tex_res) - 1)), dtype=np.float64)
    g = g[:max(1, int(np.max(np.nonzero(g > 0)[0])) + 1)] if (g > 0).any() else g
    g = np.concatenate([g, np.full(int(np.max(pt.scaffold_npix)) + 4 * p.w, g[-1])])
    k = np.arange(p.w)
    calls = []
    for s, name in enumerate(pt.names):
        a, L = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        f = max(1, -(-L // p.max_bins))          # bin factor
        w, min_n = max(4, p.w // f), max(2, p.min_px // f)
        if L // f < min_n + 2 * w:
            continue
        M = sg.block(a, a + L - 1, a, a + L - 1).astype(np.float64)
        if f > 1:
            n = L // f
            M = M[:n * f, :n * f].reshape(n, f, n, f).sum(axis=(1, 3))
            gb = np.array([g[min(d * f, len(g) - 1)] for d in range(len(g) // f + 2 * w + 2)]) * f * f
        else:
            gb = g
        adj = float(gb[np.minimum(k[:w, None] + k[None, :w] + 1, len(gb) - 1)].sum())
        hits = sorted(_scan(M, w, min_n, gb, adj, p), reverse=True)
        taken = []
        for sc, i, j, gain, drop, e in hits:
            lo, hi = i * f, min(L, j * f) - 1
            if any(not (hi < t0 or lo > t1) for t0, t1 in taken):
                continue
            taken.append((lo, hi))
            calls.append(Inversion(name, a + lo, a + hi, i == 0, j * f >= L - f + 1 or j == M.shape[0],
                                   gain, drop, e))
    return calls
