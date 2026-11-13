"""Joining and orienting chromosome-scale pieces from contact between their ends.

Two ends are joined when each is the other's best partner, the contact between their
terminal windows is at the level the map shows for truly adjacent sequence, and the next
best partner of either end is clearly weaker. An end carrying a telomere signal is a
chromosome end and is never joined (telomeres also cluster and share repeats, which gives
strong end-to-end contact between different chromosomes). The outermost pixels of each end
are skipped for the same reason.
"""
from dataclasses import dataclass
import itertools

import numpy as np


@dataclass
class Params:
    min_px: int = 150          # only pieces longer than this are joined here
    window: int = 100          # px of each end compared
    skip: int = 10             # outermost px ignored
    score_min: float = 0.6     # contact / expectation for adjacent ends
    margin: float = 0.8        # next best partner must be below this x the score
    tel_factor: float = 5.0    # telomeric end: track above this x its 99th percentile ...
    tel_floor: float = 50.0    # ... and above this value
    tel_px: int = 10


@dataclass
class Join:
    a: int                     # unit index in the (segmented) map
    a_end: str                 # 'h' (first pixel) or 't'
    b: int
    b_end: str
    score: float
    second: float
    accepted: bool
    reason: str


def _window(pt, s, end, W, K):
    a, L = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
    W = max(1, min(W, L // 2 - K))
    return (a + K, a + K + W - 1, False) if end == "h" else (a + L - K - W, a + L - K - 1, True)


def end_score(pt, sg, x, y, W, K, g):
    a0, a1, r1 = _window(pt, *x, W, K)
    b0, b1, r2 = _window(pt, *y, W, K)
    B = sg.block(a0, a1, b0, b1).astype(np.float64)
    i = np.arange(a1 - a0 + 1)
    j = np.arange(b1 - b0 + 1)
    i, j = (i[::-1] if r1 else i), (j[::-1] if r2 else j)
    expect = g[np.minimum(i[:, None] + j[None, :] + 1 + 2 * K, len(g) - 1)]
    return float(B.sum() / expect.sum())


def telomeric_ends(pt, units, p):
    tel = pt.graphs.get("telomere")
    if tel is None:
        return set()
    thr = max(p.tel_floor, p.tel_factor * float(np.percentile(tel, 99)))
    out = set()
    for s in units:
        a, L = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        if tel[a:a + p.tel_px].max() >= thr:
            out.add((s, "h"))
        if tel[a + L - p.tel_px:a + L].max() >= thr:
            out.add((s, "t"))
    return out


def find_joins(pt, sg, params=None):
    """Returns (joins, chains). A chain is a list of (unit, reversed) in chromosome order."""
    p = params or Params()
    units = [int(s) for s in np.nonzero(pt.scaffold_npix > p.min_px)[0]]
    ends = [(s, e) for s in units for e in "ht"]
    g = sg.decay_long(2 * (p.window + p.skip) + 2)
    sc = {}
    for x, y in itertools.combinations(ends, 2):
        if x[0] != y[0]:
            sc[(x, y)] = sc[(y, x)] = end_score(pt, sg, x, y, p.window, p.skip, g)
    best = {x: sorted(((sc[(x, y)], y) for y in ends if y[0] != x[0]), reverse=True)[:2] for x in ends}
    tel = telomeric_ends(pt, units, p)
    joins, link = [], {}
    for x in ends:
        if len(best[x]) < 2:
            continue
        (s1, y), (s2, _) = best[x]
        if not (best[y][0][1] == x and x < y) or s1 < 0.5 * p.score_min:
            continue
        second = max(s2, best[y][1][0])
        if x in tel or y in tel:
            ok, why = False, "an end carries a telomere signal: chromosome end, not joined"
        elif s1 < p.score_min:
            ok, why = False, "end-to-end contact below the level of adjacent sequence"
        elif second > p.margin * s1:
            ok, why = False, "another end fits almost as well: order unresolved"
        else:
            ok, why = True, "mutual best ends at the contact level of adjacent sequence"
        joins.append(Join(x[0], x[1], y[0], y[1], s1, second, ok, why))
        if ok:
            link[x], link[y] = y, x
    # walk the paths
    chains, seen = [], set()
    other = {"h": "t", "t": "h"}

    def walk(s, enter):
        chain = []
        while s not in seen:
            seen.add(s)
            chain.append((s, enter == "t"))            # entered at its tail: reversed
            nxt = link.get((s, other[enter]))
            if nxt is None:
                break
            s, enter = nxt
        return chain

    for s in units:                                    # start at free ends
        for e in "ht":
            if s not in seen and (s, e) not in link:
                chains.append(walk(s, e))
    for s in units:                                    # rings (all ends joined): open anywhere
        if s not in seen:
            chains.append(walk(s, "h"))
    out = []
    for ch in chains:                                  # keep the longest unit forward
        big = max(ch, key=lambda u: pt.scaffold_npix[u[0]])
        if big[1]:
            ch = [(s, not r) for s, r in reversed(ch)]
        out.append(ch)
    return joins, out
