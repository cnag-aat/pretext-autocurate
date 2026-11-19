"""Placement of small scaffolds: which chromosome-scale scaffold, and where if it is clear.

For each candidate piece the mean contact of its rows with every column of larger
("anchor") scaffolds is scanned:
  - the strongest window gives the anchor and a position; it must reach a fraction of what
    a truly adjacent piece would show (map's own distance decay) and beat every other anchor
  - the piece is *localised* when that position sits at a junction (gap or anchor end), its
    orientation is clear from which of its ends touches which flank, it is not repeat-like
    (coverage) and not too small; otherwise it is an *unloc* of that anchor
  - with no cis-level contact anywhere it is left alone (own chromosome, or contaminant)
"""
from dataclasses import dataclass, asdict

import numpy as np

from .signal import Signal


@dataclass
class Params:
    """Placement thresholds; sizes in full-resolution pixels."""
    max_frac: float = 0.2       # candidates are scaffolds shorter than this x the longest
    max_px: int = 150            # ... and no longer than this: the adjacency expectation is
                                 # only valid well inside the measured decay range (200 px)
    anchor_ratio: float = 5.0    # anchors are at least this many times longer than the piece
    min_anchor_px: int = 100
    w_chrom: int = 20            # half-window (px) for choosing the anchor
    w_loc: int = 5               # half-window (px) for the position
    cis_min: float = 0.5         # window contact / adjacent expectation needed to belong
    ratio_min: float = 1.3       # best anchor / next best anchor
    loc_min_px: int = 4          # smaller pieces are never localised
    loc_max_cov: float = 2.0     # above this the piece is treated as repeat: unloc only
    junction_tol: int = 5        # px from the peak to a gap or anchor end
    orient_min: float = 1.15     # contact ratio between the two orientations
    flank: int = 15              # px of anchor on each side of a junction used for orientation
    junction_max: float = 0.6    # a junction has room for an insertion if contact across it is below this ...
    room_w: int = 10             # ... or below what a gap holding the piece would show (`junction_contact`
    room_slack: float = 0.1      # is insulation(room_w): a w x w square across the boundary) plus this
    contact_share: float = 0.6   # junctions where the piece's contact is below this x its best are not candidates



@dataclass
class Placement:
    """Placement call for one small scaffold; one row of the placements TSV."""
    scaffold: str
    start_px: int
    length_px: int
    length_bp: int
    call: str            # localised | unloc | unplaced | too_small
    anchor: str
    anchor_bp: int       # 1-based position of the contact peak in the anchor (pixel-derived)
    after_px: int        # localised: insert after this original pixel (-1 otherwise)
    strand: str
    cis: float           # peak contact / adjacent expectation
    ratio: float         # best anchor / next best anchor
    orient: float        # stronger / weaker orientation score
    cov: float           # median coverage ratio of the piece
    reason: str


def _box(p, w):
    """Moving mean of `p` over [i-w, i+w], shrinking the window at the ends."""
    cs = np.concatenate(([0.0], np.cumsum(p, dtype=np.float64)))
    n = len(p)
    lo = np.clip(np.arange(n) - w, 0, n)
    hi = np.clip(np.arange(n) + w + 1, 0, n)
    return (cs[hi] - cs[lo]) / (hi - lo)


def place_scaffolds(pt, params=None, signal=None, only=None, junction_contact=None):
    """`junction_contact`: optional per-pixel contact across each boundary (insulation of
    the original map). With it, a piece is localised only at a junction with reduced
    contact (or a piece end): a junction with normal contact has no room for an insertion,
    and a piece that seems to belong there is more likely a copy of the sequence beside it.
    Among such junctions near the contact peak the one with the clearest orientation wins.

    Candidates are scaffolds up to `max_px` px and `max_frac` of the longest; `only`
    restricts them to a set of scaffold indices. Returns a list of Placement."""
    p = params or Params()
    sg = signal or Signal(pt)
    cov = sg.coverage_ratio()
    gap = pt.graphs.get("gap")
    g = sg.decay
    lens = pt.scaffold_len_bp()
    npx = pt.scaffold_npix
    scaf = pt.orig_scaffold
    longest = int(npx.max())
    out = []
    for s in range(pt.n_scaffolds):
        L = int(npx[s])
        if L == 0 or L > min(p.max_frac * longest, p.max_px) or (only is not None and s not in only):
            continue
        a = int(pt.scaffold_start[s])
        cv = float(np.median(cov[a:a + L])) if cov is not None else float("nan")
        base = dict(scaffold=pt.names[s], start_px=a, length_px=L, length_bp=int(lens[s]), cov=cv)
        anchor_px = npx[scaf] >= max(p.anchor_ratio * L, p.min_anchor_px)
        if not anchor_px.any():
            continue
        S = sg.strip(a, a + L - 1).astype(np.float32)     # piece rows x all columns
        prof = np.where(anchor_px, S.mean(0), 0.0)         # mean contact per column
        b = _box(prof, p.w_chrom)
        b[~anchor_px] = 0               # the box can spill past anchor ends; mask again
        c0 = int(np.argmax(b))
        A = int(scaf[c0])
        a0, a1 = int(pt.scaffold_start[A]), int(pt.scaffold_start[A] + npx[A] - 1)
        other = b.copy()
        other[a0:a1 + 1] = 0            # best window on any *other* anchor
        ratio = float(b[c0] / max(other.max(), 1e-3))
        # position: sharpest small window near the broad peak
        lo, hi = max(a0, c0 - 2 * p.w_chrom), min(a1, c0 + 2 * p.w_chrom)
        b5 = _box(prof[a0:a1 + 1], p.w_loc)
        c = a0 + lo - a0 + int(np.argmax(b5[lo - a0:hi - a0 + 1]))
        # expectation if the piece abutted the window: row i vs window column k is i+k px away
        d = np.arange(L)[:, None] + np.arange(1, 2 * p.w_loc + 2)[None, :]
        expect = float(g[np.minimum(d, len(g) - 1)].mean())
        cis = float(b5[c - a0] / expect)
        info = dict(anchor=pt.names[A], anchor_bp=int(round(pt.rel_coord[c] * pt.bp_per_pix)) + 1,
                    cis=cis, ratio=ratio)
        if cis < p.cis_min or ratio < p.ratio_min:
            why = "no cis-level contact with any larger scaffold" if cis < p.cis_min else \
                "contact shared between scaffolds"
            out.append(Placement(**base, call="unplaced", after_px=-1, strand="", orient=np.nan, reason=why, **info))
            continue
        # junctions of the anchor: before its first pixel, after its last, after each gap pixel
        junc = [a0 - 1, a1]
        if gap is not None:
            junc += (a0 + np.nonzero(gap[a0:a1] > 0)[0]).tolist()
        junc = np.array(sorted(set(junc)))
        t = max(1, L // 3)              # rows taken as the piece's first / last third

        def orient_at(q):
            """(plus, minus) contact scores for inserting the piece after pixel q."""
            left = S[:, max(a0, q - p.flank + 1):q + 1] if q >= a0 else None
            right = S[:, q + 1:min(a1, q + p.flank) + 1] if q < a1 else None
            plus = minus = 0.0     # '+': first pixels next to the left flank, last next to the right
            if left is not None and left.size:
                plus += left[:t].mean()
                minus += left[L - t:].mean()
            if right is not None and right.size:
                plus += right[L - t:].mean()
                minus += right[:t].mean()
            return plus, minus

        no_room = False
        if junction_contact is None:
            # junction q is the boundary q|q+1, hence the +0.5 when measuring to pixel c
            q = int(junc[np.argmin(np.abs(junc + 0.5 - c))])
            dist = abs(q + 0.5 - c)
        else:
            reach = 2 * p.w_chrom + L
            # contact per column without the piece's own surroundings: its former neighbours
            # (a trial-cut piece) touch it at saturated near-diagonal levels; NaN there
            far = prof.astype(np.float64).copy()
            far[max(0, a - 3):a + L + 3] = np.nan
            # A gap that holds a piece of L px pushes its two sides L px apart: contact across
            # it drops to `room` of normal (9 px: 0.74; 20 px: 0.58 on xgTylRafi). A small piece
            # leaves only a small dip (scaffold_119 at scaffold_13 11.04 Mb: 0.74, neighbouring
            # gaps 0.88-1.26), so the fixed junction_max alone missed most curated insertions.
            k = np.arange(p.room_w)
            gl = sg.decay_long(2 * p.room_w + L + 2)
            room = float(gl[k[:, None] + k[None, :] + 1 + L].sum() / gl[k[:, None] + k[None, :] + 1].sum())
            room_max = max(p.junction_max, room + p.room_slack)
            cands = []
            for qq in junc[(junc >= c0 - reach) & (junc <= c0 + reach)].tolist():
                w = junction_contact[max(0, qq - 2):qq + 3]
                if qq in (a0 - 1, a1) or (np.isfinite(w).any() and np.nanmin(w) < room_max):
                    pl_, mi_ = orient_at(qq)
                    # w_loc px either side: pixels right at a sequence gap read low (scaffold_2
                    # 40.5-41.5 Mb at scaffold_5 49.30 Mb: 54 at the gap, 100 five px on)
                    v = far[max(a0, qq - p.w_loc + 1):min(a1, qq + p.w_loc) + 1]
                    here = float(np.nanmean(v)) if np.isfinite(v).any() else np.nan
                    cands.append((max(pl_, mi_) / max(min(pl_, mi_), 1e-3), qq, here))   # (orient, q, contact)
            # junctions where the piece's own contact is well below its best are out (choosing by
            # orientation alone took xgTylRafi scaffold_137 to scaffold_9 41.25 Mb, contact 93,
            # past the curated 40.33 Mb gap where its contact peaks at 192); among the rest the
            # clearest orientation wins (contact alone would always send a trial-cut piece back
            # to the neighbour it was cut from: adjacent pixels are saturated)
            here_all = [c[2] for c in cands if np.isfinite(c[2])]
            if here_all:
                top = max(here_all)
                cands = [c for c in cands if not np.isfinite(c[2]) or c[2] >= p.contact_share * top]
            if cands:
                q, dist = max(cands)[1], 0.0
            else:
                q = int(junc[np.argmin(np.abs(junc + 0.5 - c))])
                dist, no_room = abs(q + 0.5 - c), True
        plus, minus = orient_at(q)
        orient = float(max(plus, minus) / max(min(plus, minus), 1e-3))
        strand = "+" if plus >= minus else "-"
        why = []
        if L < p.loc_min_px:
            why.append(f"under {p.loc_min_px} px")
        if cv > p.loc_max_cov:
            why.append("repeat-like coverage")
        if no_room:
            why.append("no junction with reduced contact near the peak")
        elif dist > p.junction_tol + L / 2:     # slack grows with piece length
            why.append("peak not at a gap or scaffold end")
        if orient < p.orient_min:
            why.append("orientation unclear")
        if why:
            out.append(Placement(**base, call="unloc", after_px=q if not no_room else -1, strand="", orient=orient,
                                 reason="; ".join(why), **info))
        else:
            out.append(Placement(**base, call="localised", after_px=q, strand=strand, orient=orient,
                                 reason="", **info))
    return out


def write_placements(rows, path):
    """Write Placements as TSV (dataclass fields as columns, floats to 3 decimals)."""
    cols = list(asdict(rows[0])) if rows else []
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in asdict(r).values()) + "\n")
