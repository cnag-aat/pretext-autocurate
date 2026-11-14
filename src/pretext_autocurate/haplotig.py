"""Haplotig (haplotypic duplication) detection for whole small scaffolds, single mode.

Evidence, all relative to the map itself:
  line     an off-diagonal line of contact parallel (or anti-parallel) to the diagonal,
           brighter than a piece that is merely adjacent to that locus would give
  cov      long-read coverage of the piece below the genome norm
  pcov     coverage of the partner locus also below the norm (reads are shared)

Limits: near-diagonal values are saturated in typical maps, so for pieces of a few pixels
the line cannot be told from adjacency and coverage carries the call. Hemizygous sex
chromosome sequence also has low coverage; it is never called on coverage alone.
"""
from dataclasses import dataclass, asdict

import numpy as np

from .signal import Signal


@dataclass
class Params:
    """Haplotig thresholds; sizes in full-resolution pixels, coverage as ratio to the norm."""
    min_px: int = 3           # smaller pieces are reported as too small to assess
    max_px: int = 400          # larger scaffolds are not haplotig candidates here
    cov_max: float = 0.8       # coverage ratio below this counts as reduced
    cov_min: float = 0.3       # ... but below this it is unmappable/repeat, not a haplotig signal
    low_line: float = 0.75     # borderline tier: line/own at least this ...
    low_cov_max: float = 1.2   # ... with coverage below this
    line_excess: float = 0.05  # line/own minus the adjacency expectation must exceed this
    n_peaks: int = 4           # line candidates examined per strand and window length
    bg_window: int = 8         # px each side for the local background (median) of a line
    bg_exclude: int = 2        # px each side of the line left out of that background


@dataclass
class Call:
    """Haplotig call for one small scaffold; one row of the haplotigs TSV."""
    scaffold: str
    start_px: int
    length_px: int
    length_bp: int
    call: str                  # haplotig_high | haplotig_medium | haplotig_low | review | none | too_small
    cov: float                 # median coverage ratio of the piece
    line_over_own: float       # line brightness / the piece's own diagonal
    line_excess: float         # line_over_own minus the adjacency expectation
    partner: str
    partner_bp: int            # 1-based start of the partner locus (pixel-derived)
    partner_strand: str        # '-': anti-parallel line (inverted copy)
    partner_cov: float
    partner_near_gap: bool
    reason: str


def _contrast(d, w, ex):
    """`d` minus its local background: median over +-w px, leaving out the central +-ex.

    A line must stand out from its surroundings, not just be bright (e.g. near the diagonal).
    """
    pad = np.pad(d, w, mode="edge")
    win = np.lib.stride_tricks.sliding_window_view(pad, 2 * w + 1)
    keep = [k for k in range(2 * w + 1) if abs(k - w) > ex]
    return d - np.median(win[:, keep], axis=1)


def _peaks(d, k, excl):
    """Indices of up to k positive maxima of `d`, each masking +-excl around it."""
    d = d.copy()
    out = []
    for _ in range(k):
        c = int(np.argmax(d))
        if d[c] <= 0:
            break
        out.append(c)
        d[max(0, c - excl):c + excl + 1] = -1
    return out


def detect_haplotigs(pt, params=None, signal=None):
    """Call every scaffold of up to `max_px` px as haplotig or not; returns a list of Call.

    For each piece the best off-diagonal line (either strand, whole or half length) is
    compared with the piece's own diagonal and with the adjacency expectation, then
    combined with the coverage of the piece and of the partner locus into a tier.
    """
    p = params or Params()
    sg = signal or Signal(pt)
    cov = sg.coverage_ratio()
    gap = pt.graphs.get("gap")
    lens = pt.scaffold_len_bp()
    n = pt.npix
    calls = []
    for s in range(pt.n_scaffolds):
        L = int(pt.scaffold_npix[s])
        if L == 0 or L > p.max_px:
            continue
        a = int(pt.scaffold_start[s])
        cv = float(np.median(cov[a:a + L])) if cov is not None else float("nan")
        if L < p.min_px:
            calls.append(Call(pt.names[s], a, L, int(lens[s]), "too_small", cv, np.nan, np.nan,
                              "", 0, "", np.nan, False, f"under {p.min_px} px"))
            continue
        S = sg.strip(a, a + L - 1).astype(np.float32)
        m = n - L + 1                   # possible partner start columns
        own = float(np.mean(S[np.arange(L), a + np.arange(L)]))   # distance-0 brightness
        w = max(3, (L + 1) // 2)
        best = None
        for strand in "+-":
            # V[i, c] = contact of piece pixel i with column c+i ('+') or of pixel L-1-i ('-'),
            # so a line starting at column c is column c of V
            V = np.stack([S[i, i:i + m] if strand == "+" else S[L - 1 - i, i:i + m] for i in range(L)])
            cs = np.cumsum(np.vstack([np.zeros((1, m), np.float32), V]), axis=0)
            # whole-length line, and the best half-length run (partial duplications)
            for d, ll in ((V.mean(0), L), (((cs[w:] - cs[:-w]) / w).max(0), w)):
                dd = _contrast(d, p.bg_window, p.bg_exclude)
                dd[max(0, a - 2 * L):a + 2 * L] = -1          # own diagonal
                for c in _peaks(dd, p.n_peaks, L):
                    # excess over what mere adjacency to the partner locus would give
                    ex = float(d[c] / max(own, 1) - sg.adjacency_level(ll))
                    pc = float(np.median(cov[c:c + L])) if cov is not None else float("nan")
                    # prefer candidates with line excess, then low partner coverage, then excess
                    key = (ex > p.line_excess, pc < p.cov_max, ex)
                    if best is None or key > best[0]:
                        best = (key, c, strand, float(d[c] / max(own, 1)), ex, pc)
        if best is None:
            calls.append(Call(pt.names[s], a, L, int(lens[s]), "none", cv, 0.0, np.nan,
                              "", 0, "", np.nan, False, "no off-diagonal line"))
            continue
        _, c, strand, ratio, ex, pc = best
        has_line, low, plow = ex > p.line_excess, p.cov_min <= cv < p.cov_max, p.cov_min <= pc < p.cov_max
        if cov is not None and (cv < p.cov_min or (has_line and pc < p.cov_min)):
            call, why = "review", "coverage near zero on the piece or its partner: unmappable repeat, not called"
        elif has_line and low and plow:
            call, why = "haplotig_high", "line, low coverage, partner low coverage"
        elif has_line and low:
            call, why = "haplotig_medium", "line and low coverage"
        elif low and cov is not None:
            call, why = "haplotig_low", "borderline: low coverage, no clear line (small piece, diverged copy, or hemizygous)"
        elif ratio >= p.low_line and cv < p.low_cov_max:
            call, why = "haplotig_low", "borderline: line as bright as the diagonal, coverage not reduced"
        elif has_line:
            call, why = "review", "line but high coverage (repeat, or collapsed partner)"
        else:
            call, why = "none", ""
        ps = int(pt.orig_scaffold[c])
        # partner_near_gap: a gap within 3 px of the partner locus (a likely cut point there)
        calls.append(Call(pt.names[s], a, L, int(lens[s]), call, cv, ratio, ex, pt.names[ps],
                          int(round(pt.rel_coord[c] * pt.bp_per_pix)) + 1, strand, pc,
                          bool(gap is not None and (gap[max(0, c - 3):c + L + 3] > 0).any()), why))
    return calls


@dataclass
class ContigHaplotig:
    """A gap-bounded contig inside a large scaffold called as a haplotig of nearby sequence."""
    scaffold: str
    start_px: int        # first pixel of the contig (original map)
    end_px: int          # last pixel
    cov: float
    partner_cov: float
    offset_px: int       # the duplicate lies this many pixels along the same scaffold
    line: float          # line value / value at distance 0
    line_over_expected: float


def detect_contig_haplotigs(pt, signal=None, min_px=8, max_px=150, cov_max=0.8, partner_cov_max=0.85,
                            line_min=0.75, line_ratio_min=2.0, big_px=400):
    """Haplotigs that are whole contigs (gap to gap) inside a chromosome-scale scaffold.

    A contig is called when it has reduced coverage, and a line of contact parallel to the
    diagonal joins it to nearby sequence of the same scaffold that also has reduced
    coverage; the line must be nearly as bright as the diagonal itself and at least
    `line_ratio_min` times what the distance alone explains. Deliberately strict: on a
    saturated map most contig-level duplications are not separable from compartment noise.
    """
    sg = signal or Signal(pt)
    cov, gap = sg.coverage_ratio(), pt.graphs.get("gap")
    if cov is None or gap is None:
        return []
    g, n = sg.decay, pt.npix
    out = []
    for s in np.nonzero(pt.scaffold_npix > big_px)[0]:
        a0 = int(pt.scaffold_start[s])
        a1 = a0 + int(pt.scaffold_npix[s]) - 1
        bounds = [a0 - 1] + (a0 + np.nonzero(gap[a0:a1] > 0)[0]).tolist() + [a1]
        for lo, hi in zip(bounds, bounds[1:]):
            a, L = lo + 1, hi - lo
            if not min_px <= L <= max_px:
                continue
            cv = float(np.median(cov[a:a + L]))
            if cv >= cov_max:
                continue
            reach = 2 * L + 50          # px: partners searched only near the contig
            c0, c1 = max(a0, a - reach), min(a1, a + L - 1 + reach)
            S = sg.block(a, a + L - 1, c0, c1).astype(np.float32)
            best = None
            for d in range(-reach, reach + 1):      # parallel lines only, offset d px
                # skip offsets so small the "line" is the saturated diagonal itself
                if abs(d) < max(4, L // 2) or a + d < c0 or a + L - 1 + d > c1:
                    continue
                val = float(S[np.arange(L), np.arange(L) + a + d - c0].mean())
                exp = float(g[min(abs(d), len(g) - 1)])
                if best is None or val - exp > best[0]:
                    best = (val - exp, d, val, exp)
            if best is None:
                continue
            _, d, val, exp = best
            pc = float(np.median(cov[a + d:a + d + L]))
            if val / g[0] >= line_min and val / max(exp, 1) >= line_ratio_min and pc < partner_cov_max:
                out.append(ContigHaplotig(pt.names[s], a, a + L - 1, cv, pc, d, float(val / g[0]), float(val / max(exp, 1))))
    return out


def write_calls(calls, path):
    """Write Calls (or ContigHaplotigs) as TSV, floats to 3 decimals."""
    cols = list(asdict(calls[0])) if calls else []
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for c in calls:
            f.write("\t".join(f"{v:.3f}" if isinstance(v, float) else str(v) for v in asdict(c).values()) + "\n")
