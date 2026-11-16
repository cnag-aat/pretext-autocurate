"""Score a savestate against a reference curation of the same map (`compare`).

Everything is measured on map pixels (one pixel = the map's bp per pixel), on the layout each
savestate produces in PretextView:

  junctions   two pieces (>= `min_px`) laid next to each other inside one painted
              chromosome, both in the chromosome proper (not Unloc / Haplotig /
              FalseDuplicate / Contaminant). A junction of one layout is supported by the other
              when the same two pixels are in one chromosome there, close together
              (`near_px`). precision = tool junctions supported by the reference; recall =
              reference junctions the tool reproduces. As in tools/compare_layouts.py.
  breaks      places where two neighbouring pixels of one input scaffold are no longer next to
              each other the same way round (a cut, or the end of a moved or inverted piece).
              A tool break matches a reference break in the same scaffold within `tol_px`.
  grouping    for each reference chromosome, the share of its sequence (chromosome proper)
              that sits in the single tool chromosome holding most of it; `grouped` is the
              share of reference chromosomes with at least 90% in one tool chromosome.
  tags        per tag name, bp tagged in the tool, in the reference and in both.
"""
import collections

import numpy as np

from .mapstate import MapState

ASIDE = {"Unloc", "Haplotig", "FalseDuplicate", "Contaminant"}
TAGS = ["HAP1", "HAP2", "Unloc", "Haplotig", "FalseDuplicate", "Contaminant", "X", "Y", "Z", "W"]


class Layout:
    """Per-pixel view of a savestate on its map."""

    def __init__(self, pt, st):
        n = pt.npix
        self.ms = ms = MapState(pt).replay(st.edits)
        self.pos = np.empty(n, dtype=np.int64)
        self.pos[ms.pix] = np.arange(n)
        paint, flags = dict(st.scaffs), dict(st.meta_flags)
        self.chrom = np.zeros(n, dtype=np.int64)
        self.aside = np.zeros(n, dtype=bool)
        self.tags = collections.defaultdict(lambda: np.zeros(n, dtype=bool))
        self.pieces = []
        for i, c in enumerate(ms.contigs):
            px = ms.pix[c.start:c.start + c.length]
            names = set(st.tag_names(flags.get(i, 0)))
            self.chrom[px] = paint.get(i, 0)
            self.aside[px] = bool(names & ASIDE)
            for t in names:
                self.tags[t][px] = True
            self.pieces.append(c)
        self.n_edits = len(st.edits)


def _junctions(L, min_px):
    """Pairs of pixels (end of one piece, start of the next) inside one painted chromosome."""
    ms, out = L.ms, []
    big = [c for c in ms.contigs if c.length >= min_px and L.chrom[ms.pix[c.start]] > 0
           and not L.aside[ms.pix[c.start]]]
    for a, b in zip(big, big[1:]):
        k = min(6, a.length // 2, b.length // 2)      # a few pixels inside: ends differ by a pixel or two
        pa, pb = int(ms.pix[a.start + a.length - 1 - k]), int(ms.pix[b.start + k])
        if L.chrom[pa] == L.chrom[pb]:
            out.append((pa, pb))
    return out


def _supported(j, L, near_px):
    pa, pb = j
    return L.chrom[pa] == L.chrom[pb] and L.chrom[pa] > 0 and abs(int(L.pos[pa]) - int(L.pos[pb])) <= near_px


def _breaks(pt, L):
    """Pixels p (of an input scaffold) after which the layout no longer continues with p+1."""
    same = pt.orig_scaffold[1:] == pt.orig_scaffold[:-1]
    step = np.abs(L.pos[1:] - L.pos[:-1])
    return np.nonzero(same & (step != 1))[0]


def _match(a, b, scaf, tol):
    """How many of breaks `a` have a break of `b` in the same scaffold within tol pixels."""
    if not len(a) or not len(b):
        return 0
    b = np.sort(b)
    k = 0
    for p in a:
        i = np.searchsorted(b, p)
        for q in b[max(0, i - 1):i + 1]:
            if abs(int(q) - int(p)) <= tol and scaf[q] == scaf[p]:
                k += 1
                break
    return k


def compare(pt, tool_st, ref_st, min_px=15, near_px=80, tol_px=3):
    """Metrics of `tool_st` against `ref_st` (both savestates of map `pt`): a flat dict."""
    T, R = Layout(pt, tool_st), Layout(pt, ref_st)
    bpp = pt.bp_per_pix
    out = {}
    jt, jr = _junctions(T, min_px), _junctions(R, min_px)
    out["junctions_tool"] = len(jt)
    out["junctions_ref"] = len(jr)
    out["junctions_supported"] = sum(_supported(j, R, near_px) for j in jt)
    out["junctions_reproduced"] = sum(_supported(j, T, near_px) for j in jr)
    out["junction_precision"] = out["junctions_supported"] / len(jt) if jt else float("nan")
    out["junction_recall"] = out["junctions_reproduced"] / len(jr) if jr else float("nan")
    bt, br = _breaks(pt, T), _breaks(pt, R)
    out["breaks_tool"], out["breaks_ref"] = len(bt), len(br)
    out["breaks_matched_tool"] = _match(bt, br, pt.orig_scaffold, tol_px)
    out["breaks_matched_ref"] = _match(br, bt, pt.orig_scaffold, tol_px)
    out["break_precision"] = out["breaks_matched_tool"] / len(bt) if len(bt) else float("nan")
    out["break_recall"] = out["breaks_matched_ref"] / len(br) if len(br) else float("nan")
    # grouping: reference chromosomes (proper) -> tool chromosomes
    proper = (R.chrom > 0) & ~R.aside
    shares, total, inmain, n90, nchr = [], 0, 0, 0, 0
    for c in np.unique(R.chrom[proper]):
        m = proper & (R.chrom == c)
        tot = int(m.sum())
        best = collections.Counter(T.chrom[m].tolist()).most_common()
        top = next((v for k, v in best if k > 0), 0)        # unpainted in the tool never counts
        nchr += 1
        n90 += top >= 0.9 * tot
        total += tot
        inmain += top
    out["ref_chromosomes"] = nchr
    out["grouped_90"] = n90 / nchr if nchr else float("nan")
    out["grouped_bp"] = inmain / total if total else float("nan")
    for t in TAGS:
        a, b = T.tags.get(t), R.tags.get(t)
        if a is None and b is None:
            continue
        a = a if a is not None else np.zeros(pt.npix, dtype=bool)
        b = b if b is not None else np.zeros(pt.npix, dtype=bool)
        out[f"tag_{t}_tool_mb"] = round(float(a.sum()) * bpp / 1e6, 2)
        out[f"tag_{t}_ref_mb"] = round(float(b.sum()) * bpp / 1e6, 2)
        out[f"tag_{t}_both_mb"] = round(float((a & b).sum()) * bpp / 1e6, 2)
    out["edits_tool"], out["edits_ref"] = T.n_edits, R.n_edits
    return out
