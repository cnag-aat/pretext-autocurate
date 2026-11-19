"""Inversions inside a scaffold in allhaps (inversion.py, as in the single rules pass).

A call is turned in place (`_layout` reverses its pixels wherever its chromosome copy ends up)
when the corners pass the level's thresholds (common.THRESHOLDS inv_*) and the other haplotype
agrees: with an alignment, the segment must align to the homolog in the opposite direction to
its flanks; without one, the Hi-C call alone decides. It is never moved: making its ends split
points let the arrangement move the pieces too (bAegFun H2_scaffold_8 23.4-25.8 Mb went to the
chromosome start at end contact 0.04; junction precision 3/6 -> 5/12).

A call whose homolog shows the same pattern at the same place is folding of the chromosome, not
a misassembly (bAegFun H2_scaffold_2 start and H1_scaffold_2 end: one chromosome end in both
haplotypes), and is only reported.
"""
from ..inversion import detect_inversions
from ..pretextio import SegmentedMap


def _project(blocks, name, a_bp, b_bp):
    """[(partner scaffold, lo bp, hi bp)] where [a_bp, b_bp) of `name` aligns (blocks of >= 50 kb)."""
    out = []
    for b in blocks or ():
        if b.qe - b.qs < 50_000:
            continue
        if b.q == name:
            s0, s1, o0, o1, other = b.qs, b.qe, b.ts, b.te, b.t
        elif b.t == name:
            s0, s1, o0, o1, other = b.ts, b.te, b.qs, b.qe, b.q
        else:
            continue
        lo, hi = max(a_bp, s0), min(b_bp, s1)
        if hi <= lo:
            continue
        f = (o1 - o0) / max(s1 - s0, 1)
        x0, x1 = o0 + (lo - s0) * f, o0 + (hi - s0) * f
        if not b.same_strand:
            x0, x1 = o1 - (hi - s0) * f, o1 - (lo - s0) * f
        out.append((other, x0, x1))
    return out


def _strand_bp(blocks, name, a_bp, b_bp):
    """{partner: aligned bp of [a_bp, b_bp), + on the same strand, - on the other}."""
    out = {}
    for b in blocks or ():
        if b.q == name:
            s0, s1, other = b.qs, b.qe, b.t
        elif b.t == name:
            s0, s1, other = b.ts, b.te, b.q
        else:
            continue
        ov = min(b_bp, s1) - max(a_bp, s0)
        if ov > 0:
            out[other] = out.get(other, 0) + (ov if b.same_strand else -ov)
    return out


def _homolog_says(blocks, name, a_bp, b_bp, flank_bp=2_000_000, min_bp=100_000):
    """True: the segment aligns to its homolog in the other direction from its flanks (an
    inversion relative to the other haplotype); False: the same direction; None: not known."""
    seg = _strand_bp(blocks, name, a_bp, b_bp)
    if not seg:
        return None
    partner = max(seg, key=lambda k: abs(seg[k]))
    s = seg[partner]
    f = _strand_bp(blocks, name, max(0, a_bp - flank_bp), a_bp).get(partner, 0) + \
        _strand_bp(blocks, name, b_bp, b_bp + flank_bp).get(partner, 0)
    if abs(s) < min_bp or abs(f) < min_bp:
        return None
    return (s > 0) != (f > 0)


def find(raw, cut_qs, blocks, th):
    """Inversion calls on the pieces after `cut_qs`; `th` the run's Thresholds. Returns (pixel
    ranges to turn in place, decision rows)."""
    vt = SegmentedMap(raw, sorted(cut_qs))
    calls = detect_inversions(vt)
    bpp = raw.bp_per_pix
    span = lambda v: (raw.rel_coord[v.start_px] * bpp, (raw.rel_coord[v.end_px] + 1) * bpp)
    name_of = lambda v: raw.names[int(raw.orig_scaffold[v.start_px])]
    by_name = {}
    for v in calls:
        by_name.setdefault(name_of(v), []).append(v)
    flips, rows = [], []
    for v in calls:
        nm = name_of(v)
        a_bp, b_bp = span(v)
        end = v.at_start or v.at_end
        ev = (f"{a_bp / 1e6:.2f}-{b_bp / 1e6:.2f} Mb looks assembled backwards: contact across the reversed ends "
              f"{v.gain:.2f}, across the assembled ones {v.drop:.2f} (1 = adjacent sequence, 0 = what sequence that far "
              f"apart shows)" + ("; at a scaffold end: half a bow-tie" if end else ""))
        mirrored = []
        for other, x0, x1 in _project(blocks, nm, a_bp, b_bp):
            for w in by_name.get(other, []):
                w0, w1 = span(w)
                if min(w1, x1) - max(w0, x0) > 0.3 * min(w1 - w0, x1 - x0):
                    mirrored.append(f"{other} {w0 / 1e6:.2f}-{w1 / 1e6:.2f} Mb")
        hom = _homolog_says(blocks, nm, a_bp, b_bp)
        strong = v.gain >= th.inv_end_gain if end else (v.gain >= th.inv_gain and v.drop <= th.inv_drop)
        if mirrored:
            act, why = "flag_inversion", ("the homolog shows the same pattern at the same place (" + ", ".join(mirrored)
                                          + "): folding of the chromosome, not a misassembly; left as it is")
        elif hom is False:
            act, why = "flag_inversion", ("the segment aligns to the other haplotype in the same direction as its "
                                          "flanks: left as assembled")
        elif not strong:
            act, why = "flag_inversion", "below this level's thresholds: proposed, not turned"
        elif end and hom is not True:
            # half a bow-tie alone is not enough here: at level 1 it turned bAegFun H2_scaffold_2's start
            # and qqTroRaco H2_scaffold_65's end (folding; the homolog's own call was just under the threshold)
            act, why = "flag_inversion", ("at a scaffold end, and no alignment to the other haplotype says it runs "
                                          "the other way: proposed, not turned")
        else:
            act = "inversion"
            why = ("turned in place; " + ("it aligns to the other haplotype in the other direction from its flanks"
                                          if hom else "no alignment to the other haplotype there: the contact decides"))
            flips.append((int(v.start_px), int(v.end_px)))
        rows.append(dict(action=act, scaffold=nm, length_mb=round((b_bp - a_bp) / 1e6, 2), partner="",
                         strand="-" if act == "inversion" else "", tag="", score=round(v.score, 2), next_best="",
                         evidence=ev + "; " + why))
    return flips, rows
