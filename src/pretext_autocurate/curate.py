"""Single-haplotype first pass: cuts, joins, haplotig calls, placement -> edits, paint, tags.

Order: misjoin cuts (preferably at gaps) and cuts around haplotig contigs; tentative cuts
at weak junctions kept only when they enable a rearrangement; inversions inside a piece (the
bow-tie: cut at both ends, rejoined reversed); joins and orientation of
large pieces; then haplotig calls and placement of small scaffolds on the cut map. The
result is a target pixel order, turned into a PretextView edit list by `plan_edits`.
"""
import collections
from dataclasses import dataclass

import numpy as np

from .haplotig import detect_haplotigs, detect_contig_haplotigs, Params as HapParams
from .misjoin import detect_misjoins, weak_junctions, Cut
from .inversion import detect_inversions, Params as InvParams
from .overlap import classify_candidates, gap_overlaps, overlap_evidence, resolve_overlaps
from .join import find_joins
from .pretextio import SegmentedMap
from .savestate import Edit, BREAK_DELTA
from .mapstate import MapState, plan_edits
from .place import place_scaffolds, Params as PlaceParams
from .savestate import SaveState, Waypoint
from .signal import Signal


@dataclass
class Decision:
    """One row of the decisions TSV, named by original scaffolds and bp (not pieces)."""
    scaffold: str
    length_px: int
    length_bp: int
    action: str          # cut | flag_cut | join | flag_join | invert | flag_invert | insert |
                         # attach_unloc | attach_haplotig | tag_haplotig | none
    target: str          # anchor scaffold (for cuts: the cut scaffold; joins: the other one)
    target_bp: int
    strand: str
    confidence: str
    evidence: str


def _vt_end_px(vt, u, e):
    """Original pixel at end e ('h' / 't') of piece u of a SegmentedMap."""
    return int(vt.scaffold_start[u]) if e == "h" else int(vt.scaffold_start[u] + vt.scaffold_npix[u] - 1)


def _find_chain(chains, u):
    return next((k for k, ch in enumerate(chains) if any(s == u for s, _ in ch)), None)


def _unlink(chains, u, end):
    """Break the chain link at unit u's end ('h' or 't'), if there is one."""
    k = _find_chain(chains, u)
    if k is None:
        return
    ch = chains[k]
    i = next(i for i, (s, _) in enumerate(ch) if s == u)
    right = (end == "t") != ch[i][1]                # does that end face the chain's right?
    parts = (ch[:i + 1], ch[i + 1:]) if right else (ch[:i], ch[i:])
    chains[k:k + 1] = [p for p in parts if len(p) > 1]


def _link(chains, a, b, a_end="t", b_end="h"):
    """Join unit a's `a_end` to unit b's `b_end` in `chains` (lists of (unit, reversed)) if
    both ends are free chain ends; units in no chain count as chains of one. Returns whether
    it did. The default joins a's tail to b's head."""
    flip = lambda ch: [(s, not r) for s, r in reversed(ch)]
    ka, kb = _find_chain(chains, a), _find_chain(chains, b)
    if ka is not None and ka == kb:                 # already in one chain
        return False
    ca = chains[ka] if ka is not None else [(a, a_end == "h")]
    cb = chains[kb] if kb is not None else [(b, b_end == "t")]
    # a's end last in its chain (facing right), b's end first (facing left)
    ra, rb = a_end == "h", b_end == "t"
    ca = ca if ca[-1] == (a, ra) else flip(ca) if ca[0] == (a, not ra) else None
    cb = cb if cb[0] == (b, rb) else flip(cb) if cb[-1] == (b, not rb) else None
    if ca is None or cb is None:                    # an end is joined already
        return False
    for k in sorted((k for k in (ka, kb) if k is not None), reverse=True):
        chains.pop(k)
    chains.append(ca + cb)
    return True


def small_scaffold_haplotigs(raw, cut_after, hap_params=None, min_px=8, big_px=400, cov_ok=0.8):
    """Cuts that take a haplotig contig off the end of a small scaffold.

    `detect_contig_haplotigs` looks only inside scaffolds over `big_px`; a smaller scaffold is
    judged whole, so one whose end contig is a haplotig and the rest not averages the two
    (xgTylRafi scaffold_30: contigs at 1.08 / 0.91 / 0.68 coverage, called a borderline
    haplotig as a whole; the curator cut at the gap, placed the rest and tagged the last
    contig Haplotig). Here each scaffold of up to `big_px` with gaps is split at them on trial
    and its contigs called one by one: a run of contigs at one end, together at least `min_px`
    and each called haplotig (line and low coverage), is cut off at its gap when the rest has
    no such call and coverage of at least `cov_ok`. Returns (cuts, haplotig piece pixel ranges).
    """
    gap = raw.graphs.get("gap")
    if gap is None:
        return [], []
    cand = {}
    for s in np.nonzero((raw.scaffold_npix <= big_px) & (raw.scaffold_npix >= 2 * min_px))[0]:
        a, L = int(raw.scaffold_start[s]), int(raw.scaffold_npix[s])
        g = [a + int(x) for x in np.nonzero(gap[a:a + L - 1] > 0)[0]]       # cut after a gap pixel
        if g and not any(a <= q < a + L - 1 for q in cut_after):
            cand[int(s)] = g
    if not cand:
        return [], []
    vt = SegmentedMap(raw, list(cut_after) + [q for g in cand.values() for q in g])
    sg = Signal(vt)
    cov = sg.coverage_ratio()
    if cov is None:
        return [], []
    calls = {c.start_px: c for c in detect_haplotigs(vt, hap_params or HapParams(), sg)}
    cuts, pieces = [], []
    for s, g in cand.items():
        a, L = int(raw.scaffold_start[s]), int(raw.scaffold_npix[s])
        bounds = [a] + [q + 1 for q in g] + [a + L]
        parts = [(lo, hi) for lo, hi in zip(bounds, bounds[1:]) if hi > lo]       # [lo, hi)
        hap = [calls.get(lo) is not None and calls[lo].call in ("haplotig_high", "haplotig_medium")
               for lo, _ in parts]
        for run in (list(range(len(parts))), list(range(len(parts) - 1, -1, -1))):   # from either end
            k = 0
            while k < len(run) and hap[run[k]]:
                k += 1
            if not 0 < k < len(parts):
                continue
            hp, rest = [parts[i] for i in run[:k]], [parts[i] for i in run[k:]]
            lo, hi = min(x for x, _ in hp), max(y for _, y in hp)
            r0, r1 = min(x for x, _ in rest), max(y for _, y in rest)
            rc = float(np.median(cov[r0:r1]))
            if hi - lo < min_px or any(hap[i] for i in run[k:]) or rc < cov_ok:
                continue
            c = calls[hp[0][0]]
            q = hi - 1 if lo == a else lo - 1            # cut after this pixel
            cuts.append(Cut(raw.names[s], q, int(round((raw.rel_coord[q] + 1) * raw.bp_per_pix)), True,
                            float("nan"), float("nan"), True,
                            f"small scaffold: the contig{'s' if k > 1 else ''} at its "
                            f"{'start' if lo == a else 'end'} ({hi - lo} px) "
                            f"{'are haplotigs' if k > 1 else 'is a haplotig'} (line {c.line_over_own:.2f} at "
                            f"{c.partner}, coverage {float(np.median(cov[lo:hi])):.2f}, partner coverage "
                            f"{c.partner_cov:.2f}), the rest not (coverage {rc:.2f}): cut at the gap"))
            pieces.append((lo, hi))
            break
    return cuts, pieces


# Level-dependent thresholds of the single-haplotype rules pass: name -> (value at level 0, 0.5,
# 1, meaning), read through allhaps' `Level` like the allhaps table. Curator, 2026-11-19: 1.2 is an
# acceptable orientation ratio (scaffold_119 into scaffold_13 at 11.04 Mb, '+', 1.26); more
# examples will tune these.
SINGLE_THRESHOLDS = {
    "orient_strict": (1.3, 1.2, 1.15, "orientation ratio to propose an insertion with a waypoint at the gap"),
    "insert_orient": (1.4, 1.2, 1.15, "orientation ratio to apply an insertion"),
    "insert_min_px": (20, 15, 8, "smallest piece (px) whose insertion is applied, not only proposed"),
    # Inversions (inversion.py): weaker calls are proposed with waypoints. Over four curated
    # genomes only icLarUrsu scaffold_13 (gain 1.25, drop 0.29) passes 1.0 / 0.3; at 0.8 / 0.4
    # bAegFun H2_scaffold_8 (1.08 / 0.32, not curated) does too. Half bow-ties at scaffold ends
    # were right about half the time (folding at chromosome ends looks the same): applied at level 1 only.
    "inv_gain": (1.2, 1.0, 0.8, "inversion with both junctions: reversed corners (1 = adjacent) needed to apply"),
    "inv_drop": (0.2, 0.3, 0.4, "inversion with both junctions: assembled corners at most this to apply"),
    "inv_end_gain": (9.9, 9.9, 1.0, "inversion at a scaffold end (half a bow-tie): reversed corner needed to apply"),
    # Displaced stretches (linescan.py). xgTylRafi at 1.5: the four the curator moved, nothing else
    # (line over the next best place 2.9, 2.7, 1.7, 1.65; the curated haplotig scaffold_2
    # 109.69-111.11 Mb 1.05).
    "move_line": (0.35, 0.25, 0.2, "displaced stretch: line contact (over adjacent sequence) at least"),
    "move_ratio": (4.0, 3.0, 2.0, "displaced stretch: line contact over contact with its own sides"),
    "move_spec": (2.0, 1.5, 1.3, "displaced stretch: line over the next best place"),
}


def single_thresholds(level=0.5):
    """{name: value} of SINGLE_THRESHOLDS at `level` (insert_min_px rounded to whole pixels)."""
    from .allhaps.common import Level
    lv = Level(level)
    out = {k: lv.at(c, d, a) for k, (c, d, a, _) in SINGLE_THRESHOLDS.items()}
    out["insert_min_px"] = int(round(out["insert_min_px"]))
    return out


def _displaced(raw, sg, cuts, th, fasta, minimap2, exact, fills, inversions):
    """Displaced stretches of `raw` (linescan.scan), each cut out (applied cuts appended to
    `cuts`, or flagged ones there switched on) and decided with the sequence when there is
    some. Copies at an end are cut off (5 px or more) and switch rows appended to `exact`.

    Returns records (dicts): d (linescan.Displaced), name (scaffold), kind ('move': into the
    place, 'haplotig': a copy of the place, tagged next to it, 'unloc': strand unclear,
    attached), t (after this pixel), strand, target (scaffold name), copies [(first px, last
    px)] cut off as Haplotigs, mid (a pixel of the moved part), al (StretchAlignment or None),
    why (evidence)."""
    from . import linescan
    lp = linescan.Params(line_min=th["move_line"], move_ratio=th["move_ratio"], spec=th["move_spec"])
    applied = sorted(c.after_px for c in cuts if c.apply)
    found = linescan.scan(raw, sg, lp, peaks=linescan.cached_peaks(raw, lp), junctions=applied)
    bpp, so, gap = raw.bp_per_pix, raw.orig_scaffold, raw.graphs.get("gap")
    start, npx = raw.scaffold_start, raw.scaffold_npix
    busy = {f["s"] for f in fills}
    inv_spans = [(lq if lq is not None else v.start_px, rq if rq is not None else v.end_px)
                 for v, apply, lq, rq in inversions if apply]
    gone = [(d.a, d.b) for d in found]
    bp_of = lambda q: int(round((raw.rel_coord[q] + 1) * bpp))          # bp after pixel q (1-based end)

    def out_of_gone(t):
        """A place inside another stretch that leaves: just before that stretch instead."""
        for _ in range(len(gone)):
            hit = [(ga, gb) for ga, gb in gone if ga <= t <= gb]
            if not hit:
                break
            t = hit[0][0] - 1
        return t

    def junctions(ts):
        """Pixels of scaffold ts after which a piece can go: gaps, applied cuts, its last pixel."""
        a0, a1 = int(start[ts]), int(start[ts] + npx[ts] - 1)
        out = {a1} | {q for q in applied if a0 <= q <= a1} | {c.after_px for c in cuts if c.apply and a0 <= c.after_px <= a1}
        if gap is not None:
            out |= {int(q) for q in a0 + np.nonzero(gap[a0:a1 + 1] > 0)[0]}
        return sorted(out)

    def cut_at(q, why):
        """Switch on a cut after q (an existing one within 2 px is used): returns its pixel."""
        near = [c for c in cuts if c.scaffold == raw.names[int(so[q])] and abs(c.after_px - q) <= 2]
        if near:
            c = min(near, key=lambda c: abs(c.after_px - q))
            if not c.apply:
                c.apply, c.reason = True, c.reason.split(": not cut")[0] + "; " + why
            return c.after_px
        cuts.append(Cut(raw.names[int(so[q])], q, bp_of(q), bool(gap is not None and gap[q] > 0), float("nan"),
                        float("nan"), True, why))
        return q

    recs = []
    for d in found:
        name, s = raw.names[d.s], d.s
        if s in busy or any(not (d.b < x0 or d.a > x1) for x0, x1 in inv_spans):
            continue                    # already moved or turned for another reason
        if any(d.a + 1 < c.after_px < d.b - 2 for c in cuts if c.apply and so[c.after_px] == s):
            continue                    # cut inside for another reason: not one stretch
        ts = int(so[d.t])
        target = raw.names[ts]
        L = d.b - d.a + 1
        a_bp, b_bp = int(round(raw.rel_coord[d.a] * bpp)), bp_of(d.b)
        where = lambda t: f"{raw.names[int(so[t])]} {bp_of(t) / 1e6:.2f} Mb"
        why = (f"{a_bp / 1e6:.2f}-{b_bp / 1e6:.2f} Mb: contact with its own sides {d.left:.2f} / {d.right:.2f} of "
               f"adjacent sequence, a line to {target} at {line_txt(d, raw)} ({d.line:.2f}; next best place "
               f"{d.extra['row_out'] / max(d.extra['row_in'], 1e-9) * d.line:.2f})")
        al = None
        if fasta:
            w0, w1 = d.extra["window"]
            t0 = int(round(raw.rel_coord[max(int(start[ts]), w0)] * bpp)) - L * bpp
            t1 = bp_of(min(int(start[ts] + npx[ts] - 1), w1)) + L * bpp
            al = linescan.align_stretch(fasta, name, a_bp, b_bp, target, t0, t1, minimap2, end_tol=max(20_000, 3 * bpp))
        # the ends first: a cut already there (within 2 px) is used, and the copies are bounded by it
        qa = cut_at(d.a - 1, "start of a displaced stretch (see the move row)")
        qb = cut_at(d.b, "end of a displaced stretch (see the move row)")
        rec = dict(d=d, name=name, kind="move", t=d.t, strand=d.strand, target=target, copies=[], al=al,
                   mid=(qa + 1 + qb) // 2, qa=qa, qb=qb)
        if al is not None and al.strand and al.aligned >= 0.7 and al.divergence <= 0.05:
            rec["kind"] = "haplotig"
            why += (f"; {al.aligned:.0%} of it aligns to {target} there ({al.divergence:.1%} divergence): a haplotypic "
                    f"copy of that place, tagged Haplotig next to it")
        else:
            if al is not None and al.strand and al.weight >= 10_000:
                rec["strand"] = al.strand
                why += f"; strand {al.strand} from the sequence ({al.weight / 1e3:.0f} kb aligned to {target} there)"
            elif d.orient >= th["insert_orient"]:
                why += f"; strand {d.strand} from contact (orientation ratio {d.orient:.2f})"
            else:
                rec["kind"] = "unloc"
                why += (f"; strand unclear (orientation ratio {d.orient:.2f}, nothing aligns): attached to {target} as "
                        f"an Unloc, waypoint at the place")
            # a copy at an end decides the place: the stretch goes right before (or after) it
            for end, cp, sw in (("head", al.head if al else None, al.head_switch if al else None),
                                ("tail", al.tail if al else None, al.tail_switch if al else None)):
                if cp is None or rec["kind"] != "move":
                    continue
                q0, q1, c0, c1 = cp
                right = (end == "tail") == (rec["strand"] == "+")       # the copy meets what follows the place
                js = junctions(ts)
                if right:
                    tp = int(start[ts]) + int(c0 / bpp)
                    cand = [q for q in js if tp - L - 2 * lp.W <= q <= tp + 1]
                    t = max(cand) if cand else None
                else:
                    tp = int(start[ts]) + int(c1 / bpp)
                    cand = [q for q in js if tp - 1 <= q <= tp + L + 2 * lp.W]
                    t = min(cand) if cand else None
                if t is not None:
                    rec["t"] = t
                why += (f"; its {end} ({(q1 - q0) / 1e3:.0f} kb) is a copy of {target} {c0 / 1e6:.2f}-{c1 / 1e6:.2f} Mb")
                if (q1 - q0) / bpp >= 5:
                    # cut off at the base the copy stops; the rest goes in, the copy is a Haplotig
                    inner = q1 if end == "head" else q0
                    qc = int(start[s]) + int(round(inner / bpp)) - 1
                    if qa + 1 <= qc < qb:
                        cut_at(qc, f"end of the copy of {target} at the {end} of a displaced stretch (see the move row)")
                        exact.append(dict(kind="cut", scaffold=name, px_bp=bp_of(qc), bp=inner, partner="", partner_bp=0,
                                          method="alignment", src="fill",
                                          evidence=f"copy of {target} {c0 / 1e6:.2f}-{c1 / 1e6:.2f} Mb ends here"))
                        rec["copies"].append((qa + 1, qc) if end == "head" else (qc + 1, qb))
                        rec["mid"] = ((qc + 1 + qb) // 2) if end == "head" else ((qa + 1 + qc) // 2)
                        why += ", cut off as a Haplotig"
                        cut_len = True
                    else:
                        cut_len = False
                else:
                    cut_len = False
                if sw is not None:
                    q_sw, t_sw = sw
                    if right:                   # the stretch, then the place's sequence after it
                        if rec["strand"] == "+":
                            face = q0 if cut_len else b_bp
                        else:
                            face = q1 + 1 if cut_len else a_bp + 1
                        row = dict(kind="switch", scaffold=name, px_bp=face, bp=q_sw + 1, partner=target,
                                   partner_bp=t_sw + 1)
                    else:                       # the place's sequence before it, then the stretch
                        row = dict(kind="switch", scaffold=target, px_bp=bp_of(rec["t"]), bp=t_sw + 1, partner=name,
                                   partner_bp=q_sw + 1)
                    exact.append(dict(row, method="alignment", src="fill",
                                      evidence=f"displaced stretch of {name} into {target}: identical around the switch"))
                    why += "; joined without a gap where the two are identical"
        rec["t"] = out_of_gone(rec["t"])
        rec["why"] = why + f"; goes in after {where(rec['t'])}, strand {rec['strand']}"
        recs.append(rec)
    return recs


def line_txt(d, raw):
    """Where a displaced stretch's line points, in Mb of the target scaffold."""
    s = int(raw.orig_scaffold[d.t])
    return f"{(d.t - raw.scaffold_start[s] + 1) * raw.bp_per_pix / 1e6:.2f} Mb"


def curate_single(pt, hap_params=None, place_params=None, unloc_min_px=3, cov_low=0.8,
                  orient_strict=None, apply_inserts=False, do_cuts=True, do_joins=True, insert_min_px=None,
                  insert_orient=None, level=0.5, fasta=None, threads=4, minimap2=None):
    """Returns (decisions, savestate). The savestate holds edits, paint, tags, waypoints, and
    `st.exact_breaks`: cuts and gapless joins placed to the base by the sequence (rows for
    `breaks.write`; not part of the savestate file).

    `fasta`: the assembly; contig-end overlaps (overlap.py) are applied only when its sequence
    confirms them (minimap2 from `minimap2` or PATH / the conda env), otherwise proposed.

    Waypoints mark positions inside contigs only (cuts, flagged cuts, proposed insertion
    loci). Standalone pieces that need a look (borderline haplotigs) carry the `vertPaint`
    tag instead, which PretextView draws as a vertical highlight band.

    An insertion into a chromosome is applied only for a piece of at least `insert_min_px`
    whose orientation ratio reaches `insert_orient`. Smaller or less clear ones are only
    proposed (attached as an unloc, waypoint at the proposed locus): on a saturated
    MQ0-only map a small piece that belongs at a gap cannot be told from a haplotig of the
    sequence next to it. `apply_inserts` applies every proposal.

    `unloc_min_px`: smaller unlocalised pieces are left alone; `cov_low`: coverage ratio
    below which a piece is not inserted (possible haplotig); `orient_strict`: orientation
    ratio needed to consider an insertion at all. `do_cuts` / `do_joins` switch off stages.
    `level` (0 conservative .. 1 aggressive) sets `orient_strict`, `insert_orient` and
    `insert_min_px` from SINGLE_THRESHOLDS unless they are given.
    """
    th = single_thresholds(level)
    orient_strict = th["orient_strict"] if orient_strict is None else orient_strict
    insert_orient = th["insert_orient"] if insert_orient is None else insert_orient
    insert_min_px = th["insert_min_px"] if insert_min_px is None else insert_min_px
    raw = pt
    sg_raw = Signal(raw)
    cuts = detect_misjoins(raw, signal=sg_raw) if do_cuts else []
    contig_haps = detect_contig_haplotigs(raw, sg_raw) if do_cuts else []
    for h in contig_haps:           # cut at the gaps on both sides of a haplotig contig
        s0 = int(raw.scaffold_start[raw.names.index(h.scaffold)])
        for q in (h.start_px - 1, h.end_px):    # cut after q; skip scaffold ends (no cut)
            if q >= s0 and raw.orig_scaffold[q] == raw.orig_scaffold[min(q + 1, raw.npix - 1)] \
                    and not any(c.after_px == q and c.apply for c in cuts):
                cuts.append(Cut(h.scaffold, q, int(round((raw.rel_coord[q] + 1) * raw.bp_per_pix)), True,
                                float("nan"), float("nan"), True,
                                f"end of a haplotig contig: copy of sequence {abs(h.offset_px)} px away, "
                                f"coverage {h.cov:.2f}, partner coverage {h.partner_cov:.2f}, line {h.line:.2f}"))
    small_haps = []
    if do_cuts:
        sc, small_haps = small_scaffold_haplotigs(raw, [c.after_px for c in cuts if c.apply], hap_params)
        cuts += sc
    # Rearrangements inside a piece: cut tentatively at weak junctions (gaps), and keep a
    # cut only if a piece then fits better somewhere else or the other way round.
    jc = sg_raw.insulation(10, 0)       # contact across every boundary of the original map
    weak = weak_junctions(raw, sg_raw, [c.after_px for c in cuts]) if do_cuts and do_joins else []
    rearranged = set()
    if weak:
        applied = [c.after_px for c in cuts if c.apply]
        # trial map with every weak junction cut; sg_raw still serves it (pixels not moved)
        vt = SegmentedMap(raw, applied + [q for q, _, _ in weak])
        jt, _ = find_joins(vt, sg_raw)
        partner = {}
        for j in jt:
            if j.accepted:
                partner[(j.a, j.a_end)], partner[(j.b, j.b_end)] = (j.b, j.b_end), (j.a, j.a_end)
        small = {}
        for q, _, _ in weak:
            A, B = int(vt.orig_scaffold[q]), int(vt.orig_scaffold[q + 1])   # pieces either side
            # 150 px = join.Params.min_px: both pieces large enough to be joined. Keep the cut
            # if either end found a different accepted partner than the piece it was cut from
            if vt.scaffold_npix[A] > 150 and vt.scaffold_npix[B] > 150:
                pa, pb = partner.get((A, "t")), partner.get((B, "h"))
                if (pa is not None and pa != (B, "h")) or (pb is not None and pb != (A, "t")):
                    rearranged.add(q)
            for u in (A, B):                    # small pieces: try placing them instead
                if insert_min_px <= vt.scaffold_npix[u] <= 150:
                    small.setdefault(u, []).append(q)
        if small:
            for pl in place_scaffolds(vt, place_params or PlaceParams(), Signal(vt), only=set(small),
                                      junction_contact=jc):
                u = vt.names.index(pl.scaffold)
                a, b = int(vt.scaffold_start[u]), int(vt.scaffold_start[u] + vt.scaffold_npix[u] - 1)
                # placed back where it came from (within 8 px): only worth a cut if reversed
                home = min(abs(pl.after_px - (a - 1)), abs(pl.after_px - b)) <= 8
                if pl.call == "localised" and pl.orient >= insert_orient and not pl.cov < cov_low \
                        and (not home or pl.strand == "-"):
                    rearranged.update(small[u])
        for q, nv, fv in weak:
            ok = q in rearranged
            cuts.append(Cut(raw.names[int(raw.orig_scaffold[q])], q,
                            int(round((raw.rel_coord[q] + 1) * raw.bp_per_pix)), True, nv, fv, ok,
                            "weak junction; a piece next to it fits better elsewhere or reversed" if ok else
                            "weak junction, but the pieces fit best as assembled (possible insertion site): not cut"))
    # Inversions inside the pieces made so far (a segment crossing a cut is no inversion). An
    # applied one is cut at both ends (left_q / right_q: cut after; None at a scaffold end) and
    # rejoined reversed after the joins; weaker ones are proposed with waypoints.
    inversions = []                     # (Inversion, apply, left_q, right_q)
    if do_cuts and do_joins:
        vt = SegmentedMap(raw, sorted(c.after_px for c in cuts if c.apply))
        gap = raw.graphs.get("gap")
        so = raw.orig_scaffold

        def snap(q):                    # to a gap within 5 px of the same scaffold, if any
            if gap is None:
                return q, False
            near = [x for x in range(q - 5, q + 6) if 0 <= x < raw.npix - 1 and gap[x] > 0 and so[x] == so[q]
                    and so[x + 1] == so[q]]
            return (min(near, key=lambda x: abs(x - q)), True) if near else (q, False)
        for v in detect_inversions(vt, InvParams()):
            end = v.at_start or v.at_end
            apply = v.gain >= th["inv_end_gain"] if end else (v.gain >= th["inv_gain"] and v.drop <= th["inv_drop"])
            lq = None if v.start_px == 0 or vt.orig_scaffold[v.start_px - 1] != vt.orig_scaffold[v.start_px] \
                else v.start_px - 1
            rq = None if v.end_px + 1 >= raw.npix or vt.orig_scaffold[v.end_px + 1] != vt.orig_scaffold[v.end_px] \
                else v.end_px
            gaps = []
            if lq is not None:
                lq, g = snap(lq)
                gaps.append(g)
            if rq is not None:
                rq, g = snap(rq)
                gaps.append(g)
            inversions.append((v, apply, lq, rq))
            if apply:
                a_, b_ = (lq + 1 if lq is not None else v.start_px), (rq if rq is not None else v.end_px)
                span = "%.2f-%.2f Mb" % (raw.rel_coord[a_] * raw.bp_per_pix / 1e6, (raw.rel_coord[b_] + 1) * raw.bp_per_pix / 1e6)
                for q, g in zip([x for x in (lq, rq) if x is not None], gaps):
                    cuts.append(Cut(raw.names[int(so[q])], q, int(round((raw.rel_coord[q] + 1) * raw.bp_per_pix)), g,
                                    float("nan"), float("nan"), True,
                                    f"end of an inversion ({span} reversed; see the invert row)"))
    # Contig-end overlaps (overlap.py): one piece's end is the other haplotype of another's end.
    # The copy on the longer piece is cut off (the curator, icLarUrsu scaffold_16 and 14) at the
    # base where the sequence says it stops, tagged Haplotig, and the other piece's end is joined
    # to the cut. Applied only when the sequence confirms it; otherwise proposed.
    overlaps = []                       # dicts (overlap.resolve_overlaps)
    exact = []                          # rows for breaks.write
    if do_cuts and do_joins:
        vt = SegmentedMap(raw, sorted(c.after_px for c in cuts if c.apply))
        overlaps = resolve_overlaps(raw, vt, fasta, minimap2, threads)
        gp = raw.graphs.get("gap")
        for rec in overlaps:
            if rec["apply"]:
                ts, ks, q = rec["tside"], rec["kside"], rec["q"]
                cuts.append(Cut(ts.name, q, rec["cut_bp"], bool(gp is not None and gp[q] > 0), float("nan"), float("nan"),
                                True, f"end of a haplotypic copy of {ks.name} {'start' if ks.head else 'end'}, "
                                      f"placed by the sequence at {rec['cut_bp']:,} bp (see the overlap row)"))
                exact += rec["exact"]
    # Gaps whose two sides are haplotypic copies of each other (sequence only; nothing changes on
    # the map): the gap is closed in the output where the copies are identical (a switch row for
    # asm), the other copy of each side becoming a haplotig. icLarUrsu scaffold_11 at 54.11 Mb:
    # 849 kb on both sides at 0.77% divergence, both kept in the curated assembly.
    gap_rows = []
    if fasta and do_cuts and raw.graphs.get("gap") is not None:
        gap_rows = gap_overlaps(raw, fasta, minimap2, threads)
        exact += [r["switch"] for r in gap_rows if r["switch"]]
    # Candidate haplotigs aligned to where their Hi-C copy line points (curator, 2026-11-20:
    # always, not only at gaps): a confirmed copy, nothing, or a gap fill whose two ends copy the
    # gap's flanks (xgTylRafi scaffold_29 into scaffold_1's 18.23 Mb gap, '-', its first 236 kb a
    # copy of the flank after it). A fill's copies of 5 px or more are cut off at the base.
    copy_calls, fills = {}, []
    if fasta and do_cuts and do_joins:
        vtc = SegmentedMap(raw, sorted(c.after_px for c in cuts if c.apply))
        taken = {r["tside"].name for r in overlaps if r["apply"]} | {r["kside"].name for r in overlaps if r["apply"]}
        cands = []
        for c in detect_haplotigs(vtc, hap_params or HapParams(), Signal(vtc)):
            if not c.partner or c.call in ("too_small", "none"):
                continue
            par, off = vtc.parent_of(c.scaffold)
            s = raw.names.index(par)
            if off != 0 or vtc.scaffold_npix[vtc.names.index(c.scaffold)] != raw.scaffold_npix[s] or par in taken:
                continue                        # whole input scaffolds only
            ppar, poff = vtc.parent_of(c.partner)
            if ppar != par:
                cands.append((par, ppar, max(0, int(c.partner_bp - 1 + poff * raw.bp_per_pix))))
        copy_calls = classify_candidates(fasta, cands, minimap2, threads)
        N_, G_ = float(raw.npix), float(raw.total_len)
        for name, cc in copy_calls.items():
            if cc.kind != "gap_fill":
                continue
            if (cc.fill[1] - cc.fill[0]) / raw.bp_per_pix < 5:
                # mostly copy, with a sliver of its own over the gap: a haplotig of that region (the
                # curator tagged such pieces whole; cutting them cost xgTylRafi break precision 95% -> 66%)
                cc.kind = "copy"
                continue
            if min(cc.left[1] - cc.left[0], cc.right[1] - cc.right[0]) / raw.bp_per_pix < 2:
                # a copy under 2 px at one junction is thin evidence for the place (xgTylRafi scaffold_23:
                # 41 and 18 kb at 4.7%; the curator kept it as an Unloc of that chromosome): proposed
                cc.kind = "fill_weak"
                continue
            s, sp = raw.names.index(name), raw.names.index(cc.partner)
            a, n = int(raw.scaffold_start[s]), int(raw.scaffold_npix[s])
            rec = dict(cc=cc, s=s, q_left=int(raw.scaffold_start[sp]) + int(max(cc.gap[0] - 1, 0) / raw.bp_per_pix))
            L = int(raw.scaffold_len_bp()[s])
            # Gapless junctions (asm): the flank runs to its switch base, then the piece goes on past
            # the identical base; at the other junction the piece runs to its switch and the flank
            # goes on. '+': the piece's start copy meets the left flank; '-': its end copy does.
            for copy, sw, at_start in ((cc.left, cc.left_switch, True), (cc.right, cc.right_switch, False)):
                if sw is None:
                    continue
                q_sw, t_sw = sw
                cut = (copy[1] - copy[0]) / raw.bp_per_pix >= 5
                left_junction = at_start == (cc.strand == "+")
                if left_junction:               # flank before the gap, then the piece
                    row = dict(kind="switch", scaffold=cc.partner, px_bp=cc.gap[0], bp=t_sw + 1, partner=name,
                               partner_bp=q_sw + 1)
                else:                           # the piece, then the flank after the gap
                    face = (copy[1] + 1 if at_start else copy[0]) if cut else (1 if at_start else L)
                    row = dict(kind="switch", scaffold=name, px_bp=face, bp=q_sw + 1, partner=cc.partner,
                               partner_bp=t_sw + 1)
                exact.append(dict(row, method="alignment", src="fill",
                                  evidence=f"gap fill of {cc.partner}'s gap at {cc.gap[0] + 1:,}: {name} and {cc.partner} "
                                           f"identical around the switch"))
            for (x0, x1), at_start in ((cc.left, True), (cc.right, False)):
                if (x1 - x0) / raw.bp_per_pix < 5:
                    continue                    # a few px: left in the piece (the switch above trims it)
                bp = x1 if at_start else x0     # the copy ends here (cut after this base)
                q = a + int(round(bp / raw.bp_per_pix)) - 1
                if a <= q < a + n - 1:
                    cuts.append(Cut(name, q, bp, False, float("nan"), float("nan"), True,
                                    f"end of its copy of {cc.partner} by the gap at {cc.gap[0] + 1:,} bp, placed by "
                                    f"the sequence; the rest fills that gap (see the fill_gap row)"))
                    exact.append(dict(kind="cut", scaffold=name, px_bp=int((raw.rel_coord[q] + 1) * G_ / N_), bp=bp,
                                      partner="", partner_bp=0, method="alignment", src="fill",
                                      evidence=f"copy of {cc.partner} by its gap at {cc.gap[0] + 1:,} ends here"))
            fills.append(rec)
    # Displaced stretches (linescan.py; curator, 2026-11-21): little contact with their own sides
    # and a line to one place elsewhere. Cut out at both ends and moved there; the sequence, when
    # there is some, gives the strand and any copy at an end (cut off as a Haplotig when 5 px or
    # more, and joined without a gap where the two are identical).
    displaced = []
    if do_cuts and do_joins:
        displaced = _displaced(raw, sg_raw, cuts, th, fasta, minimap2, exact, fills, inversions)
    cuts.sort(key=lambda c: c.after_px)
    # from here on every applied cut is a scaffold end: pieces are "scaffolds" of pt
    pt = SegmentedMap(raw, [c.after_px for c in cuts if c.apply])
    forced = {pt.names[int(pt.orig_scaffold[h.start_px])] for h in contig_haps}   # haplotigs
    forced |= {pt.names[int(pt.orig_scaffold[lo])] for lo, _ in small_haps}
    joins, chains = find_joins(pt, sg_raw) if do_joins else ([], [])
    # Applied inversions: segment B between pieces A and C goes back reversed (A's tail to B's
    # tail, B's head to C's head). Those ends are the inversion's: a join that took one of them
    # is dropped to a flag. A, B and C stay out of placement.
    inv_units, inv_rows = set(), []
    for v, apply, lq, rq in inversions:
        if not apply:
            continue
        B = int(pt.orig_scaffold[lq + 1 if lq is not None else v.start_px])
        A = int(pt.orig_scaffold[lq]) if lq is not None else None
        C = int(pt.orig_scaffold[rq + 1]) if rq is not None else None
        own = {(A, "t"), (B, "h"), (B, "t"), (C, "h")} - {(None, "t"), (None, "h")}
        for u, e in own:
            _unlink(chains, u, e)
        for j in joins:
            if j.accepted and ((j.a, j.a_end) in own or (j.b, j.b_end) in own):
                j.accepted, j.reason = False, j.reason + "; an end of an inversion, rejoined there instead"
        ok = (A is None or _link(chains, A, B, "t", "t")) and (C is None or _link(chains, B, C, "h", "h"))
        inv_units |= {u for u in (A, B, C) if u is not None}
        inv_rows.append((v, ok))
    # Applied overlaps: the copy piece leaves (Haplotig of the other piece), whose end is joined
    # to the cut; joins that took those ends become flags.
    ov_haps = {}                        # copy piece -> (kept piece, overlap record)
    for rec in overlaps:
        if not rec["apply"]:
            continue
        q, ts, ks = rec["q"], rec["tside"], rec["kside"]
        Cp = int(pt.orig_scaffold[q if ts.head else q + 1])          # the copy
        R = int(pt.orig_scaffold[q + 1 if ts.head else q])            # the rest of that scaffold
        eR = "h" if ts.head else "t"
        kp, ke = rec["keep"]
        K = int(pt.orig_scaffold[_vt_end_px(vt, kp, ke)])
        own = {(K, ke), (R, eR), (Cp, "h"), (Cp, "t")}
        for u, e in own:
            _unlink(chains, u, e)
        made = {(K, ke), (R, eR)}
        for j in joins:
            if {(j.a, j.a_end), (j.b, j.b_end)} == made:
                j.accepted, j.reason = True, j.reason + "; made: the two sides of a haplotypic overlap (see the overlap row)"
            elif j.accepted and ((j.a, j.a_end) in own or (j.b, j.b_end) in own):
                j.accepted, j.reason = False, j.reason + "; an end of a haplotypic overlap, joined past the copy instead"
        rec["linked"] = _link(chains, K, R, ke, eR)
        rec["units"] = (Cp, R, K)
        inv_units |= {Cp, R, K}
        ov_haps[Cp] = (K, rec)
    # Gap fills: the middle piece goes into the gap in the strand the alignment gives, the copies
    # cut off its ends are Haplotigs of that scaffold; none of them is placed by contact.
    for f in fills:
        cc, s = f["cc"], f["s"]
        a, n = int(raw.scaffold_start[s]), int(raw.scaffold_npix[s])
        mid = a + min(n - 1, int((cc.fill[0] + cc.fill[1]) / 2 / raw.bp_per_pix))
        f["M"] = int(pt.orig_scaffold[mid])
        f["A"] = int(pt.orig_scaffold[f["q_left"]])
        f["copies"] = sorted({int(pt.orig_scaffold[p]) for p in range(a, a + n)} - {f["M"]})
        inv_units |= {f["M"], *f["copies"]}
    # Displaced stretches: the piece and the copies cut off it stay out of placement by contact
    for r in displaced:
        r["M"] = int(pt.orig_scaffold[r["mid"]])
        r["C"] = sorted({int(pt.orig_scaffold[x]) for x0, x1 in r["copies"] for x in range(x0, x1 + 1)} - {r["M"]})
        r["A"] = int(pt.orig_scaffold[r["t"]])
        inv_units |= {r["M"], *r["C"]}
        own = {(u, e) for u in (r["M"], *r["C"]) for e in "ht"}
        for u, e in own:
            _unlink(chains, u, e)
        for j in joins:
            if j.accepted and ((j.a, j.a_end) in own or (j.b, j.b_end) in own):
                j.accepted, j.reason = False, j.reason + "; an end of a displaced stretch, moved to its place instead"
    sg = Signal(pt)                     # new Signal: decay and insulation respect the cuts
    haps = {c.scaffold: c for c in detect_haplotigs(pt, hap_params or HapParams(), sg)}
    places = {p.scaffold: p for p in place_scaffolds(pt, place_params or PlaceParams(), sg, junction_contact=jc)
              if pt.names.index(p.scaffold) not in inv_units}
    # Several pieces may go into one junction (xgTylRafi scaffold_3 at 22.97 Mb: the curator put
    # scaffold_3 46.65-47.17 and part of scaffold_30 there, side by side); each must be clear on
    # its own, and they are ordered by their contact with the two sides (see `left_share`).
    idx = {n: i for i, n in enumerate(pt.names)}

    possible = []        # scaffold indices flagged as borderline haplotigs (vertPaint tag only)
    moved = {}           # scaffold index -> (kind, anchor index, after_px, strand)
    decisions = []
    for c in cuts:
        decisions.append(Decision(c.scaffold, 0, 0, "cut" if c.apply else "flag_cut", c.scaffold, c.bp, "",
                                  "high" if c.apply and c.at_gap else "medium" if c.apply else "low",
                                  (f"contact across: near {c.near:.2f}, far {c.far:.2f} of expectation; "
                                   if c.near == c.near else "") +     # NaN: haplotig-contig cut
                                  f"{'at gap' if c.at_gap else 'no gap'}; {c.reason}"))
    proposals = [c.after_px for c in cuts if not c.apply]   # waypoints at flagged cuts too
    linked = {id(v): ok for v, ok in inv_rows}
    for v, apply, lq, rq in inversions:
        a_, b_ = (lq + 1 if lq is not None else v.start_px), (rq if rq is not None else v.end_px)
        bp0, bp1 = int(round(raw.rel_coord[a_] * raw.bp_per_pix)), int(round((raw.rel_coord[b_] + 1) * raw.bp_per_pix))
        done = apply and linked.get(id(v), False)
        one = v.at_start or v.at_end
        ev = (f"{bp0 / 1e6:.2f}-{bp1 / 1e6:.2f} Mb reversed: contact across the reversed junctions "
              f"{v.gain:.2f}, across the assembled ones {v.drop:.2f} (1 = adjacent sequence, 0 = what sequence "
              f"{(bp1 - bp0) / 1e6:.1f} Mb apart shows, {v.expected:.2f} of adjacent here)" +
              ("; at a scaffold end: half a bow-tie, one junction" if one else "; bow-tie at both ends"))
        if apply and not done:
            ev += "; cut at its ends, but an end was already joined elsewhere: not rejoined reversed"
        elif not apply:
            ev += "; proposed (waypoints at its ends), not applied at this level"
            proposals += [q for q in (lq, rq) if q is not None] + ([v.start_px] if v.at_start else []) + \
                ([v.end_px] if v.at_end else [])
        name = raw.names[int(raw.orig_scaffold[a_])]
        decisions.append(Decision(name, b_ - a_ + 1, bp1 - bp0, "invert" if done else "flag_invert", name, bp0, "-",
                                  "high" if done and not one else "medium" if done else "low", ev))
    for rec in overlaps:
        o, ts, ks, al = rec["ov"], rec["tside"], rec["kside"], rec["al"]
        wh = lambda sd: f"{sd.name} {'start' if sd.head else 'end'}"
        ev = overlap_evidence(rec, raw.bp_per_pix)
        if rec["apply"] and rec.get("linked", False):
            ev += f"; cut at {rec['cut_bp']:,} bp, the copy tagged Haplotig, {wh(ks)} joined to the rest"
            act, conf = "overlap", "high"
        else:
            ev += "; " + (rec["why"] or "applied cut, but an end was already joined elsewhere: not joined")
            act, conf = "flag_overlap", "low"
            u, e = rec["trim"]
            p = _vt_end_px(vt, u, e)
            proposals.append(p + o.overlap_px - 1 if e == "h" else p - o.overlap_px + 1)
        decisions.append(Decision(ts.name, o.overlap_px, int(al.copy_bp if al is not None else o.overlap_px * raw.bp_per_pix),
                                  act, ks.name, int(rec["cut_bp"] or ts.edge), "", conf, ev))
    for g in gap_rows:
        al = g["al"]
        if al is None or al.copy_bp < 10_000:
            continue                    # nothing like a copy: no row
        ev = (f"gap at {g['gap'][0] + 1:,} bp: the {al.copy_bp / 1e3:.0f} kb after it aligns in order with the sequence "
              f"before it ({al.aligned_bp / 1e3:.0f} kb aligned, {al.divergence:.2%} divergence, reaching "
              f"{al.keep_reach_bp / 1e3:.0f} kb from the gap)")
        if g["ok"]:
            decisions.append(Decision(g["scaffold"], 0, al.copy_bp, "close_gap", g["scaffold"], g["gap"][0] + 1, "",
                                      "high", ev + f"; asm closes it where {al.match_bp:,} bases are identical "
                                      f"({g['scaffold']} to {al.keep_switch + 1:,}, then from {al.trim_switch + 2:,}); the "
                                      f"other copy of each side becomes a haplotig"))
        elif al.aligned_bp >= 0.5 * al.copy_bp and al.divergence <= 0.05:
            decisions.append(Decision(g["scaffold"], 0, al.copy_bp, "gap_copy", g["scaffold"], g["gap"][0] + 1, "",
                                      "", ev + "; not closed: the copy stops short of the gap or has no identical "
                                                  "stretch to switch at"))
    for row in exact:
        if row.get("src") in ("gap", "fill"):
            row["applied"] = True
        else:
            row["applied"] = any(r["apply"] and r.get("linked", False) and r["tside"].name in (row["scaffold"], row["partner"])
                                 for r in overlaps)
    st_exact = [{k: v for k, v in r.items() if k not in ("applied", "src")} for r in exact if r["applied"]]
    # "<scaffold> start|end (x Mb)" for a piece end, in original scaffold coordinates
    endname = lambda u, e: "%s %s (%.2f Mb)" % (pt.parent_of(pt.names[u])[0], "start" if e == "h" else "end", (
        pt.parent_of(pt.names[u])[1] + (0 if e == "h" else pt.scaffold_npix[u])) * pt.bp_per_pix / 1e6)
    for j in joins:
        par, off = pt.parent_of(pt.names[j.b])
        decisions.append(Decision(pt.parent_of(pt.names[j.a])[0], int(pt.scaffold_npix[j.a]), 0,
                                  "join" if j.accepted else "flag_join", par,
                                  int(round((off + (0 if j.b_end == "h" else pt.scaffold_npix[j.b])) * pt.bp_per_pix)), "",
                                  "high" if j.accepted else "low",
                                  f"{endname(j.a, j.a_end)} with {endname(j.b, j.b_end)}: end contact {j.score:.2f} "
                                  f"of adjacent level, next best {j.second:.2f}; {j.reason}"))
        if not j.accepted:                      # waypoints at both ends of a flagged join
            for u, e in ((j.a, j.a_end), (j.b, j.b_end)):
                proposals.append(int(pt.scaffold_start[u] + (0 if e == "h" else pt.scaffold_npix[u] - 1)))
    n_cut_rows = len(decisions)
    for name, cc in copy_calls.items():          # fills with a thin copy at one end: proposed only
        if cc.kind != "fill_weak":
            continue
        sp = raw.names.index(cc.partner)
        q = int(raw.scaffold_start[sp]) + int(max(cc.gap[0] - 1, 0) / raw.bp_per_pix)
        proposals.append(q)
        decisions.append(Decision(name, 0, cc.fill[1] - cc.fill[0], "flag_fill_gap", cc.partner, cc.gap[0] + 1, cc.strand,
                                  "low", f"its ends look like copies of the two sides of {cc.partner}'s gap at "
                                         f"{cc.gap[0] + 1:,} bp ({cc.left[1] / 1e3:.0f} kb and "
                                         f"{(cc.right[1] - cc.right[0]) / 1e3:.0f} kb, {cc.divergence:.1%} divergence, strand "
                                         f"{cc.strand}): it may fill that gap, but a copy under 2 px is thin evidence; "
                                         f"not moved (waypoint at the gap)"))
    for f in fills:                              # gap fills: the middle in, the copies as Haplotigs
        cc, A = f["cc"], f["A"]
        moved[f["M"]] = ("", A, f["q_left"], cc.strand)
        decisions.append(Decision(pt.names[f["M"]], int(pt.scaffold_npix[f["M"]]), cc.fill[1] - cc.fill[0], "fill_gap",
                                  pt.names[A], 0, cc.strand, "high",
                                  f"its ends are copies of the two sides of {cc.partner}'s gap at {cc.gap[0] + 1:,} bp "
                                  f"({cc.left[1] / 1e3:.0f} kb and {(cc.right[1] - cc.right[0]) / 1e3:.0f} kb, "
                                  f"{cc.divergence:.1%} divergence, strand {cc.strand}): its own "
                                  f"{(cc.fill[1] - cc.fill[0]) / 1e3:.0f} kb ({cc.fill[0] / 1e6:.2f}-{cc.fill[1] / 1e6:.2f} "
                                  f"Mb) inserted into the gap; copies of 5 px or more cut off as Haplotigs"))
        for c in f["copies"]:
            moved[c] = ("Haplotig", A, -1, "+")
            decisions.append(Decision(pt.names[c], int(pt.scaffold_npix[c]), 0, "attach_haplotig", pt.names[A], 0, "",
                                      "high", f"copy of the side of {cc.partner}'s gap at {cc.gap[0] + 1:,} bp, cut off "
                                              f"the piece that fills it (see the fill_gap row)"))
    for r in displaced:                          # displaced stretches: moved, or tagged next to their place
        d, M, A = r["d"], r["M"], r["A"]
        if r["kind"] == "move":
            moved[M] = ("", A, r["t"], r["strand"])
            act, conf = "move", "high"
        elif r["kind"] == "haplotig":
            moved[M] = ("Haplotig", A, -1, "+")
            act, conf = "attach_haplotig", "high"
        else:
            moved[M] = ("Unloc", A, -1, "+")
            proposals.append(r["t"])
            act, conf = "attach_unloc", "medium"
        decisions.append(Decision(r["name"], d.b - d.a + 1, int(round((d.b - d.a + 1) * pt.bp_per_pix)), act, r["target"],
                                  int(round((raw.rel_coord[r["t"]] + 1) * pt.bp_per_pix)),
                                  r["strand"] if act == "move" else "", conf, "displaced stretch: " + r["why"]))
        for c in r["C"]:
            moved[c] = ("Haplotig", A, -1, "+")
            decisions.append(Decision(r["name"], int(pt.scaffold_npix[c]), int(round(pt.scaffold_npix[c] * pt.bp_per_pix)),
                                      "attach_haplotig", r["target"], 0, "", "high",
                                      f"copy of {r['target']} cut off the displaced stretch moved there (see the move row)"))
    for Cp, (K, rec) in ov_haps.items():         # the copies cut off by applied overlaps
        if not rec.get("linked", False):
            continue
        moved[Cp] = ("Haplotig", K, -1, "+")
        decisions.append(Decision(pt.names[Cp], int(pt.scaffold_npix[Cp]), 0, "attach_haplotig", pt.names[K], 0, "",
                                  "high", f"haplotypic copy of the end of {rec['kside'].name}, cut off by the sequence "
                                          f"(see the overlap row)"))
    for name, pl in places.items():
        h = haps.get(name)
        is_hap = (h is not None and h.call in ("haplotig_high", "haplotig_medium")) or name in forced
        ev = f"cis {pl.cis:.2f}, best/next anchor {pl.ratio:.2f}, coverage {pl.cov:.2f}"
        if h is not None and h.partner:
            ev += f"; line/own {h.line_over_own:.2f} at {h.partner}:{h.partner_bp}, partner coverage {h.partner_cov:.2f}"
        cc = copy_calls.get(pt.parent_of(name)[0]) if pt.parent_of(name)[1] == 0 else None
        if cc is not None:
            ev += (f"; sequence: {cc.aligned:.0%} of it aligns in order to {cc.partner} there at {cc.divergence:.1%} "
                   f"divergence" if cc.strand else f"; sequence: nothing aligns to {cc.partner} there")
        conf = "high" if name in forced else h.call.split("_")[1] if is_hap else ""
        if h is not None and h.call == "haplotig_low" and not is_hap and cc is not None and cc.kind == "copy":
            is_hap, conf = True, "medium"          # a borderline call the sequence confirms
            ev += ": a haplotypic copy (borderline on contact and coverage, confirmed by the sequence)"
        if name in forced:
            ev += "; contig inside a scaffold, cut out at its gaps"
        if h is not None and h.call == "haplotig_low" and not is_hap:
            ev += "; POSSIBLE HAPLOTIG (" + h.reason + "), highlighted with vertPaint"
            possible.append(idx[name])
        row = dict(scaffold=name, length_px=pl.length_px, length_bp=pl.length_bp,
                   target=pl.anchor, target_bp=pl.anchor_bp, strand="")
        if pl.call == "unplaced":
            if is_hap:
                decisions.append(Decision(**row, action="tag_haplotig", confidence=conf,
                                          evidence=ev + "; " + pl.reason))
                moved[idx[name]] = ("Haplotig", None, -1, "+")
            else:
                decisions.append(Decision(**row, action="none", confidence="", evidence=ev + "; " + pl.reason))
            continue
        A = idx[pl.anchor]
        if is_hap:
            moved[idx[name]] = ("Haplotig", A, -1, "+")
            decisions.append(Decision(**row, action="attach_haplotig", confidence=conf, evidence=ev))
        elif (pl.call == "localised" and pl.orient >= orient_strict and not pl.cov < cov_low):
            row["strand"] = pl.strand
            ev += f"; orientation ratio {pl.orient:.2f}"
            at_end = pl.after_px in (int(pt.scaffold_start[A]) - 1,
                                     int(pt.scaffold_start[A] + pt.scaffold_npix[A]) - 1)
            # extending a chromosome end is left to the curator: as an unloc it lands there anyway
            clear = (pl.length_px >= insert_min_px and pl.orient >= insert_orient
                     and not at_end)
            if apply_inserts or clear:
                moved[idx[name]] = ("", A, pl.after_px, pl.strand)
                decisions.append(Decision(**row, action="insert", confidence="high" if clear else "medium",
                                          evidence=ev))
            else:
                moved[idx[name]] = ("Unloc", A, -1, "+")
                proposals.append(pl.after_px)
                decisions.append(Decision(**row, action="attach_unloc", confidence="medium",
                                          evidence=ev + "; proposed insertion at target (waypoint), not applied"))
        elif pl.call == "localised" and not pl.cov < cov_low:
            # placed at a gap, but the orientation is short of orient_strict: say where and which
            # way, and mark the gap (xgTylRafi scaffold_119 at scaffold_13 11.04 Mb, '+' at 1.26)
            moved[idx[name]] = ("Unloc", A, -1, "+")
            proposals.append(pl.after_px)
            row["strand"] = pl.strand
            row["target_bp"] = int(round((pt.rel_coord[pl.after_px] + 1) * pt.bp_per_pix))   # the gap
            decisions.append(Decision(**row, action="attach_unloc", confidence="medium",
                                      evidence=ev + f"; localised at a gap, strand {pl.strand} but orientation ratio "
                                      f"{pl.orient:.2f} (below {orient_strict}): proposed insertion at target "
                                      f"(waypoint), not applied"))
        elif pl.length_px >= unloc_min_px:
            why = pl.reason or ("low coverage: possible haplotig" if pl.cov < cov_low else "orientation weak")
            moved[idx[name]] = ("Unloc", A, -1, "+")
            decisions.append(Decision(**row, action="attach_unloc", confidence="high" if pl.ratio >= 2 else "medium",
                                      evidence=ev + "; not localised: " + why))
        else:
            decisions.append(Decision(**row, action="none", confidence="",
                                      evidence=ev + f"; under {unloc_min_px} px, too small to act on"))

    # placement rows name pieces of the cut map: translate to original scaffolds and bp
    # (cut/join rows were written with original names already)
    for d in decisions[n_cut_rows:]:
        for field in ("scaffold", "target"):
            name = getattr(d, field)
            if name and name not in pt.names:
                continue                        # already an input scaffold (flag_fill_gap rows)
            if name:
                par, off = pt.parent_of(name)
                if field == "target":
                    d.target_bp += int(round(off * pt.bp_per_pix))
                elif par != name:
                    d.evidence = f"piece {name} (from {int(round(off * pt.bp_per_pix)) + 1} bp); " + d.evidence
                setattr(d, field, par)

    def root(s):
        """Follow anchors that are themselves moved; such pieces become unlocs of the root."""
        seen, kind, A, q, strand = set(), *moved[s]
        while A in moved and A not in seen:
            seen.add(A)
            A, q = moved[A][1], -1
            kind = kind or "Unloc"
        return kind, A, q, strand

    start, npx = pt.scaffold_start, pt.scaffold_npix
    inserts, attach = {}, {}            # anchor -> [(after_px, piece, kind, strand)]
    for s in moved:
        kind, A, q, strand = root(s)
        if A is None:
            continue
        (attach if q == -1 and kind else inserts).setdefault(A, []).append((q, s, kind, strand))
    # target layout: original order, pieces removed and re-inserted
    # perm: original pixel per output position; unit: paint id (0 = unpainted);
    # tag: 0 / 1 Unloc / 2 Haplotig; marks: waypoint positions
    perm, unit, tag, marks = [], [], [], []

    def emit(px, u, t):
        """Append original pixels `px` to the layout with paint unit `u` and tag code `t`."""
        perm.append(px)
        unit.append(np.full(len(px), u))
        tag.append(np.full(len(px), t))

    def piece(s, strand):
        """Original pixels of piece `s`, reversed for strand '-'."""
        px = np.arange(start[s], start[s] + npx[s])
        return px[::-1] if strand == "-" else px

    tags = {"": 0, "Unloc": 1, "Haplotig": 2}

    def left_share(q, t, w=20):
        """Share of piece t's contact with the junction's left side (w px up to pixel q) over
        both sides: of pieces inserted at one junction, the most left-leaning goes first."""
        a, b = int(start[t]), int(start[t] + npx[t] - 1)
        left = float(sg.block(a, b, q - w + 1, q).mean())
        right = float(sg.block(a, b, q + 1, q + w).mean())
        return left / (left + right) if left + right > 0 else 0.5

    def unit_segments(s):
        """Forward pieces of a unit with its insertions: list of (pixels, is_inserted)."""
        a0, a1 = int(start[s]), int(start[s] + npx[s] - 1)
        segs, pos = [], a0
        ins = inserts.get(s, [])
        at = collections.Counter(q for q, *_ in ins)
        order = lambda x: (x[0], -left_share(x[0], x[1]) if at[x[0]] > 1 else 0, x[1])
        for q, t, kind, strand in sorted(ins, key=order):
            if q >= pos:
                segs.append((np.arange(pos, q + 1), False))
                pos = q + 1
            segs.append((piece(t, strand), True))
        if pos <= a1:
            segs.append((np.arange(pos, a1 + 1), False))
        return segs

    # Taking a haplotig contig out of a scaffold must not split the chromosome: when no join
    # claimed either flank's free end, the flanks stay together as assembled (xgTylRafi
    # scaffold_13 became two chromosomes when the flanks' best ends tied).
    for h in contig_haps if do_joins else []:
        L, R = int(pt.orig_scaffold[max(h.start_px - 1, 0)]), int(pt.orig_scaffold[min(h.end_px + 1, pt.npix - 1)])
        if h.start_px == 0 or h.end_px + 1 >= pt.npix or L in moved or R in moved or L == R \
                or raw.orig_scaffold[h.start_px - 1] != raw.orig_scaffold[h.end_px + 1]:
            continue
        if _link(chains, L, R):
            decisions.append(Decision(pt.parent_of(pt.names[L])[0], int(npx[L]), 0, "join",
                                      pt.parent_of(pt.names[R])[0], int(round(pt.parent_of(pt.names[R])[1] * pt.bp_per_pix)),
                                      "", "medium",
                                      f"{endname(L, 't')} with {endname(R, 'h')}: flanks of the haplotig contig "
                                      f"cut out between them, kept together as assembled (no clearer join for either end)"))
    # The same for the sides of a displaced stretch moved away (unless a join took their ends)
    for r in displaced if do_joins else []:
        L, R = int(pt.orig_scaffold[r["qa"]]), int(pt.orig_scaffold[min(r["qb"] + 1, pt.npix - 1)])
        if L in moved or R in moved or L == R or raw.orig_scaffold[r["qa"]] != raw.orig_scaffold[r["qb"] + 1]:
            continue
        if _link(chains, L, R):
            decisions.append(Decision(pt.parent_of(pt.names[L])[0], int(npx[L]), 0, "join",
                                      pt.parent_of(pt.names[R])[0], int(round(pt.parent_of(pt.names[R])[1] * pt.bp_per_pix)),
                                      "", "medium",
                                      f"{endname(L, t)} with {endname(R, h)}: sides of the displaced stretch moved "
                                      f"away from between them, kept together as assembled (no clearer join for either end)"))
    in_chain = {s for ch in chains for s, _ in ch}
    # output scaffolds: joined chains plus every other piece alone, in original map order
    groups = [ch for ch in chains] + [[(int(s), False)] for s in np.argsort(start)
                                      if npx[s] > 0 and int(s) not in in_chain]
    groups.sort(key=lambda ch: min(start[s] for s, _ in ch))
    for u, ch in enumerate(groups, 1):
        s0 = ch[0][0]
        if len(ch) == 1 and s0 in moved:                # emitted with its anchor instead
            if root(s0)[1] is None:                     # haplotig with no home: tag in place
                emit(piece(s0, "+"), 0, tags["Haplotig"])
            continue
        # paint only scaffolds that changed (joined, cut, or received pieces)
        has = len(ch) > 1 or any(s in inserts or s in attach or bool(pt.was_cut[s]) for s, _ in ch)
        for s, rev in ch:
            segs = unit_segments(s)
            if rev:
                segs = [(px[::-1], ins) for px, ins in reversed(segs)]
            for px, ins in segs:
                if ins:
                    marks.append(sum(map(len, perm)))
                emit(px, u if has else 0, 0)
        for kind in ("Unloc", "Haplotig"):              # at the end of the chromosome
            rest = [x for s, _ in ch for x in attach.get(s, []) if x[2] == kind]
            for q, t, k, strand in sorted(rest, key=lambda x: -npx[x[1]]):     # largest first
                emit(piece(t, "+"), u, tags[kind])
    perm, unit, tag = np.concatenate(perm), np.concatenate(unit), np.concatenate(tag)
    st = SaveState.blank()
    where0 = np.empty(raw.npix, dtype=np.int64)
    where0[perm] = np.arange(raw.npix)
    # an explicit break is only needed where both sides stay next to each other, in order;
    # elsewhere the moves already separate them. Breaks go first, on the unedited map.
    st.edits = [Edit(c.after_px, c.after_px, BREAK_DELTA) for c in cuts
                if c.apply and where0[c.after_px + 1] - where0[c.after_px] == 1] + plan_edits(perm)
    ms = MapState(raw).replay(st.edits)     # replay as PretextView will, to check and index
    if not np.array_equal(ms.pix, perm):
        raise AssertionError("planned edits do not reproduce the target layout")
    bits = {1: 1 << st.tag_id("Unloc"), 2: 1 << st.tag_id("Haplotig")}
    vert, highlight = 1 << st.tag_id("vertPaint"), set(possible)
    ids = {}                                # paint unit -> scaffold id numbered from 1
    # paint and tags are keyed by contig index after all edits (FORMATS.md, Savestate)
    for i, c in enumerate(ms.contigs):
        if unit[c.start]:
            st.scaffs.append((i, ids.setdefault(int(unit[c.start]), len(ids) + 1)))
        mask = bits[int(tag[c.start])] if tag[c.start] else 0
        if int(pt.orig_scaffold[ms.pix[c.start]]) in highlight:
            mask |= vert                 # highlight only: "look at this", implies nothing else
        if mask:
            st.meta_flags.append((i, mask))
    where = np.empty(pt.npix, dtype=np.int64)  # original pixel -> position in the new map
    where[perm] = np.arange(pt.npix)
    marks += [int(where[max(q, 0)]) for q in proposals]
    marks += [int(where[c.after_px]) for c in cuts if c.apply]
    for k, m in enumerate(sorted(marks)[:255]):     # waypoint count and ids are u8
        x = m / pt.npix - 0.5                       # on the diagonal (FORMATS.md)
        st.waypoints.append(Waypoint(x, x, 30.0, k))
    st.exact_breaks = st_exact
    return decisions, st


def write_decisions(rows, path):
    """Write Decisions as TSV, one row each, in the given order."""
    cols = ["scaffold", "length_px", "length_bp", "action", "target", "target_bp", "strand", "confidence", "evidence"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(getattr(r, c)) for c in cols) + "\n")
