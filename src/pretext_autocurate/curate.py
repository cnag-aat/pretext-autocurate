"""Single-haplotype first pass: cuts, joins, haplotig calls, placement -> edits, paint, tags.

Order: misjoin cuts (preferably at gaps) and cuts around haplotig contigs; tentative cuts
at weak junctions kept only when they enable a rearrangement; joins and orientation of
large pieces; then haplotig calls and placement of small scaffolds on the cut map. The
result is a target pixel order, turned into a PretextView edit list by `plan_edits`.
"""
from dataclasses import dataclass

import numpy as np

from .haplotig import detect_haplotigs, detect_contig_haplotigs, Params as HapParams
from .misjoin import detect_misjoins, weak_junctions, Cut
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
    action: str          # cut | flag_cut | join | flag_join | insert | attach_unloc |
                         # attach_haplotig | tag_haplotig | none
    target: str          # anchor scaffold (for cuts: the cut scaffold; joins: the other one)
    target_bp: int
    strand: str
    confidence: str
    evidence: str


def curate_single(pt, hap_params=None, place_params=None, unloc_min_px=3, cov_low=0.8,
                  orient_strict=1.3, apply_inserts=False, do_cuts=True, do_joins=True, insert_min_px=15, insert_orient=1.3):
    """Returns (decisions, savestate). The savestate holds edits, paint, tags, waypoints.

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
    """
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
    cuts.sort(key=lambda c: c.after_px)
    # from here on every applied cut is a scaffold end: pieces are "scaffolds" of pt
    pt = SegmentedMap(raw, [c.after_px for c in cuts if c.apply])
    forced = {pt.names[int(pt.orig_scaffold[h.start_px])] for h in contig_haps}   # haplotigs
    joins, chains = find_joins(pt, sg_raw) if do_joins else ([], [])
    sg = Signal(pt)                     # new Signal: decay and insulation respect the cuts
    haps = {c.scaffold: c for c in detect_haplotigs(pt, hap_params or HapParams(), sg)}
    places = {p.scaffold: p for p in place_scaffolds(pt, place_params or PlaceParams(), sg, junction_contact=jc)}
    # one insertion per junction: the piece with the clearest orientation keeps it
    best_at = {}
    for p_ in places.values():
        if p_.call == "localised" and (p_.after_px not in best_at or p_.orient > best_at[p_.after_px].orient):
            best_at[p_.after_px] = p_
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
    for name, pl in places.items():
        h = haps.get(name)
        is_hap = (h is not None and h.call in ("haplotig_high", "haplotig_medium")) or name in forced
        ev = f"cis {pl.cis:.2f}, best/next anchor {pl.ratio:.2f}, coverage {pl.cov:.2f}"
        if h is not None and h.partner:
            ev += f"; line/own {h.line_over_own:.2f} at {h.partner}:{h.partner_bp}, partner coverage {h.partner_cov:.2f}"
        conf = "high" if name in forced else h.call.split("_")[1] if is_hap else ""
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
                     and best_at.get(pl.after_px) is pl and not at_end)
            if apply_inserts or clear:
                moved[idx[name]] = ("", A, pl.after_px, pl.strand)
                decisions.append(Decision(**row, action="insert", confidence="high" if clear else "medium",
                                          evidence=ev))
            else:
                moved[idx[name]] = ("Unloc", A, -1, "+")
                proposals.append(pl.after_px)
                decisions.append(Decision(**row, action="attach_unloc", confidence="medium",
                                          evidence=ev + "; proposed insertion at target (waypoint), not applied"))
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

    def unit_segments(s):
        """Forward pieces of a unit with its insertions: list of (pixels, is_inserted)."""
        a0, a1 = int(start[s]), int(start[s] + npx[s] - 1)
        segs, pos = [], a0
        for q, t, kind, strand in sorted(inserts.get(s, [])):
            if q >= pos:
                segs.append((np.arange(pos, q + 1), False))
                pos = q + 1
            segs.append((piece(t, strand), True))
        if pos <= a1:
            segs.append((np.arange(pos, a1 + 1), False))
        return segs

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
    return decisions, st


def write_decisions(rows, path):
    """Write Decisions as TSV, one row each, in the given order."""
    cols = ["scaffold", "length_px", "length_bp", "action", "target", "target_bp", "strand", "confidence", "evidence"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(getattr(r, c)) for c in cols) + "\n")
