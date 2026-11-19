"""The allhaps pass (`curate_allhaps`): sites and cuts, pairing, small scaffolds, merges,
arrangement, layout, rows and savestate."""
import re
import numpy as np
from ..mapstate import MapState, plan_edits
from ..pretextio import SegmentedMap
from ..savestate import SaveState, Edit, Waypoint, BREAK_DELTA
from .arrange import _adjacency, _aligned_bp, _aligned_dir, _block_index, _contig_conflicts, _homolog_dir, _project, arrange_pieces
from .common import Thresholds, _signal, _union, split_block
from .evidence import _graph_evidence, _unitig_evidence
from .pairing import _quality, merge_pairs, pair_homologs, same_prefix_pairs
from .sites import telomere_sites
from . import chromjoin, graphjoin, insert, inversions, selfdup, seqfix
from .small import place_small
from .state import Run


def curate_allhaps(pt, phasing="hic", hap_regex=r"^(H\d+)_", hi=None, do_cuts=True, blocks=None, busco_hits=None,
                   small_blocks=None, do_place=True, level=0.5, do_arrange=True, same_hap_blocks=None, sex_system=None,
                   graph_junctions=None, unitig_paths=None, end_joins=None, stop_after=None, fasta=None, threads=4,
                   minimap2=None):
    """allhaps curation: cuts (`telomere_sites`), homolog pairing, small scaffolds
    (`place_small`), merged pairs, order and orientation within each copy (end contact, then
    the other haplotype's alignment), parallel layout, HAP tags, paint, savestate.
    `pt` is the MQ0 map; `hi` the high-MQ map of the same assembly, if there is one; `blocks`
    the H2-against-H1 alignment and `small_blocks` small scaffolds against everything (both
    optional); `busco_hits` per-gene loci for HAP1 choice and false duplicates. `do_cuts`,
    `do_place`, `do_arrange` switch steps off. `level` (0 conservative .. 1 aggressive, see
    `common.THRESHOLDS`) moves every act-or-flag threshold together. `same_hap_blocks`: alignments of
    sex-linked scaffolds within one haplotype (see sexchrom), for their PARs; `sex_system`
    'XY' or 'ZW' (species knowledge) adds X/Y or Z/W tags to the sex chromosome copies.
    `graph_junctions`: hifiasm graph evidence at YaHS gaps (`graph.load`): those joins are
    right, so never cut, inverted or rearranged there. `unitig_paths`: unitig-graph paths
    across YaHS gaps (`graph.unitig_paths`), weaker evidence (see `_unitig_evidence`).
    `end_joins`: joins between scaffold ends only the graph sees (`graph.scaffold_end_joins`):
    made (`graphjoin`), and evidence for merging pairs.
    `fasta`: the assembly; contig-end overlaps within a haplotype and gaps with a copy on both
    sides are then decided by the sequence (`seqfix`), and info["exact_breaks"] holds the rows
    for breaks.tsv (cuts placed to the base, gapless joins).
    Returns (rows for the decisions table, savestate, info). `stop_after` (a step name in
    STEPS, for tests) returns the run's state (`state.Run`) after that step instead."""
    raw = pt
    th = Thresholds(level)                    # every level-dependent threshold (common.THRESHOLDS)
    run_ = Run(pt=pt, phasing=phasing, hap_regex=hap_regex, hi=hi, do_cuts=do_cuts, blocks=blocks,
               busco_hits=busco_hits, small_blocks=small_blocks, do_place=do_place, level=level,
               do_arrange=do_arrange, same_hap_blocks=same_hap_blocks, sex_system=sex_system,
               graph_junctions=graph_junctions, unitig_paths=unitig_paths, end_joins=end_joins or [], raw=raw,
               th=th, fasta=fasta, threads=threads, minimap2=minimap2)
    if stop_after is not None and stop_after not in STEPS:
        raise ValueError(f"stop_after: unknown step {stop_after!r} (steps: {', '.join(STEPS)})")
    for name in STEPS[:-1]:
        globals()[name](run_)
        if name == stop_after:
            return run_
    return _layout(run_)


# the steps of curate_allhaps, in order; each reads and sets attributes of the Run (state.py)
STEPS = ("_setup", "_small", "_merge", "_sex_chain", "_arrange", "_reference", "_follow", "_graph_place", "_rows",
         "_layout")


def _setup(run_):
    """Sex chromosomes, a first pairing, break sites (cut / invert_tip / keep) with the graph
    evidence, the cut map and the pairing on it."""
    from .. import sexchrom
    blocks, do_cuts, graph_junctions = run_.blocks, run_.do_cuts, run_.graph_junctions
    hap_regex, hi, level, raw = run_.hap_regex, run_.hi, run_.level, run_.raw
    same_hap_blocks, small_blocks, th = run_.same_hap_blocks, run_.small_blocks, run_.th
    unitig_paths = run_.unitig_paths
    # hemizygous stretches, PARs, sex chromosomes joined to another chromosome (sexchrom)
    sx = sexchrom.analyse(raw, blocks, same_hap_blocks, level=level)
    margin = th.pair_margin
    # 1. a first pairing on the uncut map, so that sites can be mirrored to the homolog; then
    # cut, and pair again on the pieces
    first_pairs, _ = pair_homologs(raw, hap_regex, blocks=blocks, margin=margin)
    sites = telomere_sites(raw, first_pairs, hi, hap_regex, blocks, gap_max=th.gap_max,
                           sharp_max=th.sharp_max, tip_max=th.tip_max) if do_cuts else []
    graph_rows = _graph_evidence(raw, sites, graph_junctions or []) + _unitig_evidence(raw, sites, unitig_paths or [])
    # two copies of one region in a scaffold: cut at a gap only when they are haplotypic
    dups = selfdup.self_duplications(raw, blocks) if blocks is not None and do_cuts else []
    sites = selfdup.apply_cuts(sites, dups)
    graph_rows += selfdup.rows(raw, dups)
    # loose scaffolds that fill a hole inside their own haplotype's copy: cut it at the gap there
    if do_cuts:
        sites, ins_rows = insert.insert_cuts(raw, first_pairs, blocks, hap_regex, sites)
        graph_rows += ins_rows
    cuts = [(x["q"], x["s"], x["reason"]) for x in sites if x["action"] == "cut"]
    # contig-end overlaps within a haplotype and gaps with a copy on both sides (seqfix)
    run_.seq_overlaps, run_.exact_breaks = [], []
    if run_.fasta and do_cuts:
        more, seq_rows, exact, applied = seqfix.find(raw, [q for q, _, _ in cuts], run_.fasta, hap_regex,
                                                     run_.threads, run_.minimap2)
        cuts += more
        graph_rows += seq_rows
        run_.seq_overlaps, run_.exact_breaks = applied, exact
    # inversions inside a scaffold: turned in place in the layout (inversions.py)
    run_.inv_flips = []
    if do_cuts:
        run_.inv_flips, inv_rows = inversions.find(raw, [q for q, _, _ in cuts], blocks, th)
        graph_rows += inv_rows
    flips = [x["tip"] for x in sites if x["action"] == "invert_tip"]
    pt = SegmentedMap(raw, [q for q, _, _ in cuts])          # from here on, units are cut pieces
    pairs, info = pair_homologs(pt, hap_regex, blocks=blocks, margin=margin)
    if small_blocks is not None:
        pairs += same_prefix_pairs(pt, info, small_blocks)
    run_.cuts, run_.flips, run_.graph_rows, run_.info, run_.pairs = cuts, flips, graph_rows, info, pairs
    run_.pt, run_.sites, run_.sx = pt, sites, sx
    run_.cut_out = sorted(q for d in dups if d["cut"] for q in d["cut"])       # gaps a copy was cut out at


def _small(run_):
    """Small scaffolds (place_small) and the state the arrangement steps share."""
    blocks, do_place, hap_regex, hi, info = run_.blocks, run_.do_place, run_.hap_regex, run_.hi, run_.info
    pairs, phasing, pt, raw, sites = run_.pairs, run_.phasing, run_.pt, run_.raw, run_.sites
    small_blocks, th = run_.small_blocks, run_.th
    # 2. small scaffolds: which chromosome and copy (unloc / candidate / extends a copy)
    small_calls = place_small(pt, pairs, info, small_blocks, hi, phasing=phasing, seq_min=th.seq_min,
                              hic_agree=th.hic_agree, hic_alone=th.hic_alone,
                              hi_unique=th.hi_unique, mq0_unique=th.mq0_unique) if do_place else []
    attach = {}
    for c in small_calls:
        if c["action"] == "unloc":
            attach.setdefault(c["copy"], []).append(c)
    hap = info["hap"]
    # Order and orientation within each chromosome copy: its pieces between the sites that
    # were kept or tip-inverted, plus attached pieces large enough to carry an end signal.
    npx_, start_ = pt.scaffold_npix, pt.scaffold_start
    # end contact on the high-MQ map: at MQ0 the homolog's reads blur which ends touch
    sg_arr = _signal(hi) if hi is not None else _signal(raw)
    # kept / inverted: split points only (not where the hifiasm graph links the two sides)
    soft = sorted(x["q"] for x in sites if x["action"] != "cut" and not x.get("graph") and not x.get("hold"))
    tips = {tuple(x["tip"]) for x in sites if x["action"] == "invert_tip"}
    chain_px, gain_min = th.chain_px, th.gain_min
    attach_link, attach_margin, orient_gain_ = th.attach_link, th.attach_margin, th.orient_gain
    mate_of = {}
    for p in pairs:
        if p.kind != "extra":
            mate_of[p.a], mate_of[p.b] = p.b, p.a
    extra_of = {}
    for p in pairs:
        if p.kind == "extra" and p.a in mate_of:
            # an extra piece homologous to copy p.a belongs in the chain of p.a's mate
            extra_of.setdefault(mate_of[p.a], []).append(p.b)
    # arranged: copy -> layout [(first px, last px, reversed)]; joined: unit -> copy whose chain
    # it is in; arr_rows: decision rows; sol: per copy, its pieces and end-contact solution
    arranged, joined, arr_rows, sol = {}, {}, [], {}
    block_idx = _block_index(blocks) if blocks is not None else None
    if block_idx is not None and small_blocks is not None:
        # unique alignments of small scaffolds to the other haplotype, which the H2-against-H1
        # alignment can miss (it reports one best target per stretch of the query)
        pre = lambda n: (re.match(hap_regex, n).group(1) if re.match(hap_regex, n) else "")
        have = {}
        for b in blocks:                                    # what the main alignment already says, per pair
            have.setdefault((b.q, b.t), []).append((b.qs, b.qe))
            have.setdefault((b.t, b.q), []).append((b.ts, b.te))
        # seen: at least half of this block is already covered by the main alignment
        seen = lambda b: _union([(max(b.qs, x), min(b.qe, y)) for x, y in have.get((b.q, b.t), [])
                                 if min(b.qe, y) > max(b.qs, x)]) >= 0.5 * (b.qe - b.qs)
        more = [b for b in small_blocks if b.primary and b.mapq >= 30 and pre(b.q) != pre(b.t) and not seen(b)]
        for name, v in _block_index(more).items():
            block_idx.setdefault(name, []).extend(v)
    run_.arr_rows, run_.arranged, run_.attach, run_.attach_link = arr_rows, arranged, attach, attach_link
    run_.attach_margin, run_.block_idx, run_.chain_px = attach_margin, block_idx, chain_px
    run_.extra_of, run_.gain_min, run_.hap, run_.joined = extra_of, gain_min, hap, joined
    run_.mate_of, run_.npx_, run_.orient_gain_, run_.sg_arr = mate_of, npx_, orient_gain_, sg_arr
    run_.small_calls, run_.soft, run_.sol, run_.start_, run_.tips = small_calls, soft, sol, start_, tips


def _merge(run_):
    """Pairs that are parts of one chromosome merged into one (merge_pairs)."""
    arr_rows, attach, do_arrange, extra_of = run_.arr_rows, run_.attach, run_.do_arrange, run_.extra_of
    hap, npx_, pairs, pt, sg_arr, th = run_.hap, run_.npx_, run_.pairs, run_.pt, run_.sg_arr, run_.th
    glinks = graphjoin.links(pt, run_.raw, run_.end_joins)
    # Pairs that are parts of one chromosome become one: the smaller pair's two units go into
    # the chains of the larger pair's two copies (same haplotype prefix on each side).
    host_of = {}                                           # merged unit -> host copy
    if do_arrange:
        size = lambda P: npx_[P.a] + npx_[P.b]
        root = {}
        find = lambda P: find(root[id(P)]) if id(P) in root else P    # union-find on pairs

        def merge(P, Q, score, text):
            """Q's two units go into the chains of P's copies (the larger pair hosts)."""
            P, Q = find(P), find(Q)
            if P is Q:
                return
            if size(Q) > size(P):
                P, Q = Q, P
            root[id(Q)] = P
            Q.kind = "merged"
            side = {hap[P.a]: P.a, hap[P.b]: P.b}
            for u_ in (Q.a, Q.b):
                host_of[u_] = side[hap[u_]]
            for old_u, new_u in list(host_of.items()):       # units merged into Q earlier follow it
                if new_u in (Q.a, Q.b):
                    host_of[old_u] = side[hap[new_u]]
            arr_rows.append(dict(
                action="merge", scaffold=pt.names[Q.a], length_mb=round(npx_[Q.a] * pt.bp_per_pix / 1e6, 2),
                partner=pt.names[Q.b], strand="", tag="", score=score, next_best="",
                evidence=f"this pair and the pair {pt.names[P.a]} / {pt.names[P.b]} are parts of one chromosome: {text}; "
                         "laid out and painted as one chromosome"))
        for v, P, Q, x, y, how in merge_pairs(pt, sg_arr, pairs, hap, link_min=th.merge_link_min,
                                              specific=th.merge_specific, link_weak=th.merge_link_weak,
                                              block_min=th.merge_block_min, graph_links=glinks):
            merge(P, Q, round(v, 2), f"their ends touch at {v:.2f} of adjacent sequence in both haplotypes (the weaker "
                                     f"of the two){how}")
        # chromosome-level contact alone (arms apart across centromeric repeats): chromjoin
        extra = {}
        for u_, h_ in host_of.items():
            extra.setdefault(h_, []).append(u_)
        for c_, lst in attach.items():
            extra.setdefault(c_, []).extend(x["u"] for x in lst)
        for c_, lst in extra_of.items():                    # further pieces bound for the copy's chain
            extra.setdefault(c_, []).extend(lst)
        for d in chromjoin.candidates(pt, pairs, hap, extra, margin=th.chrom_margin, spread_min=th.chrom_spread,
                                      min_mb=th.chrom_min_mb):
            text = (f"chromosome-level contact (repeat bins masked) between {pt.names[d['x']]} and {pt.names[d['y']]} is "
                    f"{d['ratio']:.1f}x their background, {d['ratio'] / max(d['second'], 1e-9):.1f}x either one's next "
                    f"partner, spread over {100 * d['spread']:.0f}% of the block; between their homologs "
                    f"{d['ratio_other']:.1f}x, {100 * d['spread_other']:.0f}%; their ends need not touch (chromosome arms "
                    "apart across centromeric repeats)")
            if th.chrom_act >= 0.5:
                merge(d["P"], d["Q"], round(d["ratio"], 1), text)
            else:
                arr_rows.append(dict(
                    action="merge_candidate", scaffold=pt.names[d["y"]], length_mb=round(npx_[d["y"]] * pt.bp_per_pix / 1e6, 2),
                    partner=pt.names[d["x"]], strand="", tag="", score=round(d["ratio"], 1), next_best="",
                    evidence=text + "; not merged at this level (review)"))
        # unlocs and extra pieces of a merged unit go with its host copy
        for u_, h_ in host_of.items():
            if attach.get(u_):
                attach.setdefault(h_, []).extend(attach.pop(u_))
            if extra_of.get(u_):
                extra_of.setdefault(h_, []).extend(extra_of.pop(u_))
    members = {}
    for u_, h_ in host_of.items():
        members.setdefault(h_, []).append(u_)
    run_.host_of, run_.members = host_of, members


def _sex_chain(run_):
    """Hemizygous pieces in no pair chained to the sex chromosome copy they contact."""
    from .. import sexchrom
    attach, chain_px, do_arrange, extra_of = run_.attach, run_.chain_px, run_.do_arrange, run_.extra_of
    hap, host_of, info, mate_of, members = run_.hap, run_.host_of, run_.info, run_.mate_of, run_.members
    npx_, phasing, pt, raw, start_, sx = run_.npx_, run_.phasing, run_.pt, run_.raw, run_.start_, run_.sx
    th = run_.th
    # Hemizygous pieces in no pair (a sex chromosome in several scaffolds has no homolog to
    # guide it): each goes into the chain of the sex chromosome copy it has chromosome-level
    # contact with, if that copy clearly beats every other chromosome copy; phase 1 then
    # places it by end contact (the copy's arrangement is forced). In Troglobisium the X is
    # H1_scaffold_1 25-119.5 Mb with H2_scaffold_46 (curator). Under trio phasing, only
    # within a haplotype.
    sexfrac = lambda u: float((sx["mask"][int(start_[u]):int(start_[u] + npx_[u])] > 0).mean()) if npx_[u] else 0.0
    sex_forced, sex_rows = set(), []
    if do_arrange and sx["stretches"]:
        attached = {c["u"]: k for k, lst in attach.items() for c in lst}
        groups = {c: [c] + members.get(c, []) + [x["u"] for x in attach.get(c, [])] for c in mate_of if c not in host_of}
        rng_ = lambda v: (int(start_[v]), int(start_[v] + npx_[v] - 1))
        for v in sorted(info["unpaired"], key=lambda v: -npx_[v]):
            if v in mate_of or v in host_of or v in attached or npx_[v] < chain_px or sexfrac(v) < 0.5:
                continue
            score = {}
            for c, g in groups.items():
                if c == v or (phasing == "trio" and hap[c] != hap[v]):
                    continue
                vals = [sexchrom._contact(raw, *rng_(v), *rng_(x)) for x in g if x != v]
                vals = [x for x in vals if not np.isnan(x)]
                if vals:
                    score[c] = max(vals)
            if not score:
                continue
            ranked = sorted(score.items(), key=lambda x: -x[1])
            c, best = ranked[0]
            second = ranked[1][1] if len(ranked) > 1 else 0.0
            if best < th.sex_chain_min or best < th.sex_chain_ratio * second \
                    or max(sexfrac(x) for x in groups[c]) < 0.5:
                continue
            extra_of.setdefault(c, []).append(v)
            sex_forced.add(c)
            groups[c].append(v)
            sex_rows.append(dict(
                action="sex_chromosome_join", scaffold=pt.names[v], length_mb=round(npx_[v] * pt.bp_per_pix / 1e6, 2),
                partner=pt.names[c], strand="", tag="", score=round(best, 3), next_best=round(second, 3),
                evidence=f"hemizygous (coverage {float(np.median(sx['norm'][rng_(v)[0]:rng_(v)[1] + 1])):.2f} of a normal "
                         f"copy, {100 * sexfrac(v):.0f}% sex-linked) and in no pair; chromosome-level contact {best:.3f} "
                         f"with the sex chromosome copy {pt.names[c]}, against {second:.3f} with any other: put in its "
                         "chain by end contact (see the rearrange row)"))
    span = lambda a, b: (f"{raw.names[raw.orig_scaffold[a]]} {raw.rel_coord[a] * raw.bp_per_pix / 1e6:.1f}-"
                         f"{(raw.rel_coord[b] + 1) * raw.bp_per_pix / 1e6:.1f} Mb")
    run_.sex_forced, run_.sex_rows, run_.sexfrac, run_.span = sex_forced, sex_rows, sexfrac, span


def _arrange(run_):
    """Phase 1: each copy arranged by end contact; large attached pieces placed in it."""
    arr_rows, arranged, attach, attach_link = run_.arr_rows, run_.arranged, run_.attach, run_.attach_link
    attach_margin, block_idx, chain_px = run_.attach_margin, run_.block_idx, run_.chain_px
    do_arrange, extra_of, gain_min, host_of = run_.do_arrange, run_.extra_of, run_.gain_min, run_.host_of
    joined, mate_of, members, npx_ = run_.joined, run_.mate_of, run_.members, run_.npx_
    orient_gain_, pt, raw, sex_forced = run_.orient_gain_, run_.pt, run_.raw, run_.sex_forced
    sg_arr, soft, sol, span, start_, th = run_.sg_arr, run_.soft, run_.sol, run_.span, run_.start_, run_.th
    tips, cut_out = run_.tips, run_.cut_out
    # 3. phase 1: each copy (with the units merged into it) arranged by end contact alone
    for u in (sorted(mate_of, key=lambda u: start_[u]) if do_arrange else []):
        if u in host_of:
            continue
        units = [u] + sorted(members.get(u, []), key=lambda v: -npx_[v])

        def split(v):
            """Unit v as pieces split at its kept / tip-inverted sites (gaps where the sides
            belong to one chromosome but the order is not certain)."""
            a, b = int(start_[v]), int(start_[v] + npx_[v] - 1)
            out, lo = [], a
            for q in [q for q in soft if a <= q < b]:
                if q - lo >= 4 and b - q > 4:              # no pieces of a few pixels
                    out.append((lo, q))
                    lo = q + 1
            return out + [(lo, b)]
        pieces = [p for v in units for p in split(v)]
        others = sorted([c["u"] for c in attach.get(u, []) if npx_[c["u"]] >= chain_px] +
                        [v for v in extra_of.get(u, []) if npx_[v] >= chain_px], key=lambda v: -npx_[v])
        # Extra pieces are split at their kept sites as well, so that one can go into a gap of
        # another (Troglobisium: the curator put H2_scaffold_46 0-8.59 Mb, hemizygous, into the
        # X's gap at H1_scaffold_1 109.78 Mb). The search is exact, so at most 9 pieces; an
        # extra that does not fit split goes in whole, or not at all.
        kept, n_own, block_of = [], len(pieces), {}
        for v in others:
            sp = split(v)
            if len(pieces) + len(sp) > 9:
                sp = [(int(start_[v]), int(start_[v] + npx_[v] - 1))]
            if len(pieces) + len(sp) > 9:
                break
            block_of[v] = list(range(len(pieces), len(pieces) + len(sp)))
            pieces += sp
            kept.append(v)
        others = kept
        # a merged chromosome is laid out as one in any case, and so is a sex chromosome with
        # hemizygous pieces put in its chain
        forced = len(units) > 1 or (u in sex_forced and any(v in others for v in extra_of.get(u, [])))
        sol[u] = dict(pieces=pieces, init=[(i, pieces[i] in tips) for i in range(len(pieces))], units=units,
                      others=others, forced=forced, adj=None, row=None, best=0.0, was=0.0)
        if len(pieces) < 2:
            continue
        init = [(i, pieces[i] in tips) for i in range(len(pieces))]    # as assembled, tips inverted
        # the other copy's pixel ranges (its units, unlocs and extras): where guide() looks
        mates = [mate_of[v] for v in units]
        cols = [(int(start_[v]), int(start_[v] + npx_[v] - 1)) for m_ in mates for v in
                [m_] + [c["u"] for c in attach.get(m_, [])] + extra_of.get(m_, [])]
        memo = {}

        def guide(i, pieces=pieces, cols=cols, memo=memo):
            """Direction of piece i in the other copy: from the alignments if there are any
            (and they are clear), else from the Hi-C homology line."""
            if i not in memo:
                d = _aligned_dir(raw, block_idx, pieces[i], cols) if block_idx is not None else 0
                if d:
                    by_seq.add(i)
                memo[i] = d or _homolog_dir(_signal(raw), pieces[i][0], pieces[i][1], cols)
            return memo[i]
        by_seq = set()
        arr, best, was, adj, unclear, guided = arrange_pieces(sg_arr, pieces, init, orient_gain=th.orient_gain,
                                                      guide=guide)
        sol[u].update(adj=adj, best=best, was=was)
        placed_only = None
        # Pieces of one scaffold that are apart only because a haplotypic copy was cut out
        # between them (selfdup) are one chromosome copy: joined again as assembled.
        same = [v for v in others if v in extra_of.get(u, []) and raw.orig_scaffold[start_[v]] == raw.orig_scaffold[start_[u]]
                and any(min(start_[u], start_[v]) <= q < max(start_[u], start_[v]) for q in cut_out)]
        if same and not forced and (arr == init or best - was < gain_min):
            own = sorted([p for i, p in enumerate(pieces) if i < n_own or any(i in block_of[v] for v in same)])
            arranged[u] = [(a, b, False) for a, b in own]
            for v in same:
                joined[v] = u
            sol[u]["others"], sol[u]["row"] = same, len(arr_rows)
            arr_rows.append(dict(
                action="rearrange", scaffold=pt.names[u], length_mb=round(npx_[u] * pt.bp_per_pix / 1e6, 2), partner="",
                strand="", tag="", score=round(best, 2), next_best=round(was, 2),
                evidence="pieces of one scaffold joined again as assembled (" + ", ".join(span(a, b) for a, b in own)
                         + "): a haplotypic copy was cut out between them (see the self_duplication row)"))
            continue
        if not forced and (arr == init or best - was < gain_min):
            # Not enough gain to rearrange the copy itself. Its attached pieces large enough
            # to carry an end signal still belong in it, not at its end as unlocs (curator's
            # Troglobisium curation joins all four of them, 2-10 Mb): each goes, largest first,
            # where its ends touch the copy's (left as assembled) best, if one of its two new
            # neighbours touches it at `attach_link` of adjacent sequence or more.
            ins = [v for v in others if v in block_of and v not in extra_of.get(u, [])]
            if not ins:
                continue                           # nothing to place: left as assembled
            tot = lambda a_: sum(adj[(x, rx, y, ry)] for (x, rx), (y, ry) in zip(a_, a_[1:]))
            cur, notes = list(init[:n_own]), []
            for v in ins:
                blk = [init[i] for i in block_of[v]]
                cands = []                         # (score, position, flipped, arrangement, best new link)
                for pos in range(len(cur) + 1):
                    for flip in (False, True):
                        b = [(i, not r) for i, r in reversed(blk)] if flip else blk
                        cand = cur[:pos] + b + cur[pos:]
                        left = adj[(cur[pos - 1][0], cur[pos - 1][1], b[0][0], b[0][1])] if pos > 0 else 0.0
                        right = adj[(b[-1][0], b[-1][1], cur[pos][0], cur[pos][1])] if pos < len(cur) else 0.0
                        cands.append((tot(cand), pos, flip, cand, max(left, right)))
                cands.sort(key=lambda c: -c[0])
                sc, pos, flip, cand, link = cands[0]
                alt = next((c[0] for c in cands if c[1] != pos), 0.0)        # best elsewhere
                same = next(c for c in cands if c[1] == pos and c[2] != flip)  # same place, turned
                if link < attach_link:
                    notes.append(f"{pt.names[v]} left as an unloc (best end contact {link:.2f})")
                    continue
                if sc - alt < attach_margin:
                    notes.append(f"{pt.names[v]} left as an unloc (where it goes is too close to call: "
                                 f"{sc:.2f} against {alt:.2f} elsewhere)")
                    continue
                turned = ""
                if flip and sc - same[0] < orient_gain_:
                    cand = same[3]                 # orientation too close to call: as assembled
                    turned = f", orientation too close to call ({sc:.2f} / {same[0]:.2f}), as assembled"
                cur = cand
                notes.append(f"{pt.names[v]} placed (end contact {link:.2f}, {sc - alt:.2f} better than "
                             f"elsewhere{turned})")
            if len(cur) == n_own:
                continue
            arr, placed_only = cur, notes
            others = [v for v in others if any(i in dict(arr) for i in block_of.get(v, []))]
            sol[u]["others"] = others
            best = tot(arr)
        sol[u]["row"] = len(arr_rows)
        arranged[u] = [(pieces[i][0], pieces[i][1], r) for i, r in arr]
        for v in others + units[1:]:
            joined[v] = u
        links = [adj[(x, rx, y, ry)] for (x, rx), (y, ry) in zip(arr, arr[1:])]
        arr_rows.append(dict(
            action="rearrange", scaffold=pt.names[u], length_mb=round(npx_[u] * pt.bp_per_pix / 1e6, 2), partner="",
            strand="", tag="", score=round(best, 2), next_best=round(was, 2),
            evidence="new order: " + ", ".join(span(pieces[i][0], pieces[i][1]) + (" reversed" if r else "") for i, r in arr)
                     + "; contact between neighbouring ends " + "/".join(f"{x:.2f}" for x in links)
                     + (f" (sum {best:.2f}) against {was:.2f} as assembled" if placed_only is None else
                        "; the copy's own pieces as assembled, its attached pieces put where their ends touch it "
                        "best: " + ", ".join(placed_only))
                     + ("; oriented by the other haplotype ("
                        + ("sequence alignment" if any(i in by_seq for i in guided) else "Hi-C homology line") + "): "
                        + ", ".join(span(*pieces[i]) for i in guided) if guided and placed_only is None else "")
                     + ("; orientation too close to call, left as it was, for: "
                        + ", ".join(span(*pieces[i]) for i in unclear) if unclear and placed_only is None else "")))


def _reference(run_):
    """With alignments: thresholds and helpers for following the reference copy; loose
    scaffolds put into chains."""
    arr_rows, block_idx, blocks, host_of = run_.arr_rows, run_.block_idx, run_.blocks, run_.host_of
    joined, mate_of, npx_, pairs, pt = run_.joined, run_.mate_of, run_.npx_, run_.pairs, run_.pt
    raw, sg_arr, small_calls, sol, start_ = run_.raw, run_.sg_arr, run_.small_calls, run_.sol, run_.start_
    th = run_.th
    # With sequence alignments the other haplotype is the guide for order and orientation.
    # Of the two copies of a chromosome one is the reference (arranged by end contact above);
    # the pieces of the other are put in the order and orientation in which they align along
    # it. The reference is the copy whose arrangement the other copy's contigs contradict
    # less (curator: contigs outrank scaffold joins; `_contig_conflicts`), else the copy in
    # fewer pieces. End contact is the check: if it clearly prefers another arrangement, that
    # may be a real difference between the haplotypes, so the contact arrangement stays and
    # the row says so.
    pair_of = {}
    for p in pairs:
        if p.kind not in ("extra", "merged"):
            pair_of[p.a] = pair_of[p.b] = p
    # how much more end contact the contact arrangement must score than the colinear one for
    # it to stay (higher = follow the alignment more readily)
    conflict_min = th.conflict_min
    # contig-led reference: bp of one copy's contigs contradicting the other's arrangement
    # needed to make that copy the reference, and by how much it must exceed the reverse
    contig_min, contig_ratio = th.contig_min, 2.0
    # Scaffolds that are in no pair but align, over at least half their length, to one copy of
    # a chromosome: they belong in the chain of the other copy (their own haplotype), at the
    # place where they align.
    loose, false_dup, inside = {}, {}, {}
    if block_idx is not None:
        ridx = {n: i for i, n in enumerate(raw.names)}

        def unit_at(name, bp):
            """Unit (cut piece) holding position bp of an original scaffold; None if absent."""
            i = ridx.get(name)
            if i is None or raw.scaffold_npix[i] == 0:
                return None
            return int(pt.orig_scaffold[int(raw.scaffold_start[i]) + min(int(bp / raw.bp_per_pix), int(raw.scaffold_npix[i]) - 1)])
        acc = {}
        for bk in blocks:
            if bk.mapq < 30:
                continue
            r_ = (bk.te - bk.ts) / max(bk.qe - bk.qs, 1)         # target bp per query bp
            for x, y, bp in split_block(bk, unit_at, raw.bp_per_pix):
                if x is None or y is None or x == y:
                    continue
                acc.setdefault(x, {}).setdefault(y, 0)
                acc[x][y] += bp
                acc.setdefault(y, {}).setdefault(x, 0)
                acc[y][x] += bp * r_
        min_loose = th.min_loose                         # smallest loose scaffold, px
        for v, d in acc.items():
            if v in mate_of or v in host_of or npx_[v] < min_loose:
                continue
            y = max(d, key=d.get)                          # its best target, mapped to the host copy
            c = host_of.get(y, y)
            if d[y] < 0.5 * npx_[v] * raw.bp_per_pix or c not in sol or c not in pair_of:
                continue
            loose.setdefault(pair_of[c].b if pair_of[c].a == c else pair_of[c].a, []).append(v)
    g230 = sg_arr.decay_long(230)
    seen_small = {c["u"] for c in small_calls}
    rng_of = lambda v: (int(start_[v]), int(start_[v] + npx_[v] - 1))

    def along(copy, ref_layout, p, base, base_arr):
        """Pieces `base` of a copy plus its loose scaffolds, in the order and orientation in
        which they align along the other copy's layout. Returns (pieces, arrangement, units of
        the added pieces) or None if a base piece does not align."""
        extra_u = [v for v in loose.get(copy, []) if rng_of(v) not in base and joined.get(v, copy) == copy]
        P = list(base) + [rng_of(v) for v in extra_u]
        proj = [_project(raw, block_idx, x, ref_layout) for x in P]
        if any(proj[i] is None for i in range(len(base))):
            return None
        # A loose scaffold that aligns where its own copy already has sequence is a second copy of
        # that sequence (haplotig or false duplication), not a missing piece: it stays out.
        have = sorted(x for i in range(len(base)) for x in proj[i][2])

        def covered(segs):
            """Share of the layout stretches segs already covered by the base pieces' alignments."""
            tot = cov = 0.0
            for l, h in segs:
                tot += h - l
                cov += _union([(max(l, a), min(h, b)) for a, b in have if min(h, b) > max(l, a)])
            return cov / max(tot, 1e-9)
        keep = [i for i in range(len(P)) if proj[i] is not None and (i < len(base) or covered(proj[i][2]) < 0.5)]
        for i in range(len(base), len(P)):
            if proj[i] is not None and i not in keep:
                false_dup[extra_u[i - len(base)]] = (copy, covered(proj[i][2]))
        # A loose scaffold whose stretch lies inside the span of one base piece cannot go where it
        # aligns without cutting that piece (insert.insert_cuts cuts at a gap there when it can);
        # put next to an end it would be far from its place (icLarUrsu SUPER_8): it stays out.
        def extent(segs, q=0.05):
            """Layout span of the central part (1 - 2q of the aligned length) of the stretches segs."""
            s = sorted(segs)
            w = np.array([h - l for l, h in s], dtype=float)
            if not len(s) or w.sum() <= 0:
                return None
            c = np.cumsum(w) / w.sum()
            return s[int(np.searchsorted(c, q))][0], s[min(int(np.searchsorted(c, 1 - q)), len(s) - 1)][1]
        spans = [extent(proj[b][2]) for b in range(len(base))]
        for i in [i for i in keep if i >= len(base)]:
            e = extent(proj[i][2])
            hit = next((b for b, sp in enumerate(spans) if e and sp and sp[0] + 2 < e[0] and e[1] < sp[1] - 2), None)
            if hit is not None:
                keep.remove(i)
                inside[extra_u[i - len(base)]] = (copy, P[hit])
        # sorted by position along the other copy, reversed where the alignment runs backwards
        arr = [(i, proj[i][1] < 0) for i in sorted(keep, key=lambda i: proj[i][0])]
        if not p.same_strand and p.kind == "pair":         # the second copy is shown reversed as a whole
            arr = [(i, not r) for i, r in reversed(arr)]
        return P, arr, {i: extra_u[i - len(base)] for i in keep if i >= len(base)}

    def join_rows(added, copy, ref):
        """Record loose scaffolds `added` as joined into `copy`, with a 'join' row each (small
        scaffolds have their own row from place_small)."""
        for v in added.values():
            joined[v] = copy
            if v in seen_small:
                # its own row (place_small) may say Hi-C alone left it out; this is sequence evidence
                for c in small_calls:
                    if c["u"] == v:
                        c["join_note"] = (f"but most of it aligns to {pt.names[ref]} (the other haplotype): placed in the "
                                          f"chain of {pt.names[copy]} where it aligns (sequence alignment; see the "
                                          "rearrange row)")
            else:
                arr_rows.append(dict(action="join", scaffold=pt.names[v], length_mb=round(npx_[v] * pt.bp_per_pix / 1e6, 2),
                                     partner=pt.names[copy], strand="", tag="", score="", next_best="",
                                     evidence=f"in no pair, but most of it aligns to {pt.names[ref]} (the other haplotype): placed in "
                                              f"the chain of {pt.names[copy]} where it aligns (see the rearrange row)"))

    run_.along, run_.conflict_min, run_.contig_min = along, conflict_min, contig_min
    run_.contig_ratio, run_.false_dup, run_.g230, run_.join_rows = contig_ratio, false_dup, g230, join_rows
    run_.loose, run_.pair_of, run_.inside = loose, pair_of, inside


def _follow(run_):
    """Phase 2, once per pair: the other copy follows the reference's layout (_project)."""
    along, arr_rows, arranged, attach = run_.along, run_.arr_rows, run_.arranged, run_.attach
    block_idx, conflict_min, contig_min = run_.block_idx, run_.conflict_min, run_.contig_min
    contig_ratio, false_dup, g230, hi = run_.contig_ratio, run_.false_dup, run_.g230, run_.hi
    host_of, join_rows, joined, loose, npx_ = run_.host_of, run_.join_rows, run_.joined, run_.loose, run_.npx_
    pair_of, pt, raw, sg_arr, small_calls = run_.pair_of, run_.pt, run_.raw, run_.sg_arr, run_.small_calls
    sol, span, start_, tips = run_.sol, run_.span, run_.start_, run_.tips
    # 4. phase 2, once per pair: the other copy follows the reference's layout (`_project`)
    for u in (sorted(sol, key=lambda u: start_[u]) if block_idx is not None else []):
        p = pair_of.get(u)
        if p is None or u != p.a or p.b not in sol:
            continue
        size = lambda v: sum(b - a + 1 for a, b in sol[v]["pieces"])
        # reference: the copy whose arrangement the other copy's contigs contradict less
        # (contigs outrank scaffold joins, curator 2026-11-14); when neither or both are
        # contradicted, the copy in fewer pieces, then the longer one
        ref, oth = sorted((p.a, p.b), key=lambda v: (len(sol[v]["pieces"]), -size(v)))
        lay = lambda v: arranged.get(v) or [(a, b, r) for (a, b), (_, r) in zip(sol[v]["pieces"], sol[v]["init"])]
        # x_ref: bp of oth's contigs that contradict ref's layout; x_oth the other way round
        x_ref, w_ref = _contig_conflicts(raw, block_idx, sol[oth]["pieces"], lay(ref))
        x_oth, w_oth = _contig_conflicts(raw, block_idx, sol[ref]["pieces"], lay(oth))
        ref_note, row0 = "", len(arr_rows)
        fmt = lambda w: "; ".join(f"{c} across {', '.join(s)}" for c, s in w[:3]) + (" ..." if len(w) > 3 else "")
        if x_ref >= contig_min and x_ref > contig_ratio * x_oth:
            # oth's contigs say ref's pieces are out of order or turned: oth leads
            ref_note = (f"; reference copy {pt.names[oth]}, not {pt.names[ref]} (fewer pieces): contigs of "
                        f"{pt.names[oth]} contradict the arrangement of {pt.names[ref]} over {x_ref / 1e6:.2f} Mb "
                        f"({fmt(w_ref)}), the other way round over {x_oth / 1e6:.2f} Mb; contigs outrank scaffold joins")
            ref, oth = oth, ref
        elif x_ref >= contig_min and x_oth >= contig_min:
            # contigs of each copy contradict the other: a real difference or errors in both
            arr_rows.append(dict(action="arrangement_differs", scaffold=pt.names[ref],
                                 length_mb=round(npx_[ref] * pt.bp_per_pix / 1e6, 2), partner=pt.names[oth], strand="",
                                 tag="", score=round(x_oth / 1e6, 2), next_best=round(x_ref / 1e6, 2),
                                 evidence=f"contigs of each copy contradict the other's arrangement ({pt.names[oth]} "
                                          f"against {pt.names[ref]}: {x_ref / 1e6:.2f} Mb, {fmt(w_ref)}; {pt.names[ref]} "
                                          f"against {pt.names[oth]}: {x_oth / 1e6:.2f} Mb, {fmt(w_oth)}); "
                                          f"{pt.names[ref]} (fewer pieces) kept as reference, to check"))
            row0 += 1
        elif x_ref >= contig_min:
            ref_note = (f"; contigs of {pt.names[oth]} contradict the arrangement of {pt.names[ref]} over "
                        f"{x_ref / 1e6:.2f} Mb ({fmt(w_ref)}) but its own contigs contradict {pt.names[oth]} as much: "
                        f"{pt.names[ref]} (fewer pieces) kept as reference, to check")
        elif x_oth >= contig_min:
            # the usual case: the reference's contigs show where the other copy's joins are wrong
            ref_note = (f"; contigs of {pt.names[ref]} contradict the arrangement of {pt.names[oth]} over "
                        f"{x_oth / 1e6:.2f} Mb ({fmt(w_oth)}), and not the other way round: {pt.names[ref]} leads")
        so, sr = sol[oth], sol[ref]
        n0 = len(so["pieces"])                 # oth's own pieces; loose scaffolds are added after
        layout = arranged.get(ref) or [(a, b, r) for (a, b), (_, r) in zip(sr["pieces"], sr["init"])]
        got = along(oth, layout, p, so["pieces"], None)
        if got is not None and len(got[0]) >= 2:
            P, aln, added = got
            # link: end contact of two pieces, from phase 1's table when both are oth's own
            link = lambda a, ra, b, rb: so["adj"][(a, ra, b, rb)] if so["adj"] is not None and a < n0 and b < n0 else \
                _adjacency(sg_arr, g230, P[a], ra, P[b], rb)
            # t_aln: end contact of oth's own pieces in the colinear order (comparable to so["best"]);
            # cur: oth's arrangement after phase 1
            base_aln = [x for x in aln if x[0] < n0]
            t_aln = sum(so["adj"][(a, ra, b, rb)] for (a, ra), (b, rb) in zip(base_aln, base_aln[1:])) if so["adj"] is not None else 0.0
            cur = [(so["pieces"].index((a, b)), r) for a, b, r in arranged[oth]] if oth in arranged else so["init"]
            order_txt = ", ".join(span(P[i][0], P[i][1]) + (" reversed" if r else "") for i, r in aln)
            tip_undone = hi is not None and any(P[i] in tips and not r for i, r in aln)
            # outcomes: colinear already (note it) / would undo a tip inversion (flag, keep) /
            # contact clearly prefers another arrangement (flag, keep) / follow the alignment
            if aln == cur:
                if so["row"] is not None:
                    arr_rows[so["row"]]["evidence"] += "; the same order and orientation as in the other haplotype (sequence alignment)"
            elif tip_undone and so["row"] is None:
                # an end piece was inverted on the evidence of the high-MQ map and the telomere position;
                # the alignment would put it back. Both copies may differ there (or carry the same
                # error): the inversion stays, flagged.
                arr_rows.append(dict(action="arrangement_differs", scaffold=pt.names[oth],
                                     length_mb=round(npx_[oth] * pt.bp_per_pix / 1e6, 2), partner=pt.names[ref], strand="",
                                     tag="", score=round(so["best"], 2), next_best=round(t_aln, 2),
                                     evidence=f"an end piece was inverted (see the invert_tip row), but as assembled it is "
                                              f"colinear with {pt.names[ref]} (sequence alignment): the haplotypes may differ "
                                              "here, or both carry the same arrangement; inversion kept, to check"))
            elif so["best"] - t_aln > conflict_min:
                note = (f"; the other haplotype (sequence alignment along {pt.names[ref]}) has a different arrangement: {order_txt} "
                        f"(end contact {t_aln:.2f}); end contact clearly prefers the one "
                        + ("applied" if so["row"] is not None else "assembled")
                        + ": possibly a real difference between the haplotypes, to check")
                if so["row"] is not None:
                    arr_rows[so["row"]]["evidence"] += note
                else:
                    arr_rows.append(dict(action="arrangement_differs", scaffold=pt.names[oth],
                                         length_mb=round(npx_[oth] * pt.bp_per_pix / 1e6, 2), partner=pt.names[ref], strand="",
                                         tag="", score=round(so["best"], 2), next_best=round(t_aln, 2), evidence=note[2:]))
            else:
                if so["row"] is not None:              # phase 1's row is replaced (None, dropped below)
                    arr_rows[so["row"]] = None
                    so["row"] = None
                if aln == so["init"] and not so["forced"]:
                    arranged.pop(oth, None)                # the assembled arrangement is the colinear one
                    for v in so["others"]:
                        joined.pop(v, None)
                else:
                    arranged[oth] = [(P[i][0], P[i][1], r) for i, r in aln]
                    for v in so["others"] + so["units"][1:]:
                        joined[v] = oth
                    join_rows(added, oth, ref)
                    links = "/".join(f"{link(a, ra, b, rb):.2f}" for (a, ra), (b, rb) in zip(aln, aln[1:]))
                    so["row"] = len(arr_rows)
                    arr_rows.append(dict(
                        action="rearrange", scaffold=pt.names[oth], length_mb=round(npx_[oth] * pt.bp_per_pix / 1e6, 2),
                        partner=pt.names[ref], strand="", tag="", score=round(t_aln, 2), next_best=round(so["was"], 2),
                        evidence="new order: " + order_txt + f"; order and orientation follow the other haplotype (sequence "
                                 f"alignment along {pt.names[ref]}); contact between neighbouring ends {links}; its own pieces "
                                 f"score {t_aln:.2f} this way, {so['was']:.2f} as assembled, {so['best']:.2f} at best by end "
                                 "contact alone"))
        # the reference copy's own loose scaffolds go in where they align along the other copy,
        # provided that leaves the reference's pieces as they are
        if loose.get(ref):
            lay_o = arranged.get(oth) or [(a, b, r) for (a, b), (_, r) in zip(so["pieces"], so["init"])]
            base = [(a, b) for a, b, _ in layout]
            got = along(ref, lay_o, p, base, None)
            if got is not None and got[2]:
                P, aln, added = got
                # applied only if the reference's own pieces keep their order and orientation
                if [x for x in aln if x[0] < len(base)] == [(i, r) for i, (_, _, r) in enumerate(layout)]:
                    arranged[ref] = [(P[i][0], P[i][1], r) for i, r in aln]
                    for v in sr["others"] + sr["units"][1:]:
                        joined.setdefault(v, ref)
                    join_rows(added, ref, oth)
                    arr_rows.append(dict(
                        action="rearrange", scaffold=pt.names[ref], length_mb=round(npx_[ref] * pt.bp_per_pix / 1e6, 2),
                        partner=pt.names[oth], strand="", tag="", score="", next_best="",
                        evidence="new order: " + ", ".join(span(P[i][0], P[i][1]) + (" reversed" if r else "") for i, r in aln)
                                 + f"; scaffolds that were in no pair are placed where they align along {pt.names[oth]} "
                                   "(sequence alignment); the pieces that were already there are unchanged"))
        if ref_note:                       # say why this copy led, on every row it produced
            idxs = set(range(row0, len(arr_rows))) | ({so["row"]} if so["row"] is not None else set())
            for k in sorted(idxs):
                if arr_rows[k] is not None:
                    arr_rows[k]["evidence"] += ref_note
    arr_rows = [r for r in arr_rows if r is not None]
    # Second copies found above: not joined, kept at the end of their copy, tagged FalseDuplicate
    # (curator, 2026-11-14: sequence that is not in the contig is kept out, whatever made it).
    by_unit = {c["u"]: c for c in small_calls}
    for v, (copy, frac) in false_dup.items():
        if v in joined:
            continue
        why = (f"{100 * frac:.0f}% of the stretch of the other haplotype it aligns to is already covered by "
               f"{pt.names[copy]}, its own haplotype's copy: a second copy of that sequence, kept out of the chain")
        if v in by_unit:
            by_unit[v]["third_copy"] = True
            by_unit[v]["evidence"] += "; " + why
        else:
            c = dict(u=v, copy=copy, action="unloc", third_copy=True, evidence=why, score="", ratio="")
            small_calls.append(c)
            attach.setdefault(copy, []).append(c)
    if run_.seq_overlaps:                    # copies cut off at a contig-end overlap (seqfix)
        arr_rows += seqfix.tag_copies(run_)
    # Loose scaffolds whose place is inside one piece of their copy (no gap there to cut at): left
    # out, with a row saying where they belong.
    for v, (copy, piece) in run_.inside.items():
        if v in joined:
            continue
        arr_rows.append(dict(action="insert_candidate", scaffold=pt.names[v], length_mb=round(npx_[v] * pt.bp_per_pix / 1e6, 2),
                             partner=pt.names[copy], strand="", tag="", score="", next_best="",
                             evidence=f"in no pair; it aligns to the other haplotype inside the stretch that {span(*piece)} "
                                      f"covers, where {pt.names[copy]} has no sequence of its own: a missing piece of it, "
                                      "but that piece would have to be cut to put it in (no gap there): not moved"))
    # Pieces that continue a copy beyond one of its ends (see place_small): put next to that end,
    # turned so that the overlapping ends meet.
    for c in small_calls:
        if not c.get("extend") or c["action"] != "unloc":
            continue
        v, (X, x_end, u_end) = c["u"], c["extend"]
        H = host_of.get(X, X)
        if H in arranged:
            cur = list(arranged[H])
        elif H in sol:
            cur = [(a, b, r) for (a, b), (_, r) in zip(sol[H]["pieces"], sol[H]["init"])]
        else:
            cur = [(int(start_[H]), int(start_[H] + npx_[H] - 1), False)]
        cur = [t for t in cur if t[0] != int(start_[v])]
        end_px = int(start_[X]) if x_end == "h" else int(start_[X] + npx_[X] - 1)
        k_ = next((i for i, (a, b, _) in enumerate(cur) if a <= end_px <= b), None)
        if k_ is None:
            continue
        end_left = (x_end == "h") != cur[k_][2]            # that end of X is on the left of its piece as laid out
        vr = (u_end == "h") if end_left else (u_end == "t")
        cur.insert(k_ if end_left else k_ + 1, (int(start_[v]), int(start_[v] + npx_[v] - 1), vr))
        arranged[H] = cur
        joined[v] = H
        c["action"] = "join"
        c["join_note"] = f"placed next to that end of {pt.names[X]}" + (" (reversed)" if vr else "")
        row = dict(action="rearrange", scaffold=pt.names[H], length_mb=round(npx_[H] * pt.bp_per_pix / 1e6, 2),
                   partner="", strand="", tag="", score="", next_best="",
                   evidence="new order: " + ", ".join(span(a, b) + (" reversed" if r else "") for a, b, r in cur)
                            + f"; {pt.names[v]} continues {pt.names[X]} beyond its {'start' if x_end == 'h' else 'end'} "
                              "(its non-repetitive end overlaps it; sequence alignment)")
        old_rows = [i for i, r in enumerate(arr_rows) if r["action"] == "rearrange" and r["scaffold"] == pt.names[H]]
        if old_rows:
            row["evidence"] += "; replaces: " + arr_rows[old_rows[0]]["evidence"]
            arr_rows[old_rows[0]] = row
        else:
            arr_rows.append(row)
    # a small scaffold placed in a chain is a join, no longer an unloc at the end of its copy
    for c in small_calls:
        if c["u"] in joined:
            c["action"] = "join"
    for u in attach:
        attach[u] = [c for c in attach[u] if c["u"] not in joined]
    run_.arr_rows = arr_rows


def _graph_place(run_):
    """Joins only hifiasm's graph sees, between scaffold ends: the smaller scaffold goes at that
    end of the larger one's chromosome (`graphjoin.place`)."""
    run_.arr_rows += graphjoin.place(run_)


def _rows(run_):
    """Decision rows that do not depend on the layout (small scaffolds, BUSCO, cuts, sites, sex
    chromosomes) and the BUSCO genes per unit."""
    arr_rows, blocks, busco_hits, cuts = run_.arr_rows, run_.blocks, run_.busco_hits, run_.cuts
    graph_rows, hap_regex, pt, raw = run_.graph_rows, run_.hap_regex, run_.pt, run_.raw
    sex_rows, sites, small_calls, sx = run_.sex_rows, run_.sites, run_.small_calls, run_.sx
    # 5. decision rows that do not depend on the layout
    npx, start = pt.scaffold_npix, pt.scaffold_start
    rows = []
    for c in small_calls:
        act = "join" if c["action"] == "join" else "false_duplicate" if c["third_copy"] else \
            "unloc" if c["action"] == "unloc" else "unloc_candidate"
        rows.append(dict(action=act, scaffold=pt.names[c["u"]],
                         length_mb=round(npx[c["u"]] * pt.bp_per_pix / 1e6, 2), partner=pt.names[c["copy"]], strand="",
                         tag={"false_duplicate": "FalseDuplicate", "unloc": "Unloc"}.get(act, ""), score=c["score"], next_best="",
                         evidence=c["evidence"] + ("; " + c["join_note"] if c.get("join_note") else
                                                   "; large enough to place by end contact: joined into the chromosome "
                                                   "(see the rearrange row)" if c["action"] == "join" else "")))
    rows += arr_rows + sex_rows
    mb = lambda s: npx[s] * pt.bp_per_pix / 1e6
    # BUSCO genes per unit (complete copies, and how many of them are duplicated in the unit)
    genes, dup_of = {}, {}
    if busco_hits is not None:
        from ..busco import summarise as busco_summary
        ridx = {n: i for i, n in enumerate(raw.names)}
        seen = {}
        for gene, loci in busco_hits.items():
            for sc, s0, e0, _, _ in loci:
                i = ridx.get(sc)
                if i is None or raw.scaffold_npix[i] == 0:
                    continue
                px = int(raw.scaffold_start[i]) + min(int((s0 + e0) / 2 / raw.bp_per_pix), int(raw.scaffold_npix[i]) - 1)
                u = int(pt.orig_scaffold[px])
                seen.setdefault(u, {}).setdefault(gene, 0)
                seen[u][gene] += 1
        genes = {u: (len(d), sum(v > 1 for v in d.values())) for u, d in seen.items()}
        hap_name = lambda n: (re.match(hap_regex, n).group(1) if re.match(hap_regex, n) else "")
        per_hap, _, _ = busco_summary(busco_hits, hap_name)
        from ..busco import false_duplicates
        for sc, sc2, a0, b0, shared, n_genes in false_duplicates(busco_hits, hap_name):
            if sc in ridx and raw.scaffold_npix[ridx[sc]] > 0:
                dup_of[ridx[sc]] = sc2
            rows.append(dict(action="false_duplicate_candidate", scaffold=sc,
                             length_mb=round(raw.scaffold_len_bp()[ridx[sc]] / 1e6, 2) if sc in ridx else "",
                             partner=sc2, strand="", tag="", score=shared, next_best="",
                             evidence=f"{shared} of its {n_genes} BUSCO genes have a second copy in the same haplotype on "
                                      f"{sc2} at {a0 / 1e6:.2f}-{b0 / 1e6:.2f} Mb: likely a false duplication (or haplotig) "
                                      "of that region; highlighted with vertPaint, not tagged"))
        total = len(busco_hits)
        for h, v in sorted(per_hap.items()):
            rows.append(dict(action="busco_summary", scaffold=h, length_mb="", partner="", strand="", tag="", score="",
                             next_best="", evidence=f"of {total} BUSCO genes found anywhere: complete {v['complete']} "
                                                    f"({100 * v['complete'] / max(total, 1):.1f}%), single {v['single']}, "
                                                    f"duplicated {v['duplicated']}, missing {total - v['complete']}"))
    # Sex chromosomes (sexchrom): where a sex-linked stretch meets other sequence, the coverage
    # steps; a cut there gets the reason, an uncut boundary that is not a PAR is flagged.
    bpp = raw.bp_per_pix
    mbr = lambda a, b: f"{raw.rel_coord[a] * bpp / 1e6:.2f}-{(raw.rel_coord[b] + 1) * bpp / 1e6:.2f} Mb"
    near_px = max(5, int(1e6 / bpp))                       # a cut this close to a coverage step is at it
    end_at = lambda q: next((e for e in sx["ends"] if e["s"] == raw.orig_scaffold[q] and abs(e["boundary"] - q) <= near_px), None)
    step_txt = lambda e: (f"coverage steps from {e['cov_sex']:.2f} (sex-linked, hemizygous) to {e['cov_other']:.2f} of a "
                          f"normal copy at {e['edge_bp'] / 1e6:.2f} Mb")
    for q, s, why in cuts:
        e = end_at(q)
        if e is not None and e["kind"] in ("joined", "fused"):
            why += f"; {step_txt(e)}: a sex chromosome joined to another chromosome"
        rows.append(dict(action="cut", scaffold=raw.names[s], length_mb=round(raw.rel_coord[q] * raw.bp_per_pix / 1e6, 2),
                         partner="", strand="", tag="", score="", next_best="",
                         evidence=f"cut at {raw.rel_coord[q] * raw.bp_per_pix / 1e6:.2f} Mb; " + why))
    rows += graph_rows
    for x in sites:
        if x["action"] != "cut" and not x.get("graph") and not x.get("quiet"):   # graph-backed ones: the graph_* row says it
            rows.append(dict(action=x["action"], scaffold=raw.names[x["s"]],
                             length_mb=round(raw.rel_coord[x["q"]] * raw.bp_per_pix / 1e6, 2), partner="", strand="",
                             tag="", score="", next_best="",
                             evidence=f"end piece {(x['tip'][1] - x['tip'][0] + 1) * raw.bp_per_pix / 1e6:.1f} Mb; " + x["reason"]))
    for d in sx["stretches"]:
        L = (d["b"] - d["a"] + 1) * bpp
        if d["via"]:
            o = d["via"]
            what = (f"its sequence aligns into the sex-linked stretch of {raw.names[o[0]]} ({mbr(o[1], o[2])}): the "
                    "homologous (gametologous) part of the other sex chromosome; a PAR would be at normal coverage")
        else:
            what = (f"{d['hemi_bp'] / 1e6:.1f} Mb of it at half coverage with no alignment to the other haplotype"
                    + ("" if blocks is not None else " (no alignments given: from coverage alone)")
                    + f", {d['aln_bp'] / 1e6:.1f} Mb aligned: the sex-specific part of a sex chromosome")
        rows.append(dict(action="sex_linked_region", scaffold=raw.names[d["s"]], length_mb=round(L / 1e6, 2), partner="",
                         strand="", tag="", score=round(d["cov"], 2), next_best="",
                         evidence=f"{mbr(d['a'], d['b'])}: coverage {d['cov']:.2f} of a normal copy; {what}"))
    for e in sx["ends"]:
        cut_here = any(abs(q - e["boundary"]) <= near_px for q, _, _ in cuts)
        where = f"{mbr(e['a'], e['b'])}, next to the sex-linked stretch, which ends at {e['edge_bp'] / 1e6:.2f} Mb"
        ct = "unknown" if np.isnan(e["contact"]) else f"{e['contact']:.3f}"
        if e["kind"] == "par":
            p0, p1 = e["par"]
            ev = (f"pseudoautosomal region {p0 / 1e6:.2f}-{p1 / 1e6:.2f} Mb ({where}): normal coverage "
                  f"({e['cov_other']:.2f}), chromosome-level contact with the sex-linked stretch ({ct}), and "
                  + (f"it aligns to the normal-coverage end of the other sex chromosome, {e['partner']}"
                     if e["by"] == "alignment" else
                     f"MQ0 contact with the normal-coverage end of {e['partner']} as strong as within each "
                     f"({e.get('homolog_contact', float('nan')):.2f}; reads of a PAR map to both copies)")
                  + "; recombines like an autosome: not a haplotig or false duplication")
            if cut_here:
                ev += "; a cut was made at its boundary (see the cut row): to check"
            rows.append(dict(action="par", scaffold=raw.names[e["s"]], length_mb=round((p1 - p0) / 1e6, 2),
                             partner=e["partner"], strand="", tag="", score=ct, next_best="", evidence=ev))
        elif e["kind"] == "par_candidate":
            rows.append(dict(action="par_candidate", scaffold=raw.names[e["s"]],
                             length_mb=round((e["b"] - e["a"] + 1) * bpp / 1e6, 2), partner="", strand="", tag="",
                             score=ct, next_best="",
                             evidence=f"{where}: normal coverage ({e['cov_other']:.2f}) with chromosome-level contact "
                                      f"({ct}): a PAR, or an autosome fused to the sex chromosome; "
                                      + ("not confirmed: it does not align to another sex chromosome (if both are in one "
                                         "haplotype, give --fasta so they are aligned to each other)" if e["have_aln"] else
                                         "no alignments to tell (give --fasta or --paf)")))
        elif e["kind"] == "fused" or not cut_here:
            ev = (f"{where}: {step_txt(e)}; "
                  + (f"chromosome-level contact {ct}, and it aligns to {e['other']} ({e['other_bp'] / 1e6:.1f} Mb), not to "
                     "another sex chromosome: an autosome fused to the sex chromosome (a neo-sex chromosome), or a misjoin"
                     if e["kind"] == "fused" else
                     f"no chromosome-level contact across ({ct}): another chromosome joined to the sex chromosome; not cut "
                     "(no site found here)"))
            rows.append(dict(action="sex_boundary", scaffold=raw.names[e["s"]],
                             length_mb=round(raw.rel_coord[e["boundary"]] * bpp / 1e6, 2), partner=e["other"] or "",
                             strand="", tag="", score=ct, next_best="", evidence=ev))

    run_.bpp, run_.dup_of, run_.genes, run_.mb, run_.npx, run_.rows = bpp, dup_of, genes, mb, npx, rows
    run_.start = start


def _layout(run_):
    """Layout order (pairs side by side, HAP tags, sex chromosomes, unpaired), paint, tags and
    the savestate; returns (rows, savestate, info)."""
    from .. import sexchrom
    arranged, attach, blocks, bpp, cuts = run_.arranged, run_.attach, run_.blocks, run_.bpp, run_.cuts
    dup_of, flips, genes, hap, host_of = run_.dup_of, run_.flips, run_.genes, run_.hap, run_.host_of
    info, joined, mb, npx, pairs = run_.info, run_.joined, run_.mb, run_.npx, run_.pairs
    phasing, pt, raw, rows, same_hap_blocks = run_.phasing, run_.pt, run_.raw, run_.rows, run_.same_hap_blocks
    sex_system, sexfrac, sites, small_calls = run_.sex_system, run_.sexfrac, run_.sites, run_.small_calls
    start, sx = run_.start, run_.sx
    # 6. layout: pairs side by side, HAP tags, sex chromosomes, unpaired; paint, tags, savestate
    unit_h = "% of the shorter aligned" if blocks is not None else ""
    haps = sorted({h for h in hap if h})
    cov = pt.graphs.get("coverage")
    main = [p for p in pairs if p.kind not in ("extra", "merged")]
    extra = {}
    for p in pairs:
        if p.kind == "extra":
            extra.setdefault(host_of.get(p.a, p.a), []).append(p)
    paired = np.concatenate([np.arange(start[s], start[s] + npx[s]) for p in main for s in (p.a, p.b)]) \
        if main else np.array([], dtype=int)
    # coverage of one copy: the median over paired units
    norm = float(np.median(cov[paired])) if cov is not None and len(paired) else 0.0
    # order: (unit, reversed) in final layout order; tag_of / paint_of: HAP tag and paint id
    order, tag_of, paint_of = [], {}, {}
    unloc_units, sex_like, in_chain = set(), set(), set()
    # sex chromosome copies laid out (main unit of each chain), how they relate (paired, or
    # aligned bp), and the X/Y/Z/W tag per unit (only with --sex-system)
    sex_copies, sex_links, sex_tag_of = [], {}, {}
    look_units = {c["u"] for c in small_calls if c["action"] == "candidate" and not c["third_copy"]}
    fd_units = {c["u"] for c in small_calls if c["third_copy"] and c["action"] != "join"}
    inv_flips = getattr(run_, "inv_flips", [])

    def turn_inversions(px):
        """Each inversion (inversions.py) reversed in place where its pixels lie together."""
        for a, b in inv_flips:
            idx = np.nonzero((px >= a) & (px <= b))[0]
            if len(idx) == b - a + 1 and idx[-1] - idx[0] == b - a:
                px[idx[0]:idx[-1] + 1] = px[idx[0]:idx[-1] + 1][::-1].copy()
        return px

    def pixels(s, rev=False):
        """Pixels of a unit in layout order, with any inverted end piece turned in place."""
        if s in arranged:
            px = np.concatenate([np.arange(a, b + 1)[::-1] if r else np.arange(a, b + 1) for a, b, r in arranged[s]])
            px = turn_inversions(px)
            return px[::-1] if rev else px
        px = np.arange(start[s], start[s] + npx[s])
        for a, b in flips:
            if start[s] <= a and b < start[s] + npx[s]:
                px[a - start[s]:b - start[s] + 1] = px[a - start[s]:b - start[s] + 1][::-1].copy()
        px = turn_inversions(px)
        return px[::-1] if rev else px
    # pairs largest first, the two copies side by side; paint ids 2k+1 / 2k+2
    for k, p in enumerate(sorted(main, key=lambda p: -(npx[p.a] + npx[p.b]))):
        qa, ga, ta = _quality(pt, p.a, info["tel_thr"])
        qb, gb, tb = _quality(pt, p.b, info["tel_thr"])
        first, second = p.a, p.b
        # both copies at half coverage, or both mostly in sex-linked stretches (sexchrom)
        hemi = (norm and all(float(np.median(cov[start[u]:start[u] + npx[u]])) / norm < 0.7 for u in (p.a, p.b))) \
            or all(sexfrac(u) >= 0.5 for u in (p.a, p.b))
        if p.kind == "unequal" and not hemi:
            p.kind, p.note = "pair", "lengths differ: one copy is probably fragmented or incomplete"
        if p.kind == "pair" and hemi and min(mb(p.a), mb(p.b)) >= 5:
            # both copies at about half coverage: hemizygous, i.e. X with Y (or Z with W) of similar size
            p.kind, p.note = "unequal", "both copies have about half the coverage of other pairs (hemizygous)"
            sex_like.add(id(p))
        if p.kind == "unequal":
            # convention: both sex chromosomes go to HAP1; the larger one (with the pieces in its
            # chain) first
            chain = lambda c: npx[c] + sum(npx[v] for v, u_ in joined.items() if u_ == c)
            first, second = (p.a, p.b) if chain(p.a) >= chain(p.b) else (p.b, p.a)
            tags = ("HAP1", "HAP1")
            why = "possible sex chromosome pair: both tagged HAP1 by convention"
            if id(p) in sex_like or all(sexfrac(u) >= 0.5 for u in (p.a, p.b)):
                sex_copies += [first, second]
                sex_links[frozenset((first, second))] = 1.0
        elif phasing == "trio":
            first, second = sorted((p.a, p.b), key=lambda s: haps.index(hap[s]))
            tags = ("HAP" + str(haps.index(hap[first]) + 1), "HAP" + str(haps.index(hap[second]) + 1))
            why = "trio phasing: tags follow the name prefix"
        else:
            # hic / none: whole superscaffolds may swap; the better copy becomes HAP1.
            # BUSCO decides only when the copies clearly differ (more single-copy genes, fewer
            # duplicated); otherwise fewer gaps, telomeric ends, length.
            ba, bb = genes.get(p.a, (0, 0)), genes.get(p.b, (0, 0))
            sa, sb = ba[0] - 2 * ba[1], bb[0] - 2 * bb[1]      # complete genes, duplicated ones penalised
            if genes and abs(sa - sb) > max(2, 0.02 * max(ba[0], bb[0])):
                if sb > sa:
                    first, second = p.b, p.a
                why = "copy with more complete, non-duplicated BUSCO genes is HAP1"
            else:
                if qb > qa:
                    first, second = p.b, p.a
                why = "better copy (fewer gaps, telomeres, length) is HAP1"
            if genes:
                why += f"; BUSCO genes {genes.get(first, (0, 0))[0]}/{genes.get(second, (0, 0))[0]}"
            tags = ("HAP1", "HAP2")
        # an unequal pair shares only a small region (e.g. a PAR): too little to orient by
        flip = (not p.same_strand) and (p.kind == "pair" or id(p) in sex_like)
        paint_of[first], paint_of[second] = 2 * k + 1, 2 * k + 2
        tag_of[first], tag_of[second] = tags
        for copy, rev in ((first, False), (second, flip)):
            order.append((copy, rev))
            for v, u_ in joined.items():
                if u_ == copy:
                    paint_of[v], tag_of[v] = paint_of[copy], tag_of[copy]
                    in_chain.add(v)
            # unlocs: at the end of their copy of the chromosome, painted and tagged with it
            for c in sorted(attach.get(copy, []), key=lambda c: -npx[c["u"]]):
                order.append((c["u"], False))
                paint_of[c["u"]], tag_of[c["u"]] = paint_of[copy], tag_of[copy]
                unloc_units.add(c["u"])
        for e in sorted(extra.get(p.a, []) + extra.get(p.b, []), key=lambda e: -npx[e.b]):
            if e.b in joined or e.b in fd_units:           # in a chain, or laid out as a false duplicate
                tag_of.setdefault(e.b, tag_of[second] if hap[e.b] == hap[second] else tag_of[first])
                continue
            order.append((e.b, not e.same_strand))
            tag_of[e.b] = tag_of[second] if hap[e.b] == hap[second] else tag_of[first]
            rows.append(dict(action="place_next_to_homolog", scaffold=pt.names[e.b], length_mb=round(mb(e.b), 2),
                             partner=pt.names[e.a], strand="+" if e.same_strand else "-", tag=tag_of[e.b],
                             score=round(e.score, 1), next_best=round(e.second, 1),
                             evidence=e.note + "; not painted, position along the homolog not resolved"))
        rows.append(dict(action="pair" if p.kind == "pair" else "pair_unequal", scaffold=pt.names[first],
                         length_mb=round(mb(first), 2), partner=pt.names[second],
                         strand="-" if flip else "+", tag="/".join(tags), score=round(p.score, 1),
                         next_best=round(p.second, 1),
                         evidence=f"homology ({info['source']}) {p.score:.1f}{unit_h} against next best {p.second:.1f}; "
                                  f"partner {mb(second):.1f} Mb; gaps {ga if first == p.a else gb}/{gb if first == p.a else ga}, "
                                  f"telomeric ends {ta if first == p.a else tb}/{tb if first == p.a else ta}; {why}"
                                  + ("; " + p.note if p.note else "")))
    # Hemizygous sex chromosomes: no homolog of their own size, and about half the coverage
    # of a paired copy. Candidates are put together, larger first, both HAP1 (convention).
    covr = lambda s: float(np.median(cov[start[s]:start[s] + npx[s]])) / norm if norm else float("nan")
    done = {x for x, _ in order} | set(joined)
    sex = [s for s in info["unpaired"] if s not in done and mb(s) >= 5 and (covr(s) < 0.7 or sexfrac(s) >= 0.5)]
    k = len(main)
    for s in sorted(sex, key=lambda s: -npx[s]):
        partner = [t for t in sex if t != s and hap[t] != hap[s] and (s, t) in info["line"]]
        partner = max(partner, key=lambda t: info["line"][(s, t)]) if partner else None
        rev = False
        if partner is not None and partner in {x for x, _ in order}:
            rev = not info["strand"][(s, partner)]
        order.append((s, rev))
        tag_of[s] = "HAP1"
        sex_copies.append(s)
        if partner is not None:
            sex_links[frozenset((s, partner))] = 1.0
        paint_of[s] = 2 * k + 1
        k += 1
        rows.append(dict(action="sex_chromosome_candidate", scaffold=pt.names[s], length_mb=round(mb(s), 2),
                         partner=pt.names[partner] if partner is not None else "", strand="-" if rev else "+",
                         tag="HAP1", score=round(info["line"][(s, partner)], 1) if partner is not None else "",
                         next_best="", evidence=f"no homolog of similar size; coverage {covr(s):.2f} of a paired copy "
                                                "(hemizygous); tagged HAP1 by convention"))
    # Sex chromosome labels (only with --sex-system: which copy is which needs species
    # knowledge). Copies are related when paired or when they align to each other (summed
    # over the pieces in their chains); sexchrom.sex_labels applies the system.
    chain_of = lambda c: [c] + [v for v, u_ in joined.items() if u_ == c] + [x["u"] for x in attach.get(c, [])]
    sex_idx = _block_index(list(blocks or []) + list(same_hap_blocks or [])) if (blocks or same_hap_blocks) else None
    if sex_idx and len(sex_copies) > 1:
        for i, c1 in enumerate(sex_copies):
            for c2 in sex_copies[i + 1:]:
                bp = sum(_aligned_bp(raw, sex_idx, (int(start[x]), int(start[x] + npx[x] - 1)),
                                     (int(start[y]), int(start[y] + npx[y] - 1))) for x in chain_of(c1) for y in chain_of(c2))
                if bp >= 1e6:
                    sex_links[frozenset((c1, c2))] = sex_links.get(frozenset((c1, c2)), 0.0) + bp
    sizes = {c: sum(int(npx[x]) for x in chain_of(c)) * raw.bp_per_pix for c in sex_copies}
    labels, label_note = sexchrom.sex_labels(sex_system, sex_copies, sizes, sex_links)
    for c, t in labels.items():
        for x in chain_of(c):
            sex_tag_of[x] = t
    if sx["norm"] is not None:
        par_mb = sum((e["par"][1] - e["par"][0]) / 1e6 for e in sx["ends"] if e["kind"] == "par")
        if sx["stretches"]:
            ev = (f"{len(sx['stretches'])} sex-linked stretches ({sum((d['b'] - d['a'] + 1) for d in sx['stretches']) * bpp / 1e6:.1f} Mb), "
                  f"{sum(e['kind'] == 'par' for e in sx['ends'])} PARs ({par_mb:.1f} Mb); sex chromosome copies: "
                  + (", ".join(f"{pt.names[c]} ({sizes[c] / 1e6:.1f} Mb{', ' + labels[c] if c in labels else ''})"
                               for c in sex_copies) or "none laid out as such"))
        else:
            ev = ("no hemizygous stretch found. That does not mean there are no sex chromosomes: the sample may be of "
                  "the homogametic sex (XX, ZZ), or the sex chromosomes are homomorphic; either way they look like "
                  "autosomes here and are paired as such")
        rows.append(dict(action="sex_summary", scaffold="", length_mb="", partner="", strand="", tag="", score="",
                         next_best="", evidence=ev + (f"; {label_note}" if label_note else
                                                      "" if sex_system else "; not labelled X/Y or Z/W (give --sex-system)")))
    placed = {s for s, _ in order} | set(joined)
    for s in list(info["unpaired"]) + [u for u in info.get("repetitive", {}) if u not in info["unpaired"]]:
        if s in placed:
            continue
        c = float(np.median(cov[start[s]:start[s] + npx[s]])) / norm if norm else float("nan")
        guess = ("coverage about twice a single copy: possibly a collapsed chromosome (both haplotypes in one)"
                 if c > 1.6 else "single-copy coverage: possibly a sex chromosome, or its homolog is fragmented")
        if mb(s) < 5:
            guess = "small; no single clear homolog at this resolution"
        if s in info.get("repetitive", {}):
            guess = (f"too repetitive to assign: coverage is well above a normal copy (collapsed repeats); its best "
                     f"match is {pt.names[info['repetitive'][s]]}; left unplaced")
        rows.append(dict(action="unpaired", scaffold=pt.names[s], length_mb=round(mb(s), 2), partner="", strand="",
                         tag="", score="", next_best="", evidence=f"coverage {c:.2f} of paired scaffolds; {guess}"))
    # everything else after the pairs, in assembly order
    rest = [int(s) for s in np.argsort(start) if npx[s] > 0 and int(s) not in placed and int(s) not in joined]
    perm = np.concatenate([pixels(s, rev) for s, rev in order] + [pixels(s) for s in rest])
    st = SaveState.blank()
    where = np.empty(raw.npix, dtype=np.int64)                 # original pixel -> layout position
    where[perm] = np.arange(raw.npix)
    # a cut whose two sides stay adjacent in the layout needs an explicit break edit; the
    # moves that realise the permutation come from plan_edits
    st.edits = [Edit(q, q, BREAK_DELTA) for q, _, _ in cuts if where[q + 1] - where[q] == 1] + plan_edits(perm)
    ms = MapState(raw).replay(st.edits)
    if not np.array_equal(ms.pix, perm):
        raise AssertionError("planned edits do not reproduce the target layout")
    # a waypoint at every site (PretextView holds at most 255)
    for k, q in enumerate([x["q"] for x in sites][:255]):
        x = where[q] / raw.npix - 0.5
        st.waypoints.append(Waypoint(x, x, 30.0, k))
    # paint and tags per contig of the final map (a unit's tag, else its prefix's HAP number)
    for i, c in enumerate(ms.contigs):
        s = int(pt.orig_scaffold[ms.pix[c.start]])
        if s in paint_of:
            st.scaffs.append((i, paint_of[s]))
        t = tag_of.get(s) or (("HAP" + str(haps.index(hap[s]) + 1)) if hap[s] and haps.index(hap[s]) < 2 else "")
        mask = (1 << st.tag_id(t)) if t else 0
        if s in fd_units:
            mask |= 1 << st.tag_id("FalseDuplicate")
        elif s in unloc_units:
            mask |= 1 << st.tag_id("Unloc")
        if c.scaffold in dup_of or s in look_units:
            mask |= 1 << st.tag_id("vertPaint")
        if s in sex_tag_of:
            mask |= 1 << st.tag_id(sex_tag_of[s])
        if mask:
            st.meta_flags.append((i, mask))
    info["sex"] = sx
    info["exact_breaks"] = run_.exact_breaks
    return rows, st, info
