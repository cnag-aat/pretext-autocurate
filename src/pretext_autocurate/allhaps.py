"""allhaps mode, step 1: pair homologous scaffolds of different haplotypes and lay each
pair side by side in the same orientation, tagged HAP1 / HAP2.

Homology is read from the MQ0 map: reads from sequence shared by two haplotypes are split
between the copies, so a scaffold shows a line of contact with its homolog, along the
diagonal of their block if both are assembled in the same direction and along the
anti-diagonal if not. Unrelated scaffolds show only background.

Not done here yet: cuts, joins, unlocs, haplotigs (see docs/HANDOFF.md).
"""
import re
from dataclasses import dataclass

import numpy as np

from .mapstate import MapState, plan_edits
from .pretextio import SegmentedMap
from .savestate import SaveState, Edit, Waypoint, BREAK_DELTA
from .signal import Signal


@dataclass
class Pair:
    a: int                # scaffold index, placed first (HAP1)
    b: int                # its homolog (HAP2)
    score: float          # mean of the homology line (map units at the working level)
    second: float         # best line of either scaffold with any other scaffold
    same_strand: bool
    kind: str             # pair | unequal (possible sex chromosomes) | extra
    note: str


def _level(pt, max_side=8192):
    lv = 0
    while (pt.npix >> lv) > max_side and lv < pt.n_mips - 1:
        lv += 1
    return lv


def _quality(pt, s, tel_thr):
    """Higher is better: fewer gaps, telomeres at both ends, then length."""
    a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
    gap, tel = pt.graphs.get("gap"), pt.graphs.get("telomere")
    gaps = int((gap[a:a + n] > 0).sum()) if gap is not None else 0
    k = max(1, min(10, n // 4))
    tels = 0 if tel is None else int(tel[a:a + k].max() >= tel_thr) + int(tel[a + n - k:a + n].max() >= tel_thr)
    return (-gaps, tels, n), gaps, tels


def pair_homologs(pt, hap_regex=r"^(H\d+)_", min_bins=4, score_min=None, margin=1.3, unequal=0.5, blocks=None,
                  min_seq_bp=500000, repeat_cov=1.3):
    """Returns (pairs, info). `pt` must be the MQ0 map.

    With `blocks` (sequence alignments between haplotypes, see homology.py) homology is
    taken from them: the score of two scaffolds is the percentage of the shorter one that
    aligns to the other, ranking is by aligned bases, strand by the majority of aligned
    bases. Without blocks it is estimated from the MQ0 Hi-C signal."""
    lv = _level(pt)
    M = pt.matrix(lv).astype(np.float32)
    f = pt.npix // M.shape[0]
    hap = [(re.match(hap_regex, n).group(1) if re.match(hap_regex, n) else "") for n in pt.names]
    if len({h for h in hap if h}) < 2:
        raise ValueError("fewer than two haplotypes found in scaffold names with " + hap_regex)

    def rng(s):
        a, b = int(pt.scaffold_start[s]), int(pt.scaffold_start[s] + pt.scaffold_npix[s])
        return slice(-(-a // f), b // f)

    bins = lambda s: max(1, rng(s).stop - rng(s).start)
    line, strand, total = {}, {}, {}
    if blocks is not None:
        raw = getattr(pt, "pt", pt)
        idx = {n: i for i, n in enumerate(raw.names)}
        bpp = raw.bp_per_pix

        def unit(name, bp):
            i = idx.get(name)
            if i is None or raw.scaffold_npix[i] == 0:
                return None
            px = int(raw.scaffold_start[i]) + min(int(bp / bpp), int(raw.scaffold_npix[i]) - 1)
            return int(pt.orig_scaffold[px])

        acc, any_bp = {}, {}
        for bk in blocks:
            x, y = unit(bk.q, (bk.qs + bk.qe) / 2), unit(bk.t, (bk.ts + bk.te) / 2)
            if x is None or y is None or x == y:
                continue
            key = (min(x, y), max(x, y))
            any_bp[key] = any_bp.get(key, 0) + bk.qe - bk.qs          # every mapping quality
            if bk.mapq < 30:
                continue
            a_ = acc.setdefault(key, [0, 0])
            a_[0] += bk.qe - bk.qs
            a_[1] += (bk.qe - bk.qs) if bk.same_strand else 0
        loose = {}
        for (x, y), bp in any_bp.items():
            loose.setdefault(x, []).append((bp, y))
            loose.setdefault(y, []).append((bp, x))
        length = lambda s: pt.scaffold_npix[s] * bpp
        big = [s for s in range(pt.n_scaffolds) if hap[s] and length(s) >= min_seq_bp]
        keep = set(big)
        for (x, y), (bp, same_bp) in acc.items():
            if x in keep and y in keep and hap[x] != hap[y]:
                line[(x, y)] = line[(y, x)] = 100.0 * bp / min(length(x), length(y))
                total[(x, y)] = total[(y, x)] = float(bp)
                strand[(x, y)] = strand[(y, x)] = same_bp >= 0.5 * bp
        if score_min is None:
            score_min = 20.0          # percent of the shorter scaffold aligned
    else:
        big = [s for s in range(pt.n_scaffolds) if bins(s) >= min_bins and hap[s]]
        for i, x in enumerate(big):
            for y in big[i + 1:]:
                if hap[x] == hap[y]:
                    continue
                B = M[rng(x), rng(y)]
                if B.shape[0] > B.shape[1]:
                    B = B.T
                # the block's own diagonal and anti-diagonal. (Sliding the shorter scaffold along
                # the longer one was tried and made pairing of cut pieces worse; with sequence
                # homology from the assembly this Hi-C estimate becomes a cross-check only.)
                nx = B.shape[0]
                ii = np.arange(nx)
                jj = np.round(ii * (B.shape[1] - 1) / max(nx - 1, 1)).astype(int)
                d, a = float(B[ii, jj].mean()), float(B[ii, B.shape[1] - 1 - jj].mean())
                line[(x, y)] = line[(y, x)] = max(d, a)
                # ranking uses the total along the line: a few noisy bins of a tiny scaffold must
                # not outrank a line that runs the whole length of a large one
                total[(x, y)] = total[(y, x)] = max(d, a) * min(B.shape)
                strand[(x, y)] = strand[(y, x)] = d >= a
    vals = np.array(list(line.values()))
    background = float(np.median(vals)) if len(vals) else 0.0
    best, floor, ranked = {}, {}, {}
    for x in big:
        cand = sorted(((total[(x, y)], y) for y in big if (x, y) in line), reverse=True)
        ranked[x] = cand
        best[x] = cand[:2] + [(0.0, -1)] * (2 - len(cand[:2]))
        cand = [(line[(x, y)], y) for _, y in cand]
        # a homolog must stand well above what this scaffold shows with unrelated ones
        # (large scaffolds have lower values overall, so the level is set per scaffold)
        floor[x] = score_min if score_min is not None else max(2.0, 5 * float(np.median([c[0] for c in cand] or [0])))
    tel = pt.graphs.get("telomere")
    tel_thr = max(50.0, 5 * float(np.percentile(tel, 99))) if tel is not None else np.inf
    covg = pt.graphs.get("coverage")
    covn = None
    if covg is not None:
        largest = pt.scaffold_npix[pt.orig_scaffold] >= 0.2 * pt.scaffold_npix.max()
        covn = covg / max(float(np.median(covg[largest])), 1e-9)
    pairs, used, repetitive = [], set(), {}
    for x in sorted(big, key=lambda s: -pt.scaffold_npix[s]):
        if x in used:
            continue
        (t1, y), (t2, _) = best[x]
        if y >= 0 and best[y][0][1] == x and covn is not None:
            # too repetitive to assign: long-read coverage well above a normal copy means
            # collapsed repeats; such a scaffold matches several places and is left unplaced
            rep = [u for u in (x, y) if float(np.median(covn[pt.scaffold_start[u]:pt.scaffold_start[u] + pt.scaffold_npix[u]])) >= repeat_cov]
            if rep:
                for u in rep:
                    repetitive[u] = (y if u == x else x)
                continue
        if y < 0 or y in used or best[y][0][1] != x or line[(x, y)] < max(floor[x], floor[y]):
            continue
        # the runner-up that matters is one that is not itself a piece of the same homolog
        # (a scaffold whose own best partner is x or y is a co-partner, added below)
        def rival(u, v):
            co = [z for _, z in ranked[u][1:] if best[z][0][1] == u]
            # pieces of one homolog cannot add up to much more than the scaffold they match;
            # if they do (microchromosomes contact each other a lot), they are rivals
            fits = bins(v) + sum(bins(z) for z in co) <= 1.3 * bins(u)
            return max([t for t, z in ranked[u][1:] if not (fits and z in co)] + [0.0])
        if max(rival(x, y), rival(y, x)) * margin > t1:
            continue
        s1 = line[(x, y)]
        second = max([line[(x, z)] for _, z in best[x][1:] if z >= 0] + [line[(y, z)] for _, z in best[y][1:] if z >= 0] + [0.0])
        ratio = min(pt.scaffold_npix[x], pt.scaffold_npix[y]) / max(pt.scaffold_npix[x], pt.scaffold_npix[y])
        kind = "pair" if ratio >= unequal else "unequal"
        note = "" if kind == "pair" else (
            f"lengths differ {1 / ratio:.1f}-fold: possible sex chromosomes (homology in a shared region only), "
            "or one copy is fragmented; orientation left as assembled")
        pairs.append(Pair(x, y, s1, second, strand[(x, y)], kind, note))
        used.update((x, y))
    for x in sorted(big, key=lambda s: -pt.scaffold_npix[s]):     # pieces of a homolog split over scaffolds
        if x in used:
            continue
        t1, y = best[x][0]
        s1 = line[(x, y)] if y >= 0 else 0.0
        main = [p for p in pairs if y in (p.a, p.b) and p.kind != "extra"]
        if not main:
            continue
        mate = main[0].b if main[0].a == y else main[0].a
        taken = bins(mate) + sum(bins(p.b) for p in pairs if p.kind == "extra" and p.a == y)
        # a further piece must be clear (2x the next best), comparable to the main homology,
        # and fit in the length of the scaffold it matches
        sized = sorted((abs(np.log(min(bins(p.a), bins(p.b)) / bins(x))), p.score) for p in pairs if p.kind == "pair")
        typical = float(np.median([sc for _, sc in sized[:5]])) if sized else 0.0
        if (s1 >= floor[x] and best[x][1][0] * 2.0 <= t1 and s1 >= 0.4 * typical
                and taken + bins(x) <= 1.3 * bins(y)):
            pairs.append(Pair(y, x, s1, line[(x, best[x][1][1])] if best[x][1][1] >= 0 else 0.0, strand[(x, y)], "extra",
                              "further piece homologous to the same scaffold"))
            used.add(x)
    info = dict(level=lv, bin_px=f, background=background, big=big, hap=hap, line=line, strand=strand,
                source="sequence alignment" if blocks is not None else "MQ0 Hi-C",
                unpaired=[s for s in big if s not in used], tel_thr=tel_thr,
                repetitive=repetitive)
    return pairs, info


def telomere_sites(pt, pairs, hi=None, hap_regex=r"^(H\d+)_", blocks=None, contact_max=0.5, near_px=3, end_px=15,
                   tol_frac=0.01, broad_min=0.6, invert_min=1.3):
    """Gaps inside a scaffold that carry a telomere signal with reduced contact across, plus
    the homologous gap in the paired scaffold (dual scaffolding). Each is decided as:

      cut         the two sides have different homologs in the other haplotype, or share
                  no chromosome-level contact: two chromosomes joined end to end
      invert_tip  the sides belong together and the end piece fits better the other way
                  round (the telomere then ends up at the scaffold end)
      keep        the sides belong together as assembled; flagged for the curator

    Returns a list of dicts: q (site is after this pixel), s (scaffold), action, reason, tip
    (pixel range of the end piece).
    """
    tel, gap = pt.graphs.get("telomere"), pt.graphs.get("gap")
    if tel is None or gap is None:
        return []
    thr = max(50.0, 5 * float(np.percentile(tel, 99)))
    sg0, sgh = Signal(pt), Signal(hi if hi is not None else pt)
    near = sgh.insulation(10, 0)
    which = "high-MQ map" if hi is not None else "MQ0 map"

    def contact(q):
        w = near[max(0, q - 2):q + 3]
        return float(np.nanmin(w)) if np.isfinite(w).any() else float("nan")

    sites, carried = {}, {}
    for q in np.nonzero(gap > 0)[0]:
        q = int(q)
        s = int(pt.orig_scaffold[q])
        a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        if q - a < end_px or a + n - 1 - q < end_px:
            continue
        t = float(tel[max(a, q - near_px):q + near_px + 1].max())
        c = contact(q)
        if t >= thr and c < contact_max:
            sites[q] = (s, f"telomere signal ({t:.0f}) at a gap inside the scaffold, contact across {c:.2f} on the {which}")
    partner = {}
    for p in pairs:
        if p.kind == "pair":
            partner[p.a], partner[p.b] = (p.b, p.same_strand), (p.a, p.same_strand)
    for q, (s, _) in list(sites.items()):
        if s not in partner:
            continue
        y, same = partner[s]
        a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        b, m = int(pt.scaffold_start[y]), int(pt.scaffold_npix[y])
        cand = b + np.nonzero(gap[b:b + m] > 0)[0]
        cand = [int(c) for c in cand if c - b >= end_px and b + m - 1 - c >= end_px
                and contact(int(c)) < contact_max and not any(abs(int(c) - k) <= 5 for k in sites)]
        hit = None
        if blocks is not None:
            # Dual scaffolding with alignments: a gap in the partner whose two flanks are
            # homologous to sequence on opposite sides of this site. Cutting both then gives
            # pieces that match one to one, whatever way round each scaffold joined them.
            bpp, sname, yname = pt.bp_per_pix, pt.names[s], pt.names[y]

            def to_s(g, side):
                """Position in s (pixels) homologous to the partner sequence just left (-1)
                or right (+1) of gap pixel g, from the nearest alignment end within 2 Mb."""
                gbp, best_d, best = (g - b) * bpp, 2e6, None
                for bk in blocks:
                    if bk.q == yname and bk.t == sname:
                        own, other = (bk.qs, bk.qe), (bk.ts, bk.te)
                    elif bk.t == yname and bk.q == sname:
                        own, other = (bk.ts, bk.te), (bk.qs, bk.qe)
                    else:
                        continue
                    k = 1 if side < 0 else 0                 # block end facing the gap
                    d = (gbp - own[k]) * (1 if side < 0 else -1)
                    if -2 * bpp <= d < best_d:
                        best_d, best = d, a + int((other[k] if bk.same_strand else other[1 - k]) / bpp)
                return best
            for c in cand:
                lft, rgt = to_s(c, -1), to_s(c, +1)
                if lft is not None and rgt is not None and min(lft, rgt) < q < max(lft, rgt) \
                        and abs(lft - rgt) > 10 * max(30, int(tol_frac * n)):
                    hit = c
                    break
        else:
            fr = (q - a) / n
            target = b + int(round((fr if same else 1 - fr) * m))
            tol = max(30, int(tol_frac * m))
            near_c = [c for c in cand if abs(c - target) <= tol]
            hit = min(near_c, key=lambda c: abs(c - target)) if near_c else None
        if hit is not None:
            carried[hit] = q
            sites[hit] = (y, f"gap matching the telomere-bearing gap in {pt.names[s]} (other haplotype): the two sides "
                             f"are homologous to the two sides of that gap; contact across {contact(hit):.2f} on the {which}")
    if not sites:
        return []
    # 1. homologs of the two sides, after cutting every site tentatively
    vt = SegmentedMap(pt, list(sites))
    pt_pairs, _ = pair_homologs(vt, hap_regex, blocks=blocks)
    mate = {}
    for p in pt_pairs:
        mate.setdefault(p.a, p.b)
        mate.setdefault(p.b, p.a)
    parent = lambda u: vt.parent_of(vt.names[u])[0]
    # 2. chromosome-level contact between the sides (MQ0, binned, background removed)
    f = 8
    M = sg0.binned(f)
    bg = float(M[::37, ::41].mean())
    nb = M.shape[0]
    scb = pt.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, pt.npix - 1)]
    D = 40
    g = np.array([np.diagonal(M, d)[scb[:nb - d] == scb[d:]].mean() for d in range(D)])
    gh = sgh.decay_long(400)
    out = []
    for q, (s, why) in sorted(sites.items()):
        a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        left, right = q - a + 1, n - (q - a + 1)
        tip_right = right <= left
        tip = (q + 1, a + n - 1) if tip_right else (a, q)
        L, Rr = int(vt.orig_scaffold[q]), int(vt.orig_scaffold[q + 1])
        ml, mr = mate.get(L), mate.get(Rr)
        b = q // f
        T = min((q - a + 1) // f, (a + n - 1 - q) // f, 15)
        broad = float("nan")
        if T - 2 >= 2:
            rows, cols = np.arange(b - T + 1, b - 1), np.arange(b + 3, b + T + 1)
            B = M[np.ix_(rows, cols)]
            ex = g[np.minimum(cols[None, :] - rows[:, None], D - 1)]
            broad = float((B.mean() - bg) / max(ex.mean() - bg, 1e-6))
        if ml is not None and mr is not None and parent(ml) != parent(mr):
            action = "cut"
            why += f"; the two sides have different homologs ({vt.names[ml]} and {vt.names[mr]}): two chromosomes joined"
        elif not broad >= broad_min:
            action = "cut"
            why += (f"; chromosome-level contact between the sides is {broad:.2f} of what one chromosome shows: "
                    "two chromosomes joined end to end")
        else:
            # same chromosome. Does the end piece fit better the other way round? (high-MQ map)
            W = max(3, min((tip[1] - tip[0] + 1) // 2, 50))
            body = (q - W + 1, q) if tip_right else (q + 1, q + W)
            near_w = (q + 1, q + W) if tip_right else (q - W + 1, q)
            far_w = (tip[1] - W + 1, tip[1]) if tip_right else (tip[0], tip[0] + W - 1)
            i = np.arange(W)
            ex = gh[np.minimum(i[:, None] + i[None, :] + 1, len(gh) - 1)].mean()
            as_is = float(sgh.block(body[0], body[1], near_w[0], near_w[1]).mean() / ex)
            inv = float(sgh.block(body[0], body[1], far_w[0], far_w[1]).mean() / ex)
            why += f"; the sides share chromosome-level contact ({broad:.2f}): same chromosome"
            if inv >= invert_min * max(as_is, 1e-3) and (tip[1] - tip[0] + 1) <= 0.4 * n:
                action = "invert_tip"
                why += (f"; the end piece fits better inverted ({inv:.2f} against {as_is:.2f} as assembled, {which}), "
                        "which also puts the telomere at the scaffold end")
            else:
                action = "keep"
                why += f"; end piece as assembled {as_is:.2f}, inverted {inv:.2f} ({which}): left as it is"
        out.append(dict(q=q, s=s, action=action, reason=why, tip=tip))
    # dual scaffolding: a site found only as the homologous point of a misjoin is a misjoin too
    # (at MQ0 its own contact test is unreliable: reads from the other haplotype bridge it)
    by_q = {x["q"]: x for x in out}
    for c, q in carried.items():
        if by_q[q]["action"] == "cut" and by_q[c]["action"] != "cut":
            by_q[c]["action"] = "cut"
            by_q[c]["reason"] += f"; cut because the homologous site in {pt.names[by_q[q]['s']]} is a misjoin"
    return out


def curate_allhaps(pt, phasing="hic", hap_regex=r"^(H\d+)_", hi=None, do_cuts=True, blocks=None, busco_hits=None):
    """allhaps curation so far: telomere-supported cuts, homolog pairing, parallel layout,
    HAP tags. `hi` is the high-MQ map of the same assembly, if there is one.
    Returns (rows for the decisions table, savestate, info)."""
    raw = pt
    first_pairs, _ = pair_homologs(raw, hap_regex, blocks=blocks)
    sites = telomere_sites(raw, first_pairs, hi, hap_regex, blocks) if do_cuts else []
    cuts = [(x["q"], x["s"], x["reason"]) for x in sites if x["action"] == "cut"]
    flips = [x["tip"] for x in sites if x["action"] == "invert_tip"]
    pt = SegmentedMap(raw, [q for q, _, _ in cuts])
    pairs, info = pair_homologs(pt, hap_regex, blocks=blocks)
    hap = info["hap"]
    unit_h = "% of the shorter aligned" if blocks is not None else ""
    haps = sorted({h for h in hap if h})
    cov = pt.graphs.get("coverage")
    npx, start = pt.scaffold_npix, pt.scaffold_start
    main = [p for p in pairs if p.kind != "extra"]
    extra = {}
    for p in pairs:
        if p.kind == "extra":
            extra.setdefault(p.a, []).append(p)
    paired = np.concatenate([np.arange(start[s], start[s] + npx[s]) for p in main for s in (p.a, p.b)]) \
        if main else np.array([], dtype=int)
    norm = float(np.median(cov[paired])) if cov is not None and len(paired) else 0.0
    rows, order, tag_of, paint_of = [], [], {}, {}
    mb = lambda s: npx[s] * pt.bp_per_pix / 1e6
    # BUSCO genes per unit (complete copies, and how many of them are duplicated in the unit)
    genes, dup_of = {}, {}
    if busco_hits is not None:
        from .busco import summarise as busco_summary
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
        from .busco import false_duplicates
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
    for q, s, why in cuts:
        rows.append(dict(action="cut", scaffold=raw.names[s], length_mb=round(raw.rel_coord[q] * raw.bp_per_pix / 1e6, 2),
                         partner="", strand="", tag="", score="", next_best="",
                         evidence=f"cut at {raw.rel_coord[q] * raw.bp_per_pix / 1e6:.2f} Mb; " + why))
    for x in sites:
        if x["action"] != "cut":
            rows.append(dict(action=x["action"], scaffold=raw.names[x["s"]],
                             length_mb=round(raw.rel_coord[x["q"]] * raw.bp_per_pix / 1e6, 2), partner="", strand="",
                             tag="", score="", next_best="",
                             evidence=f"end piece {(x['tip'][1] - x['tip'][0] + 1) * raw.bp_per_pix / 1e6:.1f} Mb; " + x["reason"]))

    def pixels(s, rev=False):
        """Pixels of a unit in layout order, with any inverted end piece turned in place."""
        px = np.arange(start[s], start[s] + npx[s])
        for a, b in flips:
            if start[s] <= a and b < start[s] + npx[s]:
                px[a - start[s]:b - start[s] + 1] = px[a - start[s]:b - start[s] + 1][::-1].copy()
        return px[::-1] if rev else px
    for k, p in enumerate(sorted(main, key=lambda p: -(npx[p.a] + npx[p.b]))):
        qa, ga, ta = _quality(pt, p.a, info["tel_thr"])
        qb, gb, tb = _quality(pt, p.b, info["tel_thr"])
        first, second = p.a, p.b
        hemi = norm and all(float(np.median(cov[start[u]:start[u] + npx[u]])) / norm < 0.7 for u in (p.a, p.b))
        if p.kind == "unequal" and not hemi:
            p.kind, p.note = "pair", "lengths differ: one copy is probably fragmented or incomplete"
        if p.kind == "unequal":
            # convention: both sex chromosomes go to HAP1; the larger one first
            first, second = (p.a, p.b) if npx[p.a] >= npx[p.b] else (p.b, p.a)
            tags = ("HAP1", "HAP1")
            why = "possible sex chromosome pair: both tagged HAP1 by convention"
        elif phasing == "trio":
            first, second = sorted((p.a, p.b), key=lambda s: haps.index(hap[s]))
            tags = ("HAP" + str(haps.index(hap[first]) + 1), "HAP" + str(haps.index(hap[second]) + 1))
            why = "trio phasing: tags follow the name prefix"
        else:
            # hic / none: whole superscaffolds may swap; the better copy becomes HAP1.
            # BUSCO decides only when the copies clearly differ (more single-copy genes, fewer
            # duplicated); otherwise fewer gaps, telomeric ends, length.
            ba, bb = genes.get(p.a, (0, 0)), genes.get(p.b, (0, 0))
            sa, sb = ba[0] - 2 * ba[1], bb[0] - 2 * bb[1]
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
        order.append((first, False))
        # an unequal pair shares only a small region (e.g. a PAR): too little to orient by
        flip = (not p.same_strand) and p.kind == "pair"
        order.append((second, flip))
        paint_of[first], paint_of[second] = 2 * k + 1, 2 * k + 2
        tag_of[first], tag_of[second] = tags
        for e in sorted(extra.get(p.a, []) + extra.get(p.b, []), key=lambda e: -npx[e.b]):
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
    sex = [s for s in info["unpaired"] if mb(s) >= 5 and covr(s) < 0.7]
    k = len(main)
    for s in sorted(sex, key=lambda s: -npx[s]):
        partner = [t for t in sex if t != s and hap[t] != hap[s] and (s, t) in info["line"]]
        partner = max(partner, key=lambda t: info["line"][(s, t)]) if partner else None
        rev = False
        if partner is not None and partner in {x for x, _ in order}:
            rev = not info["strand"][(s, partner)]
        order.append((s, rev))
        tag_of[s] = "HAP1"
        paint_of[s] = 2 * k + 1
        k += 1
        rows.append(dict(action="sex_chromosome_candidate", scaffold=pt.names[s], length_mb=round(mb(s), 2),
                         partner=pt.names[partner] if partner is not None else "", strand="-" if rev else "+",
                         tag="HAP1", score=round(info["line"][(s, partner)], 1) if partner is not None else "",
                         next_best="", evidence=f"no homolog of similar size; coverage {covr(s):.2f} of a paired copy "
                                                "(hemizygous); tagged HAP1 by convention; which is X/Y/Z/W is not determined"))
    placed = {s for s, _ in order}
    for s in info["unpaired"]:
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
    rest = [int(s) for s in np.argsort(start) if npx[s] > 0 and int(s) not in placed]
    perm = np.concatenate([pixels(s, rev) for s, rev in order] + [pixels(s) for s in rest])
    st = SaveState.blank()
    where = np.empty(raw.npix, dtype=np.int64)
    where[perm] = np.arange(raw.npix)
    st.edits = [Edit(q, q, BREAK_DELTA) for q, _, _ in cuts if where[q + 1] - where[q] == 1] + plan_edits(perm)
    ms = MapState(raw).replay(st.edits)
    if not np.array_equal(ms.pix, perm):
        raise AssertionError("planned edits do not reproduce the target layout")
    for k, q in enumerate([x["q"] for x in sites][:255]):
        x = where[q] / raw.npix - 0.5
        st.waypoints.append(Waypoint(x, x, 30.0, k))
    for i, c in enumerate(ms.contigs):
        s = int(pt.orig_scaffold[ms.pix[c.start]])
        if s in paint_of:
            st.scaffs.append((i, paint_of[s]))
        t = tag_of.get(s) or (("HAP" + str(haps.index(hap[s]) + 1)) if hap[s] and haps.index(hap[s]) < 2 else "")
        mask = (1 << st.tag_id(t)) if t else 0
        if c.scaffold in dup_of:
            mask |= 1 << st.tag_id("vertPaint")
        if mask:
            st.meta_flags.append((i, mask))
    return rows, st, info


def write_rows(rows, path):
    cols = ["action", "scaffold", "length_mb", "partner", "strand", "tag", "score", "next_best", "evidence"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
