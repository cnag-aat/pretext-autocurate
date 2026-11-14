"""allhaps mode, step 1: pair homologous scaffolds of different haplotypes and lay each
pair side by side in the same orientation, tagged HAP1 / HAP2.

Homology is read from the MQ0 map: reads from sequence shared by two haplotypes are split
between the copies, so a scaffold shows a line of contact with its homolog, along the
diagonal of their block if both are assembled in the same direction and along the
anti-diagonal if not. Unrelated scaffolds show only background.

Also here: telomere-supported cuts, small scaffolds attached as unlocs (`place_small`),
homologs within one haplotype file (`same_prefix_pairs`). Not done yet: joins, haplotig
tags (see docs/HANDOFF.md).
"""
import re
from dataclasses import dataclass

import numpy as np

from .mapstate import MapState, plan_edits
from .pretextio import SegmentedMap
from .savestate import SaveState, Edit, Waypoint, BREAK_DELTA
from .signal import Signal


@dataclass
class Level:
    """One setting for how readily the tool acts: 0 = conservative (act only on the clearest
    evidence, flag the rest), 0.5 = default, 1 = aggressive (act on weaker evidence too).
    Every threshold that decides between acting and flagging is read through `at`."""
    a: float = 0.5

    def at(self, conservative, default, aggressive):
        a = min(max(self.a, 0.0), 1.0)
        return conservative + (default - conservative) * a / 0.5 if a <= 0.5 else \
            default + (aggressive - default) * (a - 0.5) / 0.5


@dataclass
class Pair:
    a: int                # scaffold index, placed first (HAP1)
    b: int                # its homolog (HAP2)
    score: float          # mean of the homology line (map units at the working level)
    second: float         # best line of either scaffold with any other scaffold
    same_strand: bool
    kind: str             # pair | unequal (possible sex chromosomes) | extra
    note: str


def _signal(m):
    """One Signal per map, kept on the map object: its decay curve, insulation profiles,
    binned matrix and decoded tiles are reused between steps and between sweep levels."""
    m = getattr(m, "pt", m)
    if "_signal" not in m.__dict__:
        m.__dict__["_signal"] = Signal(m)
    return m.__dict__["_signal"]


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
    base = getattr(pt, "pt", pt)
    if ("_mip", lv) not in base.__dict__:
        base.__dict__[("_mip", lv)] = base.matrix(lv).astype(np.float32)
    M = base.__dict__[("_mip", lv)]
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
                   tol_frac=0.01, broad_min=0.6, invert_min=1.3, gap_max=0.6, sharp_max=0.1, merge_px=5, tip_max=0.5, use_breaks=True):
    """Candidate break sites, plus the homologous gap in the paired scaffold (dual scaffolding):
      - a gap inside a scaffold with contact across below `gap_max` of expectation (a telomere
        signal there is supporting evidence, not a requirement: telomere tracks are not
        reliable in every species);
      - a point with no gap where contact drops to background on both scales (`sharp_max`),
        with normal contact on each side: cut outright.
    Each gap site is decided as:

      cut         the two sides have different homologs in the other haplotype, or share
                  no chromosome-level contact: two chromosomes joined end to end
      invert_tip  the sides belong together and the end piece fits better the other way
                  round (the telomere then ends up at the scaffold end)
      keep        the sides belong together as assembled; flagged for the curator

    Returns a list of dicts: q (site is after this pixel), s (scaffold), action, reason, tip
    (pixel range of the end piece).
    """
    tel, gap = pt.graphs.get("telomere"), pt.graphs.get("gap")
    if gap is None:
        return []
    if tel is None:
        tel = np.zeros(pt.npix)
    thr = max(50.0, 5 * float(np.percentile(tel, 99)))
    sg0 = _signal(pt)
    sgh = _signal(hi) if hi is not None else sg0
    near = sgh.insulation(10, 0)
    which = "high-MQ map" if hi is not None else "MQ0 map"

    near0 = sg0.insulation(10, 0) if hi is not None else near

    def contact(q, near=near):
        w = near[max(0, q - 2):q + 3]
        return float(np.nanmin(w)) if np.isfinite(w).any() else float("nan")

    sites, carried, site_c = {}, {}, {}
    for q in np.nonzero(gap > 0)[0]:
        q = int(q)
        s = int(pt.orig_scaffold[q])
        a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        if q - a < end_px or a + n - 1 - q < end_px:
            continue
        t = float(tel[max(a, q - near_px):q + near_px + 1].max())
        c = contact(q)
        if t >= thr and c < contact_max:
            pass
        elif not (c < gap_max and (near0 is near or contact(q, near0) < gap_max)):
            continue                           # without a telomere: reduced on the MQ0 map as well as the high-MQ one
        close = [k for k in sites if abs(k - q) <= merge_px]
        if close:                              # several gap pixels at one place: keep the weakest point
            if c >= site_c[close[0]]:
                continue
            del sites[close[0]]
        site_c[q] = c
        sites[q] = (s, (f"telomere signal ({t:.0f}) at a gap inside the scaffold, " if t >= thr else "gap inside the scaffold, ")
                    + f"contact across {c:.2f} of expectation on the {which}")
    # sharp breaks away from any gap: contact at background on both scales (dense maps only:
    # on a sparse high-MQ map a 10 px window is too noisy for this)
    sharp = {}
    from .misjoin import detect_misjoins, Params as MisjoinParams
    dense = sg0.decay[min(20, len(sg0.decay) - 1)] >= 0.25 * sg0.decay[0]
    if sharp_max > 0 and dense:
        for cdt in detect_misjoins(pt, MisjoinParams(cut_max=sharp_max), sg0):
            q = int(cdt.after_px)
            s = int(pt.orig_scaffold[q])
            a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
            if (not cdt.apply or cdt.at_gap or q - a < end_px or a + n - 1 - q < end_px
                    or any(abs(k - q) <= 2 * merge_px for k in sites) or not cdt.far < sharp_max):
                continue
            sharp[q] = (s, f"no gap here, but contact across drops to background on the MQ0 map ({cdt.near:.2f} within "
                           f"10 px, {cdt.far:.2f} at 10-40 px) with normal contact on each side")
    # gaps where the alignment to the other haplotype breaks: candidates even with normal contact
    if blocks is not None and use_breaks:
        for bq, kind, lt, rt in alignment_breaks(pt, blocks):
            near_gaps = [x for x in range(max(0, bq - 2), min(pt.npix, bq + 3)) if gap[x] > 0]
            if not near_gaps:
                continue
            q = min(near_gaps, key=lambda x: abs(x - bq))
            s = int(pt.orig_scaffold[q])
            a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
            if q - a < 5 or a + n - 1 - q < 5 or any(abs(k - q) <= 2 * merge_px for k in sharp):
                continue
            note = f"the alignment to the other haplotype breaks here ({kind}: {lt} | {rt})"
            close = [k for k in sites if abs(k - q) <= merge_px]
            if close:
                sites[close[0]] = (sites[close[0]][0], sites[close[0]][1] + "; " + note)
                continue
            site_c[q] = contact(q)
            sites[q] = (s, f"gap inside the scaffold; {note}; contact across {site_c[q]:.2f} of expectation on the {which}")
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
                and contact(int(c)) < max(contact_max, gap_max) and (int(c) in sites or not any(abs(int(c) - k) <= 5 for k in sites))]
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
        if hit is not None and hit in sites:
            carried[hit] = q                   # already a site on its own evidence; the homology adds to it
        elif hit is not None:
            carried[hit] = q
            sites[hit] = (y, f"gap matching the telomere-bearing gap in {pt.names[s]} (other haplotype): the two sides "
                             f"are homologous to the two sides of that gap; contact across {contact(hit):.2f} on the {which}")
    if not sites and not sharp:
        return []
    # 1. homologs of the two sides, after cutting every site tentatively
    vt = SegmentedMap(pt, list(sites) + list(sharp))
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
    t5, t3 = pt.graphs.get("telomere_5p"), pt.graphs.get("telomere_3p")
    from .tracks import telomere_threshold
    thr_s = telomere_threshold(t5 + t3) if t5 is not None and t3 is not None else None

    def strand_at(q, tip_right):
        """With telomere strand tracks: which side of the gap carries a telomere facing it.
        A chromosome's right end reads the motif forward (3p), its left end reverse (5p), so
        a 3p peak just left of the gap and/or a 5p peak just right of it face the gap.
        Returns 'both', 'tip' (only the end piece), 'body' (only the long side) or None."""
        if thr_s is None:
            return None
        # a pixel spans tens of kb: telomeres either side of a gap often fall in the gap's pixel
        # or its neighbour, so each side's window includes the pixel across the gap
        l3 = float(t3[max(0, q - near_px):q + 2].max()) >= thr_s
        r5 = float(t5[q:q + 2 + near_px].max()) >= thr_s
        if l3 and r5:
            return "both"
        if not (l3 or r5):
            return None
        left_side = l3                          # the facing telomere is on the left of the gap
        return "tip" if left_side != tip_right else "body"
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
        strand = strand_at(q, tip_right)
        if ml is not None and mr is not None and parent(ml) != parent(mr):
            action = "cut"
            why += f"; the two sides have different homologs ({vt.names[ml]} and {vt.names[mr]}): two chromosomes joined"
        elif strand in ("both", "body"):
            action = "cut"
            why += ("; telomeres on both sides face the gap: two chromosome ends joined" if strand == "both" else
                    "; the telomere on the long side faces the gap: that side ends a chromosome, the end piece "
                    "belongs elsewhere") + " (telomere strand tracks)"
        elif not np.isfinite(broad) and site_c.get(q, contact(q)) < tip_max:
            action = "cut"
            why += ("; one side is too short to measure chromosome-level contact with the other, and contact across "
                    "the gap is clearly reduced: cut")
        elif broad < broad_min:
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
            if strand == "tip" and (tip[1] - tip[0] + 1) <= 0.4 * n:
                action = "invert_tip"
                why += ("; the end piece's telomere faces the gap (telomere strand tracks): inverted, it ends the "
                        f"chromosome (contact: inverted {inv:.2f}, as assembled {as_is:.2f}, {which})")
            elif inv >= invert_min * max(as_is, 1e-3) and (tip[1] - tip[0] + 1) <= 0.4 * n:
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
    for q, (s, why) in sharp.items():
        out.append(dict(q=q, s=s, action="cut", reason=why, tip=(q, q)))
    out.sort(key=lambda x: x["q"])
    return out


def _union(iv):
    tot, end = 0, -1
    for a, b in sorted(iv):
        if b > end:
            tot += b - max(a, end)
            end = b
    return tot


def _small_alignments(pt, small_blocks):
    """{query unit: {target unit: [intervals on the query, intervals with mapq >= 30, bases on
    the same strand, bases]}} from the alignment of small scaffolds against everything.
    A unit's alignments to itself (internal repeats) are left out."""
    raw = getattr(pt, "pt", pt)
    idx = {n: i for i, n in enumerate(raw.names)}
    bpp = raw.bp_per_pix

    def unit(name, bp):
        i = idx.get(name)
        if i is None or raw.scaffold_npix[i] == 0:
            return None
        return int(pt.orig_scaffold[int(raw.scaffold_start[i]) + min(int(bp / bpp), int(raw.scaffold_npix[i]) - 1)])

    out = {}
    for bk in small_blocks:
        x, y = unit(bk.q, (bk.qs + bk.qe) / 2), unit(bk.t, (bk.ts + bk.te) / 2)
        if x is None or y is None or x == y:
            continue
        e = out.setdefault(x, {}).setdefault(y, [[], [], 0, 0, []])
        # with the density of matches, and where it lands in the target scaffold
        e[4].append((bk.qs, bk.qe, bk.matches / max(bk.qe - bk.qs, 1), bk.t, bk.ts, bk.te, bk.same_strand))
        e[0].append((bk.qs, bk.qe))
        if bk.mapq >= 30:
            e[1].append((bk.qs, bk.qe))
        e[2] += (bk.qe - bk.qs) if bk.same_strand else 0
        e[3] += bk.qe - bk.qs
    return out


def same_prefix_pairs(pt, info, small_blocks, min_bp=500000, frac_min=0.5):
    """Homologs that ended up in the same haplotype file (both `H1_`, say): two unpaired
    scaffolds of one prefix that are each other's best match, with at least half of the
    shorter one aligned uniquely to the other. The H2-against-H1 alignment cannot show these."""
    aln = _small_alignments(pt, small_blocks)
    bpp, hap, npx = pt.bp_per_pix, info["hap"], pt.scaffold_npix
    free = set(info["unpaired"]) | set(info.get("repetitive", {}))
    free = {u for u in free if npx[u] * bpp >= min_bp}

    def best(u):
        t = {v: _union(e[0]) for v, e in aln.get(u, {}).items()}
        return max(t, key=t.get) if t else None

    out = []
    for x in sorted(free, key=lambda u: -npx[u]):
        y = best(x)
        if y is None or y not in free or hap[y] != hap[x] or best(y) != x or (npx[y], y) > (npx[x], x):
            continue
        e = aln[y][x]
        uniq = _union(e[1]) / (npx[y] * bpp)
        if uniq < frac_min:
            continue
        out.append(Pair(x, y, 100.0 * uniq, 0.0, e[2] >= 0.5 * e[3], "pair",
                        f"both copies carry the prefix {hap[x]}_: homologs that ended up in the same haplotype file"))
    return out


def place_small(pt, pairs, info, small_blocks=None, hi=None, seq_min=0.5, seq_margin=2.0, hic_agree=1.5,
                hic_alone=2.5, hic_alone_px=20, hi_margin=2.0, repeat_cov=1.5, dense_min=0.5, repeat_frac=0.9,
                phasing="hic", copy_margin=1.5, cis_frac=0.2):
    """Which chromosome each remaining small scaffold belongs to.

    Two independent lines of evidence: sequence (share of the scaffold that aligns to a
    chromosome's copies, from the small-against-all alignment) and Hi-C (mean contact with
    each chromosome on the MQ0 map, best copy; confirmed on the high-MQ map when it is the
    only evidence). Most small scaffolds of an allhaps assembly are repeats that contact a
    repeat-rich chromosome at MQ0 and nothing at high MQ; those are left alone.

      unloc      both lines agree, or Hi-C alone is strong on both maps: attach
      candidate  sequence says one chromosome but Hi-C does not confirm, or long-read
                 coverage is repeat-like (collapsed repeats match several places): report,
                 highlight, do not move

    Which copy (curator, 2026-11-13: a small piece may go to the other prefix's copy, also
    under hic phasing; never under trio): sequence found in one copy only -> the copy that
    lacks it; found in both -> the copy it contacts on the high-MQ map, else the shorter one.

    Returns a list of dicts (u, copy, action, third_copy, evidence)."""
    npx, start, hap = pt.scaffold_npix, pt.scaffold_start, info["hap"]
    main = [p for p in pairs if p.kind != "extra"]
    chrom, copies = {}, {}
    for k, p in enumerate(main):
        chrom[p.a] = chrom[p.b] = k
        copies[k] = (p.a, p.b)
    taken = set(chrom) | {p.b for p in pairs if p.kind == "extra"}
    small = [u for u in range(pt.n_scaffolds) if npx[u] > 0 and u not in taken and hap[u]]
    if not small or not main:
        return []
    aln = _small_alignments(pt, small_blocks) if small_blocks is not None else {}
    bpp = pt.bp_per_pix
    seq, nonrep = {}, {}
    for u in small:
        per, L = {}, npx[u] * bpp
        # Non-repetitive sequence: stretches with a dense alignment somewhere (most bases matching).
        # A repeat array aligns loosely to every copy of the array, and how much of it each
        # haplotype assembled says nothing; where there is non-repetitive sequence the decision
        # rests on that alone (curator, 2026-11-14).
        dense = [(a, b, d) for e in aln.get(u, {}).values() for a, b, d, *_ in e[4] if d >= dense_min]
        D = sorted((a, b) for a, b, _ in dense)
        Ld = _union(D)
        if Ld >= max(50000, 0.1 * L) and Ld < 0.9 * L:
            top = max(d for _, _, d in dense)

            def within(iv):
                return _union([(max(a, x), min(b, y)) for a, b in iv for x, y in D if min(b, y) > max(a, x)])
            for v, e in aln.get(u, {}).items():
                if v in chrom:
                    per.setdefault(chrom[v], {})[v] = [(a, b) for a, b, d, *_ in e[4] if d >= 0.7 * top]
            f = sorted(((within([i for iv in d.values() for i in iv]) / Ld, k) for k, d in per.items()), reverse=True)
            if f and f[0][0] > 0:
                k = f[0][1]
                seq[u] = (f[0][0], f[1][0] if len(f) > 1 else 0.0, k, {v: within(iv) / Ld for v, iv in per[k].items()})
                nonrep[u] = (Ld, L)
            continue
        for v, e in aln.get(u, {}).items():
            if v in chrom:
                per.setdefault(chrom[v], {})[v] = e[0]
        f = sorted(((_union([i for iv in d.values() for i in iv]) / L, k) for k, d in per.items()), reverse=True)
        if f:
            k = f[0][1]
            seq[u] = (f[0][0], f[1][0] if len(f) > 1 else 0.0, k,
                      {v: _union(iv) / L for v, iv in per[k].items()})
    raw_ = getattr(pt, "pt", pt)
    ridx_ = {n: i for i, n in enumerate(raw_.names)}

    def _extends(u, X):
        """(X, end of X, end of u) if u's dense alignments to copy X lie at an end of X and at an
        end of u, with the strand that makes u continue X beyond that end; else None."""
        ents = [x for x in aln.get(u, {}).get(X, [[], [], 0, 0, []])[4]]
        if not ents:
            return None
        top = max(d for e in aln.get(u, {}).values() for _, _, d, *_ in e[4])
        ents = [x for x in ents if x[2] >= 0.7 * top and x[2] >= dense_min]
        if not ents:
            return None
        L = npx[u] * bpp
        q0, q1 = min(x[0] for x in ents), max(x[1] for x in ents)
        tol_q = max(3 * bpp, 0.05 * L)
        if q1 - q0 > 0.9 * L:
            return None                                    # nothing of it reaches beyond the overlap
        u_end = "h" if q0 <= tol_q else "t" if L - q1 <= tol_q else None
        px = lambda name, bp: int(raw_.scaffold_start[ridx_[name]] + bp / bpp)
        t0, t1 = min(px(x[3], x[4]) for x in ents), max(px(x[3], x[5]) for x in ents)
        sX, eX = int(start[X]), int(start[X] + npx[X] - 1)
        tol_x = max(3, int(0.03 * npx[X]))
        x_end = "h" if t0 - sX <= tol_x else "t" if eX - t1 <= tol_x else None
        if u_end is None or x_end is None:
            return None
        same = sum(x[1] - x[0] for x in ents if x[6]) >= 0.5 * sum(x[1] - x[0] for x in ents)
        if same != (u_end != x_end):                       # the strand must let u run on past X's end
            return None
        return (X, x_end, u_end)

    need = [u for u in small if npx[u] >= hic_alone_px or (u in seq and seq[u][0] >= seq_min)]
    name = lambda k: "/".join(pt.names[v] for v in copies[k])
    cs = np.concatenate(([0.0], np.zeros(pt.npix)))
    units = sorted(chrom)

    def contact(sg, u):
        a, L = int(start[u]), int(npx[u])
        prof = sg.strip(a, a + L - 1).astype(np.float32).mean(0)
        cs[1:] = np.cumsum(prof, dtype=np.float64)
        per = {v: (cs[start[v] + npx[v]] - cs[start[v]]) / npx[v] for v in units}
        sc = sorted(((max(per[v] for v in copies[k]), k) for k in copies), reverse=True)
        return sc, per

    sg0 = _signal(pt)
    sgh = _signal(hi) if hi is not None else None
    rept = pt.graphs.get("repeat")
    covg, covn = pt.graphs.get("coverage"), None
    if covg is not None:
        largest = npx[pt.orig_scaffold] >= 0.2 * npx.max()
        covn = covg / max(float(np.median(covg[largest])), 1e-9)
    out = []
    for u in sorted(need, key=lambda u: start[u]):
        sc0, per0 = contact(sg0, u)
        r0 = sc0[0][0] / max(sc0[1][0], 1e-6)
        s = seq.get(u)
        seq_ok = s is not None and s[0] >= seq_min and s[0] >= seq_margin * s[1]
        k, action, why = None, None, ""
        if seq_ok:
            k = s[2]
            what = "it" if u not in nonrep else (f"its non-repetitive sequence ({nonrep[u][0] / 1e3:.0f} kb of "
                                                 f"{nonrep[u][1] / 1e3:.0f} kb; the rest aligns only loosely, as repeats do)")
            why = f"{100 * min(s[0], 1):.0f}% of {what} aligns to the chromosome {name(k)} (next chromosome {100 * s[1]:.0f}%)"
            if sc0[0][1] == k and r0 >= hic_agree:
                action = "unloc"
                why += f"; Hi-C agrees: MQ0 contact with it is {r0:.1f}x the next chromosome"
            else:
                action = "candidate"
                why += (f"; Hi-C does not confirm it (strongest MQ0 contact is with {name(sc0[0][1])}, "
                        f"{r0:.1f}x the next): not moved")
        elif npx[u] >= hic_alone_px and r0 >= hic_alone and not (s is not None and s[0] >= seq_min and s[2] != sc0[0][1]):
            k = sc0[0][1]
            ok, extra = True, ""
            if sgh is not None:
                sch, perh = contact(sgh, u)
                rh = sch[0][0] / max(sch[1][0], 1e-6)
                ok = sch[0][1] == k and rh >= hi_margin
                extra = f", and {rh:.1f}x on the high-MQ map"
            if ok:
                action = "unloc"
                why = f"Hi-C alone: MQ0 contact with the chromosome {name(k)} is {r0:.1f}x the next chromosome{extra}"
        if action is None:
            continue
        if action == "unloc" and rept is not None and not seq_ok:
            # Hi-C alone on a scaffold that is nearly all repeat: its contacts may come from the
            # repeat's other copies. With sequence evidence the decision already rests on its unique
            # (non-repetitive) alignments, and repeat content does not count (curator, 2026-11-14).
            rf = float(rept[start[u]:start[u] + npx[u]].mean())
            if rf >= repeat_frac:
                action = "candidate"
                why += f"; but {100 * rf:.0f}% of it is repeat (repeat track) and only Hi-C places it: not moved"
        if action == "unloc" and covn is not None:
            cv = float(np.median(covn[start[u]:start[u] + npx[u]]))
            if cv >= repeat_cov:
                action = "candidate"
                why += f"; but its coverage is {cv:.1f}x a normal copy (collapsed repeats): not moved"
        a_, b_ = copies[k]
        own = [v for v in (a_, b_) if hap[v] == hap[u]]
        fr = s[3] if seq_ok else {}
        fa, fb = fr.get(a_, 0.0), fr.get(b_, 0.0)
        how, ext = "", None
        one_copy = seq_ok and max(fa, fb) >= seq_min and min(fa, fb) < 0.5 * seq_min
        if one_copy and phasing != "trio" and main[k].kind != "unequal":
            ext = _extends(u, a_ if fa > fb else b_)
        if ext is not None:
            # Its non-repetitive end overlaps an end of the one copy that has that sequence, and the
            # rest of it reaches beyond: it continues that copy, whatever its prefix (curator,
            # 2026-11-14: Hi-C phasing is not perfect; partitioning homologous sequence between the
            # haplotypes is more parsimonious than a segmental duplication against a deletion).
            copy = ext[0]
            how = (f"its non-repetitive end overlaps the {'start' if ext[1] == 'h' else 'end'} of {pt.names[copy]}, the "
                   "only copy with that sequence, and the rest of it reaches beyond: it continues that copy")
        elif main[k].kind == "unequal" or len(own) == 0:
            # sex chromosomes, or neither copy carries its prefix: the copy with the evidence
            copy = max((a_, b_), key=lambda v: (fr.get(v, 0.0), per0[v]))
        elif len(own) == 2:
            # both copies carry its prefix (homologs in one haplotype file): the shorter one, to
            # balance the haplotypes (curator, 2026-11-14)
            copy = a_ if npx[a_] < npx[b_] else b_
            how = f"both copies carry its prefix; put with the shorter copy, {pt.names[copy]}, to balance the haplotypes"
        elif phasing == "trio" or not seq_ok:
            copy = own[0]
        elif max(fa, fb) >= seq_min and min(fa, fb) < 0.5 * seq_min:
            # in one copy only: it belongs to the copy that lacks it
            copy = a_ if fa < fb else b_
            how = f"only {pt.names[b_ if copy == a_ else a_]} has this sequence, so it goes to {pt.names[copy]}"
        else:
            # in both copies already. The high-MQ map separates haplotypes; failing that, the
            # shorter copy is the one more likely to lack sequence
            copy = own[0]
            if sgh is not None:
                _, perh = contact(sgh, u)
                hi_, lo_ = max(perh[a_], perh[b_]), min(perh[a_], perh[b_])
                if hi_ >= copy_margin * max(lo_, 1e-6) and hi_ >= 0.5:
                    copy = a_ if perh[a_] > perh[b_] else b_
                    how = f"high-MQ contact is {hi_ / max(lo_, 1e-6):.1f}x stronger with {pt.names[copy]}"
            if not how and npx[a_] != npx[b_]:
                copy = a_ if npx[a_] < npx[b_] else b_
                how = f"no haplotype-specific contact; put with the shorter copy, {pt.names[copy]}"
        if how:
            why += "; " + how
        if hap[copy] != hap[u]:
            why += f" (prefix {hap[u]}_ in a {hap[copy]}_ superscaffold)"
        if action == "unloc" and len(own) == 1 and main[k].kind != "unequal":
            # A piece of this chromosome also touches the copy of its own haplotype (cis contact).
            # Contact with the other haplotype's copy only is the homolog line alone: the piece is
            # then a chromosome of its own whose homolog is joined to this one in the other
            # haplotype, and the place to look is a misjoin there.
            other = b_ if own[0] == a_ else a_
            if per0[own[0]] < cis_frac * per0[other]:
                action = "candidate"
                why += (f"; but it hardly touches {pt.names[own[0]]}, the copy of its own haplotype "
                        f"({per0[own[0]] / max(per0[other], 1e-6):.2f} of its contact with {pt.names[other]}): probably "
                        f"the homolog of part of {pt.names[other]}, which would then be a misjoin; not moved")
        third = bool(seq_ok and all(s[3].get(v, 0.0) >= seq_min for v in (a_, b_)))
        if third:
            why += ("; the sequence is present in both copies of the chromosome already: a third copy, "
                    "tagged FalseDuplicate")
        out.append(dict(u=u, copy=copy, action=action, third_copy=third, evidence=why, extend=ext,
                        score=round(100 * min(s[0], 1), 1) if seq_ok else round(r0, 1), ratio=round(r0, 1)))
    return out


def _adjacency(sg, g, P, rP, Q, rQ, window=100):
    """Contact between the right end of piece P and the left end of piece Q as laid out
    (pieces are inclusive pixel ranges; r = reversed), relative to what adjacent sequence
    shows at the same distances. The outermost pixels are skipped."""
    def win(piece, tail):
        L = piece[1] - piece[0] + 1
        K = min(10, L // 10)
        W = max(1, min(window, L // 2 - K))
        return (piece[1] - K - W + 1, piece[1] - K, K, True) if tail else (piece[0] + K, piece[0] + K + W - 1, K, False)
    a0, a1, k1, t1 = win(P, not rP)
    b0, b1, k2, t2 = win(Q, rQ)
    B = sg.block(a0, a1, b0, b1).astype(np.float64)
    i, j = np.arange(a1 - a0 + 1), np.arange(b1 - b0 + 1)
    i, j = (i[::-1] if t1 else i), (j[::-1] if t2 else j)
    expect = g[np.minimum(i[:, None] + j[None, :] + 1 + k1 + k2, len(g) - 1)]
    return float(B.sum() / max(expect.sum(), 1e-9))


def _homolog_dir(sg, a, b, cols, max_rows=120, k=6):
    """Direction of the homology line of rows a..b within the column ranges `cols` (the other
    copy of the chromosome): +1 if the homolog runs the same way, -1 if reversed, 0 if no
    line stands out. Sums the map along diagonals and anti-diagonals and takes the sharpest
    one (a line is a few pixels wide; cis contact is broad and cancels in the contrast)."""
    L = b - a + 1
    if L > max_rows:
        a += (L - max_rows) // 2
        b = a + max_rows - 1
    W = b - a + 1
    St = sg.strip(a, b).astype(np.float32)
    best = {1: 0.0, -1: 0.0}
    noise = False
    for lo, hi in cols:
        M = hi - lo + 1 - W + 1
        if M <= 4 * k:
            continue
        F, R = np.zeros(M, np.float32), np.zeros(M, np.float32)
        for i in range(W):
            F += St[i, lo + i:lo + i + M]
            R += St[i, lo + W - 1 - i:lo + W - 1 - i + M]
        for d, D in ((1, F / W), (-1, R / W)):
            c = D[k:-k] - 0.5 * (D[:-2 * k] + D[2 * k:])
            best[d] = max(best[d], float(c.max()))
            noise = True
    if not noise:
        return 0
    d = 1 if best[1] >= best[-1] else -1
    if best[d] < 2.0 * best[-d]:                # the other direction shows what "no line" gives here
        return 0
    return d


def alignment_breaks(raw, blocks, min_len=100000, min_px=30, end_bp=300000, repeat_max=0.8):
    """Places where the alignment to the other haplotype stops being colinear: between two
    consecutive long alignments of a scaffold the other side changes scaffold, strand, or
    jumps. Returns [(pixel, kind, text for the left side, text for the right side)].

    At a gap this means the scaffold was joined differently from the other haplotype (a
    misjoin or misorientation in one of the two). Inside a contig it usually marks where
    the other haplotype is in several scaffolds, which says where those belong."""
    bpp = raw.bp_per_pix
    ix = {n: i for i, n in enumerate(raw.names)}
    by = {}
    for bk in blocks:
        if bk.mapq >= 30 and bk.qe - bk.qs >= min_len:
            by.setdefault(bk.q, []).append((bk.qs, bk.qe, bk.t, bk.ts, bk.te, bk.same_strand))
            by.setdefault(bk.t, []).append((bk.ts, bk.te, bk.q, bk.qs, bk.qe, bk.same_strand))
    txt = lambda x: f"{x[2]} {x[3] / 1e6:.2f}-{x[4] / 1e6:.2f} Mb {'+' if x[5] else '-'}"
    lens = dict(zip(raw.names, raw.scaffold_len_bp().tolist()))
    rep = raw.graphs.get("repeat")
    out = []
    for name, L in by.items():
        i = ix.get(name)
        if i is None or raw.scaffold_npix[i] < min_px:
            continue
        chain = []
        for x in sorted(L):
            if not chain or x[1] > chain[-1][1]:            # skip alignments contained in the previous one
                chain.append(x)
        for A, B in zip(chain, chain[1:]):
            if A[2] != B[2]:
                # If each side's alignment runs to the end of its scaffold, the other haplotype is
                # simply in two scaffolds there: nothing disagrees (it is where those two join).
                la, lb = lens.get(A[2], 0), lens.get(B[2], 0)
                left_rest = (la - A[4]) if A[5] else A[3]
                right_rest = B[3] if B[5] else (lb - B[4])
                if left_rest < end_bp and right_rest < end_bp:
                    continue
                kind = "the two sides align to different scaffolds"
            elif A[5] != B[5]:
                kind = "the two sides align in opposite orientations"
            else:
                d = (B[3] - A[4]) if A[5] else (A[3] - B[4])
                if -3e5 <= d <= max(1e6, 3 * abs(B[0] - A[1]) + 3e5):
                    continue
                kind = "the two sides align to distant places"
            q = int(raw.scaffold_start[i] + min((A[1] + B[0]) / 2 / bpp, raw.scaffold_npix[i] - 1))
            if rep is not None and float(rep[max(0, q - 2):q + 3].mean()) >= repeat_max:
                continue                                    # inside a repeat: alignments end there anyway
            out.append((q, kind, txt(A), txt(B)))
    return out


def _block_index(blocks):
    """{scaffold name: [(start, end, other name, other start, other end, same strand)]} with
    every alignment entered from both of its sides."""
    idx = {}
    for bk in blocks:
        if bk.mapq < 30:
            continue
        idx.setdefault(bk.q, []).append((bk.qs, bk.qe, bk.t, bk.ts, bk.te, bk.same_strand))
        idx.setdefault(bk.t, []).append((bk.ts, bk.te, bk.q, bk.qs, bk.qe, bk.same_strand))
    return idx


def _aligned_dir(raw, idx, piece, cols, min_bp=100000, min_frac=0.1, clear=0.7):
    """Direction of a piece relative to the other copy from sequence alignments: +1 if at
    least 70% of its bases aligned to the pixel ranges `cols` are on the same strand, -1 if
    at most 30%, 0 if too little aligns or the strands are mixed."""
    bpp = raw.bp_per_pix
    name = raw.names[raw.orig_scaffold[piece[0]]]
    lo, hi = raw.rel_coord[piece[0]] * bpp, (raw.rel_coord[piece[1]] + 1) * bpp
    targets = {}
    for a, b in cols:
        targets.setdefault(raw.names[raw.orig_scaffold[a]], []).append((raw.rel_coord[a] * bpp, (raw.rel_coord[b] + 1) * bpp))
    tot = same = 0.0
    for s0, e0, other, os_, oe_, st in idx.get(name, []):
        ov = min(e0, hi) - max(s0, lo)
        if ov <= 0 or other not in targets:
            continue
        # the part of the other side that corresponds to the overlap (alignments are long)
        f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
        if not st:
            f0, f1 = 1 - f1, 1 - f0
        o0, o1 = os_ + f0 * (oe_ - os_), os_ + f1 * (oe_ - os_)
        if any(min(o1, t1) - max(o0, t0) > 0 for t0, t1 in targets[other]):
            tot += ov
            same += ov if st else 0.0
    if tot < max(min_bp, min_frac * (hi - lo)):
        return 0
    return 1 if same >= clear * tot else -1 if same <= (1 - clear) * tot else 0


def _project(raw, idx, piece, layout, min_bp=100000, min_frac=0.1, clear=0.7):
    """Where a piece lies along the other copy as it is laid out, from sequence alignments.
    `layout`: the other copy's pieces in order, (first pixel, last pixel, reversed).
    Returns (position in layout pixels, direction +1 / -1, the aligned stretches of the layout)
    or None when too little of the piece aligns there or the strands are mixed."""
    bpp = raw.bp_per_pix
    name = raw.names[raw.orig_scaffold[piece[0]]]
    lo, hi = raw.rel_coord[piece[0]] * bpp, (raw.rel_coord[piece[1]] + 1) * bpp
    ref, off = {}, 0
    for a, b, rev in layout:
        ref.setdefault(raw.names[raw.orig_scaffold[a]], []).append((raw.rel_coord[a] * bpp, (raw.rel_coord[b] + 1) * bpp, off, rev))
        off += b - a + 1
    tot = same = pos = 0.0
    segs = []
    for s0, e0, other, os_, oe_, st in idx.get(name, []):
        ov = min(e0, hi) - max(s0, lo)
        if ov <= 0 or other not in ref:
            continue
        f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
        if not st:
            f0, f1 = 1 - f1, 1 - f0
        o0, o1 = os_ + f0 * (oe_ - os_), os_ + f1 * (oe_ - os_)
        for t0, t1, off_k, rev_k in ref[other]:
            w = min(o1, t1) - max(o0, t0)
            if w <= 0:
                continue
            to_px = lambda x: off_k + ((t1 - x) if rev_k else (x - t0)) / bpp
            mid = (max(o0, t0) + min(o1, t1)) / 2
            pos += w * to_px(mid)
            tot += w
            same += w if st != rev_k else 0.0
            segs.append(tuple(sorted((to_px(max(o0, t0)), to_px(min(o1, t1))))))
    if tot < max(min_bp, min_frac * (hi - lo)):
        return None
    d = 1 if same >= clear * tot else -1 if same <= (1 - clear) * tot else 0
    return (pos / tot, d, segs) if d else None


def arrange_pieces(sg, pieces, init, max_n=9, orient_gain=0.12, guide=None):
    """Order and orientation of the pieces of one chromosome copy that give the most contact
    between neighbouring ends (exact search over all orders and orientations).
    `pieces`: inclusive pixel ranges; `init`: the arrangement to compare with, a list of
    (index, reversed). Returns (best arrangement, its score, score of init, link scores,
    pieces whose orientation was too close to call and was left as it was)."""
    n = len(pieces)
    g = sg.decay_long(230)
    adj = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                for ri in (False, True):
                    for rj in (False, True):
                        adj[(i, ri, j, rj)] = _adjacency(sg, g, pieces[i], ri, pieces[j], rj)
    total = lambda arr: sum(adj[(a, ra, b, rb)] for (a, ra), (b, rb) in zip(arr, arr[1:]))
    if n > max_n or n < 2:
        return list(init), total(init), total(init), adj, [], []
    best = {}                                   # (mask, last, reversed) -> (score, previous state)
    for i in range(n):
        for r in (False, True):
            best[(1 << i, i, r)] = (0.0, None)
    for mask in range(1, 1 << n):
        for i in range(n):
            if not mask >> i & 1:
                continue
            for r in (False, True):
                cur = best.get((mask, i, r))
                if cur is None:
                    continue
                for j in range(n):
                    if mask >> j & 1:
                        continue
                    for rj in (False, True):
                        key = (mask | 1 << j, j, rj)
                        v = cur[0] + adj[(i, r, j, rj)]
                        if key not in best or v > best[key][0]:
                            best[key] = (v, (mask, i, r))
    full = (1 << n) - 1
    end = max(((best[(full, i, r)][0], i, r) for i in range(n) for r in (False, True)))
    arr, state = [], (full, end[1], end[2])
    while state is not None:
        arr.append((state[1], state[2]))
        state = best[state][1]
    arr.reverse()
    big = max(range(n), key=lambda i: pieces[i][1] - pieces[i][0])
    if dict(arr)[big] != dict(init)[big]:       # keep the largest piece the way it was
        arr = [(i, not r) for i, r in reversed(arr)]
    # A piece whose two orientations score about the same: end contact cannot orient it.
    # Then the other haplotype is the guide (`guide(i)`: +1 / -1 / 0, direction of piece i's
    # homology line in the other copy): the piece is turned so that its homolog runs the same
    # way as that of the nearest piece whose orientation is clear. Without a guide it keeps
    # the orientation it had.
    r_of = dict(arr)
    flip_cost, guided = {}, []
    for k_, (i, r) in enumerate(arr):
        alt = arr[:k_] + [(i, not r)] + arr[k_ + 1:]
        flip_cost[i] = total(arr) - total(alt)
    unclear = [i for i, _ in arr if flip_cost[i] < orient_gain]
    if unclear and len(unclear) == n:           # nothing is clear: the largest piece is the reference
        unclear.remove(big)
    order_ = [i for i, _ in arr]
    for i in list(unclear):
        ref = min((j for j in order_ if j not in unclear), key=lambda j: abs(order_.index(j) - order_.index(i)), default=None)
        di, dj = (guide(i), guide(ref)) if guide is not None and ref is not None else (0, 0)
        if di and dj:
            r_of[i] = r_of[ref] ^ (di != dj)
            unclear.remove(i)                   # settled by the homolog
            guided.append(i)
        else:
            r_of[i] = dict(init)[i]
    arr = [(i, r_of[i]) for i in order_]
    if dict(arr)[big] != dict(init)[big]:
        arr = [(i, not r) for i, r in reversed(arr)]
    return arr, total(arr), total(init), adj, unclear, guided


def merge_pairs(pt, sg, pairs, hap, min_px=60, link_min=0.5, margin=1.2, f=8, specific=3.0, debug=None):
    """Homolog pairs that are parts of one chromosome: an end of one pair's copy touches an
    end of the other pair's copy at the level of adjacent sequence, and the corresponding
    ends of their two homologs do the same (the join is seen in both haplotypes, which
    telomere clustering or repeats shared between chromosomes do not produce). No end may carry
    a telomere signal, and the contact must be `specific` times what either pair shows with its
    third-best partner (a piece of a chromosome has at most two neighbours). Works on the
    map binned by f pixels. Returns [(score, host pair, member pair, x, y)]: the member's
    unit y joins the host's unit x; each pair end is used once."""
    M = sg.binned(f)
    nb = M.shape[0]
    raw = getattr(pt, "pt", pt)
    scb = raw.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, raw.npix - 1)]
    D = 40
    g = np.array([np.diagonal(M, d)[scb[:nb - d] == scb[d:]].mean() for d in range(D)])
    start, npx = pt.scaffold_start, pt.scaffold_npix

    def win(u, tail):
        a, b = -(-int(start[u]) // f), int(start[u] + npx[u]) // f - 1
        n = b - a + 1
        K = 1 if n >= 6 else 0
        W = max(1, min(12, n // 2 - K))
        return (b - K - W + 1, b - K, K, True) if tail else (a + K, a + K + W - 1, K, False)

    def adj(x, ex, y, ey):
        a0, a1, k1, t1 = win(x, ex == "t")
        b0, b1, k2, t2 = win(y, ey == "t")
        i, j = np.arange(a1 - a0 + 1), np.arange(b1 - b0 + 1)
        i, j = (i[::-1] if t1 else i), (j[::-1] if t2 else j)
        expect = g[np.minimum(i[:, None] + j[None, :] + 1 + k1 + k2, D - 1)]
        return float(M[a0:a1 + 1, b0:b1 + 1].sum() / max(expect.sum(), 1e-9))

    other = {"h": "t", "t": "h"}
    main = [p for p in pairs if p.kind == "pair" and min(npx[p.a], npx[p.b]) >= min_px and hap[p.a] != hap[p.b]]
    tel, telomeric = pt.graphs.get("telomere"), set()
    if tel is not None:
        thr = max(50.0, 5 * float(np.percentile(tel, 99)))
        for P in main:
            for u in (P.a, P.b):
                a0, n = int(start[u]), int(npx[u])
                if tel[a0:a0 + 10].max() >= thr:
                    telomeric.add((u, "h"))
                if tel[a0 + n - 10:a0 + n].max() >= thr:
                    telomeric.add((u, "t"))
    allc = []
    for i, P in enumerate(main):
        for Q in main[i + 1:]:
            a, b = P.a, P.b
            c, d = (Q.a, Q.b) if hap[Q.a] == hap[a] else (Q.b, Q.a)
            combos = []
            for ea in "ht":
                for ec in "ht":
                    eb = ea if P.same_strand else other[ea]
                    ed = ec if Q.same_strand else other[ec]
                    v = min(adj(a, ea, c, ec), adj(b, eb, d, ed))
                    combos.append((v, ea, ec, bool({(a, ea), (c, ec), (b, eb), (d, ed)} & telomeric)))
            combos.sort(key=lambda x: -x[0])
            allc.append((combos[0][0], id(P), combos[0][1], id(Q), combos[0][2], P, Q, a, c, combos[0][3], combos[1][0]))
    # what each pair shows with its third-best partner: a piece of a chromosome has at most two
    # neighbours, so anything beyond the second is the level of unspecific contact for that pair
    by_pair = {}
    for x in allc:
        by_pair.setdefault(x[1], []).append(x[0])
        by_pair.setdefault(x[3], []).append(x[0])
    third = {k: (sorted(v, reverse=True)[2] if len(v) > 2 else 0.0) for k, v in by_pair.items()}
    if debug is not None:
        debug.extend((x[0], x[10], third[x[1]], third[x[3]], x[9], pt.names[x[7]], pt.names[x[8]]) for x in allc)
    cand = []
    for x in allc:
        v, ip, ea, iq, ec, P, Q, a, c, at_tel, second = x
        # Adjacent ends: no telomere signal, and far more contact than the pairs show with
        # partners they cannot be next to. Small chromosomes that cluster in the nucleus (bird
        # microchromosomes) touch many others at a similar level.
        if at_tel or v < specific * max(third[ip], third[iq]):
            continue
        cand.append((v, ip, ea, iq, ec, P, Q, a, c))
    cand.sort(key=lambda x: -x[0])
    used, group, out = set(), {}, []
    for k, (v, ip, ea, iq, ec, P, Q, a, c) in enumerate(cand):
        if v < link_min:
            break
        if (ip, ea) in used or (iq, ec) in used or group.get(ip, ip) == group.get(iq, iq):
            continue
        rival = max([w for w, jp, fa, jq, fc, *_ in cand[k + 1:] if (jp, fa) in ((ip, ea), (iq, ec)) or (jq, fc) in ((ip, ea), (iq, ec))]
                    + [0.0])
        if rival * margin > v:
            continue
        used.update(((ip, ea), (iq, ec)))
        ga, gb = group.get(ip, ip), group.get(iq, iq)
        for key in [x for x, gid in group.items() if gid == gb] + [iq]:
            group[key] = ga
        group[ip] = ga
        out.append((v, P, Q, a, c))
    return out


def curate_allhaps(pt, phasing="hic", hap_regex=r"^(H\d+)_", hi=None, do_cuts=True, blocks=None, busco_hits=None,
                   small_blocks=None, do_place=True, level=0.5, do_arrange=True):
    """allhaps curation so far: cuts at gaps with reduced contact and at sharp breaks, homolog
    pairing, parallel layout, HAP tags, unlocs. `hi` is the high-MQ map of the same assembly,
    if there is one. `level` (0 conservative .. 1 aggressive, see `Level`) moves every
    act-or-flag threshold together.
    Returns (rows for the decisions table, savestate, info)."""
    raw = pt
    lv = Level(level)
    margin = lv.at(1.5, 1.3, 1.15)
    first_pairs, _ = pair_homologs(raw, hap_regex, blocks=blocks, margin=margin)
    sites = telomere_sites(raw, first_pairs, hi, hap_regex, blocks, gap_max=lv.at(0.45, 0.6, 0.75),
                           sharp_max=lv.at(0.05, 0.1, 0.2), tip_max=lv.at(0.35, 0.5, 0.65)) if do_cuts else []
    cuts = [(x["q"], x["s"], x["reason"]) for x in sites if x["action"] == "cut"]
    flips = [x["tip"] for x in sites if x["action"] == "invert_tip"]
    pt = SegmentedMap(raw, [q for q, _, _ in cuts])
    pairs, info = pair_homologs(pt, hap_regex, blocks=blocks, margin=margin)
    if small_blocks is not None:
        pairs += same_prefix_pairs(pt, info, small_blocks)
    small_calls = place_small(pt, pairs, info, small_blocks, hi, phasing=phasing, seq_min=lv.at(0.6, 0.5, 0.4),
                              hic_agree=lv.at(2.0, 1.5, 1.2), hic_alone=lv.at(3.5, 2.5, 2.0)) if do_place else []
    attach = {}
    for c in small_calls:
        if c["action"] == "unloc":
            attach.setdefault(c["copy"], []).append(c)
    hap = info["hap"]
    # Order and orientation within each chromosome copy: its pieces between the sites that
    # were kept or tip-inverted, plus attached pieces large enough to carry an end signal.
    npx_, start_ = pt.scaffold_npix, pt.scaffold_start
    sg_arr = _signal(hi) if hi is not None else _signal(raw)
    soft = sorted(x["q"] for x in sites if x["action"] != "cut")
    tips = {tuple(x["tip"]) for x in sites if x["action"] == "invert_tip"}
    chain_px, gain_min = lv.at(60, 40, 20), lv.at(0.3, 0.15, 0.05)
    mate_of = {}
    for p in pairs:
        if p.kind != "extra":
            mate_of[p.a], mate_of[p.b] = p.b, p.a
    extra_of = {}
    for p in pairs:
        if p.kind == "extra" and p.a in mate_of:
            extra_of.setdefault(mate_of[p.a], []).append(p.b)
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
        seen = lambda b: _union([(max(b.qs, x), min(b.qe, y)) for x, y in have.get((b.q, b.t), [])
                                 if min(b.qe, y) > max(b.qs, x)]) >= 0.5 * (b.qe - b.qs)
        more = [b for b in small_blocks if b.primary and b.mapq >= 30 and pre(b.q) != pre(b.t) and not seen(b)]
        for name, v in _block_index(more).items():
            block_idx.setdefault(name, []).extend(v)
    # Pairs that are parts of one chromosome become one: the smaller pair's two units go into
    # the chains of the larger pair's two copies (same haplotype prefix on each side).
    host_of = {}
    if do_arrange:
        size = lambda P: npx_[P.a] + npx_[P.b]
        root = {}
        find = lambda P: find(root[id(P)]) if id(P) in root else P
        for v, P, Q, x, y in merge_pairs(pt, sg_arr, pairs, hap, link_min=lv.at(0.7, 0.5, 0.3),
                                         specific=lv.at(4.0, 3.0, 2.0)):
            P, Q = find(P), find(Q)
            if P is Q:
                continue
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
                partner=pt.names[Q.b], strand="", tag="", score=round(v, 2), next_best="",
                evidence=f"this pair and the pair {pt.names[P.a]} / {pt.names[P.b]} are parts of one chromosome: their ends "
                         f"touch at {v:.2f} of adjacent sequence in both haplotypes (the weaker of the two); "
                         "laid out and painted as one chromosome"))
        for u_, h_ in host_of.items():
            if attach.get(u_):
                attach.setdefault(h_, []).extend(attach.pop(u_))
            if extra_of.get(u_):
                extra_of.setdefault(h_, []).extend(extra_of.pop(u_))
    members = {}
    for u_, h_ in host_of.items():
        members.setdefault(h_, []).append(u_)
    span = lambda a, b: (f"{raw.names[raw.orig_scaffold[a]]} {raw.rel_coord[a] * raw.bp_per_pix / 1e6:.1f}-"
                         f"{(raw.rel_coord[b] + 1) * raw.bp_per_pix / 1e6:.1f} Mb")
    for u in (sorted(mate_of, key=lambda u: start_[u]) if do_arrange else []):
        if u in host_of:
            continue
        units = [u] + sorted(members.get(u, []), key=lambda v: -npx_[v])
        pieces = []
        for v in units:
            a, b = int(start_[v]), int(start_[v] + npx_[v] - 1)
            lo = a
            for q in [q for q in soft if a <= q < b]:
                if q - lo >= 4 and b - q > 4:              # no pieces of a few pixels
                    pieces.append((lo, q))
                    lo = q + 1
            pieces.append((lo, b))
        others = sorted([c["u"] for c in attach.get(u, []) if npx_[c["u"]] >= chain_px] +
                        [v for v in extra_of.get(u, []) if npx_[v] >= chain_px], key=lambda v: -npx_[v])
        others = others[:max(0, 9 - len(pieces))]              # the search is exact, so at most 9 pieces
        pieces += [(int(start_[v]), int(start_[v] + npx_[v] - 1)) for v in others]
        forced = len(units) > 1                                # a merged chromosome is laid out as one in any case
        sol[u] = dict(pieces=pieces, init=[(i, pieces[i] in tips) for i in range(len(pieces))], units=units,
                      others=others, forced=forced, adj=None, row=None, best=0.0, was=0.0)
        if len(pieces) < 2:
            continue
        init = [(i, pieces[i] in tips) for i in range(len(pieces))]
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
        arr, best, was, adj, unclear, guided = arrange_pieces(sg_arr, pieces, init, orient_gain=lv.at(0.2, 0.12, 0.05),
                                                      guide=guide)
        sol[u].update(adj=adj, best=best, was=was)
        if not forced and (arr == init or best - was < gain_min):
            continue
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
                     + f" (sum {best:.2f}) against {was:.2f} as assembled"
                     + ("; oriented by the other haplotype ("
                        + ("sequence alignment" if any(i in by_seq for i in guided) else "Hi-C homology line") + "): "
                        + ", ".join(span(*pieces[i]) for i in guided) if guided else "")
                     + ("; orientation too close to call, left as it was, for: "
                        + ", ".join(span(*pieces[i]) for i in unclear) if unclear else "")))
    # With sequence alignments the other haplotype is the guide for order and orientation.
    # Of the two copies of a chromosome the one in fewer pieces is the reference (arranged
    # by end contact above); the pieces of the other are put in the order and orientation in
    # which they align along it. End contact is the check: if it clearly prefers another
    # arrangement, that may be a real difference between the haplotypes, so the contact
    # arrangement stays and the row says so.
    pair_of = {}
    for p in pairs:
        if p.kind not in ("extra", "merged"):
            pair_of[p.a] = pair_of[p.b] = p
    conflict_min = lv.at(0.3, 0.5, 0.8)
    # Scaffolds that are in no pair but align, over at least half their length, to one copy of
    # a chromosome: they belong in the chain of the other copy (their own haplotype), at the
    # place where they align.
    loose, false_dup = {}, {}
    if block_idx is not None:
        ridx = {n: i for i, n in enumerate(raw.names)}

        def unit_at(name, bp):
            i = ridx.get(name)
            if i is None or raw.scaffold_npix[i] == 0:
                return None
            return int(pt.orig_scaffold[int(raw.scaffold_start[i]) + min(int(bp / raw.bp_per_pix), int(raw.scaffold_npix[i]) - 1)])
        acc = {}
        for bk in blocks:
            if bk.mapq < 30:
                continue
            x, y = unit_at(bk.q, (bk.qs + bk.qe) / 2), unit_at(bk.t, (bk.ts + bk.te) / 2)
            if x is None or y is None or x == y:
                continue
            acc.setdefault(x, {}).setdefault(y, 0)
            acc[x][y] += bk.qe - bk.qs
            acc.setdefault(y, {}).setdefault(x, 0)
            acc[y][x] += bk.te - bk.ts
        min_loose = lv.at(8, 4, 2)
        for v, d in acc.items():
            if v in mate_of or v in host_of or npx_[v] < min_loose:
                continue
            y = max(d, key=d.get)
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
            tot = cov = 0.0
            for l, h in segs:
                tot += h - l
                cov += _union([(max(l, a), min(h, b)) for a, b in have if min(h, b) > max(l, a)])
            return cov / max(tot, 1e-9)
        keep = [i for i in range(len(P)) if proj[i] is not None and (i < len(base) or covered(proj[i][2]) < 0.5)]
        for i in range(len(base), len(P)):
            if proj[i] is not None and i not in keep:
                false_dup[extra_u[i - len(base)]] = (copy, covered(proj[i][2]))
        arr = [(i, proj[i][1] < 0) for i in sorted(keep, key=lambda i: proj[i][0])]
        if not p.same_strand and p.kind == "pair":         # the second copy is shown reversed as a whole
            arr = [(i, not r) for i, r in reversed(arr)]
        return P, arr, {i: extra_u[i - len(base)] for i in keep if i >= len(base)}

    def join_rows(added, copy, ref):
        for v in added.values():
            joined[v] = copy
            if v not in seen_small:
                arr_rows.append(dict(action="join", scaffold=pt.names[v], length_mb=round(npx_[v] * pt.bp_per_pix / 1e6, 2),
                                     partner=pt.names[copy], strand="", tag="", score="", next_best="",
                                     evidence=f"in no pair, but most of it aligns to {pt.names[ref]} (the other haplotype): placed in "
                                              f"the chain of {pt.names[copy]} where it aligns (see the rearrange row)"))

    for u in (sorted(sol, key=lambda u: start_[u]) if block_idx is not None else []):
        p = pair_of.get(u)
        if p is None or u != p.a or p.b not in sol:
            continue
        size = lambda v: sum(b - a + 1 for a, b in sol[v]["pieces"])
        ref, oth = sorted((p.a, p.b), key=lambda v: (len(sol[v]["pieces"]), -size(v)))
        so, sr = sol[oth], sol[ref]
        n0 = len(so["pieces"])
        layout = arranged.get(ref) or [(a, b, r) for (a, b), (_, r) in zip(sr["pieces"], sr["init"])]
        got = along(oth, layout, p, so["pieces"], None)
        if got is not None and len(got[0]) >= 2:
            P, aln, added = got
            link = lambda a, ra, b, rb: so["adj"][(a, ra, b, rb)] if so["adj"] is not None and a < n0 and b < n0 else \
                _adjacency(sg_arr, g230, P[a], ra, P[b], rb)
            base_aln = [x for x in aln if x[0] < n0]
            t_aln = sum(so["adj"][(a, ra, b, rb)] for (a, ra), (b, rb) in zip(base_aln, base_aln[1:])) if so["adj"] is not None else 0.0
            cur = [(so["pieces"].index((a, b)), r) for a, b, r in arranged[oth]] if oth in arranged else so["init"]
            order_txt = ", ".join(span(P[i][0], P[i][1]) + (" reversed" if r else "") for i, r in aln)
            tip_undone = hi is not None and any(P[i] in tips and not r for i, r in aln)
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
                if so["row"] is not None:
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
    for c in small_calls:
        if c["u"] in joined:
            c["action"] = "join"
    for u in attach:
        attach[u] = [c for c in attach[u] if c["u"] not in joined]
    unit_h = "% of the shorter aligned" if blocks is not None else ""
    haps = sorted({h for h in hap if h})
    cov = pt.graphs.get("coverage")
    npx, start = pt.scaffold_npix, pt.scaffold_start
    main = [p for p in pairs if p.kind not in ("extra", "merged")]
    extra = {}
    for p in pairs:
        if p.kind == "extra":
            extra.setdefault(host_of.get(p.a, p.a), []).append(p)
    paired = np.concatenate([np.arange(start[s], start[s] + npx[s]) for p in main for s in (p.a, p.b)]) \
        if main else np.array([], dtype=int)
    norm = float(np.median(cov[paired])) if cov is not None and len(paired) else 0.0
    rows, order, tag_of, paint_of = [], [], {}, {}
    unloc_units, sex_like, in_chain = set(), set(), set()
    look_units = {c["u"] for c in small_calls if c["action"] == "candidate" and not c["third_copy"]}
    fd_units = {c["u"] for c in small_calls if c["third_copy"] and c["action"] != "join"}
    for c in small_calls:
        act = "join" if c["action"] == "join" else "false_duplicate" if c["third_copy"] else \
            "unloc" if c["action"] == "unloc" else "unloc_candidate"
        rows.append(dict(action=act, scaffold=pt.names[c["u"]],
                         length_mb=round(npx[c["u"]] * pt.bp_per_pix / 1e6, 2), partner=pt.names[c["copy"]], strand="",
                         tag={"false_duplicate": "FalseDuplicate", "unloc": "Unloc"}.get(act, ""), score=c["score"], next_best="",
                         evidence=c["evidence"] + ("; " + c["join_note"] if c.get("join_note") else
                                                   "; large enough to place by end contact: joined into the chromosome "
                                                   "(see the rearrange row)" if c["action"] == "join" else "")))
    rows += arr_rows
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
        if s in arranged:
            px = np.concatenate([np.arange(a, b + 1)[::-1] if r else np.arange(a, b + 1) for a, b, r in arranged[s]])
            return px[::-1] if rev else px
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
        if p.kind == "pair" and hemi and min(mb(p.a), mb(p.b)) >= 5:
            # both copies at about half coverage: hemizygous, i.e. X with Y (or Z with W) of similar size
            p.kind, p.note = "unequal", "both copies have about half the coverage of other pairs (hemizygous)"
            sex_like.add(id(p))
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
    sex = [s for s in info["unpaired"] if s not in done and mb(s) >= 5 and covr(s) < 0.7]
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
    rest = [int(s) for s in np.argsort(start) if npx[s] > 0 and int(s) not in placed and int(s) not in joined]
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
        if s in fd_units:
            mask |= 1 << st.tag_id("FalseDuplicate")
        elif s in unloc_units:
            mask |= 1 << st.tag_id("Unloc")
        if c.scaffold in dup_of or s in look_units:
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
