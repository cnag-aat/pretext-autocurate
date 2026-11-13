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
                   tol_frac=0.01, broad_min=0.6, invert_min=1.3, gap_max=0.6, sharp_max=0.1, merge_px=5, tip_max=0.5):
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
        e = out.setdefault(x, {}).setdefault(y, [[], [], 0, 0])
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
                hic_alone=2.5, hic_alone_px=20, hi_margin=2.0, repeat_cov=1.5,
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
    seq = {}
    for u in small:
        per, L = {}, npx[u] * bpp
        for v, e in aln.get(u, {}).items():
            if v in chrom:
                per.setdefault(chrom[v], {})[v] = e[0]
        f = sorted(((_union([i for iv in d.values() for i in iv]) / L, k) for k, d in per.items()), reverse=True)
        if f:
            k = f[0][1]
            seq[u] = (f[0][0], f[1][0] if len(f) > 1 else 0.0, k,
                      {v: _union(iv) / L for v, iv in per[k].items()})
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
            why = f"{100 * min(s[0], 1):.0f}% of it aligns to the chromosome {name(k)} (next chromosome {100 * s[1]:.0f}%)"
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
        if action == "unloc" and covn is not None:
            cv = float(np.median(covn[start[u]:start[u] + npx[u]]))
            if cv >= repeat_cov:
                action = "candidate"
                why += f"; but its coverage is {cv:.1f}x a normal copy (collapsed repeats): not moved"
        a_, b_ = copies[k]
        own = [v for v in (a_, b_) if hap[v] == hap[u]]
        fr = s[3] if seq_ok else {}
        fa, fb = fr.get(a_, 0.0), fr.get(b_, 0.0)
        how = ""
        if main[k].kind == "unequal" or len(own) != 1:
            # sex chromosomes, same-prefix pair: the copy with the evidence
            copy = max((a_, b_), key=lambda v: (fr.get(v, 0.0), per0[v]))
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
            why += ("; the sequence is present in both copies of the chromosome already: a third copy "
                    "(haplotig or false duplication), highlighted")
        out.append(dict(u=u, copy=copy, action=action, third_copy=third, evidence=why,
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
        return list(init), total(init), total(init), adj, []
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
    flip_cost = {}
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
        else:
            r_of[i] = dict(init)[i]
    arr = [(i, r_of[i]) for i in order_]
    if dict(arr)[big] != dict(init)[big]:
        arr = [(i, not r) for i, r in reversed(arr)]
    return arr, total(arr), total(init), adj, unclear


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
    arranged, joined, arr_rows = {}, {}, []
    span = lambda a, b: (f"{raw.names[raw.orig_scaffold[a]]} {raw.rel_coord[a] * raw.bp_per_pix / 1e6:.1f}-"
                         f"{(raw.rel_coord[b] + 1) * raw.bp_per_pix / 1e6:.1f} Mb")
    for u in (sorted(mate_of, key=lambda u: start_[u]) if do_arrange else []):
        a, b = int(start_[u]), int(start_[u] + npx_[u] - 1)
        pieces, lo = [], a
        for q in [q for q in soft if a <= q < b]:
            pieces.append((lo, q))
            lo = q + 1
        pieces.append((lo, b))
        own = len(pieces)
        others = [c["u"] for c in attach.get(u, []) if npx_[c["u"]] >= chain_px] + \
                 [v for v in extra_of.get(u, []) if npx_[v] >= chain_px]
        pieces += [(int(start_[v]), int(start_[v] + npx_[v] - 1)) for v in others]
        if len(pieces) < 2 or len(pieces) > 9 or min(q - p_ for p_, q in pieces) < 4:
            continue
        init = [(i, pieces[i] in tips) for i in range(len(pieces))]
        m_ = mate_of[u]
        cols = [(int(start_[v]), int(start_[v] + npx_[v] - 1)) for v in
                [m_] + [c["u"] for c in attach.get(m_, [])] + extra_of.get(m_, [])]
        memo = {}

        def guide(i, pieces=pieces, cols=cols, memo=memo):
            if i not in memo:
                memo[i] = _homolog_dir(_signal(raw), pieces[i][0], pieces[i][1], cols)
            return memo[i]
        arr, best, was, adj, unclear = arrange_pieces(sg_arr, pieces, init, orient_gain=lv.at(0.2, 0.12, 0.05),
                                                      guide=guide)
        if arr == init or best - was < gain_min:
            continue
        arranged[u] = [(pieces[i][0], pieces[i][1], r) for i, r in arr]
        for v in others:
            joined[v] = u
        links = [adj[(x, rx, y, ry)] for (x, rx), (y, ry) in zip(arr, arr[1:])]
        arr_rows.append(dict(
            action="rearrange", scaffold=pt.names[u], length_mb=round(npx_[u] * pt.bp_per_pix / 1e6, 2), partner="",
            strand="", tag="", score=round(best, 2), next_best=round(was, 2),
            evidence="new order: " + ", ".join(span(pieces[i][0], pieces[i][1]) + (" reversed" if r else "") for i, r in arr)
                     + "; contact between neighbouring ends " + "/".join(f"{x:.2f}" for x in links)
                     + f" (sum {best:.2f}) against {was:.2f} as assembled"
                     + ("; orientation too close to call, left as it was, for: "
                        + ", ".join(span(*pieces[i]) for i in unclear) if unclear else "")))
    for c in small_calls:
        if c["u"] in joined:
            c["action"] = "join"
    for u in attach:
        attach[u] = [c for c in attach[u] if c["u"] not in joined]
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
    unloc_units, sex_like, in_chain = set(), set(), set()
    look_units = {c["u"] for c in small_calls if c["action"] == "candidate" or c["third_copy"]}
    for c in small_calls:
        rows.append(dict(action={"unloc": "unloc", "join": "join"}.get(c["action"], "unloc_candidate"), scaffold=pt.names[c["u"]],
                         length_mb=round(npx[c["u"]] * pt.bp_per_pix / 1e6, 2), partner=pt.names[c["copy"]], strand="",
                         tag="Unloc" if c["action"] == "unloc" else "", score=c["score"], next_best="",
                         evidence=c["evidence"] + ("; large enough to place by end contact: joined into the chromosome "
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
            if e.b in joined:
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
        if s in unloc_units:
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
