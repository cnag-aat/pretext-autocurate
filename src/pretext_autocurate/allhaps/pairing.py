"""Homolog pairs: `pair_homologs` (alignment blocks or the MQ0 homology line),
`same_prefix_pairs` (homologs in one haplotype file), `merge_pairs` (pairs that are
pieces of one chromosome)."""
import re
import numpy as np
from .common import Pair, _level, _mean0, _union, split_block
from .small import _small_alignments


def _quality(pt, s, tel_thr):
    """Higher is better: fewer gaps, telomeres at both ends, then length.
    Returns (sort key, gap pixels, telomeric ends 0-2); `tel_thr` is in telomere-track units."""
    a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
    gap, tel = pt.graphs.get("gap"), pt.graphs.get("telomere")
    gaps = int((gap[a:a + n] > 0).sum()) if gap is not None else 0
    k = max(1, min(10, n // 4))                     # end window: up to 10 px, a quarter of short ones
    tels = 0 if tel is None else int(tel[a:a + k].max() >= tel_thr) + int(tel[a + n - k:a + n].max() >= tel_thr)
    return (-gaps, tels, n), gaps, tels


def pair_homologs(pt, hap_regex=r"^(H\d+)_", min_bins=4, score_min=None, margin=1.3, unequal=0.5, blocks=None,
                  min_seq_bp=500000, repeat_cov=1.3):
    """Returns (pairs, info). `pt` must be the MQ0 map.

    With `blocks` (sequence alignments between haplotypes, see homology.py) homology is
    taken from them: the score of two scaffolds is the percentage of the shorter one that
    aligns to the other, ranking is by aligned bases, strand by the majority of aligned
    bases. Without blocks it is estimated from the MQ0 Hi-C signal.

    `pt` may be a SegmentedMap (cut pieces count as units). Thresholds: `min_bins` smallest
    unit in matrix bins (Hi-C only); `min_seq_bp` smallest unit in bp (alignments only);
    `score_min` homology floor (percent aligned, or map units; per-scaffold if None);
    `margin` how far the best partner must beat a rival; `unequal` length ratio below which
    a pair is 'unequal' (possible sex chromosomes); `repeat_cov` normalised long-read
    coverage above which a scaffold is too repetitive to pair."""
    # one whole-map matrix at a coarse mip level, cached on the map (reused across sweep levels)
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
        """Matrix bins lying wholly inside unit s (pixels rounded inwards to bins of f px)."""
        a, b = int(pt.scaffold_start[s]), int(pt.scaffold_start[s] + pt.scaffold_npix[s])
        return slice(-(-a // f), b // f)

    bins = lambda s: max(1, rng(s).stop - rng(s).start)
    line, strand, total = {}, {}, {}
    if blocks is not None:
        raw = getattr(pt, "pt", pt)
        idx = {n: i for i, n in enumerate(raw.names)}
        bpp = raw.bp_per_pix

        def unit(name, bp):
            """Unit of `pt` (a cut piece if `pt` is segmented) holding position bp of a scaffold
            of the original assembly; None if the name is not in the map."""
            i = idx.get(name)
            if i is None or raw.scaffold_npix[i] == 0:
                return None
            px = int(raw.scaffold_start[i]) + min(int(bp / bpp), int(raw.scaffold_npix[i]) - 1)
            return int(pt.orig_scaffold[px])

        # per unit pair: aligned bp at mapq >= 30 (and how much on the same strand) for the
        # score, and bp at any mapq for `loose`; a block is shared out over the units it covers
        # (split_block: one long block can run across cuts)
        acc, any_bp = {}, {}
        for bk in blocks:
            for x, y, bp in split_block(bk, unit, bpp):
                if x is None or y is None or x == y:
                    continue
                key = (min(x, y), max(x, y))
                any_bp[key] = any_bp.get(key, 0) + bp          # every mapping quality
                if bk.mapq < 30:
                    continue
                a_ = acc.setdefault(key, [0, 0])
                a_[0] += bp
                a_[1] += bp if bk.same_strand else 0
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
    # greedy, largest first: a pair is two units that are each other's best partner, above
    # both floors, with no rival within `margin`
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
            """Total homology of u's strongest competitor to partner v, co-partners excluded."""
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


def same_prefix_pairs(pt, info, small_blocks, min_bp=500000, frac_min=0.5):
    """Homologs that ended up in the same haplotype file (both `H1_`, say): two unpaired
    scaffolds of one prefix that are each other's best match, with at least half of the
    shorter one aligned uniquely to the other. The H2-against-H1 alignment cannot show these.
    `min_bp`: smallest scaffold considered; `frac_min`: share of the shorter one that must
    align at mapq >= 30. Returns Pair objects (score: that share, in percent)."""
    aln = _small_alignments(pt, small_blocks)
    bpp, hap, npx = pt.bp_per_pix, info["hap"], pt.scaffold_npix
    free = set(info["unpaired"]) | set(info.get("repetitive", {}))
    free = {u for u in free if npx[u] * bpp >= min_bp}

    def best(u):
        """The unit u has most aligned bases to (any mapq), or None."""
        t = {v: _union(e[0]) for v, e in aln.get(u, {}).items()}
        return max(t, key=t.get) if t else None

    out = []
    for x in sorted(free, key=lambda u: -npx[u]):
        y = best(x)
        # mutual best, same prefix, and each pair entered once (x the larger)
        if y is None or y not in free or hap[y] != hap[x] or best(y) != x or (npx[y], y) > (npx[x], x):
            continue
        e = aln[y][x]                          # the shorter one (y) as query: unique share of it
        uniq = _union(e[1]) / (npx[y] * bpp)
        if uniq < frac_min:
            continue
        out.append(Pair(x, y, 100.0 * uniq, 0.0, e[2] >= 0.5 * e[3], "pair",
                        f"both copies carry the prefix {hap[x]}_: homologs that ended up in the same haplotype file"))
    return out


def merge_pairs(pt, sg, pairs, hap, min_px=60, link_min=0.5, margin=1.2, f=8, specific=3.0, debug=None,
                link_weak=None, block_min=None, graph_links=()):
    """Homolog pairs that are parts of one chromosome: an end of one pair's copy touches an
    end of the other pair's copy at the level of adjacent sequence, and the corresponding
    ends of their two homologs do the same (the join is seen in both haplotypes, which
    telomere clustering or repeats shared between chromosomes do not produce). No end may carry
    a telomere signal, and the contact must be `specific` times what either pair shows with its
    third-best partner (a piece of a chromosome has at most two neighbours). Works on the
    map binned by f pixels. Returns [(score, host pair, member pair, x, y)]: the member's
    unit y joins the host's unit x; each pair end is used once.

    `min_px`: smallest copy considered (pixels); `link_min`: least end contact (fraction of
    adjacent sequence, the weaker haplotype); `margin`: how far a link must beat any other
    candidate link at the same ends.

    Chromosome-level contact (Troglobisium: joins the curator found obvious had end contact
    0.25-0.5, while their whole blocks showed 7-36x the trans background, against 3.5x for a
    false candidate): with `block_min`, the mean contact over the whole block of two copies,
    over the mean of each copy's blocks with the other copies of its haplotype, must reach
    `block_min` in both haplotypes for a link that is under `link_min` but at least
    `link_weak` (still `specific`). `graph_links`: {frozenset({(unit, end), (unit, end)})},
    copy ends that hifiasm's graph joins (directly or through one small scaffold); a pair of
    pairs whose ends it joins in either haplotype is merged at those ends when the
    chromosome-level test passes, whatever the end contact (the graph link is unique).
    Returns [(score, host pair, member pair, x, y, how)]."""
    M = sg.binned(f)
    nb = M.shape[0]
    raw = getattr(pt, "pt", pt)
    scb = raw.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, raw.npix - 1)]
    D = 40
    # expected contact at each distance (bins) within one scaffold of the original assembly
    g = np.array([_mean0(np.diagonal(M, d)[scb[:nb - d] == scb[d:]]) for d in range(D)])
    start, npx = pt.scaffold_start, pt.scaffold_npix

    def win(u, tail):
        """End window of unit u in bins: (first, last, skipped bins K, at the tail); up to 12
        bins, skipping the outermost bin when the unit is long enough."""
        a, b = -(-int(start[u]) // f), int(start[u] + npx[u]) // f - 1
        n = b - a + 1
        K = 1 if n >= 6 else 0
        W = max(1, min(12, n // 2 - K))
        return (b - K - W + 1, b - K, K, True) if tail else (a + K, a + K + W - 1, K, False)

    def adj(x, ex, y, ey):
        """Contact between end ex of unit x and end ey of unit y ('h' head, 't' tail), relative
        to what adjacent sequence shows at the same distances (as `_adjacency`, binned)."""
        a0, a1, k1, t1 = win(x, ex == "t")
        b0, b1, k2, t2 = win(y, ey == "t")
        i, j = np.arange(a1 - a0 + 1), np.arange(b1 - b0 + 1)
        i, j = (i[::-1] if t1 else i), (j[::-1] if t2 else j)
        expect = g[np.minimum(i[:, None] + j[None, :] + 1 + k1 + k2, D - 1)]
        return float(M[a0:a1 + 1, b0:b1 + 1].sum() / max(expect.sum(), 1e-9))

    other = {"h": "t", "t": "h"}
    glinks = set(graph_links)
    main = [p for p in pairs if p.kind == "pair" and min(npx[p.a], npx[p.b]) >= min_px and hap[p.a] != hap[p.b]]
    # ends with a telomere signal in their last 10 px: a chromosome ends there, nothing joins
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
    # allc: per pair of pairs, the best of the four end combinations: (score, id P, end of a,
    # id Q, end of c, P, Q, a, c, an end is telomeric, second-best combination's score)
    allc = []
    for i, P in enumerate(main):
        for Q in main[i + 1:]:
            a, b = P.a, P.b                     # a and c in one haplotype, b and d in the other
            c, d = (Q.a, Q.b) if hap[Q.a] == hap[a] else (Q.b, Q.a)
            combos = []
            for ea in "ht":
                for ec in "ht":
                    # the homologous ends in the other haplotype (swapped if a pair is on
                    # opposite strands); the join must show in both, so the weaker one counts
                    eb = ea if P.same_strand else other[ea]
                    ed = ec if Q.same_strand else other[ec]
                    v1, v2 = adj(a, ea, c, ec), adj(b, eb, d, ed)
                    gl = frozenset({(a, ea), (c, ec)}) in glinks or frozenset({(b, eb), (d, ed)}) in glinks
                    combos.append((min(v1, v2), ea, ec, bool({(a, ea), (c, ec), (b, eb), (d, ed)} & telomeric), v1, v2, gl))
            combos.sort(key=lambda x: -x[0])
            best = next((x for x in combos if x[6]), combos[0])          # the ends the graph joins, if any
            second = next(x[0] for x in combos if x is not best)
            allc.append((best[0], id(P), best[1], id(Q), best[2], P, Q, a, c, best[3], second, best[4], best[5], best[6],
                         b, d))
    # what each pair shows with its third-best partner: a piece of a chromosome has at most two
    # neighbours, so anything beyond the second is the level of unspecific contact for that pair
    by_pair = {}
    for x in allc:
        by_pair.setdefault(x[1], []).append(x[0])
        by_pair.setdefault(x[3], []).append(x[0])
    third = {k: (sorted(v, reverse=True)[2] if len(v) > 2 else 0.0) for k, v in by_pair.items()}
    if debug is not None:
        # (weaker contact, second-best end combination, third-best partners, telomeric, the two
        # units of one haplotype, contact in that haplotype and in the other)
        debug.extend((x[0], x[10], third[x[1]], third[x[3]], x[9], pt.names[x[7]], pt.names[x[8]], x[11], x[12], x[13])
                     for x in allc)
    # chromosome-level contact: whole-block mean of two copies over each one's mean block
    # with the other copies of its haplotype (its trans background)
    units = [u for P in main for u in (P.a, P.b)]
    brng = lambda u: slice(-(-int(start[u]) // f), max(int(start[u] + npx[u]) // f, -(-int(start[u]) // f) + 1))
    bmean, bg_ = {}, {}

    def block(x, y):
        k = (min(x, y), max(x, y))
        if k not in bmean:
            bmean[k] = float(M[brng(x), brng(y)].mean())
        return bmean[k]

    def bg(x):
        if x not in bg_:
            bg_[x] = max(float(np.mean([block(x, y) for y in units if y != x and hap[y] == hap[x]] or [0.0])), 1e-9)
        return bg_[x]
    ratio = lambda x, y: block(x, y) / max(bg(x), bg(y))
    cand = []
    for x in allc:
        v, ip, ea, iq, ec, P, Q, a, c, at_tel, second, _, _, gl, b, d = x
        # Adjacent ends: no telomere signal, and far more contact than the pairs show with
        # partners they cannot be next to. Small chromosomes that cluster in the nucleus (bird
        # microchromosomes) touch many others at a similar level.
        spec = v >= specific * max(third[ip], third[iq])
        how = None
        if spec and not at_tel and v >= link_min:
            how = "end"
        elif block_min is not None and min(ratio(a, c), ratio(b, d)) >= block_min:
            if gl:
                how = "graph"
            elif spec and not at_tel and link_weak is not None and v >= link_weak:
                how = "block"
        if how is None and not (spec and not at_tel):
            continue                    # (a specific link under the thresholds still counts as a rival)
        cand.append((v, ip, ea, iq, ec, P, Q, a, c, how, min(ratio(a, c), ratio(b, d)) if how not in (None, "end") else None))
    # greedy, strongest first (graph links first): each pair end joins at most once, no cycles
    # (group: union of pairs already linked), and a link must beat by `margin` the next
    # candidate at its ends (a graph link is unique by itself)
    cand.sort(key=lambda x: (x[9] != "graph", -x[0]))
    used, group, out = set(), {}, []
    for k, (v, ip, ea, iq, ec, P, Q, a, c, how, br) in enumerate(cand):
        if how is None:
            continue
        if (ip, ea) in used or (iq, ec) in used or group.get(ip, ip) == group.get(iq, iq):
            continue
        rival = max([w for w, jp, fa, jq, fc, *_ in cand[k + 1:] if (jp, fa) in ((ip, ea), (iq, ec)) or (jq, fc) in ((ip, ea), (iq, ec))]
                    + [0.0])
        if how != "graph" and rival * margin > v:
            continue
        used.update(((ip, ea), (iq, ec)))
        ga, gb = group.get(ip, ip), group.get(iq, iq)
        for key in [x for x, gid in group.items() if gid == gb] + [iq]:
            group[key] = ga
        group[ip] = ga
        text = ("" if how == "end" else
                f"; whole-block contact {br:.1f}x the trans background in both haplotypes"
                + ("; hifiasm's graph joins these ends" if how == "graph" else ""))
        out.append((v, P, Q, a, c, text))
    return out
