"""Small scaffolds: which chromosome and copy (`place_small`): unloc, candidate, join when
it extends a copy, FalseDuplicate when it is a second copy."""
import numpy as np
from .common import _signal, _union


def _small_alignments(pt, small_blocks):
    """{query unit: {target unit: [intervals on the query, intervals with mapq >= 30, bases on
    the same strand, bases]}} from the alignment of small scaffolds against everything.
    A unit's alignments to itself (internal repeats) are left out. Intervals are in bp of
    the original query scaffold. Each entry also has, last, a list of (qs, qe, match density,
    target name, ts, te, same strand) per block."""
    raw = getattr(pt, "pt", pt)
    idx = {n: i for i, n in enumerate(raw.names)}
    bpp = raw.bp_per_pix

    def unit(name, bp):
        """Unit of `pt` holding position bp of an original scaffold; None if not in the map."""
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


def place_small(pt, pairs, info, small_blocks=None, hi=None, seq_min=0.5, seq_margin=2.0, hic_agree=1.5,
                hic_alone=2.5, hic_alone_px=20, hi_margin=2.0, repeat_cov=1.5, dense_min=0.5, repeat_frac=0.9,
                phasing="hic", copy_margin=1.5, cis_frac=0.2, hi_unique=5.0, mq0_unique=6.0, unique_floor=0.5):
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
    When its non-repetitive end overlaps an end of the only copy with that sequence and the
    rest reaches beyond, it continues that copy (`extend`; joined there in `curate_allhaps`).
    Sequence already in both copies makes it a third copy (`third_copy`, FalseDuplicate).

    Thresholds: `seq_min` share of the scaffold aligned to the chromosome and `seq_margin`
    over the next chromosome; `hic_agree` / `hic_alone` MQ0 contact ratio best / next
    chromosome with / without sequence (`hic_alone_px` smallest scaffold for Hi-C alone,
    `hi_margin` the same ratio on the high-MQ map); `repeat_cov` normalised coverage and
    `repeat_frac` repeat-track share that make a Hi-C-only call a candidate; `dense_min`
    match density of a non-repetitive alignment; `copy_margin` high-MQ contact ratio to pick
    a copy; `cis_frac` least contact with its own haplotype's copy, relative to the other.

    Pieces smaller than `hic_alone_px` (curator, 2026-11-16): unlocalised when their signal
    on the high-MQ map is unique to one homolog pair (best pair `hi_unique` times the next,
    and at least `unique_floor` mean contact, so that a handful of reads does not count);
    unique on the MQ0 map as well is stronger evidence, and without a high-MQ map it is the
    evidence (`mq0_unique` times the next pair, same floor).

    Returns a list of dicts (u, copy, action, third_copy, evidence, extend, score, ratio)."""
    npx, start, hap = pt.scaffold_npix, pt.scaffold_start, info["hap"]
    main = [p for p in pairs if p.kind != "extra"]
    # chromosome k = main pair k; chrom: unit -> k, copies: k -> its two units
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
    # seq[u] = (best chromosome's share, next chromosome's share, best k, {copy: share});
    # nonrep[u] = (non-repetitive bp, length bp) when the shares are of non-repetitive sequence
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
        # enough dense sequence to judge on (50 kb or 10%), but not nearly all of it: a scaffold
        # that is all dense is judged on all its alignments below
        if Ld >= max(50000, 0.1 * L) and Ld < 0.9 * L:
            top = max(d for _, _, d in dense)

            def within(iv):
                """Bases of intervals iv that fall in the dense (non-repetitive) stretches D."""
                return _union([(max(a, x), min(b, y)) for a, b in iv for x, y in D if min(b, y) > max(a, x)])
            # per copy, alignments nearly as dense as u's best (70%): its true homolog, not a
            # looser paralog
            for v, e in aln.get(u, {}).items():
                if v in chrom:
                    per.setdefault(chrom[v], {})[v] = [(a, b) for a, b, d, *_ in e[4] if d >= 0.7 * top]
            f = sorted(((within([i for iv in d.values() for i in iv]) / Ld, k) for k, d in per.items()), reverse=True)
            if f and f[0][0] > 0:
                k = f[0][1]
                seq[u] = (f[0][0], f[1][0] if len(f) > 1 else 0.0, k, {v: within(iv) / Ld for v, iv in per[k].items()})
                nonrep[u] = (Ld, L)
            continue
        # otherwise: share of the whole scaffold aligned to each chromosome (any density)
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
        # "at an end": within 3 px or 5% of u (bp), within 3 px or 3% of X (pixels)
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

    # every small scaffold: sequence, Hi-C alone (large enough), or a unique signal (small)
    need = list(small)
    name = lambda k: "/".join(pt.names[v] for v in copies[k])
    cs = np.concatenate(([0.0], np.zeros(pt.npix)))
    units = sorted(chrom)

    def contact(sg, u):
        """Mean contact of u with each paired unit on map sg (per pixel of that unit).
        Returns (chromosomes sorted by their better copy, [(value, k)]; {unit: value})."""
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
        r0 = sc0[0][0] / max(sc0[1][0] if len(sc0) > 1 else 0.0, 1e-6)      # MQ0 contact, best chromosome over the next
        s = seq.get(u)
        seq_ok = s is not None and s[0] >= seq_min and s[0] >= seq_margin * s[1]
        k, action, why, small_route = None, None, "", False
        # which chromosome: sequence first (unloc if Hi-C agrees, else candidate); Hi-C alone
        # only when strong, not contradicted by sequence, and confirmed on the high-MQ map
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
                rh = sch[0][0] / max(sch[1][0] if len(sch) > 1 else 0.0, 1e-6)
                ok = sch[0][1] == k and rh >= hi_margin
                extra = f", and {rh:.1f}x on the high-MQ map"
            if ok:
                action = "unloc"
                why = f"Hi-C alone: MQ0 contact with the chromosome {name(k)} is {r0:.1f}x the next chromosome{extra}"
        elif npx[u] < hic_alone_px and not (s is not None and s[0] >= seq_min):
            # too small to weigh contact ratios: unlocalised only where its signal is unique to
            # one homolog pair, on the high-MQ map (curator), or on the MQ0 map without one
            small_route = True
            u0 = r0 >= mq0_unique and sc0[0][0] >= unique_floor
            if sgh is not None:
                sch, _ = contact(sgh, u)
                rh = sch[0][0] / max(sch[1][0] if len(sch) > 1 else 0.0, 1e-6)
                if rh >= hi_unique and sch[0][0] >= unique_floor and (sch[0][1] == sc0[0][1] or not u0):
                    k = sch[0][1]
                    action = "unloc"
                    why = (f"small ({npx[u]} px): its high-MQ signal is unique to the chromosome {name(k)} "
                           f"({min(rh, 99):.0f}x the next pair)"
                           + (f", and so is its MQ0 signal ({min(r0, 99):.0f}x)" if u0 and sc0[0][1] == k else ""))
            elif u0:
                k = sc0[0][1]
                action = "unloc"
                why = (f"small ({npx[u]} px): its MQ0 signal is unique to the chromosome {name(k)} "
                       f"({min(r0, 99):.0f}x the next pair; no high-MQ map)")
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
        if small_route and action != "unloc":
            # the small-piece rule either places it or says nothing: no candidate rows for
            # repeat-like pieces, and the other routes (alignment into a copy) stay open
            continue
        # which copy of chromosome k (rules in the docstring); own: the copies with u's prefix
        a_, b_ = copies[k]
        own = [v for v in (a_, b_) if hap[v] == hap[u]]
        fr = s[3] if seq_ok else {}
        fa, fb = fr.get(a_, 0.0), fr.get(b_, 0.0)
        how, ext = "", None
        # in one copy only: aligned well to one, under half of seq_min to the other
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
            copy = own[0]                          # trio phasing is trusted; Hi-C alone cannot pick
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
        # a second copy of sequence in the same haplotype is FalseDuplicate (curator rule)
        third = bool(seq_ok and all(s[3].get(v, 0.0) >= seq_min for v in (a_, b_)))
        if third:
            why += ("; the sequence is present in both copies of the chromosome already: a third copy, "
                    "tagged FalseDuplicate")
        out.append(dict(u=u, copy=copy, action=action, third_copy=third, evidence=why, extend=ext,
                        score=round(100 * min(s[0], 1), 1) if seq_ok else round(r0, 1), ratio=round(r0, 1)))
    return out
