"""Two copies of one region inside one scaffold (HANDOFF 5: two haplotypes in one scaffold).

From the alignment between the haplotypes (primary alignments): two stretches of the same
scaffold that align to overlapping stretches of one scaffold of the other haplotype. The
scaffold then holds two copies where the other haplotype holds one. That is either

  haplotypic   both alleles in one scaffold: the copy that is the other haplotype's allele
               shares its reads with that allele, so it and the other haplotype's copy are
               each at about half the normal read depth, while the scaffold's own allele is
               at normal depth (Flagger's "false duplication", here across haplotypes);
  repeat       the copies (or the other haplotype's copy) at about twice normal depth or
               more: collapsed repeats or a paralog family, not an assembly error;
  duplication  every copy at normal depth: a real (segmental or tandem) duplication this
               haplotype carries, left as it is.

Curator's rule: split only at a gap between the copies, and only with evidence that they
are haplotypic and not paralogous. In a `haplotypic` pair with a gap between the copies the
half-depth copy is cut out at the nearest gap on each side (a piece of its own, which the
small-scaffold step then judges, e.g. as a false duplicate); everything else is a review row. Depth is the map's
coverage track (long reads), relative to the median over the large scaffolds; without it
nothing is cut."""
import numpy as np


def _project(bk, t):
    """Query position of target position t along block bk (linear, strand-aware)."""
    f = (t - bk.ts) / max(bk.te - bk.ts, 1)
    return bk.qs + f * (bk.qe - bk.qs) if bk.same_strand else bk.qe - f * (bk.qe - bk.qs)


def self_duplications(raw, blocks, min_bp=50000, min_shared=50000, min_ident=0.5, repeat_depth=1.5,
                      half_max=0.7, normal_min=0.8):
    """Returns [dict] per pair of copies: s (scaffold index), c1, c2 (bp ranges of the copies),
    t, t0, t1 (other haplotype's scaffold and the shared stretch), inverted, d1, d2, dt (depth
    of each copy and of the other haplotype's copy, nan without a coverage track), gaps (gap
    pixels between the copies), kind, cut (gap pixels to cut after, or None), evidence."""
    idx = {n: i for i, n in enumerate(raw.names)}
    bpp = raw.bp_per_pix
    cov, gap = raw.graphs.get("coverage"), raw.graphs.get("gap")
    norm = None
    if cov is not None:
        largest = raw.scaffold_npix[raw.orig_scaffold] >= 0.2 * raw.scaffold_npix.max()
        norm = max(float(np.median(cov[largest])), 1e-9)

    def px(s, bp):
        return int(raw.scaffold_start[s]) + min(int(bp / bpp), int(raw.scaffold_npix[s]) - 1)

    def depth(s, a, b):
        if norm is None:
            return float("nan")
        p0, p1 = px(s, a), px(s, b)
        return float(np.median(cov[p0:max(p1, p0 + 1)])) / norm

    by = {}
    for bk in blocks:
        if bk.qe - bk.qs >= min_bp and bk.matches >= min_ident * (bk.qe - bk.qs) \
                and bk.q in idx and bk.t in idx and raw.scaffold_npix[idx[bk.q]] > 0:
            by.setdefault((bk.q, bk.t), []).append(bk)
    out = []
    for (q, t), xs in sorted(by.items()):
        xs.sort(key=lambda b: b.qs)
        s, ts_ = idx[q], idx[t]
        for i, x in enumerate(xs):
            for y in xs[i + 1:]:
                if y.qs < x.qe:                        # overlapping in the query: one copy
                    continue
                t0, t1 = max(x.ts, y.ts), min(x.te, y.te)
                if t1 - t0 < min_shared:
                    continue
                c1 = tuple(sorted((_project(x, t0), _project(x, t1))))
                c2 = tuple(sorted((_project(y, t0), _project(y, t1))))
                if c2[0] < c1[1]:
                    continue
                d1, d2, dt = depth(s, *c1), depth(s, *c2), depth(ts_, t0, t1)
                g0, g1 = px(s, c1[1]), px(s, c2[0])
                gaps = [int(p) for p in range(g0, g1 + 1) if gap is not None and gap[p] > 0]
                where = (f"{c1[0] / 1e6:.2f}-{c1[1] / 1e6:.2f} and {c2[0] / 1e6:.2f}-{c2[1] / 1e6:.2f} Mb "
                         f"({'inverted' if x.same_strand != y.same_strand else 'same orientation'}) both align to "
                         f"{t} {t0 / 1e6:.2f}-{t1 / 1e6:.2f} Mb ({(t1 - t0) / 1e3:.0f} kb)")
                dtxt = f"read depth {d1:.2f} / {d2:.2f}, the other haplotype's copy {dt:.2f} (1 = normal)"
                cut = None
                if norm is None:
                    kind, why = "unknown", "no coverage track: depth unknown, left as it is"
                elif max(d1, d2, dt) >= repeat_depth:
                    kind, why = "repeat", "depth about twice normal or more: collapsed repeats or paralogs, left as it is"
                elif min(d1, d2) <= half_max and dt <= half_max and max(d1, d2) >= normal_min:
                    kind = "haplotypic"
                    half = c1 if d1 <= d2 else c2
                    why = (f"the copy at {half[0] / 1e6:.2f}-{half[1] / 1e6:.2f} Mb and the other haplotype's copy are at "
                           "about half depth (they share their reads: the same allele twice), the other copy at normal "
                           "depth: both haplotypes in one scaffold")
                    if gaps:
                        # the half-depth copy is cut out at the nearest gap on each side (one of
                        # them between the copies); it is then a piece of its own (false duplicate)
                        a_, b_ = px(s, half[0]), px(s, half[1])
                        all_gaps = np.nonzero(gap[int(raw.scaffold_start[s]):
                                                  int(raw.scaffold_start[s] + raw.scaffold_npix[s])] > 0)[0] \
                            + int(raw.scaffold_start[s])
                        # (2 px of slack: bp-to-pixel rounding can put a gap on the copy's edge)
                        left = [int(p) for p in all_gaps if p <= a_ + 2]
                        right = [int(p) for p in all_gaps if p >= b_ - 2 and p > a_ + 2]
                        cut = [p for p in (max(left) if left else None, min(right) if right else None) if p is not None]
                        why += f"; cut out at {'the gaps' if len(cut) > 1 else 'the gap'} around that copy"
                    else:
                        why += "; no gap between the copies: not cut (curator's rule)"
                else:
                    kind, why = "duplication", ("every copy at about normal depth: a duplication this haplotype carries "
                                                "(segmental or tandem), not two alleles; left as it is")
                out.append(dict(s=s, c1=c1, c2=c2, t=t, t0=t0, t1=t1, inverted=x.same_strand != y.same_strand,
                                d1=d1, d2=d2, dt=dt, gaps=gaps, kind=kind, cut=cut,
                                evidence=f"{where}; {len(gaps)} gap(s) between the copies; {dtxt}; {why}"))
    return out


def apply_cuts(sites, dups):
    """Cut sites for haplotypic copies: an existing site at that gap becomes a cut (the
    reason says why), else a new cut site is added. Returns the sites, sorted."""
    by_q = {x["q"]: x for x in sites}
    why = "two copies of one region in this scaffold, haplotypic (see the self_duplication row)"
    for d in dups:
        if not d["cut"]:
            continue
        if any(by_q.get(q, {}).get("graph") for q in d["cut"]):      # the hifiasm graph backs a join there
            d["cut"] = None
            d["evidence"] += "; but the hifiasm graph backs the join at a gap around it: not cut"
            continue
        for q in d["cut"]:
            x = by_q.get(q)
            if x is None:
                by_q[q] = dict(q=q, s=d["s"], action="cut", reason="gap inside the scaffold; " + why, tip=(q, q))
            elif x["action"] != "cut":
                x["action"], x["reason"] = "cut", x["reason"] + "; overruled: " + why
    return sorted(by_q.values(), key=lambda x: x["q"])


def rows(raw, dups):
    """Decision rows `self_duplication`, one per pair of copies."""
    return [dict(action="self_duplication", scaffold=raw.names[d["s"]],
                 length_mb=round((d["c2"][1] - d["c1"][0]) / 1e6, 2), partner=d["t"], strand="-" if d["inverted"] else "+",
                 tag="", score=d["kind"], next_best="", evidence=d["evidence"]) for d in dups]
