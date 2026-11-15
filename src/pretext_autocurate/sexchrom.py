"""Sex chromosomes from coverage and alignment: hemizygous stretches, pseudoautosomal regions
(PARs), and sex chromosomes joined to another chromosome.

Coverage model (measured on both test genomes; the map's `coverage` track divided by its
median over the large scaffolds): sequence present in both haplotypes is at about 1; a
hemizygous stretch, present once in the individual (the X- or Y-specific part of an XY male,
the Z- or W-specific part of a ZW female, an X0 male's X), is at about 0.5; a collapsed
region at about 2. A PAR is in both sex chromosomes and recombines like an autosome: normal
coverage, and it aligns to the other sex chromosome's PAR.

Half coverage alone is not enough. Divergent heterozygous stretches of autosomes (1-4 Mb,
about one per chromosome in Troglobisium) are at half coverage in both copies too, but they
align to each other. So:

- hemizygous: half coverage and not aligned to the other haplotype (`aligned_fraction`);
- sex-linked stretch: a half-coverage stretch holding at least `core_min` bp of hemizygous
  sequence, plus a half-coverage stretch that aligns into one (the homologous, gametologous
  part of the other sex chromosome: in Troglobisium the Y, H2_scaffold_4 0-33 Mb, aligns to
  the X, H1_scaffold_1 25-58 Mb, both at half coverage);
- next to a sex-linked stretch, normal-coverage sequence is a PAR when it has
  chromosome-level Hi-C contact with the stretch and aligns to the normal-coverage end of
  another sex-linked scaffold (owl: Z H1_scaffold_3 0-6 Mb with W H2_scaffold_8 28-34.5 Mb);
  with no such contact it is another chromosome joined to the sex chromosome (Troglobisium:
  an autosome on H1_scaffold_1 0-25 Mb and H2_scaffold_4 33-58 Mb); with contact but
  aligned to an ordinary chromosome, an autosome fused to the sex chromosome or a misjoin.

When both sex chromosomes were assembled into one haplotype, their PARs are not in the
H2-against-H1 alignment; `same_hap_candidates` names the scaffolds to align within each
haplotype (homology.align_within) and those alignments are passed as `extra_blocks`.

Which copy is X or Y (Z or W) is not decided from the data: that needs knowledge of the
species, given as --sex-system (`sex_labels`): XY / ZW, X0 / Z0 (one sex chromosome), or
three-chromosome systems such as X1X2Y, where the Y carries parts of both X's.

Only what is hemizygous in this individual shows. A sample of the homogametic sex (XX, ZZ),
or homomorphic sex chromosomes (young, still recombining over most of their length), look
like autosomes: they are paired as such, and finding nothing is not evidence of no sex
chromosomes.
"""
import re

import numpy as np

# per-pixel states of the segmentation
HEMI, HALF_ALN, OTHER = 0, 3, 1


def _mip(raw, max_side=8192):
    """The whole-map matrix at the finest mip level of at most `max_side` bins, cached on the
    map under the same key as allhaps.pair_homologs, and its bin size in pixels."""
    lv = 0
    while (raw.npix >> lv) > max_side and lv < raw.n_mips - 1:
        lv += 1
    if ("_mip", lv) not in raw.__dict__:
        raw.__dict__[("_mip", lv)] = raw.matrix(lv).astype(np.float32)
    M = raw.__dict__[("_mip", lv)]
    return M, raw.npix // M.shape[0]


def coverage_norm(raw):
    """The coverage track divided by its median over the scaffolds of at least a fifth of the
    largest one's length (so about 1 for a normal copy), or None without a coverage track."""
    cov = raw.graphs.get("coverage")
    if cov is None:
        return None
    largest = raw.scaffold_npix[raw.orig_scaffold] >= 0.2 * raw.scaffold_npix.max()
    return cov.astype(np.float64) / max(float(np.median(cov[largest])), 1e-9)


def block_index(blocks, min_mapq=30):
    """{scaffold: [(start, end, other scaffold, other start, other end, same strand)]} with
    every alignment of mapq >= `min_mapq` entered from both sides (as allhaps._block_index)."""
    idx = {}
    for bk in blocks or []:
        if bk.mapq < min_mapq or bk.q == bk.t:
            continue
        idx.setdefault(bk.q, []).append((bk.qs, bk.qe, bk.t, bk.ts, bk.te, bk.same_strand))
        idx.setdefault(bk.t, []).append((bk.ts, bk.te, bk.q, bk.qs, bk.qe, bk.same_strand))
    return idx


def aligned_fraction(raw, idx):
    """Per pixel, the share of its bases inside an alignment to another scaffold."""
    bpp = raw.bp_per_pix
    out = np.zeros(raw.npix)
    ix = {n: i for i, n in enumerate(raw.names)}
    for name, lst in idx.items():
        i = ix.get(name)
        if i is None or raw.scaffold_npix[i] == 0:
            continue
        base, n = int(raw.scaffold_start[i]), int(raw.scaffold_npix[i])
        merged = []
        for s, e, *_ in sorted(lst):                     # union of the aligned intervals
            if merged and s <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([s, e])
        for s, e in merged:
            p0, p1 = int(s // bpp), min(int((e - 1) // bpp), n - 1)
            if p0 >= n:
                continue
            if p0 == p1:
                out[base + p0] += (e - s) / bpp
            else:
                out[base + p0] += ((p0 + 1) * bpp - s) / bpp
                out[base + p0 + 1:base + p1] += 1.0
                out[base + p1] += (e - p1 * bpp) / bpp
    return np.minimum(out, 1.0)


def _runs(state, min_len):
    """Runs [state, first, last] (inclusive) of an integer array, with runs shorter than
    `min_len` absorbed into a neighbour (the longer one) until none is left."""
    runs = []
    for i, v in enumerate(state):
        if runs and runs[-1][0] == v:
            runs[-1][2] = i
        else:
            runs.append([int(v), i, i])
    while len(runs) > 1:
        k = min(range(len(runs)), key=lambda k: runs[k][2] - runs[k][1])
        if runs[k][2] - runs[k][1] + 1 >= min_len:
            break
        if k == 0:
            j = 1
        elif k == len(runs) - 1:
            j = k - 1
        else:
            j = k - 1 if runs[k - 1][2] - runs[k - 1][1] >= runs[k + 1][2] - runs[k + 1][1] else k + 1
        runs[j][1], runs[j][2] = min(runs[j][1], runs[k][1]), max(runs[j][2], runs[k][2])
        del runs[k]
        merged = [runs[0]]                               # neighbours of one state become one run
        for r in runs[1:]:
            if r[0] == merged[-1][0]:
                merged[-1][2] = r[2]
            else:
                merged.append(r)
        runs = merged
    return runs


def _contact(raw, a, b, c, d):
    """Chromosome-level contact between pixel ranges a..b and c..d (inclusive) on the MQ0
    map: their mean contact over the geometric mean of each one's contact with itself. One
    chromosome gives 0.1-0.4 over 10 Mb, two chromosomes under 0.03. NaN if either range
    holds fewer than two whole bins."""
    M, f = _mip(raw)
    x, y = slice(-(-a // f), (b + 1) // f), slice(-(-c // f), (d + 1) // f)
    if x.stop - x.start < 2 or y.stop - y.start < 2:
        return float("nan")
    xy, xx, yy = float(M[x, y].mean()), float(M[x, x].mean()), float(M[y, y].mean())
    return xy / max(np.sqrt(xx * yy), 1e-12)


def analyse(raw, blocks=None, extra_blocks=None, level=0.5, half_max=0.72, aligned_min=0.5):
    """Sex-linked stretches and what lies next to them. `raw` is the MQ0 map (uncut), `blocks`
    the H2-against-H1 alignment, `extra_blocks` alignments within a haplotype (optional).

    Returns dict(
      stretches=[dict(s, a, b, hemi_bp, aln_bp, cov, via)]: sex-linked stretches (pixels
          a..b inclusive of scaffold s; `via` the stretch it aligns to, for the gametologous
          part of the other sex chromosome, else None),
      ends=[dict(s, side, boundary, a, b, kind, contact, partner, par, cov_sex, cov_other)]:
          the sequence next to each stretch (pixels a..b, the stretch edge at `boundary`);
          kind 'par' / 'par_candidate' / 'joined' / 'fused' (see the module docstring),
      mask: per-pixel 2 = hemizygous, 1 = other sex-linked, 0 = not sex-linked,
      par_mask: per-pixel 2 = PAR, 1 = PAR candidate,
      norm: the normalised coverage, or None (then everything else is empty)).
    `level` (0 conservative .. 1 aggressive) sets the hemizygous bp a stretch needs."""
    out = dict(stretches=[], ends=[], mask=np.zeros(raw.npix, np.int8), par_mask=np.zeros(raw.npix, np.int8),
               norm=None)
    norm = coverage_norm(raw)
    if norm is None:
        return out
    out["norm"] = norm
    lvl = min(max(level, 0.0), 1.0)
    at = lambda c, d, a: c + (d - c) * 2 * lvl if lvl <= 0.5 else d + (a - d) * 2 * (lvl - 0.5)
    # core_min: hemizygous bp that make a stretch sex-linked. Under 4 Mb, low-coverage
    # chromosome ends (subtelomeric repeats, 0.3-0.4) and heterozygous insertions qualify.
    core_min = at(8e6, 5e6, 4e6)
    same_min = at(0.08, 0.05, 0.03)                    # chromosome-level contact, see _contact
    bpp = raw.bp_per_pix
    idx = block_index(list(blocks or []) + list(extra_blocks or []))
    have_aln = bool(idx)
    if not have_aln:
        # without alignments, divergent heterozygous stretches cannot be told from hemizygous ones
        core_min *= 2
    aln = aligned_fraction(raw, idx) if have_aln else np.zeros(raw.npix)
    w = max(2, int(round(250e3 / bpp)))                  # smoothing half-window, px (~0.25 Mb)
    min_len = max(3, int(round(500e3 / bpp)))            # shortest run kept, px (~0.5 Mb)
    blip = 1.5e6 / bpp                                   # normal coverage inside a stretch, px
    ix = {n: i for i, n in enumerate(raw.names)}
    half = {}                                            # scaffold -> [[a, b, hemi px, aligned px]] (global px)
    runs_of = {}
    for s in range(raw.n_scaffolds):
        n = int(raw.scaffold_npix[s])
        if n * bpp < 1e6:
            continue
        a0 = int(raw.scaffold_start[s])
        c = norm[a0:a0 + n]
        al = aln[a0:a0 + n]
        cs = np.array([np.median(c[max(0, i - w):i + w + 1]) for i in range(n)])
        als = np.convolve(al, np.ones(2 * w + 1), "same") / np.convolve(np.ones(n), np.ones(2 * w + 1), "same")
        st = np.where(cs < half_max, np.where(als >= aligned_min, HALF_ALN, HEMI), OTHER)
        runs = _runs(st, min_len)
        runs_of[s] = [(v, a0 + x, a0 + y) for v, x, y in runs]
        # half-coverage stretches: half runs, joined across short normal runs between them
        cur = None
        for k, (v, x, y) in enumerate(runs):
            if v != OTHER:
                if cur is None:
                    cur = [a0 + x, a0 + y, 0, 0]
                cur[1] = a0 + y
                cur[2 if v == HEMI else 3] += y - x + 1
            elif cur is not None and y - x + 1 <= blip and k + 1 < len(runs):
                cur[1] = a0 + y                              # a blip: the next run is half again
            elif cur is not None:
                half.setdefault(s, []).append(cur)
                cur = None
        if cur is not None:
            half.setdefault(s, []).append(cur)
    # trim a stretch that ended on a blip back to its last half run
    for s, lst in half.items():
        for h in lst:
            last = [y for v, x, y in runs_of[s] if v != OTHER and h[0] <= x <= h[1]]
            h[1] = max(last) if last else h[1]
    def single_copy(s, h):
        """A stretch at about half coverage (not the 0.3-0.4 of repeat-rich ends) that is a
        large share of its scaffold, or large in itself (a sex chromosome joined to another):
        a sex chromosome rather than an insertion inside an autosome."""
        size = (h[1] - h[0] + 1) * bpp
        return (0.38 <= float(np.median(norm[h[0]:h[1] + 1])) <= 0.68
                and (size >= 0.3 * raw.scaffold_npix[s] * bpp or size >= 3 * core_min))
    sex = {(s, h[0], h[1]): dict(s=s, a=h[0], b=h[1], hemi_bp=h[2] * bpp, aln_bp=h[3] * bpp, via=None)
           for s, lst in half.items() for h in lst if h[2] * bpp >= core_min and single_copy(s, h)}

    def lands(name, lo, hi):
        """{other scaffold: [(bp, other start, other end)]}: where bp lo..hi of `name` aligns."""
        got = {}
        for s0, e0, o, os_, oe_, same in idx.get(name, []):
            ov = min(e0, hi) - max(s0, lo)
            if ov <= 0:
                continue
            f0, f1 = (max(s0, lo) - s0) / (e0 - s0), (min(e0, hi) - s0) / (e0 - s0)
            if not same:
                f0, f1 = 1 - f1, 1 - f0
            got.setdefault(o, []).append((ov, os_ + f0 * (oe_ - os_), os_ + f1 * (oe_ - os_)))
        return got

    def in_sex(o, x0, x1):
        """The sex-linked stretch of scaffold `o` holding bp x0..x1 (its midpoint), or None."""
        i = ix.get(o)
        if i is None:
            return None
        px = int(raw.scaffold_start[i]) + int((x0 + x1) / 2 // bpp)
        return next((k for k in sex if k[0] == i and k[1] <= px <= k[2]), None)

    loc = lambda px: raw.rel_coord[px] * bpp            # bp of a pixel within its scaffold
    # gametologous stretches: half coverage, aligned mostly into a sex-linked stretch elsewhere
    for _ in range(3):
        added = False
        for s, lst in half.items():
            for h in lst:
                if (s, h[0], h[1]) in sex or (h[1] - h[0] + 1) * bpp < 1e6 \
                        or not 0.38 <= float(np.median(norm[h[0]:h[1] + 1])) <= 0.68:
                    continue
                into, tot = {}, 0.0
                for o, hits in lands(raw.names[s], loc(h[0]), loc(h[1]) + bpp).items():
                    for ov, x0, x1 in hits:
                        tot += ov
                        k = in_sex(o, x0, x1)
                        if k is not None and k[0] != s:
                            into[k] = into.get(k, 0.0) + ov
                if into and max(into.values()) >= 0.5 * max(tot, (h[1] - h[0] + 1) * bpp * 0.5):
                    k = max(into, key=into.get)
                    sex[(s, h[0], h[1])] = dict(s=s, a=h[0], b=h[1], hemi_bp=h[2] * bpp, aln_bp=h[3] * bpp, via=k)
                    added = True
        if not added:
            break
    for k, d in sex.items():
        d["cov"] = float(np.median(norm[d["a"]:d["b"] + 1]))
        for v, x, y in runs_of[d["s"]]:
            if d["a"] <= x and y <= d["b"]:
                out["mask"][x:y + 1] = 2 if v == HEMI else 1
        out["mask"][d["a"]:d["b"] + 1] = np.maximum(out["mask"][d["a"]:d["b"] + 1], 1)
    out["stretches"] = sorted(sex.values(), key=lambda d: d["a"])
    # what lies next to each stretch, up to the next stretch or the scaffold end
    win = int(10e6 / bpp)                                # contact windows, px from the boundary
    for d in out["stretches"]:
        s = d["s"]
        a0, n = int(raw.scaffold_start[s]), int(raw.scaffold_npix[s])
        mine = sorted((x["a"], x["b"]) for x in out["stretches"] if x["s"] == s)
        for side in ("left", "right"):
            if side == "left":
                lo = max([b + 1 for a, b in mine if b < d["a"]] + [a0])
                hi = d["a"] - 1
                bnd = d["a"]
                na, nb, sa, sb = max(lo, hi - win + 1), hi, d["a"], min(d["b"], d["a"] + win - 1)
            else:
                lo = d["b"] + 1
                hi = min([a - 1 for a, b in mine if a > d["b"]] + [a0 + n - 1])
                bnd = d["b"]
                na, nb, sa, sb = lo, min(hi, lo + win - 1), max(d["a"], d["b"] - win + 1), d["b"]
            if (hi - lo + 1) < min_len:
                continue
            if any(x["a"] <= lo and hi <= x["b"] for x in out["stretches"]):
                continue
            ct = _contact(raw, na, nb, sa, sb)
            # where the neighbouring sequence aligns: to the non-sex-linked end of another
            # sex-linked scaffold (a PAR), or to an ordinary chromosome
            to_par, to_other, span = {}, {}, {}
            for o, hits in lands(raw.names[s], loc(lo), loc(hi) + bpp).items():
                i = ix.get(o)
                for ov, x0, x1 in hits:
                    if i is not None and i != s and any(x["s"] == i for x in out["stretches"]) and in_sex(o, x0, x1) is None:
                        to_par[o] = to_par.get(o, 0.0) + ov
                    else:
                        to_other[o] = to_other.get(o, 0.0) + ov
            L = (hi - lo + 1) * bpp
            partner, par = None, None
            if to_par and max(to_par.values()) >= max(100e3, 0.5 * sum(to_par.values())):
                partner = max(to_par, key=to_par.get)
                # the PAR: the part of the neighbouring sequence that aligns to the partner
                xs = [(max(s0, loc(lo)), min(e0, loc(hi) + bpp)) for s0, e0, o, *_ in idx.get(raw.names[s], [])
                      if o == partner and min(e0, loc(hi) + bpp) > max(s0, loc(lo))]
                par = (min(x for x, _ in xs), max(y for _, y in xs))
            other = max(to_other, key=to_other.get) if to_other else None
            if not np.isnan(ct) and ct < same_min:
                kind = "joined"
            elif partner is not None:
                kind = "par"
            elif other is not None and to_other[other] >= 0.5 * L:
                kind = "fused"
            else:
                kind = "par_candidate"
            edge_bp = loc(bnd) if side == "left" else loc(bnd) + bpp     # where the stretch ends, bp
            out["ends"].append(dict(s=s, side=side, boundary=bnd, edge_bp=edge_bp, a=lo, b=hi, kind=kind, contact=ct,
                                    partner=partner,
                                    par=par, by="alignment" if partner is not None else "", win=(na, nb),
                                    other=other, other_bp=to_other.get(other, 0.0) if other else 0.0, have_aln=have_aln,
                                    cov_sex=float(np.median(norm[sa:sb + 1])), cov_other=float(np.median(norm[na:nb + 1]))))
    # Not confirmed by alignment (none given, or both copies in one haplotype and not aligned
    # to each other): two candidate ends of different scaffolds whose MQ0 contact is as strong
    # as within each (the reads of a PAR map to both copies, like a homolog's).
    for e in out["ends"]:
        if e["kind"] != "par_candidate":
            continue
        best = max(((_contact(raw, *e["win"], *f["win"]), f) for f in out["ends"]
                    if f is not e and f["s"] != e["s"] and f["kind"] in ("par_candidate", "par")),
                   key=lambda x: (-1 if np.isnan(x[0]) else x[0]), default=(float("nan"), None))
        if best[1] is not None and best[0] >= 0.5:
            e.update(kind="par", partner=raw.names[best[1]["s"]], by="MQ0 Hi-C", homolog_contact=best[0],
                     par=(loc(e["a"]), loc(e["b"]) + bpp))
    for e in out["ends"]:
        a0, n = int(raw.scaffold_start[e["s"]]), int(raw.scaffold_npix[e["s"]])
        if e["kind"] == "par":
            p0 = a0 + int(e["par"][0] // bpp)
            p1 = a0 + min(int((e["par"][1] - 1) // bpp), n - 1)
            out["par_mask"][p0:p1 + 1] = 2
        elif e["kind"] == "par_candidate":
            out["par_mask"][e["a"]:e["b"] + 1] = np.maximum(out["par_mask"][e["a"]:e["b"] + 1], 1)
    return out


# --sex-system values. Two-chromosome systems: (label of the larger, of the smaller). One
# chromosome (X0 / Z0: the heterogametic sex has a single sex chromosome). Three-chromosome
# systems: (label of the copy homologous to both others, numbered label of the two others):
# in X1X2Y the Y carries parts of both X's.
SEX_SYSTEMS = {"XY": ("X", "Y"), "ZW": ("Z", "W"), "X0": ("X",), "Z0": ("Z",),
               "X1X2Y": ("Y", "X"), "XY1Y2": ("X", "Y"), "Z1Z2W": ("W", "Z"), "ZW1W2": ("Z", "W")}


def sex_labels(system, copies, size, links):
    """Tags for the sex chromosome copies found, by the species' sex system (`SEX_SYSTEMS`).

    `copies`: the copies (any hashable ids); `size`: {copy: bp, with the pieces in its chain};
    `links`: {frozenset((a, b)): weight > 0} for copies paired or aligned to each other.
    Returns ({copy: tag}, note); copies that cannot be told apart are left without a tag and
    the note says why. Nothing found is not evidence against sex chromosomes: the sample may
    be of the homogametic sex, or the sex chromosomes homomorphic."""
    if not system or not copies:
        return {}, ""
    lab = SEX_SYSTEMS[system]
    big = sorted(dict.fromkeys(copies), key=lambda c: -size.get(c, 0))
    linked = lambda a, b: links.get(frozenset((a, b)), 0) > 0
    deg = {c: sum(linked(c, d) for d in big if d != c) for c in big}
    if len(lab) == 1:
        note = (f"--sex-system {system}: every sex chromosome copy tagged {lab[0]}")
        if any(deg.values()):
            note += ("; but some are homologous to each other, as an X with a Y (Z with W) would be: check the sex "
                     "system, or whether one of them is a second copy")
        return {c: lab[0] for c in big}, note
    if system in ("XY", "ZW"):
        pair = next(((a, b) for i, a in enumerate(big) for b in big[i + 1:] if linked(a, b)), None)
        if pair is None:
            pair = tuple(big[:2])
            why = "no two copies are homologous; by size"
        else:
            why = "the homologous pair"
        out = {pair[0]: lab[0]}
        if len(pair) > 1:
            out[pair[1]] = lab[1]
        note = (f"--sex-system {system}: {lab[0]} the larger and {lab[1]} the smaller of {why} (check: in some species "
                f"{lab[1]} is the larger)" if len(pair) > 1 else
                f"--sex-system {system}: one sex chromosome copy found, tagged {lab[0]} ({lab[1]} not found: "
                "unassembled, or homomorphic with its partner)")
        rest = [c for c in big if c not in out]
        if rest:
            note += f"; {len(rest)} further cop{'y' if len(rest) == 1 else 'ies'} left untagged"
        return out, note
    # three chromosomes: the one homologous to both others is the single one
    one, two = lab
    hub = max(big, key=lambda c: (deg[c], size.get(c, 0)))
    if len(big) < 3 or deg[hub] < 2:
        return {}, (f"--sex-system {system}: {len(big)} sex chromosome cop{'y' if len(big) == 1 else 'ies'} found, "
                    f"{'none homologous to two others' if len(big) >= 3 else 'the system has three'}: not tagged")
    others = sorted([c for c in big if c != hub and linked(hub, c)], key=lambda c: -size.get(c, 0))[:2]
    out = {hub: one, others[0]: two + "1", others[1]: two + "2"}
    note = (f"--sex-system {system}: {one} is the copy homologous to both others; {two}1 the larger of these, {two}2 "
            "the smaller")
    rest = [c for c in big if c not in out]
    if rest:
        note += f"; {len(rest)} further cop{'y' if len(rest) == 1 else 'ies'} left untagged"
    return out, note


def same_hap_candidates(raw, blocks=None, hap_regex=r"^(H\d+)_", level=0.5):
    """Scaffolds with a sex-linked stretch, grouped by haplotype prefix, for the prefixes that
    hold two or more: their PARs may only be found by aligning them to each other."""
    res = analyse(raw, blocks, level=level)
    by = {}
    for d in res["stretches"]:
        m = re.match(hap_regex, raw.names[d["s"]])
        by.setdefault(m.group(1) if m else "", set()).add(raw.names[d["s"]])
    return {h: sorted(v) for h, v in by.items() if len(v) >= 2}


def _bedgraph(raw, mask, path, values):
    """Write the runs of `mask` with a value in `values` ({mask value: track value}) as a
    bedgraph in scaffold coordinates. Returns the number of lines."""
    bpp = raw.bp_per_pix
    lens = raw.scaffold_len_bp()
    n = 0
    with open(path, "w") as f:
        for s in range(raw.n_scaffolds):
            a0, k = int(raw.scaffold_start[s]), int(raw.scaffold_npix[s])
            seg = mask[a0:a0 + k]
            i = 0
            while i < k:
                j = i
                while j + 1 < k and seg[j + 1] == seg[i]:
                    j += 1
                if int(seg[i]) in values:
                    L = int(lens[s])
                    f.write(f"{raw.names[s]}\t{int(i * bpp)}\t{min(L, int((j + 1) * bpp))}\t{values[int(seg[i])]}\n")
                    n += 1
                i = j + 1
    return n


def write_tracks(raw, res, out_dir):
    """Bedgraphs for PretextGraph: `sex_chromosome` (2 hemizygous, 1 the rest of a sex-linked
    stretch) and `par` (2 PAR, 1 PAR candidate). Returns {track name: path} for those with
    any interval."""
    import os
    os.makedirs(out_dir, exist_ok=True)
    got = {}
    for name, mask in (("sex_chromosome", res["mask"]), ("par", res["par_mask"])):
        path = os.path.join(out_dir, name + ".bedgraph")
        if _bedgraph(raw, mask, path, {2: 2, 1: 1}):
            got[name] = path
    return got
