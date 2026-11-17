"""Break sites at gaps and sharp breaks (`telomere_sites`), each decided cut / invert_tip /
keep, and breaks of colinearity with the other haplotype (`alignment_breaks`)."""
import numpy as np
from ..pretextio import SegmentedMap
from .common import _mean0, _signal
from .pairing import pair_homologs


def telomere_sites(pt, pairs, hi=None, hap_regex=r"^(H\d+)_", blocks=None, contact_max=0.5, near_px=3, end_px=15,
                   tol_frac=0.01, broad_min=0.6, invert_min=1.3, gap_max=0.6, sharp_max=0.1, merge_px=5, tip_max=0.5, use_breaks=True):
    """Candidate break sites, plus the homologous gap in the paired scaffold (dual scaffolding):
      - a gap inside a scaffold with contact across below `gap_max` of expectation (a telomere
        signal there is supporting evidence, not a requirement: telomere tracks are not
        reliable in every species);
      - a point with no gap where contact drops to background on both scales (`sharp_max`),
        with normal contact on each side: cut outright;
      - with `blocks` and `use_breaks`: a gap where the alignment to the other haplotype
        stops being colinear (`alignment_breaks`), whatever the contact.
    Each gap site is decided as:

      cut         the two sides have different homologs in the other haplotype, telomeres
                  (strand tracks) face the gap from both sides or from the long side, or
                  the sides share no chromosome-level contact: two chromosomes joined
      invert_tip  the sides belong together and the end piece fits better the other way
                  round, by contact or by its telomere facing the gap (the telomere then
                  ends up at the scaffold end)
      keep        the sides belong together as assembled; flagged for the curator
    A site that is the homologous gap of a cut in the paired scaffold is cut as well.

    Units: positions in pixels (`near_px`, `end_px`, `merge_px`); `contact_max`, `gap_max`,
    `sharp_max`, `tip_max`, `broad_min` are contact as a fraction of what adjacent sequence
    shows; `invert_min` is the ratio inverted / as assembled needed to invert; `tol_frac`
    the tolerance, as a fraction of scaffold length, for the matching gap in the partner.

    Returns a list of dicts: q (site is after this pixel), s (scaffold), action, reason, tip
    (pixel range of the end piece).
    """
    tel, gap = pt.graphs.get("telomere"), pt.graphs.get("gap")
    if gap is None:
        return []
    if tel is None:
        tel = np.zeros(pt.npix)
    thr = max(50.0, 5 * float(np.percentile(tel, 99)))     # well above the track's own noise
    sg0 = _signal(pt)
    sgh = _signal(hi) if hi is not None else sg0
    # contact across each point within 10 px, relative to expectation: on the high-MQ map if
    # there is one (MQ0 reads from the other haplotype bridge misjoins), else on MQ0
    near = sgh.insulation(10, 0)
    which = "high-MQ map" if hi is not None else "MQ0 map"

    near0 = sg0.insulation(10, 0) if hi is not None else near

    def contact(q, near=near):
        """Lowest contact across within 2 px of q (fraction of expectation; nan if unknown)."""
        w = near[max(0, q - 2):q + 3]
        return float(np.nanmin(w)) if np.isfinite(w).any() else float("nan")

    # sites: {pixel: (scaffold, reason)}; carried: {partner site: the site it mirrors};
    # site_c: contact across each site
    sites, carried, site_c = {}, {}, {}
    for q in np.nonzero(gap > 0)[0]:
        q = int(q)
        s = int(pt.orig_scaffold[q])
        a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
        if q - a < end_px or a + n - 1 - q < end_px:
            continue                           # too near a scaffold end for the contact window
        t = float(tel[max(a, q - near_px):q + near_px + 1].max())
        c = contact(q)
        if t >= thr and c < contact_max:
            pass                               # telomere at the gap: a looser contact limit
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
    from ..misjoin import detect_misjoins, Params as MisjoinParams
    dense = sg0.decay[min(20, len(sg0.decay) - 1)] >= 0.25 * sg0.decay[0]
    if sharp_max > 0 and dense:
        for cdt in detect_misjoins(pt, MisjoinParams(cut_max=sharp_max), sg0):
            q = int(cdt.after_px)
            s = int(pt.orig_scaffold[q])
            a, n = int(pt.scaffold_start[s]), int(pt.scaffold_npix[s])
            # gap sites are handled above; a sharp break near one is the same event
            if (not cdt.apply or cdt.at_gap or q - a < end_px or a + n - 1 - q < end_px
                    or any(abs(k - q) <= 2 * merge_px for k in sites) or not cdt.far < sharp_max):
                continue
            sharp[q] = (s, f"no gap here, but contact across drops to background on the MQ0 map ({cdt.near:.2f} within "
                           f"10 px, {cdt.far:.2f} at 10-40 px) with normal contact on each side")
    # gaps where the alignment to the other haplotype breaks: candidates even with normal contact
    if blocks is not None and use_breaks:
        for bq, kind, lt, rt in alignment_breaks(pt, blocks):
            # only breaks at a gap (within 2 px) are sites: inside a contig a break usually marks
            # where the other haplotype is in several scaffolds (contigs over scaffold joins)
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
            if close:                          # already a site on contact: the break adds to its reason
                sites[close[0]] = (sites[close[0]][0], sites[close[0]][1] + "; " + note)
                continue
            site_c[q] = contact(q)
            sites[q] = (s, f"gap inside the scaffold; {note}; contact across {site_c[q]:.2f} of expectation on the {which}")
    # Dual scaffolding: for each site, the matching gap in the homolog (other haplotype). A
    # misjoin made by the scaffolder in one haplotype often has its twin in the other.
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
        # partner gaps away from its ends, with somewhat reduced contact, and either a site
        # already or not next to one
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
            # the two flanks must map on either side of q and well apart (ten times the
            # positional tolerance), so that the partner gap really separates the same sequence
            for c in cand:
                lft, rgt = to_s(c, -1), to_s(c, +1)
                if lft is not None and rgt is not None and min(lft, rgt) < q < max(lft, rgt) \
                        and abs(lft - rgt) > 10 * max(30, int(tol_frac * n)):
                    hit = c
                    break
        else:
            # without alignments: the gap at the same relative position along the homolog
            # (mirrored if the pair is on opposite strands), within tol pixels
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
    bg = float(M[::37, ::41].mean())           # background: a sparse sample of the whole map
    nb = M.shape[0]
    scb = pt.orig_scaffold[np.minimum(np.arange(nb) * f + f // 2, pt.npix - 1)]
    D = 40
    # expected contact at each distance (in bins) within one scaffold: the binned decay curve
    g = np.array([_mean0(np.diagonal(M, d)[scb[:nb - d] == scb[d:]]) for d in range(D)])
    gh = sgh.decay_long(400)                   # decay in pixels on the map used for tip inversion
    # 3. telomere strand tracks (telo_scan.py 5p / 3p), when the map has them
    t5, t3 = pt.graphs.get("telomere_5p"), pt.graphs.get("telomere_3p")
    from ..tracks import telomere_threshold
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
        tip_right = right <= left                  # the end piece is the shorter side
        tip = (q + 1, a + n - 1) if tip_right else (a, q)
        L, Rr = int(vt.orig_scaffold[q]), int(vt.orig_scaffold[q + 1])
        ml, mr = mate.get(L), mate.get(Rr)
        # broad: mean contact between up to 15 bins (120 px) on each side, skipping the 2-3 bins
        # next to the gap, relative to one chromosome at those distances (background removed)
        b = q // f
        T = min((q - a + 1) // f, (a + n - 1 - q) // f, 15)
        broad = float("nan")
        if T - 2 >= 2:
            rows, cols = np.arange(b - T + 1, b - 1), np.arange(b + 3, b + T + 1)
            B = M[np.ix_(rows, cols)]
            ex = g[np.minimum(cols[None, :] - rows[:, None], D - 1)]
            broad = float((B.mean() - bg) / max(ex.mean() - bg, 1e-6))
        strand = strand_at(q, tip_right)
        # decision, first rule that applies: different homologs -> cut; telomeres facing the gap
        # from both sides or the long side -> cut; too short to measure but reduced -> cut;
        # little chromosome-level contact -> cut; else same chromosome -> invert_tip or keep
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
            # Compare the body's last W px against the tip's near end (as assembled) and its far
            # end (what would touch the body if the tip were inverted), each over expectation.
            W = max(3, min((tip[1] - tip[0] + 1) // 2, 50))
            body = (q - W + 1, q) if tip_right else (q + 1, q + W)
            near_w = (q + 1, q + W) if tip_right else (q - W + 1, q)
            far_w = (tip[1] - W + 1, tip[1]) if tip_right else (tip[0], tip[0] + W - 1)
            i = np.arange(W)
            ex = gh[np.minimum(i[:, None] + i[None, :] + 1, len(gh) - 1)].mean()
            as_is = float(sgh.block(body[0], body[1], near_w[0], near_w[1]).mean() / ex)
            inv = float(sgh.block(body[0], body[1], far_w[0], far_w[1]).mean() / ex)
            why += f"; the sides share chromosome-level contact ({broad:.2f}): same chromosome"
            # an end piece over 40% of the scaffold is not a tip: inverting it would be a large
            # rearrangement on weak evidence
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
        by_q[c]["mate"], by_q[q]["mate"] = q, c      # the homologous site in the paired scaffold
        if by_q[q]["action"] == "cut" and by_q[c]["action"] != "cut":
            by_q[c]["action"] = "cut"
            by_q[c]["reason"] += f"; cut because the homologous site in {pt.names[by_q[q]['s']]} is a misjoin"
    for q, (s, why) in sharp.items():
        out.append(dict(q=q, s=s, action="cut", reason=why, tip=(q, q)))
    out.sort(key=lambda x: x["q"])
    return out


def alignment_breaks(raw, blocks, min_len=100000, min_px=30, end_bp=300000, repeat_max=0.8):
    """Places where the alignment to the other haplotype stops being colinear: between two
    consecutive long alignments of a scaffold the other side changes scaffold, strand, or
    jumps. Returns [(pixel, kind, text for the left side, text for the right side)].

    At a gap this means the scaffold was joined differently from the other haplotype (a
    misjoin or misorientation in one of the two). Inside a contig it usually marks where
    the other haplotype is in several scaffolds, which says where those belong.

    Only alignments of mapq >= 30 and `min_len` bp count; scaffolds under `min_px` pixels are
    skipped; `end_bp`: how near its scaffold end an alignment must stop for a change of
    scaffold to be read as the other haplotype simply continuing in the next scaffold;
    `repeat_max`: breaks where the repeat track is at least this share are ignored."""
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
                # same scaffold and strand: a jump on the other side, beyond a small overlap
                # (300 kb) or a gap of 1 Mb or three times the distance on this side (indels)
                d = (B[3] - A[4]) if A[5] else (A[3] - B[4])
                if -3e5 <= d <= max(1e6, 3 * abs(B[0] - A[1]) + 3e5):
                    continue
                kind = "the two sides align to distant places"
            # the break pixel: midway between the end of A and the start of B on this scaffold
            q = int(raw.scaffold_start[i] + min((A[1] + B[0]) / 2 / bpp, raw.scaffold_npix[i] - 1))
            if rep is not None and float(rep[max(0, q - 2):q + 3].mean()) >= repeat_max:
                continue                                    # inside a repeat: alignments end there anyway
            out.append((q, kind, txt(A), txt(B)))
    return out
