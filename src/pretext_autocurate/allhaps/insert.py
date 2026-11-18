"""Cuts that let a loose scaffold go inside its own haplotype's copy of a chromosome.

A scaffold in no homolog pair that aligns (over most of its length) to a stretch of the other
haplotype's scaffold O, where O's partner S (same haplotype as the loose scaffold) has no
sequence, is a missing piece of S. When that hole is in the middle of S, the chain can only put
the piece next to one of S's ends (icLarUrsu, SUPER_8: H2_scaffold_53 and H2_scaffold_65 fill
SUPER_8 30.3-34.0 Mb, which H2_scaffold_8 skips at 29.65 Mb; both were put at its
start, 30 Mb away). If S has a gap where its alignment to O jumps over the stretch (within
slack of the flanks), and little sequence of its own between the flanks, S is cut
there, and the arrangement step (`driver.along`) then puts the piece between the two halves.
No gap near the jump: no cut (curator's rule: cut only at gaps)."""
import re

import numpy as np


def insert_cuts(raw, pairs, blocks, hap_regex, sites, min_bp=200000, min_frac=0.5, max_off=2_000_000, slack=250_000,
                min_mapq=30):
    """New cut sites (dicts as in `telomere_sites`) at gaps of S where loose scaffolds fill a hole,
    and their evidence. `pairs` are on `raw` (the first pairing). Returns (sites, rows)."""
    if blocks is None:
        return sites, []
    bpp = raw.bp_per_pix
    idx = {n: i for i, n in enumerate(raw.names)}
    hap = {}
    for i, n in enumerate(raw.names):
        m = re.match(hap_regex, n)
        hap[i] = m.group(1) if m else ""
    partner = {}
    for p in pairs:
        if p.kind in ("pair", "unequal", "extra"):
            partner.setdefault(p.a, p.b)
            partner.setdefault(p.b, p.a)
    in_pair = set(partner)
    # alignments by (scaffold, other scaffold): (own start, own end, other start, other end, same strand)
    by = {}
    for bk in blocks:
        if bk.mapq < min_mapq or bk.q not in idx or bk.t not in idx:
            continue
        q, t = idx[bk.q], idx[bk.t]
        by.setdefault((q, t), []).append((bk.qs, bk.qe, bk.ts, bk.te, bk.same_strand))
        by.setdefault((t, q), []).append((bk.ts, bk.te, bk.qs, bk.qe, bk.same_strand))
    length = {i: int(raw.scaffold_npix[i]) * bpp for i in range(len(raw.names))}
    gap = raw.graphs.get("gap")
    if gap is None:
        return sites, []
    by_q = {x["q"]: x for x in sites}
    rows = []
    for v in range(len(raw.names)):
        if v in in_pair or length[v] < min_bp:
            continue
        tot = {}
        for (a, o), bl in by.items():
            if a == v and hap[o] != hap[v]:
                tot[o] = tot.get(o, 0) + sum(e - s for s, e, _, _, _ in bl)
        if not tot:
            continue
        o = max(tot, key=tot.get)
        if tot[o] < min_frac * length[v] or o not in partner:
            continue
        S = partner[o]
        if hap[S] != hap[v] or S == v:
            continue
        t0 = min(b[2] for b in by[(v, o)])
        t1 = max(b[3] for b in by[(v, o)])
        sb = by.get((S, o), [])
        if not sb:
            continue
        # S must leave the stretch empty (else v is a second copy: haplotig or false duplication)
        cov = sum(max(0, min(oe, t1) - max(os_, t0)) for _, _, os_, oe, _ in sb)
        if cov > 0.5 * (t1 - t0):
            continue
        # S positions just outside the stretch, from the blocks ending before it and starting after it
        before = [b for b in sb if b[3] <= t0 + 0.1 * (t1 - t0)]
        after = [b for b in sb if b[2] >= t1 - 0.1 * (t1 - t0)]
        if not before or not after:
            continue                                  # the stretch is beyond an end of S: no hole
        lb, ra = max(before, key=lambda b: b[3]), min(after, key=lambda b: b[2])
        s_at = lambda b, x: b[0] + (x - b[2]) / max(b[3] - b[2], 1) * (b[1] - b[0]) if b[4] else \
            b[1] - (x - b[2]) / max(b[3] - b[2], 1) * (b[1] - b[0])
        sa, sz = s_at(lb, lb[3]), s_at(ra, ra[2])
        f0, f1 = sorted((sa, sz))
        # S's own sequence between the flanks must be short next to the hole: otherwise S has its
        # own version of the stretch and v is an alternative copy (bAegFun H1_scaffold_168 in
        # H1_scaffold_6: 1.76 Mb of S between the flanks for a 0.88 Mb hole)
        if f1 - f0 > max(slack, 0.5 * (t1 - t0)):
            continue
        j = (f0 + f1) / 2
        if j < max_off or j > length[S] - max_off:
            continue                                  # at an end of S: the chain handles it
        s0, n = int(raw.scaffold_start[S]), int(raw.scaffold_npix[S])
        gp = np.nonzero(gap[s0:s0 + n] > 0)[0] + s0
        # the gap must be where the alignment jumps (within slack of the flanks): a gap further
        # off would put S's sequence between it and the hole on the wrong side of v (icLarUrsu
        # H2_scaffold_8 jumps at 29.65 Mb, its gap is at 30.67 Mb: not cut)
        gp = [int(g) for g in gp if f0 - slack <= (g - s0 + 0.5) * bpp <= f1 + slack]
        if not gp:
            continue
        q = min(gp, key=lambda g: abs((g - s0 + 0.5) * bpp - j))
        why = (f"gap inside the scaffold; {raw.names[v]} ({length[v] / 1e6:.2f} Mb, in no pair) fills a hole here: "
               f"{raw.names[S]} skips {raw.names[o]} {t0 / 1e6:.2f}-{t1 / 1e6:.2f} Mb at {j / 1e6:.2f} Mb, and "
               f"{raw.names[v]} aligns there (sequence alignment); cut at the gap at {(q - s0) * bpp / 1e6:.2f} Mb "
               "so it can go in")
        x = by_q.get(q)
        if x is not None and x.get("graph"):
            rows.append(dict(action="insert_candidate", scaffold=raw.names[v], length_mb=round(length[v] / 1e6, 2),
                             partner=raw.names[S], strand="", tag="", score="", next_best="",
                             evidence=why.replace("; cut at", "; would cut at") + "; but the hifiasm graph backs the "
                                                                                  "join at that gap: not cut"))
            continue
        if x is None:
            by_q[q] = dict(q=q, s=S, action="cut", reason=why, tip=(q, q))
        elif x["action"] != "cut":
            x["action"], x["reason"] = "cut", x["reason"] + "; overruled: " + why
        elif raw.names[v] not in x["reason"]:
            x["reason"] += "; also: " + why.split("; ", 1)[1]
    return sorted(by_q.values(), key=lambda x: x["q"]), rows
