"""Layout fixes from the curator's review of qqTroRaco (2026-11-20).

Shrapnel that is specific to one chromosome (one homolog pair) goes with it as an Unloc, even
when nothing says where in it: then a curator who wants to localise it finds it next to the pair
it belongs to (curator, 2026-11-20). That covers pieces in no pair (a piece cut off a scaffold:
qqTroRaco H1_scaffold_62 10.05-13.51 Mb), further pieces homologous to a scaffold
(H1_scaffold_76, homologous to H2_scaffold_42), and small pairs of their own whose two copies both
belong to a larger pair (H2_scaffold_69 with H1_scaffold_72 + H1_scaffold_2 18.2-21.8 Mb, all
pieces of the H1_scaffold_2 / H2_scaffold_59 + H2_scaffold_8 chromosome). Specific means: MQ0
contact with one pair at least `ratio` times that with any other pair.

The second copy of a pair is shown parallel to the first: its direction is read from how the two
arranged chains align to each other, not from the input scaffolds (after a rearrangement the
input strand no longer says it: qqTroRaco H2_scaffold_6 pieces were laid out anti-parallel to
H1_scaffold_33).
"""
import numpy as np


def orient_pair(raw, blocks, px1, px2, min_bp=1_000_000):
    """True if the chain laid out as pixels `px2` runs against the one laid out as `px1` along
    their alignment (blocks), False if with it, None if they share less than `min_bp`."""
    if not blocks:
        return None
    n = raw.npix
    pos1, pos2 = np.full(n, -1, dtype=np.int64), np.full(n, -1, dtype=np.int64)
    pos1[px1], pos2[px2] = np.arange(len(px1)), np.arange(len(px2))
    idx = {nm: i for i, nm in enumerate(raw.names)}
    bpp = raw.bp_per_pix
    pix = lambda name, bp: int(raw.scaffold_start[idx[name]]) + min(int(bp / bpp), int(raw.scaffold_npix[idx[name]]) - 1)
    xs, ys, ws = [], [], []
    for b in blocks:
        if b.q not in idx or b.t not in idx:
            continue
        f = np.linspace(0.0, 1.0, 8)                    # 8 points along the block
        qq = b.qs + f * (b.qe - 1 - b.qs)
        tt = b.ts + f * (b.te - 1 - b.ts) if b.same_strand else b.te - 1 - f * (b.te - 1 - b.ts)
        for qbp, tbp in zip(qq, tt):
            p, q = pix(b.q, qbp), pix(b.t, tbp)
            for a, c in ((p, q), (q, p)):
                if pos1[a] >= 0 and pos2[c] >= 0:
                    xs.append(pos1[a])
                    ys.append(pos2[c])
                    ws.append((b.qe - b.qs) / 8)
    if sum(ws) < min_bp or len(xs) < 4:
        return None
    x, y, w = np.array(xs, float), np.array(ys, float), np.array(ws, float)
    cov = np.sum(w * (x - np.average(x, weights=w)) * (y - np.average(y, weights=w)))
    return bool(cov < 0)


def hic_dir(sg, px1, px2, min_rows=10, min_corr=0.8):
    """Majority direction of the Hi-C homology line over the runs of the first chain (`hic_runs`),
    weighted by rows: True if the second chain runs against the first, False if with it, None
    when unclear (no run with a clear line, or no two-thirds majority)."""
    runs = hic_runs(sg, px1, px2, min_rows, min_corr)
    votes = {True: sum(w for _, d, w in runs if d), False: sum(w for _, d, w in runs if not d)}
    tot = votes[True] + votes[False]
    if not tot:
        return None
    win = votes[True] > votes[False]
    return win if votes[win] >= 2 * tot / 3 else None


def hic_runs(sg, px1, px2, min_rows=10, min_corr=0.8):
    """The Hi-C homology line between two chains laid out as `px1` and `px2`, per run of the first
    (a piece laid out in one direction): for each row with a clear peak, the column of the
    second where it peaks. Returns [((lo px, hi px), runs backwards, rows)] for the runs whose
    peaks follow a line (|correlation| >= `min_corr` over >= `min_rows` rows). bAegFun
    H1_scaffold_40 / H1_scaffold_33 (two H1 scaffolds, not in the H2-vs-H1 alignment): forwards at
    0.96 over 36 rows, while the small-scaffold alignment said the other strand (curator: invert);
    H1_scaffold_31, joined after H1_scaffold_40, runs the other way."""
    px1, px2 = np.asarray(px1), np.asarray(px2)
    s1, s2 = np.unique(px1), np.unique(px2)
    runs = lambda s: np.split(s, np.nonzero(np.diff(s) != 1)[0] + 1)
    M = np.zeros((len(s1), len(s2)), dtype=np.float32)
    i0 = 0
    for r in runs(s1):
        j0 = 0
        for c in runs(s2):
            M[i0:i0 + len(r), j0:j0 + len(c)] = sg.block(int(r[0]), int(r[-1]), int(c[0]), int(c[-1]))
            j0 += len(c)
        i0 += len(r)
    M = M[np.searchsorted(s1, px1)][:, np.searchsorted(s2, px2)]       # rows and columns in layout order
    med = np.median(M, axis=1)
    # one vote per run of the first chain (a piece laid out in one direction): a chain can hold a
    # piece with no line to the other copy (bAegFun H1_scaffold_31 joined after H1_scaffold_40),
    # whose noisy rows would blur the correlation of the whole chain
    step = np.diff(px1)
    cuts = np.nonzero((step != 1) & (step != -1))[0] + 1
    out = []
    for rows in np.split(np.arange(len(px1)), cuts):
        pk = [(i, int(np.argmax(M[i]))) for i in rows if M[i].max() > 2 * med[i] + 1]
        if len(pk) < min_rows:
            continue
        x, y = np.array(pk, dtype=float).T
        if x.std() == 0 or y.std() == 0:
            continue
        r = float(np.corrcoef(x, y)[0, 1])
        if abs(r) >= min_corr:
            out.append(((int(px1[rows].min()), int(px1[rows].max())), r < 0, len(pk)))
    return out


def assign(run_, main, pixels, sg, sexfrac, skip, ratio=3.0, small_mb=10.0):
    """Shrapnel specific to one pair -> (attachments {copy: [unit]}, dissolved pair ids, rows).
    `main`: the pairs to be laid out; `pixels(u)`: a unit's laid-out pixels (with its chain);
    `skip`: units already placed elsewhere (in a chain, attached, false duplicates)."""
    pt, raw, hap, npx = run_.pt, run_.raw, run_.hap, run_.npx_
    joined, info = run_.joined, run_.info
    mb = lambda u: npx[u] * pt.bp_per_pix / 1e6
    chain = lambda c: [c] + [v for v, x in joined.items() if x == c]
    size = lambda c: sum(mb(v) for v in chain(c))
    small = [p for p in main if max(size(p.a), size(p.b)) < small_mb]
    in_pair = {u for p in main for u in (p.a, p.b)}
    extra_b = {p.b for p in run_.pairs if p.kind == "extra"}
    loose = [int(u) for u in range(len(pt.names)) if npx[u] > 0 and u not in in_pair and u not in joined
             and u not in skip and (u in info.get("unpaired", {}) or u in extra_b)
             and u not in info.get("repetitive", {}) and sexfrac(u) < 0.5]
    big = [p for p in main if p not in small]
    cols = {id(p): [pixels(c) for c in (p.a, p.b)] for p in big}

    def rows_of(px):
        """MQ0 map rows of pixels `px` (each contiguous run read once) against every column."""
        px = np.sort(np.asarray(px))
        runs = np.split(px, np.nonzero(np.diff(px) != 1)[0] + 1)
        return np.vstack([sg.strip(int(r[0]), int(r[-1])).astype(np.float32) for r in runs])

    def best(u, own=None):
        """(pair, its contact, next pair's contact) for unit u's chain against every larger pair."""
        S = rows_of(pixels(u))
        sc = sorted(((max(float(S[:, c].mean()) for c in cols[id(p)]), p) for p in big if p is not own),
                    key=lambda x: -x[0])
        if len(sc) < 2:
            return None
        return sc[0][1], sc[0][0], sc[1][0]

    def copy_in(p, u):
        """The copy of pair p for unit u: the one of its own haplotype prefix if there is one."""
        same = [c for c in (p.a, p.b) if hap[c] == hap[u]]
        return same[0] if len(same) == 1 else max((p.a, p.b), key=lambda c: size(c))

    out, dissolved, rows = {}, set(), []

    def row(u, c, why):
        rows.append(dict(action="unloc", scaffold=pt.names[u], length_mb=round(size(u), 2), partner=pt.names[c],
                         strand="", tag="Unloc", score="", next_best="", evidence=why))
    for u in loose:
        r = best(u)
        if r and r[1] >= ratio * max(r[2], 1e-6):
            c = copy_in(r[0], u)
            out.setdefault(c, []).append(u)
            row(u, c, f"in no pair, but its contact with this chromosome pair is {r[1] / max(r[2], 1e-6):.1f}x that "
                      f"with any other ({r[1]:.1f} against {r[2]:.1f}): kept next to it as an Unloc, where a curator can "
                      "localise it")
    for p in small:
        ra, rb = best(p.a), best(p.b)
        if ra and rb and ra[0] is rb[0] and min(ra[1] / max(ra[2], 1e-6), rb[1] / max(rb[2], 1e-6)) >= ratio:
            dissolved.add(id(p))
            for u, r in ((p.a, ra), (p.b, rb)):
                c = copy_in(r[0], u)
                out.setdefault(c, []).append(u)
                row(u, c, f"a small pair of its own (with {pt.names[p.b if u == p.a else p.a]}), but both copies "
                          f"belong to this chromosome pair: contact {r[1]:.1f} against {r[2]:.1f} for any other pair; "
                          "kept next to it as an Unloc (with the pieces in its chain)")
    return out, dissolved, rows
