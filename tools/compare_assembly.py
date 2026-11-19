"""Score a layout against a curated assembly instead of a savestate (when the curation was
done in several rounds and only the final FASTA is kept).

    python3 tools/align_chunked.py CURATED.fa(.gz) DIPLOID.fa diploid_vs_curated.paf -t 6 \\
        --chunk-mb 0.1 --mm-args "-K 20M --cap-kalloc 1g -f 0.002 -w 50"
    python3 tools/compare_assembly.py MAP PAF [SAVESTATE ...] [--plot out.png] [--min-mb 1]

Every pixel of the map is projected through the alignment onto the curated assembly
(chromosome, position, strand); both haplotypes align to a primary assembly, so each piece
of either one gets a place. For each savestate (none: the map as assembled):

  grouping   per painted chromosome, the share of its bp on its main curated chromosome
             (purity); per curated chromosome and haplotype, how many layout chromosomes hold
             it (split);
  junctions  each junction between two input scaffolds in a painted chromosome: right when
             both sides land on one curated chromosome within --near-mb, else wrong (other
             chromosome, or far apart: order);
  missed     input scaffolds of one haplotype (both >= --min-mb) that sit next to each other
             along a curated chromosome but are not next to each other in the layout;

and a D-GENIES-style dot plot (layout position against curated position) per savestate."""
import argparse
import bisect
import collections
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate import MapState, Pretext, SaveState     # noqa: E402


def project(pt, paf, min_len=20000):
    """Per pixel: (curated chromosome index or -1, position bp, strand +1/-1); and the names."""
    blocks, tnames, tlen = collections.defaultdict(list), {}, {}
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if int(c[10]) < min_len:
                continue
            t = tnames.setdefault(c[5], len(tnames))
            tlen[t] = int(c[6])
            blocks[c[0]].append((int(c[2]), int(c[3]), t, int(c[7]), int(c[8]), 1 if c[4] == "+" else -1, int(c[10])))
    for v in blocks.values():
        v.sort()
    n = pt.npix
    chrom = np.full(n, -1, dtype=np.int64)
    pos = np.zeros(n, dtype=np.float64)
    strand = np.zeros(n, dtype=np.int8)
    bpp = pt.bp_per_pix
    for s, name in enumerate(pt.names):
        bl = blocks.get(name)
        if not bl or pt.scaffold_npix[s] == 0:
            continue
        starts = [b[0] for b in bl]
        a = int(pt.scaffold_start[s])
        for k in range(int(pt.scaffold_npix[s])):
            q = (k + 0.5) * bpp
            i = bisect.bisect_right(starts, q) - 1
            best = None
            for j in range(max(0, i - 20), i + 1):        # blocks starting at or before q that cover it
                b = bl[j]
                if b[0] <= q < b[1] and (best is None or b[6] > best[6]):
                    best = b
            if best is None:
                continue
            f_ = (q - best[0]) / max(best[1] - best[0], 1)
            chrom[a + k] = best[2]
            pos[a + k] = best[3] + f_ * (best[4] - best[3]) if best[5] > 0 else best[4] - f_ * (best[4] - best[3])
            strand[a + k] = best[5]
    names = [None] * len(tnames)
    for k, v in tnames.items():
        names[v] = k
    return chrom, pos, strand, names, tlen


def score(pt, st, proj, near_mb=2.0, min_mb=1.0, chrom_prefix="SUPER"):
    chrom, pos, strand, names, _ = proj
    ms = MapState(pt).replay(st.edits) if st is not None else MapState(pt)
    paint = dict(st.scaffs) if st is not None else {}
    bpp = pt.bp_per_pix
    hap = lambda s: pt.names[s].split("_")[0]
    # painted chromosomes: consecutive contigs with one paint id (unpainted: each its own)
    groups, cur, last = [], None, None
    for i, c in enumerate(ms.contigs):
        pid = paint.get(i) if st is not None else ("asm", c.scaffold)
        if pid is None:
            pid = ("unpainted", i)
        if pid != last:
            cur = []
            groups.append(cur)
            last = pid
        cur.append(c)
    pur, tot = 0.0, 0.0
    split = collections.defaultdict(set)
    right, wrong_chr, wrong_ord, joins = 0, 0, 0, []
    absent, unloc, unplaced = 0, 0, 0
    adj = set()
    # chromosome-level sequences of the curation (SUPER_N by the ToL convention); the rest are
    # pieces it left unplaced
    is_chrom = np.array([n.startswith(chrom_prefix) for n in names])
    # contigs tagged Unloc in the savestate (they go at the end of a chromosome by convention)
    unloc_at = set()
    if st is not None:
        for i, mask in st.meta_flags:
            if "Unloc" in st.tag_names(mask) and i < len(ms.contigs):
                unloc_at.add(ms.contigs[i].start)
    # a curated unloc belongs to its chromosome: SUPER_11_unloc_1 -> SUPER_11
    base = np.array([names.index(n) if n in names else -1 for n in
                     (re.sub(r"_unloc_\d+$", "", n) for n in names)])
    base = np.where(base >= 0, base, np.arange(len(names)))
    hap_px = np.array([pt.names[int(s)].split("_")[0] for s in pt.orig_scaffold])
    end_px = max(3, int(1e6 / bpp))                   # look at the last ~1 Mb of each side

    def side_end(c, at_end):
        """(chromosome, position) of the side of contig c facing a junction: the chromosome most of
        its last ~1 Mb aligns to (unlocs counted with their chromosome), and the position of the
        pixels on it nearest the junction."""
        n = min(end_px, c.length)
        lay = np.arange(c.start + c.length - n, c.start + c.length) if at_end else np.arange(c.start, c.start + n)
        px = ms.pix[lay]
        ok = chrom[px] >= 0
        if not ok.any():
            return None
        b = base[chrom[px[ok]]]
        main = np.bincount(b).argmax()
        sel = px[ok][b == main]
        near = sel[-5:] if at_end else sel[:5]
        k = near[len(near) // 2]
        return int(main), float(pos[k]), int(chrom[k])
    for gi, g in enumerate(groups):
        px = np.concatenate([ms.pix[c.start:c.start + c.length] for c in g])
        cc = chrom[px]
        ok = cc >= 0
        if ok.sum():
            vals, cnt = np.unique(cc[ok], return_counts=True)
            pur += cnt.max() * bpp
            tot += ok.sum() * bpp
            for v, k in zip(vals, cnt):
                if k * bpp >= 1e6:
                    split[(v, hap(int(pt.orig_scaffold[px[0]])))].add(gi)
        for c1, c2 in zip(g, g[1:]):
            s1, s2 = int(c1.scaffold), int(c2.scaffold)
            adj.add(frozenset((s1, s2)))
            if c1.start in unloc_at or c2.start in unloc_at:
                continue                               # an unloc sits at the chromosome's end by convention
            e1, e2 = side_end(c1, True), side_end(c2, False)
            if s1 == s2 or e1 is None or e2 is None:
                continue
            (b1, x1, t1), (b2, x2, t2) = e1, e2
            if b1 != b2 and not (is_chrom[b1] and is_chrom[b2]):
                unplaced += 1                          # the curation left one side unplaced: it says nothing
                joins.append(("involves sequence the curation left unplaced", s1, s2, names[t1], names[t2]))
            elif b1 != b2:
                wrong_chr += 1
                joins.append(("other chromosome", s1, s2, names[t1], names[t2]))
            elif t1 != t2:
                unloc += 1                             # a curated unloc placed next to its chromosome
                joins.append(("unloc placed (same chromosome)", s1, s2, names[t1], names[t2]))
            elif abs(x1 - x2) > near_mb * 1e6:
                # far apart: wrong only if the stretch between them holds sequence of this haplotype
                # that the layout put somewhere else (else the haplotype lacks it: order right)
                lo, hi = sorted((x1, x2))
                # the two end windows that place the junction do not count; the rest of either
                # side does (a piece put at the wrong end of a scaffold skips that scaffold)
                n1, n2 = min(end_px, c1.length), min(end_px, c2.length)
                ends = set(ms.pix[c1.start + c1.length - n1:c1.start + c1.length].tolist()) | \
                    set(ms.pix[c2.start:c2.start + n2].tolist())
                h = hap_px[ms.pix[c1.start]]
                inb = np.nonzero((chrom == t1) & (pos > lo) & (pos < hi) & (hap_px == h))[0]
                skipped = sum(1 for p in inb.tolist() if p not in ends) * bpp
                if skipped >= max(1e6, 0.25 * (hi - lo)):
                    wrong_ord += 1
                    joins.append((f"far apart ({(hi - lo) / 1e6:.0f} Mb; skips {skipped / 1e6:.1f} Mb of this haplotype)",
                                  s1, s2, names[t1], names[t2]))
                else:
                    absent += 1
                    joins.append((f"order right, {(hi - lo) / 1e6:.0f} Mb between them absent from this haplotype",
                                  s1, s2, names[t1], names[t2]))
            else:
                right += 1
    # missed: neighbours along a curated chromosome (same haplotype, both input scaffolds >= min_mb)
    missed = []
    big = {s for s in range(len(pt.names)) if pt.scaffold_npix[s] * bpp >= min_mb * 1e6}
    for h in sorted({hap(s) for s in range(len(pt.names))}):
        for t in range(len(names)):
            sel = np.nonzero((chrom == t) & np.array([hap(int(s)) == h for s in pt.orig_scaffold]))[0]
            if not len(sel):
                continue
            order = sel[np.argsort(pos[sel])]
            seq = [int(pt.orig_scaffold[p]) for p in order]
            # a scaffold counts on this chromosome only with >= min_mb aligned to it (scattered
            # repeat hits of a scaffold that lies elsewhere do not make it a neighbour)
            cnt = collections.Counter(seq)
            seq = [s for s in seq if cnt[s] * bpp >= min_mb * 1e6]
            runs = [s for k, s in enumerate(seq) if k == 0 or s != seq[k - 1]]
            for s1, s2 in zip(runs, runs[1:]):
                if s1 in big and s2 in big and s1 != s2 and frozenset((s1, s2)) not in adj:
                    missed.append((h, names[t], min(s1, s2), max(s1, s2)))
    missed = sorted(set(missed))
    nsplit = sum(1 for v in split.values() if len(v) > 1)
    return dict(groups=len(groups), purity=pur / max(tot, 1), right=right, wrong_chromosome=wrong_chr, wrong_order=wrong_ord,
                absent=absent, unloc=unloc, unplaced=unplaced,
                split_chromosomes=nsplit, missed=missed, joins=joins, ms=ms)


def plot(pt, ms, proj, out, title):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    chrom, pos, _, names, tlen = proj
    offs, acc = {}, 0.0
    for t in sorted(range(len(names)), key=lambda t: -tlen[t]):
        offs[t] = acc
        acc += tlen[t]
    x = np.arange(pt.npix)
    p = ms.pix
    ok = chrom[p] >= 0
    y = np.array([offs[int(t)] for t in chrom[p][ok]]) + pos[p][ok]
    hap = np.array([pt.names[int(pt.orig_scaffold[q])].startswith("H1") for q in p[ok]])
    fig, ax = plt.subplots(figsize=(9, 9))
    ax.scatter(x[ok][hap] * pt.bp_per_pix / 1e9, y[hap] / 1e9, s=0.3, c="tab:blue", label="H1")
    ax.scatter(x[ok][~hap] * pt.bp_per_pix / 1e9, y[~hap] / 1e9, s=0.3, c="tab:orange", label="H2")
    for t, o in offs.items():
        if tlen[t] > 5e6:
            ax.axhline(o / 1e9, color="0.85", lw=0.4)
    ax.set_xlabel("layout position (Gb)")
    ax.set_ylabel("curated assembly position (Gb, chromosomes largest first)")
    ax.set_title(title)
    ax.legend(markerscale=20, loc="lower right")
    plt.tight_layout()
    plt.savefig(out, dpi=120)


def read_paf(paf, min_len=0):
    """PAF rows by query name: (qs, qe, strand, tname, tlen, ts, te, matches, blocklen, mapq)."""
    rows = collections.defaultdict(list)
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if int(c[10]) < min_len:
                continue
            rows[c[0]].append((int(c[2]), int(c[3]), c[4], c[5], int(c[6]), int(c[7]), int(c[8]),
                               int(c[9]), int(c[10]), int(c[11])))
    return rows


def group_names(pt, ms, st):
    """Painted chromosomes (as in `score`) named H1_chr01..., unpainted pieces by scaffold."""
    paint = dict(st.scaffs) if st is not None else {}
    out, last, nchr, seen = [], None, collections.Counter(), collections.Counter()
    for i, c in enumerate(ms.contigs):
        pid = paint.get(i) if st is not None else ("asm", c.scaffold)
        if pid is None:
            pid = ("unpainted", i)
        if pid != last:
            sname = pt.names[c.scaffold]
            if st is not None and not (isinstance(pid, tuple) and pid[0] == "unpainted"):
                h = sname.split("_")[0]
                nchr[h] += 1
                name = f"{h}_chr{nchr[h]:02d}"
            else:
                seen[sname] += 1
                name = sname if seen[sname] == 1 else f"{sname}_part{seen[sname]}"
            out.append([name, []])
            last = pid
        out[-1][1].append(c)
    return out


def write_dgenies(pt, ms, st, rows, prefix, label):
    """Lift every PAF block through the layout onto its chromosomes (query) and write the PAF
    and query index for D-GENIES' "plot alignment" mode (target index: write_target_idx)."""
    bpp = pt.bp_per_pix
    slen = {s: int(pt.scaffold_npix[s]) * bpp for s in range(len(pt.names))}
    groups = group_names(pt, ms, st)
    lens = {}
    with open(f"{prefix}.{label}.paf", "w") as paf:
        recs = []
        for name, contigs in groups:
            off = 0.0
            for c in contigs:
                s = c.scaffold
                lo, hi = c.lo * bpp, min((c.hi + 1) * bpp, slen[s])
                for qs, qe, strand, tn, tl, ts, te, m, bl, mq in rows.get(pt.names[s], ()):
                    a, b = max(qs, lo), min(qe, hi)
                    if b <= a:
                        continue
                    f0, f1 = (a - qs) / (qe - qs), (b - qs) / (qe - qs)
                    if strand == "+":
                        t0, t1 = ts + f0 * (te - ts), ts + f1 * (te - ts)
                    else:
                        t0, t1 = te - f1 * (te - ts), te - f0 * (te - ts)
                    if c.inverted:
                        q0, q1 = off + (hi - b), off + (hi - a)
                        strand2 = "-" if strand == "+" else "+"
                    else:
                        q0, q1 = off + (a - lo), off + (b - lo)
                        strand2 = strand
                    k = (b - a) / (qe - qs)
                    recs.append((name, int(q0), int(q1), strand2, tn, tl, int(t0), int(t1),
                                 max(1, int(m * k)), max(1, int(bl * k)), mq))
                off += max(hi - lo, 0)
            lens[name] = int(off)
        for r in recs:
            n = r[0]
            paf.write("\t".join(map(str, (n, lens[n], min(r[1], lens[n]), min(r[2], lens[n]), *r[3:]))) + "\n")
    with open(f"{prefix}.{label}.query.idx", "w") as f:
        f.write(f"{label}\n")
        for name, _ in groups:
            f.write(f"{name}\t{lens[name]}\n")


def write_target_idx(rows, prefix, name="curated"):
    tl = {}
    for v in rows.values():
        for r in v:
            tl[r[3]] = r[4]
    with open(f"{prefix}.target.idx", "w") as f:
        f.write(f"{name}\n")
        for t, ln in sorted(tl.items(), key=lambda x: -x[1]):
            f.write(f"{t}\t{ln}\n")


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("map")
    ap.add_argument("paf")
    ap.add_argument("savestates", nargs="*")
    ap.add_argument("--near-mb", type=float, default=2.0)
    ap.add_argument("--min-mb", type=float, default=1.0)
    ap.add_argument("--plot", help="PNG prefix: one dot plot per layout")
    ap.add_argument("--list", action="store_true", help="list wrong and missed joins")
    ap.add_argument("--dgenies", metavar="PREFIX", help="per layout, a PAF lifted onto its chromosomes and "
                    "query/target .idx files for D-GENIES (plot alignment mode)")
    ap.add_argument("--dgenies-min-len", type=int, default=5000, help="drop blocks shorter than this (bp)")
    a = ap.parse_args(argv)
    pt = Pretext(a.map)
    proj = project(pt, a.paf)
    rows = None
    if a.dgenies:
        rows = read_paf(a.paf, a.dgenies_min_len)
        write_target_idx(rows, a.dgenies)
    cov = (proj[0] >= 0).mean()
    print(f"{100 * cov:.1f}% of the map's pixels placed on the curated assembly")
    for k, path in enumerate([None] + a.savestates):
        st = SaveState.load(path) if path else None
        r = score(pt, st, proj, a.near_mb, a.min_mb)
        label = os.path.basename(path) if path else "as assembled"
        print(f"{label}: {r['groups']} chromosomes/pieces, purity {100 * r['purity']:.1f}% of bp, junctions right "
              f"{r['right']} (+{r['absent']} across sequence this haplotype lacks), unloc placed {r['unloc']}, with curated-unplaced {r['unplaced']}, wrong "
              f"chromosome {r['wrong_chromosome']}, wrong order {r['wrong_order']}, curated chromosomes "
              f"split {r['split_chromosomes']}, missed joins {len(r['missed'])}")
        if a.list:
            for kind, s1, s2, c1, c2 in r["joins"]:
                lab = "wrong" if kind.startswith(("other", "far")) else "note"
                print(f"   {lab}: {pt.names[s1]} + {pt.names[s2]}: {kind} ({c1} / {c2})")
            for h, t, s1, s2 in r["missed"]:
                print(f"   missed: {pt.names[s1]} + {pt.names[s2]} (adjacent on {t})")
        if a.plot:
            plot(pt, r["ms"], proj, f"{a.plot}.{'asm' if path is None else ('tool' if len(a.savestates) == 1 else f'layout{k}')}.png", label)
        if a.dgenies:
            tag = "asm" if path is None else f"layout{k}"
            write_dgenies(pt, r["ms"], st, rows, a.dgenies, tag)
            print(f"   D-GENIES: {a.dgenies}.{tag}.paf + {a.dgenies}.{tag}.query.idx (target {a.dgenies}.target.idx)")


if __name__ == "__main__":
    main()
