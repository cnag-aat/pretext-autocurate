"""Two copies of one region inside one scaffold, from a haplotype-against-haplotype PAF
(primary alignments): two stretches of the same query scaffold that align to overlapping
stretches of one target scaffold. The query haplotype then holds two copies where the
target holds one: both alleles in one scaffold (haplotypic duplication), a false
duplication of one allele, or a paralog the target haplotype has collapsed or lacks.

    python3 tools/self_dups.py PAF [--min-bp 50000] [--fai FASTA.fai] [--gaps AGP ...]

Prints one line per duplicated stretch pair: query scaffold, the two copies (Mb), their
orientation to each other, the shared target stretch, identity of each copy to the target,
the distance between the copies, and whether a gap (`N` in the AGPs, or a YaHS scaffold
join) lies between them."""
import argparse
import sys


def read(paf, min_bp):
    out = []
    with open(paf) as f:
        for line in f:
            c = line.split("\t")
            if int(c[10]) < min_bp:
                continue
            de = next((float(x[5:]) for x in c[12:] if x.startswith("de:f:")), float("nan"))
            out.append(dict(q=c[0], qs=int(c[2]), qe=int(c[3]), strand=c[4], t=c[5], ts=int(c[7]), te=int(c[8]),
                            ident=int(c[9]) / max(int(c[10]), 1), mapq=int(c[11]), qlen=int(c[1]), de=de))
    return out


def gaps_from_agp(paths):
    """{scaffold: [gap start bp, ...]} from YaHS AGPs (lines of type N/U)."""
    g = {}
    for p in paths:
        with open(p) as f:
            for line in f:
                if line.startswith("#"):
                    continue
                c = line.rstrip("\n").split("\t")
                if len(c) > 4 and c[4] in ("N", "U"):
                    g.setdefault(c[0], []).append(int(c[1]))
    return g


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("paf")
    ap.add_argument("--min-bp", type=int, default=50000)
    ap.add_argument("--min-shared", type=int, default=50000, help="target bp both copies cover")
    ap.add_argument("--gaps", nargs="*", default=[], help="AGPs of the query haplotype (gap positions)")
    ap.add_argument("--rename", default="", help="prefix to add to AGP scaffold names, e.g. H2_")
    ap.add_argument("--min-ident", type=float, default=0.5, help="matches / alignment length, both copies")
    ap.add_argument("--map", help=".pretext map with a coverage track: depth of each copy and of the target")
    a = ap.parse_args(argv)
    depth = lambda name, s, e: float("nan")
    if a.map:
        sys.path.insert(0, __file__.rsplit("/tools/", 1)[0] + "/src")
        import numpy as np
        from pretext_autocurate import Pretext
        pt = Pretext(a.map)
        cov = pt.graphs["coverage"].astype(float)
        idx = {n: i for i, n in enumerate(pt.names)}
        largest = pt.scaffold_npix[pt.orig_scaffold] >= 0.2 * pt.scaffold_npix.max()
        norm = float(np.median(cov[largest]))

        def depth(name, s, e):
            """Median coverage over [s, e) bp of a scaffold, relative to the large scaffolds."""
            i = idx.get(name)
            if i is None or pt.scaffold_npix[i] == 0:
                return float("nan")
            a0, n = int(pt.scaffold_start[i]), int(pt.scaffold_npix[i])
            p0 = a0 + min(int(s / pt.bp_per_pix), n - 1)
            p1 = a0 + min(max(int(e / pt.bp_per_pix), p0 - a0 + 1), n)
            return float(np.median(cov[p0:p1])) / norm
    al = read(a.paf, a.min_bp)
    gaps = {a.rename + k: v for k, v in gaps_from_agp(a.gaps).items()}
    by = {}
    for x in al:
        by.setdefault((x["q"], x["t"]), []).append(x)
    found = []
    for (q, t), xs in by.items():
        xs.sort(key=lambda x: x["qs"])
        for i, x in enumerate(xs):
            for y in xs[i + 1:]:
                if y["qs"] < x["qe"]:                     # overlapping in the query: same copy
                    continue
                sh = min(x["te"], y["te"]) - max(x["ts"], y["ts"])
                if sh < a.min_shared or min(x["ident"], y["ident"]) < a.min_ident:
                    continue
                between = [p for p in gaps.get(q, []) if x["qe"] <= p <= y["qs"]]
                found.append((q, x, y, sh, between))
    found.sort(key=lambda f: (f[0], f[1]["qs"]))
    def project(b, t):
        """Query position of target position t along block b (linear, strand-aware)."""
        f = (t - b["ts"]) / max(b["te"] - b["ts"], 1)
        return b["qs"] + f * (b["qe"] - b["qs"]) if b["strand"] == "+" else b["qe"] - f * (b["qe"] - b["qs"])

    print("query\tqlen_mb\tcopy1_mb\tcopy2_mb\trel_orient\ttarget\tshared_kb\tident1\tident2\tmapq\t"
          "dist_kb\tgaps_between\tdepth1\tdepth2\tdepth_target")
    for q, x, y, sh, between in found:
        t0, t1 = max(x["ts"], y["ts"]), min(x["te"], y["te"])
        c1 = sorted((project(x, t0), project(x, t1)))
        c2 = sorted((project(y, t0), project(y, t1)))
        print(f"{q}\t{x['qlen'] / 1e6:.1f}\t{c1[0] / 1e6:.2f}-{c1[1] / 1e6:.2f}\t{c2[0] / 1e6:.2f}-{c2[1] / 1e6:.2f}\t"
              f"{'same' if x['strand'] == y['strand'] else 'inverted'}\t{x['t']}:{t0 / 1e6:.2f}-{t1 / 1e6:.2f}\t"
              f"{sh / 1e3:.0f}\t{x['ident']:.3f}\t{y['ident']:.3f}\t"
              f"{x['mapq']}/{y['mapq']}\t{(c2[0] - c1[1]) / 1e3:.0f}\t{len(between)}\t"
              f"{depth(q, *c1):.2f}\t{depth(q, *c2):.2f}\t{depth(x['t'], t0, t1):.2f}")
    print(f"# {len(found)} pairs in {len({f[0] for f in found})} query scaffolds", file=sys.stderr)


if __name__ == "__main__":
    main()
