"""Joins between scaffolds that a reference curation makes and the tool does not.

    python3 tools/missed_joins.py MAP TOOL_SAVESTATE REF_SAVESTATE [--min-mb 1] [--decisions TSV]

A join is two consecutive pieces of one painted chromosome in the reference that come from
different input scaffolds. It is missed when the tool does not put the two pieces next to
each other (within `--near-px`) in one painted chromosome. For each, the two pieces (input
scaffold and Mb range) and what the tool did with each side: its chromosome (painted
neighbours) or that it is unpainted; with `--decisions`, the tool's rows for both scaffolds."""
import argparse
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate import Pretext, SaveState            # noqa: E402
from pretext_autocurate.score import Layout, _supported      # noqa: E402


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("map")
    ap.add_argument("tool")
    ap.add_argument("ref")
    ap.add_argument("--min-mb", type=float, default=1.0, help="smallest input scaffold on either side")
    ap.add_argument("--near-px", type=int, default=80)
    ap.add_argument("--decisions")
    a = ap.parse_args(argv)
    pt = Pretext(a.map)
    T, R = Layout(pt, SaveState.load(a.tool)), Layout(pt, SaveState.load(a.ref))
    bpp = pt.bp_per_pix
    mb = lambda s: pt.scaffold_npix[s] * bpp / 1e6
    rows = {}
    if a.decisions:
        with open(a.decisions) as f:
            next(f)
            for line in f:
                c = line.rstrip("\n").split("\t")
                rows.setdefault(c[1].split(":")[0], []).append(c)

    def piece(c, L):
        px = L.ms.pix[c.start:c.start + c.length]
        s = int(pt.orig_scaffold[px[0]])
        r0, r1 = sorted((int(pt.rel_coord[px[0]]), int(pt.rel_coord[px[-1]])))
        return s, f"{pt.names[s]} {r0 * bpp / 1e6:.2f}-{(r1 + 1) * bpp / 1e6:.2f}{' rev' if px[0] > px[-1] else ''}"

    def where(p, L):
        """The tool's chromosome around pixel p: names of the painted pieces in it."""
        ch = L.chrom[p]
        if ch == 0:
            return "unpainted" + (" (" + ",".join(t for t in ("Unloc", "FalseDuplicate", "Haplotig") if L.tags[t][p]) + ")"
                                  if any(L.tags[t][p] for t in ("Unloc", "FalseDuplicate", "Haplotig")) else "")
        names = []
        for c in L.ms.contigs:
            q = int(L.ms.pix[c.start])
            if L.chrom[q] == ch:
                n = pt.names[int(pt.orig_scaffold[q])]
                if not names or names[-1] != n:
                    names.append(n)
        return f"chromosome {ch}: " + " | ".join(names)

    ms = R.ms
    seq = [c for c in ms.contigs if R.chrom[ms.pix[c.start]] > 0 and not R.aside[ms.pix[c.start]]]
    n = 0
    for x, y in zip(seq, seq[1:]):
        px, py = int(ms.pix[x.start + x.length - 1]), int(ms.pix[y.start])
        if R.chrom[px] != R.chrom[py]:
            continue
        sx, tx = piece(x, R)
        sy, ty = piece(y, R)
        if sx == sy or min(mb(sx), mb(sy)) < a.min_mb:
            continue
        k = min(6, x.length // 2, y.length // 2)
        j = (int(ms.pix[x.start + x.length - 1 - k]), int(ms.pix[y.start + k]))
        if _supported(j, T, a.near_px):
            continue
        n += 1
        print(f"{n}. {tx} ({mb(sx):.1f} Mb)  +  {ty} ({mb(sy):.1f} Mb)")
        print(f"     tool, left side:  {where(px, T)}")
        print(f"     tool, right side: {where(py, T)}")
        for s in (sx, sy):
            for c in rows.get(pt.names[s], []):
                if c[0] not in ("pair", "busco_summary"):
                    print(f"       {c[0]} {c[1]} {c[3]}: {c[-1][:220]}")
    print(f"{n} joins missed (both scaffolds >= {a.min_mb} Mb)")


if __name__ == "__main__":
    main()
