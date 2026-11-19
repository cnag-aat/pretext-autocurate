"""Contig-end overlaps (haplotypic copies at scaffold ends), exact breakpoints and gapless joins."""
import os
import random
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import breaks                                    # noqa: E402
from pretext_autocurate.graph import AgpRow, Faidx, read_agp, write_agp_rows   # noqa: E402
from pretext_autocurate.overlap import _copy_chain, _match_runs          # noqa: E402

LAR = os.path.join(ROOT, "seed", "fixtures", "icLarUrsu.pri.yahs.mq0.pretext")
LAR_FA = os.path.join(ROOT, "seed", "fixtures", "icLarUrsu.pri.yahs.fa")
GAP = ("scaffold", "yes", "proximity_ligation")


def rand_seq(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def mutate(s, at):
    s = list(s)
    for k in at:
        s[k] = "A" if s[k] != "A" else "C"
    return "".join(s)


def write_fa(path, recs):
    with open(path, "w") as f:
        for n, s in recs.items():
            f.write(f">{n}\n" + "".join(s[k:k + 60] + "\n" for k in range(0, len(s), 60)))


class MatchRuns(unittest.TestCase):
    def test_plus(self):
        # 5 matches, a substitution, 12 matches, 2 inserted query bases, 4 matches
        runs = _match_runs(100, 124, 1000, "+", ":5*ag:12+tt:4", 10)
        self.assertEqual(runs, [(106, 1006, 12)])

    def test_minus(self):
        # reverse strand: the query runs down from its end while the target runs up
        runs = _match_runs(100, 130, 1000, "-", ":10-ac:15", 10)
        self.assertEqual(runs, [(129, 1000, 10), (119, 1012, 15)])


class CopyChain(unittest.TestCase):
    def test_first_meeting_with_the_kept_end(self):
        # (a0, a1, b0, b1, len, de): a = depth from the trimmed end, b = depth from the kept end.
        # The copy runs to the kept end at a = 600; a repeat of its last part meets it again at 900.
        k = 100_000                     # real sizes: the chain allows 5 kb of slack
        items = [(0, 2 * k, 4 * k, 6 * k, 2 * k, 0.01), (2 * k, 4 * k, 2 * k, 4 * k, 2 * k, 0.01),
                 (4 * k, 6 * k, 0, 2 * k, 2 * k, 0.01), (7 * k, 9 * k, 0, 2 * k, 2 * k, 0.01)]
        ch = _copy_chain(items, max_gap=10 * k)
        self.assertEqual(ch[-1][1], 6 * k)
        self.assertEqual(len(ch), 3)


class RefineCuts(unittest.TestCase):
    def test_moves_both_lines(self):
        with tempfile.TemporaryDirectory() as d:
            agp = os.path.join(d, "p.agp")
            head = ["##agp-version\t2.1\n", "# HiC MAP RESOLUTION: 1000.000000 bp/texel\n"]
            objs = {"Scaffold_1": [AgpRow("Scaffold_1", "W", "s1", 1, 5000, "+", ("Painted",))],
                    "Scaffold_2": [AgpRow("Scaffold_2", "W", "s1", 5001, 9000, "+", ())]}
            write_agp_rows(agp, head, objs)
            rep = breaks.refine_cuts(agp, agp + ".out", [dict(kind="cut", scaffold="s1", px_bp=5000, bp=5321)], 2002)
            self.assertTrue(rep[0][1])
            _, o = read_agp(agp + ".out")
            self.assertEqual((o["Scaffold_1"][0].end, o["Scaffold_2"][0].beg), (5321, 5322))
            # a cut the savestate no longer has is left alone
            rep = breaks.refine_cuts(agp, agp + ".out2", [dict(kind="cut", scaffold="s1", px_bp=7000, bp=7100)], 2002)
            self.assertFalse(rep[0][1])


class Switches(unittest.TestCase):
    def test_overlap_between_two_scaffolds(self):
        rng = random.Random(1)
        K = rand_seq(2000, rng)
        T = mutate(K[1500:], [50, 450]) + rand_seq(1000, rng)     # T's first 500 bp copy K's last 500
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "in.fa")
            write_fa(inp, {"K": K, "T": T})
            fa = Faidx(inp)
            # pretext-to-asm-like output: K, a gap, T past its copy; the copy as haplotig H_1
            prim = {"SUPER_1": [AgpRow("SUPER_1", "W", "K", 1, 2000, "+", ()), AgpRow("SUPER_1", "U", "100", 0, 0, "", GAP),
                                AgpRow("SUPER_1", "W", "T", 501, 1500, "+", ())]}
            hap = {"H_1": [AgpRow("H_1", "W", "T", 1, 500, "+", ())]}
            write_agp_rows(os.path.join(d, "a.1.primary.agp"), ["##agp-version\t2.1\n"], prim)
            write_agp_rows(os.path.join(d, "a.1.haplotigs.agp"), ["##agp-version\t2.1\n"], hap)
            write_fa(os.path.join(d, "a.1.primary.fa"), {"SUPER_1": K + "N" * 100 + T[500:]})
            write_fa(os.path.join(d, "a.1.haplotigs.fa"), {"H_1": T[:500]})
            sw = dict(kind="switch", scaffold="K", px_bp=2000, bp=1750, partner="T", partner_bp=250)
            rep = []
            self.assertEqual(breaks.apply_switches(d, fa, [sw], rep), 1)
            out = Faidx(os.path.join(d, "a.1.primary.fa"))
            self.assertEqual(out.fetch("SUPER_1", 0, 10 ** 6), K[:1750] + T[250:])
            h = Faidx(os.path.join(d, "a.1.haplotigs.fa"))
            self.assertEqual(h.fetch("H_1", 0, 10 ** 6), T[:250] + K[1750:])

    def test_gap_inside_a_scaffold(self):
        rng = random.Random(2)
        A = rand_seq(1000, rng)
        B = mutate(A[700:], [20, 280]) + rand_seq(700, rng)       # B's first 300 bp copy A's last 300
        X = A + "N" * 100 + B
        with tempfile.TemporaryDirectory() as d:
            inp = os.path.join(d, "in.fa")
            write_fa(inp, {"X": X})
            fa = Faidx(inp)
            prim = {"SUPER_1": [AgpRow("SUPER_1", "W", "X", 1, 1000, "+", ()), AgpRow("SUPER_1", "N", "100", 0, 0, "", GAP),
                                AgpRow("SUPER_1", "W", "X", 1101, 2100, "+", ())]}
            write_agp_rows(os.path.join(d, "a.1.primary.agp"), ["##agp-version\t2.1\n"], prim)
            write_fa(os.path.join(d, "a.1.primary.fa"), {"SUPER_1": X})
            # A[849] (0-based) is B[149], X[1249]
            sw = dict(kind="switch", scaffold="X", px_bp=1000, bp=850, partner="X", partner_bp=1250)
            rep = []
            self.assertEqual(breaks.apply_switches(d, fa, [sw], rep), 1)
            out = Faidx(os.path.join(d, "a.1.primary.fa"))
            self.assertEqual(out.fetch("SUPER_1", 0, 10 ** 6), X[:850] + X[1250:])
            h = Faidx(os.path.join(d, "a.1.haplotigs.fa"))       # created: the other path through the bubble
            self.assertEqual(h.fetch("H_1", 0, 10 ** 6), X[1100:1250] + X[850:1000])


@unittest.skipUnless(os.path.exists(LAR) and os.path.exists(LAR_FA), "missing fixture")
class LarinusOverlaps(unittest.TestCase):
    """icLarUrsu.pri.yahs: scaffold_16 0.26-6.32 Mb is a 0.1-3% divergent copy of scaffold_22's last
    6.3 Mb; all of scaffold_32 is a copy of scaffold_14's start. Nothing else passes."""

    def test_detected_and_placed_by_the_sequence(self):
        from pretext_autocurate import Pretext
        from pretext_autocurate.homology import find_tool
        from pretext_autocurate.overlap import Side, align_overlap, detect_end_overlaps
        pt = Pretext(LAR)
        got = {frozenset({(pt.names[o.a], o.a_end), (pt.names[o.b], o.b_end)}) for o in detect_end_overlaps(pt)}
        self.assertEqual(got, {frozenset({("scaffold_16", "h"), ("scaffold_22", "t")}),
                               frozenset({("scaffold_14", "h"), ("scaffold_32", "t")})})
        mm = find_tool("minimap2")
        if mm is None:
            self.skipTest("no minimap2")
        fa = Faidx(LAR_FA)
        al = align_overlap(fa, Side("scaffold_16", 0, True), Side("scaffold_22", fa.length("scaffold_22"), False),
                           int(118 * pt.bp_per_pix * 1.5 + 20 * pt.bp_per_pix), minimap2=mm)
        self.assertAlmostEqual(al.copy_bp, 6_316_846, delta=5000)
        self.assertLess(al.divergence, 0.02)
        k, t = al.keep_switch, al.trim_switch
        self.assertEqual(fa.fetch("scaffold_22", k - 200, k + 201), fa.fetch("scaffold_16", t - 200, t + 201))
        pt.close()


if __name__ == "__main__":
    unittest.main()
