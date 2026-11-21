"""Displaced stretches (linescan.py) on xgTylRafi: the four stretches the curator moved, each
found with its place and strand, and nothing else."""
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))

from pretext_autocurate import Pretext  # noqa: E402
from pretext_autocurate.signal import Signal  # noqa: E402
from pretext_autocurate import linescan  # noqa: E402

F = os.path.join(ROOT, "seed", "fixtures")
TYL = os.path.join(F, "xgTylRafi_asm4_run1.hfsm.asm.hic.hap2.yhs_mq10_mq0.extensions.HR.pretext")
TYL_FA = os.path.join(F, "xgTylRafi_asm4_run1.hfsm.asm.hic.hap2.yhs_mq10.fa")


@unittest.skipUnless(os.path.exists(TYL), "missing fixture")
class TestLineScan(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(TYL)
        cls.found = linescan.scan(cls.pt, Signal(cls.pt), linescan.Params())

    def mb(self, q):
        s = int(self.pt.orig_scaffold[q])
        return self.pt.names[s], (q - self.pt.scaffold_start[s] + 1) * self.pt.bp_per_pix / 1e6

    def test_curated_moves(self):
        got = []
        for d in self.found:
            name, a = self.mb(d.a - 1)
            tname, t = self.mb(d.t)
            got.append((name, a, tname, t, d.strand))
        # (scaffold, start Mb, target, place Mb, strand) as the curator has them (savestate_4)
        want = [("scaffold_9", 14.80, "scaffold_10", 56.43, "+"),
                ("scaffold_10", 56.45, "scaffold_11", 40.80, "-"),
                ("scaffold_2", 40.49, "scaffold_5", 49.26, "+"),
                ("scaffold_3", 74.93, "scaffold_1", 166.64, "-")]
        self.assertEqual(len(got), 4, got)
        for name, a, tname, t, strand in want:
            hit = [g for g in got if g[0] == name and abs(g[1] - a) <= 0.15]
            self.assertEqual(len(hit), 1, (name, a, got))
            self.assertEqual(hit[0][2], tname)
            self.assertAlmostEqual(hit[0][3], t, delta=0.2)
            self.assertEqual(hit[0][4], strand)

    @unittest.skipUnless(os.path.exists(TYL_FA), "no FASTA")
    def test_tail_copy(self):
        """scaffold_9 15.30-15.65 Mb is a haplotypic copy of scaffold_10 just after its gap at
        56.72 Mb (the curator tagged scaffold_9 15.30-15.61 Haplotig)."""
        al = linescan.align_stretch(TYL_FA, "scaffold_9", 14820000, 15700000, "scaffold_10", 55500000, 58000000,
                                    end_tol=71000)
        if al is None:
            self.skipTest("no minimap2")
        self.assertEqual(al.strand, "+")
        self.assertIsNone(al.head)
        q0, q1, t0, t1 = al.tail
        self.assertAlmostEqual(q0 / 1e6, 15.30, delta=0.02)
        self.assertAlmostEqual(t0 / 1e6, 56.72, delta=0.02)
        self.assertIsNotNone(al.tail_switch)


if __name__ == "__main__":
    unittest.main()
