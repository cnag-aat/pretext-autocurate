"""Inversions inside a scaffold (bow-tie), end-to-end chain links, and track names."""
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext                                  # noqa: E402
from pretext_autocurate.curate import _link, _unlink, curate_single     # noqa: E402
from pretext_autocurate.inversion import detect_inversions              # noqa: E402
from pretext_autocurate.mapstate import MapState                        # noqa: E402
from pretext_autocurate.pretextio import standard_graph_names           # noqa: E402
from pretext_autocurate.synth import toy_map                            # noqa: E402

LAR = os.path.join(ROOT, "seed", "fixtures", "icLarUrsu.pri.yahs.mq0.pretext")
CHROMS = {"A": 300, "B": 150}
# scaffold_1 holds chromosome A with 120-220 assembled backwards; scaffold_2 is intact
TOY = [("scaffold_1", [("A", 1, 0, 120, False), ("A", 1, 120, 220, True), ("A", 1, 220, 300, False)]),
       ("scaffold_2", [("B", 1, 0, 150, False)])]


class ToyInversion(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.tmp = tempfile.TemporaryDirectory()
        path = os.path.join(cls.tmp.name, "toy.pretext")
        cls.truth = toy_map(path, TOY, CHROMS, npix=512)
        cls.pt = Pretext(path)

    @classmethod
    def tearDownClass(cls):
        cls.pt.close()
        cls.tmp.cleanup()

    def chrom_a(self, st):
        ms = MapState(self.pt).replay(st.edits)
        return [self.truth[int(q)][2] for q in ms.pix if self.truth[int(q)] and self.truth[int(q)][0] == "A"]

    def test_detected(self):
        calls = detect_inversions(self.pt)
        self.assertEqual([c.scaffold for c in calls], ["scaffold_1"])
        c = calls[0]
        a = int(self.pt.scaffold_start[0])
        # the segment is the reversed piece, give or take the gap pixels at its ends
        self.assertLessEqual(abs((c.start_px - a) - 137), 3)
        self.assertLessEqual(abs((c.end_px - a) - 251), 3)
        self.assertFalse(c.at_start or c.at_end)
        self.assertGreater(c.gain, 0.8)
        self.assertLess(c.drop, 0.2)

    def test_applied_at_level_1(self):
        dec, st = curate_single(self.pt, level=1)
        self.assertEqual([d.action for d in dec if "invert" in d.action], ["invert"])
        pos = self.chrom_a(st)
        d = np.diff(pos)
        self.assertTrue((d >= 0).all() or (d <= 0).all(), "chromosome A in order after the inversion")

    def test_proposed_at_default_level(self):
        dec, st = curate_single(self.pt)
        inv = [d for d in dec if "invert" in d.action]
        self.assertEqual([d.action for d in inv], ["flag_invert"])
        self.assertIn("proposed", inv[0].evidence)
        self.assertFalse(any(d.action == "cut" for d in dec))
        self.assertFalse(np.all(np.diff(self.chrom_a(st)) >= 0))       # left as assembled


class EndLinks(unittest.TestCase):
    def test_tail_to_tail_and_head_to_head(self):
        chains = []
        self.assertTrue(_link(chains, 0, 1, "t", "t"))       # A's tail to B's tail: B reversed
        self.assertEqual(chains, [[(0, False), (1, True)]])
        self.assertTrue(_link(chains, 1, 2, "h", "h"))       # B's head to C's head
        self.assertEqual(chains, [[(0, False), (1, True), (2, False)]])

    def test_unlink(self):
        chains = [[(0, False), (1, False), (2, False)]]
        _unlink(chains, 1, "t")                               # 1's tail faces 2
        self.assertEqual(chains, [[(0, False), (1, False)]])
        chains = [[(0, False), (1, True), (2, False)]]
        _unlink(chains, 1, "t")                               # reversed: 1's tail faces 0
        self.assertEqual(chains, [[(1, True), (2, False)]])
        _unlink(chains, 5, "h")                               # not in a chain: nothing
        self.assertEqual(chains, [[(1, True), (2, False)]])


class AllhapsInversion(unittest.TestCase):
    """allhaps: projection onto the homolog (folding guard) and reading the arrangement."""

    def test_project_through_the_alignment(self):
        from pretext_autocurate.allhaps.inversions import _project
        from pretext_autocurate.homology import Block
        fwd = Block("H2_s1", 1_000_000, 3_000_000, "H1_s1", 5_000_000, 7_000_000, True, 1, 60)
        rev = Block("H2_s1", 1_000_000, 3_000_000, "H1_s1", 5_000_000, 7_000_000, False, 1, 60)
        self.assertEqual(_project([fwd], "H2_s1", 1_500_000, 2_000_000), [("H1_s1", 5_500_000, 6_000_000)])
        self.assertEqual(_project([rev], "H2_s1", 1_500_000, 2_000_000), [("H1_s1", 6_000_000, 6_500_000)])
        self.assertEqual(_project([fwd], "H1_s1", 6_000_000, 6_500_000), [("H2_s1", 2_000_000, 2_500_000)])
        self.assertEqual(_project([fwd], "H1_s9", 0, 10), [])

    def test_homolog_direction(self):
        from pretext_autocurate.allhaps.inversions import _homolog_says
        from pretext_autocurate.homology import Block
        B = lambda qs, qe, ts, te, same: Block("H1_s1", qs, qe, "H2_s1", ts, te, same, 1, 60)
        flanks = [B(0, 4_000_000, 0, 4_000_000, True), B(6_000_000, 10_000_000, 6_000_000, 10_000_000, True)]
        inverted = flanks + [B(4_000_000, 6_000_000, 4_000_000, 6_000_000, False)]
        same = flanks + [B(4_000_000, 6_000_000, 4_000_000, 6_000_000, True)]
        self.assertTrue(_homolog_says(inverted, "H1_s1", 4_000_000, 6_000_000))
        self.assertFalse(_homolog_says(same, "H1_s1", 4_000_000, 6_000_000))
        self.assertIsNone(_homolog_says(flanks, "H1_s1", 4_000_000, 6_000_000))     # the segment does not align


class TrackNames(unittest.TestCase):
    def test_old_names(self):
        g = {"GAPs": 1, "ONTcov": 2, "TELOMERES": 3, "5p_telomere": 4}
        names = standard_graph_names(g)
        self.assertEqual(set(g), {"gap", "coverage", "telomere", "5p_telomere"})
        self.assertEqual(names["gap"], "GAPs")

    def test_standard_names_win(self):
        g = {"gap": 1, "GAPs": 2, "coverage": 3}
        standard_graph_names(g)
        self.assertEqual(g, {"gap": 1, "GAPs": 2, "coverage": 3})


@unittest.skipUnless(os.path.exists(LAR), "missing fixture")
class LarinusInversion(unittest.TestCase):
    """icLarUrsu.pri.yahs: the curator reversed scaffold_13 around 41-46 Mb (and made 39.6-41.2 an
    Unloc); scaffold_4 68.8-78.2 Mb shows both orientations (heterozygous?) and was left alone."""

    def test_scaffold_13_called_scaffold_4_not(self):
        pt = Pretext(LAR)
        self.assertIn("gap", pt.graphs)                 # 'GAPs' in this older map
        calls = {c.scaffold: c for c in detect_inversions(pt)}
        c = calls.get("scaffold_13")
        self.assertIsNotNone(c)
        mb = lambda q: (q - int(pt.scaffold_start[pt.names.index(c.scaffold)])) * pt.bp_per_pix / 1e6
        self.assertAlmostEqual(mb(c.start_px), 39.7, delta=0.6)
        self.assertAlmostEqual(mb(c.end_px + 1), 46.1, delta=0.6)
        self.assertGreater(c.gain, 1.0)
        self.assertLess(c.drop, 0.3)
        self.assertNotIn("scaffold_4", calls)
        pt.close()


if __name__ == "__main__":
    unittest.main()
