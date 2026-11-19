"""Telomeres at cut ends and keeping haplotig flanks together (rules pass, single mode)."""
import os
import sys
import types
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate.curate import _link                     # noqa: E402
from pretext_autocurate.join import Params, telomeric_ends      # noqa: E402
from pretext_autocurate.pretextio import SegmentedMap           # noqa: E402


def raw_map(n_per=200, n_scaf=2):
    n = n_per * n_scaf
    so = np.repeat(np.arange(n_scaf), n_per)
    return types.SimpleNamespace(npix=n, n_scaffolds=n_scaf, names=[f"s{i}" for i in range(n_scaf)],
                                 orig_scaffold=so, rel_coord=np.arange(n) - so * n_per, graphs={})


class TelomereAtCut(unittest.TestCase):
    def setUp(self):
        self.raw = raw = raw_map()
        t5, t3 = np.zeros(raw.npix), np.zeros(raw.npix)
        t5[99:101] = 1000       # a left chromosome end at the cut (after pixel 99), in both pixels
        t5[199] = 1000          # wrong strand at an original scaffold end: still counts
        raw.graphs.update(telomere_5p=t5, telomere_3p=t3, telomere=t5 + t3)
        self.seg = SegmentedMap(raw, [99])      # pieces: s0:1 (0-99), s0:2 (100-199), s1
        self.p = Params(tel_px=3)

    def test_strand_decides_at_a_cut(self):
        tel = telomeric_ends(self.seg, [0, 1, 2], self.p)
        self.assertIn((1, "h"), tel)            # its own telomere
        self.assertNotIn((0, "t"), tel)         # the piece before the cut: not an end
        self.assertIn((1, "t"), tel)            # original scaffold end: strand not asked

    def test_without_strand_tracks_both_sides_count(self):
        g = self.raw.graphs
        del g["telomere_5p"], g["telomere_3p"]
        tel = telomeric_ends(self.seg, [0, 1, 2], self.p)
        self.assertIn((0, "t"), tel)
        self.assertIn((1, "h"), tel)


class Link(unittest.TestCase):
    def test_free_ends(self):
        chains = [[(0, False), (1, False)], [(3, False), (2, True)]]     # 2's head is last
        self.assertTrue(_link(chains, 1, 2))
        self.assertEqual(chains, [[(0, False), (1, False), (2, False), (3, True)]])

    def test_reversed_and_singletons(self):
        chains = [[(1, True), (0, True)]]       # 1's tail is the first end
        self.assertTrue(_link(chains, 1, 5))
        self.assertEqual(chains, [[(0, False), (1, False), (5, False)]])

    def test_joined_end_left_alone(self):
        chains = [[(0, False), (1, False)], [(2, False), (3, False)]]
        before = [list(c) for c in chains]
        self.assertFalse(_link(chains, 0, 2))   # 0's tail is joined to 1
        self.assertFalse(_link(chains, 0, 1))   # already one chain
        self.assertEqual(chains, before)


class SingleLevels(unittest.TestCase):
    def test_level_moves_insertion_thresholds(self):
        from pretext_autocurate.curate import single_thresholds
        d, c, a = single_thresholds(0.5), single_thresholds(0), single_thresholds(1)
        self.assertEqual((d["insert_orient"], d["orient_strict"], d["insert_min_px"]), (1.2, 1.2, 15))
        self.assertTrue(c["insert_orient"] > d["insert_orient"] > a["insert_orient"])
        self.assertTrue(c["insert_min_px"] > d["insert_min_px"] > a["insert_min_px"])
        self.assertIsInstance(single_thresholds(0.3)["insert_min_px"], int)


if __name__ == "__main__":
    unittest.main()
