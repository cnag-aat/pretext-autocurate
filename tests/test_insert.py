"""allhaps.insert: cuts that let a loose scaffold go inside its own haplotype's copy (no fixtures)."""
import os
import sys
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate.allhaps.common import Pair          # noqa: E402
from pretext_autocurate.allhaps.insert import insert_cuts   # noqa: E402
from pretext_autocurate.homology import Block               # noqa: E402

BPP = 100_000          # bp per pixel


class _Raw:
    """Just what insert_cuts reads from a map: names, pixel ranges, bp per pixel, gap track."""

    def __init__(self, scaffolds, gaps_bp):
        self.names = [n for n, _ in scaffolds]
        self.scaffold_npix = np.array([bp // BPP for _, bp in scaffolds])
        self.scaffold_start = np.concatenate([[0], np.cumsum(self.scaffold_npix)[:-1]])
        self.bp_per_pix = BPP
        gap = np.zeros(int(self.scaffold_npix.sum()))
        for name, bp in gaps_bp:
            i = self.names.index(name)
            gap[int(self.scaffold_start[i] + bp // BPP)] = 1
        self.graphs = {"gap": gap}


def _blk(q, qs, qe, t, ts, te):
    return Block(q, qs, qe, t, ts, te, True, qe - qs, 60)


def _case(s_jump_after, gap_bp, s_extra=0):
    """H1_O (50 Mb) and its homolog H2_S, which lacks O 20-22 Mb; H2_v (2 Mb) aligns there.
    S aligns to O 0-20 Mb up to s_jump_after, then (after s_extra bp of its own) to O 22-50 Mb."""
    s_len = 48_000_000 + s_extra
    raw = _Raw([("H1_O", 50_000_000), ("H2_S", s_len), ("H2_v", 2_000_000)], [("H2_S", gap_bp)])
    a = s_jump_after
    blocks = [_blk("H2_S", 0, a, "H1_O", 20_000_000 - a, 20_000_000),
              _blk("H2_S", a + s_extra, s_len, "H1_O", 22_000_000, 22_000_000 + s_len - a - s_extra),
              _blk("H2_v", 0, 2_000_000, "H1_O", 20_000_000, 22_000_000)]
    pairs = [Pair(0, 1, 90.0, 10.0, True, "pair", "")]
    return insert_cuts(raw, pairs, blocks, r"^(H\d+)_", [])


class TestInsertCuts(unittest.TestCase):
    def test_gap_at_the_jump_is_cut(self):
        sites, rows = _case(20_000_000, 20_050_000)
        self.assertEqual(len(sites), 1)
        self.assertEqual(sites[0]["action"], "cut")
        self.assertEqual(sites[0]["s"], 1)
        self.assertIn("H2_v", sites[0]["reason"])

    def test_gap_a_megabase_off_is_not_cut(self):
        # icLarUrsu H2_scaffold_8: the jump at 29.65 Mb, the gap at 30.67 Mb
        sites, _ = _case(20_000_000, 21_000_000)
        self.assertEqual(sites, [])

    def test_own_sequence_in_the_hole_is_not_cut(self):
        # bAegFun H1_scaffold_6: S has more of its own sequence between the flanks than half the hole
        sites, _ = _case(20_000_000, 20_500_000, s_extra=1_800_000)
        self.assertEqual(sites, [])


if __name__ == "__main__":
    unittest.main()
