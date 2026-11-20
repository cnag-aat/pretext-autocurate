"""allhaps layout fixes: which way the second copy of a pair runs (alignment of the arranged
chains, else the Hi-C homology line)."""
import os
import sys
import types
import unittest

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate.allhaps.shrapnel import hic_dir, orient_pair    # noqa: E402
from pretext_autocurate.homology import Block                           # noqa: E402


def fake_raw():
    # two scaffolds of 100 px each, 10 kb per pixel
    return types.SimpleNamespace(npix=200, names=["H1_a", "H2_b"], bp_per_pix=10_000.0,
                                 scaffold_start=np.array([0, 100]), scaffold_npix=np.array([100, 100]))


class OrientPair(unittest.TestCase):
    def test_follows_the_arranged_chains(self):
        raw = fake_raw()
        same = [Block("H2_b", 0, 1_000_000, "H1_a", 0, 1_000_000, True, 1, 60)]
        a = np.arange(0, 100)
        b = np.arange(100, 200)
        self.assertIs(orient_pair(raw, same, a, b), False)        # laid out as assembled: parallel
        self.assertTrue(orient_pair(raw, same, a, b[::-1]))       # the second chain reversed: anti-parallel
        rev = [Block("H2_b", 0, 1_000_000, "H1_a", 0, 1_000_000, False, 1, 60)]
        self.assertTrue(orient_pair(raw, rev, a, b))
        self.assertIsNone(orient_pair(raw, [], a, b))
        self.assertIsNone(orient_pair(raw, [Block("H2_b", 0, 50_000, "H1_a", 0, 50_000, True, 1, 60)], a, b))


class HicDir(unittest.TestCase):
    def sg(self, M):
        return types.SimpleNamespace(block=lambda r0, r1, c0, c1: M[r0:r1 + 1, c0:c1 + 1])

    def test_line_direction(self):
        n = 60
        M = np.ones((2 * n, 2 * n), dtype=np.float32)
        for i in range(n):                                        # a parallel copy line
            M[i, n + i] = M[n + i, i] = 50
        a, b = np.arange(n), np.arange(n, 2 * n)
        self.assertIs(hic_dir(self.sg(M), a, b), False)
        self.assertTrue(hic_dir(self.sg(M), a, b[::-1]))
        flat = np.ones((2 * n, 2 * n), dtype=np.float32)
        self.assertIsNone(hic_dir(self.sg(flat), a, b))


if __name__ == "__main__":
    unittest.main()
