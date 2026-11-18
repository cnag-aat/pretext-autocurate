"""The adjacency solver (`solver.solve`, `run --engine adjacency`) on toy genomes (no fixtures)."""
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pretext_autocurate import MapState, Pretext                 # noqa: E402
from pretext_autocurate.solver import _homology, contigs, solve   # noqa: E402
from pretext_autocurate.synth import toy_blocks, toy_map          # noqa: E402
from test_synthetic import CHROMS, DIPLOID                        # noqa: E402

# A and B in one H1 scaffold, the right way round: a YaHS join that is right
JOINED = [("H1_scaffold_1", [("A", 1, 0, 60, False), ("A", 1, 60, 100, False)]),
          ("H1_scaffold_2", [("B", 1, 0, 80, False)]),
          ("H1_scaffold_3", [("C", 1, 0, 50, False)]),
          ("H2_scaffold_1", [("A", 2, 0, 100, False)]),
          ("H2_scaffold_2", [("B", 2, 0, 80, False)]),
          ("H2_scaffold_3", [("C", 2, 0, 50, False)])]


class TestSolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.TemporaryDirectory()
        cls.maps = {}
        for name, sc in (("misjoin", DIPLOID), ("right", JOINED)):
            p = os.path.join(cls.dir.name, name + ".pretext")
            toy_map(p, sc, CHROMS)
            cls.maps[name] = (Pretext(p), toy_blocks(sc))

    @classmethod
    def tearDownClass(cls):
        cls.dir.cleanup()

    def test_contigs_split_at_gaps(self):
        pt, _ = self.maps["misjoin"]
        nodes, yahs = contigs(pt, r"^(H\d+)_")
        self.assertEqual(len(nodes), 6)                   # H1_scaffold_1 in two, the rest whole
        self.assertEqual(len(yahs), 1)
        u, v, q = yahs[0]
        self.assertEqual((nodes[u].s, nodes[v].s), (0, 0))
        self.assertEqual(nodes[u].b, q)                   # the cut goes after the gap pixel

    def test_misjoin_cut_and_pairs(self):
        pt, blocks = self.maps["misjoin"]
        rows, st, info = solve(pt, blocks=blocks)
        self.assertEqual((info["cuts"], info["joins"], info["pairs"]), (1, 0, 3))
        cut = [r for r in rows if r["action"] == "cut"]
        self.assertEqual(cut[0]["scaffold"], "H1_scaffold_1")
        self.assertIn("hic", cut[0]["evidence"])          # the row names the evidence
        ms = MapState(pt).replay(st.edits)                 # the savestate is a valid permutation
        self.assertEqual(sorted(ms.pix.tolist()), list(range(pt.npix)))
        tags = {st.tag_names(m)[0] for _, m in st.meta_flags}
        self.assertEqual(tags, {"HAP1", "HAP2"})

    def test_right_join_kept(self):
        pt, blocks = self.maps["right"]
        rows, st, info = solve(pt, blocks=blocks)
        self.assertEqual(info["cuts"], 0)
        self.assertEqual(info["pairs"], 3)
        kept = [r for r in rows if r["action"] == "kept"]
        self.assertEqual(len(kept), 1)

    def test_conservative_level_keeps_more(self):
        pt, blocks = self.maps["misjoin"]
        cuts = [solve(pt, blocks=blocks, level=lv)[2]["cuts"] for lv in (0, 0.5, 1)]
        self.assertEqual(cuts, sorted(cuts))               # never fewer cuts at a higher level

    def test_homology(self):
        olen = np.array([10_000_000, 10_000_000])
        self.assertEqual(_homology((0, 5_000_000, 1), (0, 5_200_000, -1), olen), 1.0)      # facing, near
        self.assertEqual(_homology((0, 5_000_000, 1), (0, 5_200_000, 1), olen), -1.0)      # same way
        self.assertEqual(_homology((0, 1_000_000, 1), (0, 9_000_000, -1), olen), -1.0)     # far apart
        self.assertEqual(_homology((0, 2_000_000, 1), (1, 0, -1), olen), -1.0)             # room on 0
        self.assertIsNone(_homology((0, 9_900_000, 1), (1, 100_000, -1), olen))           # both at ends
        self.assertIsNone(_homology(None, (1, 0, -1), olen))


if __name__ == "__main__":
    unittest.main()
