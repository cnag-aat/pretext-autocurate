"""The adjacency solver (`solver.solve`, `run --engine adjacency`) on toy genomes (no fixtures)."""
import os
import sys
import tempfile
import unittest
from unittest import mock

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pretext_autocurate import MapState, Pretext                 # noqa: E402
from pretext_autocurate import solver                             # noqa: E402
from pretext_autocurate.solver import HOM_ACROSS, _homology, contigs, solve   # noqa: E402
from pretext_autocurate.homology import Block                      # noqa: E402
from pretext_autocurate.synth import toy_blocks, toy_map          # noqa: E402
from test_synthetic import CHROMS, DIPLOID                        # noqa: E402

# A and B in one H1 scaffold, the right way round: a YaHS join that is right
JOINED = [("H1_scaffold_1", [("A", 1, 0, 60, False), ("A", 1, 60, 100, False)]),
          ("H1_scaffold_2", [("B", 1, 0, 80, False)]),
          ("H1_scaffold_3", [("C", 1, 0, 50, False)]),
          ("H2_scaffold_1", [("A", 2, 0, 100, False)]),
          ("H2_scaffold_2", [("B", 2, 0, 80, False)]),
          ("H2_scaffold_3", [("C", 2, 0, 50, False)])]
# a 6-px piece of B's end (under two 8-px bins) inside A: needs full-resolution contact
SMALL = [("H1_scaffold_1", [("A", 1, 0, 60, False), ("B", 1, 74, 80, False), ("A", 1, 60, 100, False)]),
         ("H1_scaffold_2", [("B", 1, 0, 74, False)]),
         ("H1_scaffold_3", [("C", 1, 0, 50, False)]),
         ("H2_scaffold_1", [("A", 2, 0, 100, False)]),
         ("H2_scaffold_2", [("B", 2, 0, 80, False)]),
         ("H2_scaffold_3", [("C", 2, 0, 50, False)])]


class TestSolver(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.TemporaryDirectory()
        cls.maps = {}
        for name, sc in (("misjoin", DIPLOID), ("right", JOINED), ("small", SMALL)):
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

    def test_small_contig_at_full_resolution(self):
        pt, blocks = self.maps["small"]
        # binned only: the 6-px piece is merged into A and cannot move
        self.assertEqual(solve(pt, blocks=blocks)[2]["cuts"], 0)
        rows, st, info = solve(pt, blocks=blocks, fine_px=2)
        self.assertEqual(info["cuts"], 2)                  # out of A on both sides
        joins = [r for r in rows if r["action"] == "join"]
        self.assertTrue(any("H1_scaffold_2" in r["partner"] and "0.61-0.68" in r["scaffold"] for r in joins), joins)
        ms = MapState(pt).replay(st.edits)
        self.assertEqual(sorted(ms.pix.tolist()), list(range(pt.npix)))

    def test_mutual_best_lifts_specificity_veto(self):
        pt, _ = self.maps["right"]
        # without YaHS or alignment, the right A-A join (contact 0.7, 1.84x the next end) is not
        # specific enough at level 0 (2x); each other's best, it is made once 0.7 counts as near
        self.assertEqual(solve(pt, level=0, no_yahs=True)[2]["cuts"], 1)
        with mock.patch.dict(solver.LEVELS, {"mutual_min": (0.6, 0.6, 0.6)}):
            self.assertEqual(solve(pt, level=0, no_yahs=True)[2]["cuts"], 0)

    def test_pairs_laid_out_in_parallel(self):
        # H2's copy of A stored reversed, one contig per path: it must be turned so the homolog
        # off-diagonal runs parallel to the diagonal, not crossing it (an X)
        sc = [("H1_scaffold_1", [("A", 1, 0, 100, False)]), ("H1_scaffold_2", [("B", 1, 0, 80, False)]),
              ("H2_scaffold_1", [("A", 2, 0, 100, True)]), ("H2_scaffold_2", [("B", 2, 0, 80, False)])]
        p = os.path.join(self.dir.name, "reversed.pretext")
        cell = toy_map(p, sc, CHROMS)
        pt = Pretext(p)
        rows, st, info = solve(pt, blocks=toy_blocks(sc))
        self.assertEqual(info["pairs"], 2)
        pos = np.empty(pt.npix, dtype=np.int64)
        pos[MapState(pt).replay(st.edits).pix] = np.arange(pt.npix)
        for ch in ("A", "B"):
            x = {cell[q][2]: pos[q] for q in cell if cell[q] and cell[q][0] == ch and cell[q][1] == 1}
            y = {cell[q][2]: pos[q] for q in cell if cell[q] and cell[q][0] == ch and cell[q][1] == 2}
            k = sorted(set(x) & set(y))
            self.assertGreater(np.corrcoef([x[i] for i in k], [y[i] for i in k])[0, 1], 0.9, ch)

    def test_rephase_missing_copy(self):
        # hap 2's copy of B 50-80 sits in an H1 scaffold (hifiasm mis-phasing); an all-vs-all PAF
        # shows it on H1's own B, whose 50-80 has no H2 homolog: it moves to H2 and joins H2's B
        sc = [("H1_scaffold_1", [("A", 1, 0, 100, False)]), ("H1_scaffold_2", [("B", 1, 0, 80, False)]),
              ("H1_scaffold_3", [("B", 2, 50, 80, False)]),
              ("H2_scaffold_1", [("A", 2, 0, 100, False)]), ("H2_scaffold_2", [("B", 2, 0, 50, False)])]
        p = os.path.join(self.dir.name, "misphased.pretext")
        toy_map(p, sc, CHROMS)
        pt = Pretext(p)
        same = Block("H1_scaffold_3", 0, 300_000, "H1_scaffold_2", 500_000, 800_000, True, 300_000, 60)
        blocks = toy_blocks(sc) + [same]
        hap_of = {n: n[:2] for n in pt.names}
        moved, _ = solver.rephase(pt, blocks, hap_of)
        self.assertEqual(moved, {"H1_scaffold_3": "H2"})
        # a third copy (H2 already has B 50-80) stays where it is
        sc3 = sc + [("H2_scaffold_3", [("B", 2, 50, 80, False)])]
        p3 = os.path.join(self.dir.name, "thirdcopy.pretext")
        toy_map(p3, sc3, CHROMS)
        pt3 = Pretext(p3)
        self.assertEqual(solver.rephase(pt3, toy_blocks(sc3) + [same], {n: n[:2] for n in pt3.names})[0], {})
        rows, st, info = solve(pt, blocks=blocks)
        self.assertEqual(info["rephased"], 1)
        joins = [r for r in rows if r["action"] == "join"]
        self.assertTrue(any("H1_scaffold_3" in r["scaffold"] + r["partner"] and "H2_scaffold_2" in r["scaffold"] + r["partner"]
                            for r in joins), joins)

    def test_homology(self):
        # landing: (scaffold, position, direction out, contig, bp of that contig beyond it)
        olen = np.array([10_000_000, 10_000_000])
        L = lambda s, p, d, c=0, room=0.0: (s, p, d, c, room)
        self.assertEqual(_homology(L(0, 5_000_000, 1), L(0, 5_200_000, -1)), 1.0)            # facing, one contig
        self.assertEqual(_homology(L(0, 5_000_000, 1), L(0, 5_200_000, 1)), -1.0)            # same way
        self.assertEqual(_homology(L(0, 1_000_000, 1), L(0, 9_000_000, -1)), -1.0)           # far apart
        self.assertEqual(_homology(L(0, 5_000_000, 1), L(0, 5_200_000, -1, 1)), HOM_ACROSS)  # across a YaHS join
        self.assertEqual(_homology(L(0, 2_000_000, 1, 0, 8e6), L(1, 0, -1), olen), -1.0)     # its contig goes on
        self.assertEqual(_homology(L(0, 2_000_000, 1), L(1, 0, -1), olen), -HOM_ACROSS)      # its scaffold goes on
        self.assertIsNone(_homology(L(0, 9_900_000, 1), L(1, 100_000, -1), olen))            # both at ends
        self.assertIsNone(_homology(None, L(1, 0, -1), olen))


if __name__ == "__main__":
    unittest.main()
