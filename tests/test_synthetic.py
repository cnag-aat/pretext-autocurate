"""End-to-end tests on small synthetic maps (`pretext_autocurate.synth`): they need no
fixtures, so they run everywhere, including GitHub Actions on every push.

The toy genome: chromosomes A (100 px), B (80 px) and C (50 px) in two haplotypes. H1 has
A and B misjoined in one scaffold (a gap between them); every other chromosome copy is a
scaffold of its own."""
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, SaveState, MapState, cli          # noqa: E402
from pretext_autocurate.allhaps import curate_allhaps                     # noqa: E402
from pretext_autocurate.curate import curate_single                       # noqa: E402
from pretext_autocurate.pretextio import decode_bc4                       # noqa: E402
from pretext_autocurate.score import compare                              # noqa: E402
from pretext_autocurate.synth import encode_bc4, toy_map, write_pretext   # noqa: E402

CHROMS = {"A": 100, "B": 80, "C": 50}
DIPLOID = [("H1_scaffold_1", [("A", 1, 0, 100, False), ("B", 1, 0, 80, False)]),    # misjoin
           ("H1_scaffold_2", [("C", 1, 0, 50, False)]),
           ("H2_scaffold_1", [("A", 2, 0, 100, False)]),
           ("H2_scaffold_2", [("B", 2, 0, 80, False)]),
           ("H2_scaffold_3", [("C", 2, 0, 50, False)])]
SINGLE = [("scaffold_1", [("A", 1, 0, 100, False), ("B", 1, 0, 80, False)]),       # misjoin
          ("scaffold_2", [("C", 1, 0, 50, False)])]


class TestFormat(unittest.TestCase):
    def test_bc4_roundtrip(self):
        rng = np.random.default_rng(1)
        smooth = np.clip(np.add.outer(np.arange(64), np.arange(64)) * 2 + rng.integers(0, 3, (64, 64)), 0, 255)
        img = smooth.astype(np.uint8)
        dec = decode_bc4(encode_bc4(img), 64)
        self.assertLessEqual(int(np.abs(dec.astype(int) - img.astype(int)).max()), 2)
        flat = np.full((8, 8), 77, dtype=np.uint8)
        self.assertTrue((decode_bc4(encode_bc4(flat), 8) == 77).all())

    def test_write_read(self):
        """Header, pixels per scaffold, every mip level, tracks and the header hash."""
        with tempfile.TemporaryDirectory() as d:
            npix = 512
            i = np.arange(npix)
            m = np.clip(200.0 / (1 + np.abs(i[:, None] - i[None, :])), 0, 255)
            p = os.path.join(d, "x.pretext")
            write_pretext(p, ["a", "b"], [3_000_000, 1_000_000], m, graphs={"gap": np.zeros(npix), "coverage": i})
            pt = Pretext(p)
            self.assertEqual((pt.npix, pt.names, pt.total_len), (512, ["a", "b"], 4_000_000))
            self.assertEqual(pt.scaffold_npix.tolist(), [384, 128])
            self.assertTrue(np.array_equal(pt.graphs["coverage"], i))
            got = pt.matrix(0).astype(int)
            self.assertLessEqual(int(np.abs(got - np.rint(m)).max()), 30)        # BC4 is lossy
            self.assertEqual(got[10, 10], 200)
            for lv in range(1, pt.n_mips):
                self.assertEqual(pt.matrix(lv).shape, (512 >> lv, 512 >> lv))
            self.assertEqual(len(pt.hash), 8)
            pt.close()


class TestTwoHaplotypes(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.d = tempfile.mkdtemp()
        cls.map = os.path.join(cls.d, "toy_diploid.pretext")
        toy_map(cls.map, DIPLOID, CHROMS)
        cls.pt = Pretext(cls.map)
        cls.rows, cls.st, cls.info = curate_allhaps(cls.pt, phasing="hic")

    def test_misjoin_cut(self):
        cuts = [r for r in self.rows if r["action"] == "cut"]
        self.assertEqual([r["scaffold"] for r in cuts], ["H1_scaffold_1"])
        self.assertAlmostEqual(float(cuts[0]["length_mb"]), 1.0, delta=0.05)     # A is 100 px of 10 kb
        self.assertIn("different homologs", cuts[0]["evidence"])

    def test_pairs(self):
        pairs = {frozenset((r["scaffold"], r["partner"])) for r in self.rows if r["action"] == "pair"}
        self.assertEqual(pairs, {frozenset(("H1_scaffold_1:1", "H2_scaffold_1")),
                                 frozenset(("H1_scaffold_1:2", "H2_scaffold_2")),
                                 frozenset(("H1_scaffold_2", "H2_scaffold_3"))})

    def test_savestate(self):
        """The edits are a permutation of the map; every piece is painted and tagged HAP1/HAP2;
        each chromosome's two copies sit next to each other."""
        ms = MapState(self.pt).replay(self.st.edits)
        self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(self.pt.npix)))
        flags = dict(self.st.meta_flags)
        tags = [set(self.st.tag_names(flags.get(i, 0))) for i in range(len(ms.contigs))]
        self.assertTrue(all(t & {"HAP1", "HAP2"} for t in tags))
        order = [self.pt.names[c.scaffold] for c in ms.contigs]
        self.assertEqual(len(ms.contigs), 6)                      # five scaffolds, one cut
        paint = dict(self.st.scaffs)
        self.assertEqual(len({paint.get(i) for i in range(len(ms.contigs))}), 6)
        i = order.index("H2_scaffold_3")
        self.assertIn("H1_scaffold_2", order[max(0, i - 1):i + 2])

    def test_score_against_itself_and_truth(self):
        """score.compare runs on a real savestate; the tool's layout against itself is
        perfect, and against the uncurated map it has edits to show."""
        m = compare(self.pt, self.st, self.st)
        self.assertEqual((m["break_precision"], m["break_recall"]), (1.0, 1.0))
        blank = SaveState.blank()
        self.assertEqual(compare(self.pt, blank, self.st)["break_recall"], 0.0)

    def test_cli_run(self):
        """The command line end to end: savestate (bound to the map), decisions, review list,
        AGP as PretextView exports it, run info."""
        out = os.path.join(self.d, "out")
        cli.main(["run", "--map", self.map, "--out-dir", out])
        name = os.path.basename(self.map) + ".savestate_autocurate"
        st = SaveState.load(os.path.join(out, name))
        self.assertEqual(st.header_hash, self.pt.hash)
        for f in ("autocurate.decisions.tsv", "autocurate.review.tsv", "autocurate.run.txt"):
            self.assertTrue(os.path.getsize(os.path.join(out, f)) > 0, f)
        with open(os.path.join(out, "autocurate.agp")) as f:
            agp = f.read().splitlines()
        self.assertEqual(agp[0], "##agp-version\t2.1")
        self.assertTrue(agp[2].startswith("# HiC MAP RESOLUTION: "))
        w = [l.split("\t") for l in agp if not l.startswith("#") and l.split("\t")[4] == "W"]
        self.assertEqual(len(w), 6)
        self.assertEqual({x[5] for x in w}, {n for n, _ in DIPLOID})
        # a sweep: one directory per level, and the summary table
        cli.main(["run", "--map", self.map, "--out-dir", out, "--sweep", "0,1"])
        self.assertTrue(os.path.exists(os.path.join(out, "level_0.00", name + "_level0.00")))
        self.assertTrue(os.path.exists(os.path.join(out, "autocurate.sweep.tsv")))


class TestSingleHaplotype(unittest.TestCase):
    def test_misjoin_cut(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "toy_single.pretext")
            toy_map(p, SINGLE, CHROMS)
            pt = Pretext(p)
            dec, st = curate_single(pt)
            cuts = [x for x in dec if x.action == "cut"]
            self.assertEqual([x.scaffold for x in cuts], ["scaffold_1"])
            ms = MapState(pt).replay(st.edits)
            self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(pt.npix)))
            pt.close()


if __name__ == "__main__":
    unittest.main()
