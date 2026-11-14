"""Evidence tracks: bedgraph/BED to pixels, telomere strands at a break site."""
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext  # noqa: E402
from pretext_autocurate import tracks  # noqa: E402
from pretext_autocurate.allhaps import pair_homologs, telomere_sites  # noqa: E402

M0 = os.path.join(ROOT, "seed", "fixtures", "bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid_mq0.extensions.HR.pretext")


@unittest.skipUnless(os.path.exists(M0), "missing fixture")
class TestTracks(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(M0)
        cls.hi = Pretext(M0.replace("diploid_mq0", "diploid_mq10"), load_graphs=False)

    def test_to_pixels(self):
        pt, bpp = self.pt, self.pt.bp_per_pix
        i = pt.names.index("H1_scaffold_10")
        a = int(pt.scaffold_start[i])
        with tempfile.TemporaryDirectory() as d:
            bg, bed = os.path.join(d, "t.bg"), os.path.join(d, "r.bed")
            with open(bg, "w") as f:
                f.write(f"H1_scaffold_10\t0\t10000\t800\nH1_scaffold_10\t{int(5 * bpp)}\t{int(5 * bpp) + 10000}\t300\n"
                        "no_such_scaffold\t0\t100\t5\n")
            with open(bed, "w") as f:
                f.write(f"H1_scaffold_10\t{int(2 * bpp)}\t{int(3 * bpp)}\nH1_scaffold_10\t{int(4 * bpp)}\t{int(4.5 * bpp)}\n")
            v, used = tracks.to_pixels(pt, bg)
            self.assertEqual(used, 2)
            self.assertEqual(v[a], 800)
            self.assertEqual(v[a + 5], 300)
            self.assertEqual(v[a + 1], 0)
            r, _ = tracks.to_pixels(pt, bed, value=False)
            self.assertAlmostEqual(r[a + 2], 1.0, places=3)
            self.assertAlmostEqual(r[a + 4], 0.5, places=2)
            self.assertAlmostEqual(r[a + 3], 0.0, places=3)

    def test_busco_tracks(self):
        """Duplicate pairs get one label per pair of scaffolds (both places coloured alike);
        a pair sharing one gene is left out; genes absent from one haplotype are listed."""
        hits = {f"g{i}": [("H1_a", 1000 * i, 1000 * i + 500, "+", 1.0), ("H1_b", 1000 * i, 1000 * i + 500, "+", 1.0),
                          ("H2_a", 1000 * i, 1000 * i + 500, "+", 1.0)] for i in range(3)}
        hits["solo"] = [("H1_c", 10, 20, "+", 1.0), ("H1_d", 10, 20, "+", 1.0), ("H2_a", 5, 9, "+", 1.0)]
        hits["only_h1"] = [("H1_a", 9000, 9500, "+", 1.0)]
        hits["tandem"] = [("H2_a", 100000, 100500, "+", 1.0), ("H2_a", 120000, 120500, "+", 1.0),
                          ("H1_a", 100000, 100500, "+", 1.0)]
        with tempfile.TemporaryDirectory() as d:
            out = tracks.busco_tracks(hits, lambda n: n.split("_")[0], d)
            self.assertEqual(set(out), {"busco_dup_pairs", "busco_missing_hap"})
            lines = [l.split() for l in open(out["busco_dup_pairs"])]
            self.assertEqual({l[3] for l in lines}, {"d1"})                       # one pair: H1_a / H1_b
            self.assertEqual({l[0] for l in lines}, {"H1_a", "H1_b"})
            self.assertEqual(len(lines), 6)
            miss = [l.split() for l in open(out["busco_missing_hap"])]
            self.assertEqual(miss, [["H1_a", "9000", "9500", "1"]])

    def test_strands_decide_a_telomere_site(self):
        """H1_scaffold_2 at 166.3 Mb: the end piece is inverted (curator). With strand tracks:
        a telomere on the end piece facing the gap keeps that call; telomeres on both sides
        facing the gap make it a cut (two chromosome ends)."""
        pt = self.pt
        pairs, _ = pair_homologs(pt)
        i = pt.names.index("H1_scaffold_2")
        q = int(pt.scaffold_start[i] + round(166.26e6 / pt.bp_per_pix)) - 1
        saved = {k: pt.graphs.get(k) for k in ("telomere", "telomere_5p", "telomere_3p")}
        try:
            for case, (l3, r5) in {"tip": (False, True), "both": (True, True)}.items():
                t5, t3 = np.zeros(pt.npix), np.zeros(pt.npix)
                if l3:
                    t3[q - 1] = 900
                if r5:
                    t5[q + 2] = 900
                pt.graphs["telomere_5p"], pt.graphs["telomere_3p"] = t5, t3
                pt.graphs["telomere"] = saved["telomere"]
                sites = telomere_sites(pt, pairs, self.hi)
                act = {x["q"]: x for x in sites}
                site = min(act, key=lambda k: abs(k - q))
                self.assertLessEqual(abs(site - q), 3)
                want = "invert_tip" if case == "tip" else "cut"
                self.assertEqual(act[site]["action"], want, act[site]["reason"])
                self.assertIn("telomere strand", act[site]["reason"])
        finally:
            for k, v in saved.items():
                if v is None:
                    pt.graphs.pop(k, None)
                else:
                    pt.graphs[k] = v


if __name__ == "__main__":
    unittest.main()
