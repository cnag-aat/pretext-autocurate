"""score.compare: a savestate against a reference curation of the same map."""
import math
import os
import sys
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, SaveState          # noqa: E402
from pretext_autocurate.score import compare               # noqa: E402

FIX = os.path.join(ROOT, "seed", "fixtures")
TYL = os.path.join(FIX, "xgTylRafi_asm4_run1.hfsm.asm.hic.hap2.yhs_mq10_mq0.extensions.HR.pretext")
TRO = os.path.join(FIX, "qqTroRaco.hfsm.asm.bp.yhs_mq10.diploid_mq0.extensions.HR.pretext")


@unittest.skipUnless(os.path.exists(TYL), "missing fixture")
class TestScore(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(TYL, load_graphs=False)
        cls.ref = SaveState.load(TYL + ".savestate_4")

    def test_reference_against_itself(self):
        m = compare(self.pt, self.ref, self.ref)
        for k in ("junction_precision", "junction_recall", "break_precision", "break_recall",
                  "grouped_90", "grouped_bp"):
            self.assertEqual(m[k], 1.0, k)
        self.assertGreater(m["junctions_ref"], 10)
        self.assertEqual(m["tag_Haplotig_both_mb"], m["tag_Haplotig_ref_mb"])

    def test_unedited_map(self):
        """The map as assembled: no junctions, no breaks, nothing painted."""
        blank = SaveState.blank()
        m = compare(self.pt, blank, self.ref)
        self.assertEqual((m["junctions_tool"], m["junctions_reproduced"], m["breaks_tool"]), (0, 0, 0))
        self.assertTrue(math.isnan(m["junction_precision"]))
        self.assertEqual((m["junction_recall"], m["break_recall"], m["grouped_90"]), (0.0, 0.0, 0.0))


@unittest.skipUnless(os.path.exists(TRO), "missing fixture")
class TestScoreTwoHaplotypes(unittest.TestCase):
    def test_tags(self):
        """Tag agreement by bp: the curator's X and Y against themselves."""
        pt = Pretext(TRO, load_graphs=False)
        ref = SaveState.load(TRO + ".savestate_final")
        m = compare(pt, ref, ref)
        self.assertAlmostEqual(m["tag_X_both_mb"], m["tag_X_ref_mb"])
        self.assertGreater(m["tag_X_ref_mb"], 100)
        self.assertEqual(m["ref_chromosomes"], 128)


if __name__ == "__main__":
    unittest.main()
