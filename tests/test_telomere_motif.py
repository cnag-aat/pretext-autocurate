"""Telomere motif detected from the scaffold ends (tracks.detect_motif)."""
import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate.tracks import _primitive, _unit_name, detect_motif     # noqa: E402

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LAR_FA = os.path.join(ROOT, "seed", "fixtures", "icLarUrsu.pri.yahs.fa")


def genome(path, motif, n=8, size=100_000, seed=1):
    rng = random.Random(seed)
    rc = motif[::-1].translate(str.maketrans("ACGT", "TGCA"))
    with open(path, "w") as f:
        for i in range(n):
            body = "".join(rng.choice("ACGT") for _ in range(size))
            s = (rc * 300 + body + motif * 300) if motif else body
            f.write(f">scaffold_{i + 1}\n" + "".join(s[k:k + 60] + "\n" for k in range(0, len(s), 60)))


class Motif(unittest.TestCase):
    def test_names_and_units(self):
        self.assertEqual(_unit_name("AGGGTT"), "TTAGGG")
        self.assertEqual(_unit_name("GGGTTAGTTT"), "TTTGGGTTAG")
        self.assertFalse(_primitive("ATATAT"))         # AT repeated
        self.assertFalse(_primitive("AAAAAC"))         # near-homopolymer
        self.assertTrue(_primitive("TTAGG"))

    def test_detected_at_both_ends(self):
        with tempfile.TemporaryDirectory() as d:
            for motif in ("TTAGGG", "TTTGGGTTAG"):
                p = os.path.join(d, motif + ".fa")
                genome(p, motif)
                self.assertEqual(detect_motif(p)[0], motif)

    def test_none_without_telomeres(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "none.fa")
            genome(p, "")
            self.assertIsNone(detect_motif(p)[0])

    @unittest.skipUnless(os.path.exists(LAR_FA), "missing fixture")
    def test_larinus(self):
        self.assertEqual(detect_motif(LAR_FA)[0], "TTTGGGTTAG")


if __name__ == "__main__":
    unittest.main()
