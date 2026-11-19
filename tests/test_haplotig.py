"""Haplotig detector regression against the curated xgTylRafi map.

The curation is a guide, not exact truth: the bounds below are deliberately loose and
only guard against the detector getting worse.
"""
import collections
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, SaveState, MapState, header_hash, cli
from pretext_autocurate.haplotig import detect_haplotigs
from pretext_autocurate.signal import Signal

TYL = os.path.join(ROOT, "seed", "fixtures", "xgTylRafi_asm4_run1.hfsm.asm.hic.hap2.yhs_mq10_mq0.extensions.HR.pretext")


@unittest.skipUnless(os.path.exists(TYL), "missing fixture")
class TestHaplotig(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(TYL)
        cls.st = SaveState.load(TYL + ".savestate_4")
        cls.calls = detect_haplotigs(cls.pt)

    def test_decay_is_monotone_overall(self):
        g = Signal(self.pt).decay
        self.assertGreater(g[0], g[10])
        self.assertGreater(g[10], g[100])

    def test_agreement_with_curation(self):
        ms = MapState(self.pt).replay(self.st.edits)
        fl = dict(self.st.meta_flags)
        truth = {}
        for i, c in enumerate(ms.contigs):
            if c.length == self.pt.scaffold_npix[c.scaffold]:
                truth[self.pt.names[c.scaffold]] = self.st.tag_names(fl.get(i, 0))
        tab = collections.Counter()
        for c in self.calls:
            tab[(c.call, "Haplotig" in truth.get(c.scaffold, []))] += 1
        called = ("haplotig_high", "haplotig_medium")
        tp = sum(tab[(k, True)] for k in called)
        fp = sum(tab[(k, False)] for k in called)
        self.assertGreaterEqual(tp, 70)                    # of 117 whole-scaffold haplotigs
        self.assertLessEqual(fp, 10)
        self.assertLessEqual(tab[("haplotig_high", False)], 5)
        self.assertLessEqual(tab[("none", True)], 12)

    def test_contaminant_not_called(self):
        by = {c.scaffold: c for c in self.calls}
        self.assertEqual(by["scaffold_18"].call if "scaffold_18" in by else "none", "none")

    def test_against_reviewed_tags(self):
        """Curator's review of the tool's tags (2026-11-13): confident calls are nearly all
        right, borderline calls catch most of the rest, zero-coverage repeats are not called."""
        rev = SaveState.load(TYL + ".savestate_4_curator_reviewed_tags")
        ms = MapState(self.pt).replay(rev.edits)
        fl = dict(rev.meta_flags)
        hap = {self.pt.names[c.scaffold] for i, c in enumerate(ms.contigs)
               if c.length == self.pt.scaffold_npix[c.scaffold] and "Haplotig" in rev.tag_names(fl.get(i, 0))}
        by = {c.scaffold: c.call for c in self.calls}
        sure = [n for n, k in by.items() if k in ("haplotig_high", "haplotig_medium")]
        low = [n for n, k in by.items() if k == "haplotig_low"]
        self.assertGreaterEqual(sum(n in hap for n in sure), 80)
        self.assertLessEqual(sum(n not in hap for n in sure), 3)
        missed = [n for n in hap if n in by and n not in sure]
        self.assertGreaterEqual(sum(n in low for n in missed), 0.8 * len(missed))
        self.assertLessEqual(len(low), 90)
        for n in ("scaffold_196", "scaffold_237"):        # repeats with near-zero coverage
            self.assertNotIn(by[n], ("haplotig_high", "haplotig_medium"))

    def test_contig_haplotigs(self):
        """Strict contig-level calls: the two clearest curated cases, and at most two more."""
        from pretext_autocurate.haplotig import detect_contig_haplotigs
        calls = detect_contig_haplotigs(self.pt)
        bpp = self.pt.bp_per_pix
        where = {(c.scaffold, round(self.pt.rel_coord[c.start_px] * bpp / 1e6, 1)) for c in calls}
        self.assertIn(("scaffold_3", 119.8), where)
        self.assertIn(("scaffold_13", 23.6), where)
        self.assertLessEqual(len(calls), 4)
        gap = self.pt.graphs["gap"] > 0
        for c in calls:                       # delimited by gaps or scaffold ends
            self.assertTrue(gap[c.end_px] or c.end_px + 1 == self.pt.npix or
                            self.pt.orig_scaffold[c.end_px + 1] != self.pt.orig_scaffold[c.end_px])

    def test_small_scaffold_haplotig_contig(self):
        """scaffold_30: the curator cut its last contig off at the gap as a haplotig; scaffold_28
        (a 6 px low-coverage end contig) was tagged whole, so it must not be split."""
        from pretext_autocurate.curate import small_scaffold_haplotigs
        cuts, pieces = small_scaffold_haplotigs(self.pt, [])
        by = {c.scaffold: c for c in cuts}
        self.assertIn("scaffold_30", by)
        self.assertNotIn("scaffold_28", by)
        self.assertLessEqual(len(cuts), 3)
        s = self.pt.names.index("scaffold_30")
        a, L = int(self.pt.scaffold_start[s]), int(self.pt.scaffold_npix[s])
        self.assertTrue(self.pt.graphs["gap"][by["scaffold_30"].after_px] > 0)    # at the gap
        self.assertIn((by["scaffold_30"].after_px + 1, a + L), pieces)          # the end contig

    def test_cli_writes_loadable_state(self):
        with tempfile.TemporaryDirectory() as d:
            cli.main(["haplotigs", TYL, "--out-dir", d])
            st = SaveState.load(os.path.join(d, os.path.basename(TYL) + ".savestate_autocurate"))
            self.assertEqual(st.header_hash, header_hash(TYL))
            self.assertEqual(st.edits, [])
            n = len(MapState(self.pt).contigs)
            self.assertTrue(all(i < n and m == 1 for i, m in st.meta_flags))
            self.assertEqual(len(st.meta_flags), sum(c.call in ("haplotig_high", "haplotig_medium") for c in self.calls))


if __name__ == "__main__":
    unittest.main()
