"""End-to-end single-mode pass on xgTylRafi, compared loosely with the curation."""
import collections
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, SaveState, MapState, header_hash, cli
from pretext_autocurate.curate import curate_single

TYL = os.path.join(ROOT, "seed", "fixtures", "xgTylRafi_asm4_run1.hfsm.asm.hic.hap2.yhs_mq10_mq0.extensions.HR.pretext")


@unittest.skipUnless(os.path.exists(TYL), "missing fixture")
class TestRun(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(TYL)
        cls.dec, cls.st = curate_single(cls.pt)
        cur = SaveState.load(TYL + ".savestate_4")
        ms = MapState(cls.pt).replay(cur.edits)
        paint, fl = dict(cur.scaffs), dict(cur.meta_flags)
        cls.chrom = np.zeros(cls.pt.npix, dtype=int)
        cls.truth = {}
        for i, c in enumerate(ms.contigs):
            cls.chrom[ms.pix[c.start:c.start + c.length]] = paint.get(i, 0)
            if c.length == cls.pt.scaffold_npix[c.scaffold]:
                cls.truth[cls.pt.names[c.scaffold]] = (cur.tag_names(fl.get(i, 0)), paint.get(i, 0))

    def test_state_is_consistent(self):
        ms = MapState(self.pt).replay(self.st.edits)
        self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(self.pt.npix)))
        n = len(ms.contigs)
        self.assertTrue(all(i < n for i, _ in self.st.scaffs + self.st.meta_flags))
        # new contigs come only from applied cuts and from insertions splitting their target
        cuts = sum(d.action == "cut" for d in self.dec)
        ins = sum(d.action == "insert" for d in self.dec)
        joins = sum(d.action == "join" for d in self.dec)
        base = int((self.pt.scaffold_npix > 0).sum()) + cuts
        # a join of two pieces that were cut apart and put back in the same order merges them again
        self.assertTrue(base - joins <= n <= base + ins, (n, base, ins, joins))
        self.assertLessEqual(sum(e.is_break for e in self.st.edits), cuts)   # breaks only where sides stay adjacent
        self.assertLessEqual(len(self.st.waypoints), 255)

    def test_chromosome_assignment_matches_curation(self):
        ok = bad = 0
        for d in self.dec:
            ch = self.truth.get(d.scaffold, ([], 0))[1]
            if d.action.startswith("attach") and ch and "piece " not in d.evidence:
                px = self.pt.scaffold_start[self.pt.names.index(d.target)] + int((d.target_bp - 1) / self.pt.bp_per_pix)
                if self.chrom[px] == ch:
                    ok += 1
                else:
                    bad += 1
        self.assertGreaterEqual(ok, 170)
        self.assertLessEqual(bad, 4)

    def test_joins_and_orientation(self):
        """scaffold_2 88.7-109.7 Mb (reversed) + scaffold_15 19.8-24.4 Mb (reversed) +
        scaffold_2 111.1 Mb onward, as the curator built chromosome 3; no join at telomeres."""
        joins = [d for d in self.dec if d.action == "join"]
        self.assertGreaterEqual(len(joins), 8)
        ms = MapState(self.pt).replay(self.st.edits)
        pos = np.empty(self.pt.npix, dtype=int)
        pos[ms.pix] = np.arange(self.pt.npix)
        bpp = self.pt.bp_per_pix

        def px(name, mb):
            return int(self.pt.scaffold_start[self.pt.names.index(name)] + mb * 1e6 / bpp)
        a0, a1 = pos[px("scaffold_2", 90)], pos[px("scaffold_2", 108)]
        b0, b1 = pos[px("scaffold_15", 20.5)], pos[px("scaffold_15", 23.5)]
        c0, c1 = pos[px("scaffold_2", 113)], pos[px("scaffold_2", 170)]
        self.assertTrue(a1 < a0 < b1 < b0 < c0 < c1, (a0, a1, b0, b1, c0, c1))
        self.assertLess(c0 - a0, (113 - 88.7 + 24.4 - 19.5) * 1e6 / bpp + 100)   # contiguous, nothing else between
        paint = dict(self.st.scaffs)
        ids = {paint.get(int(ms.cid[p])) for p in (a0, b0, c0)}
        self.assertEqual(len(ids), 1)
        self.assertNotIn(None, ids)
        flagged = [d.evidence for d in self.dec if d.action == "flag_join"]
        self.assertTrue(any("telomere" in e and "scaffold_11" in e and "scaffold_16" in e for e in flagged))

    def test_inversion_inside_a_piece(self):
        """scaffold_2 53.9-88.7 Mb: the curator's order is 82.0->88.7 then 82.0->53.9. The
        tool gets there by a tentative cut at the weak junction at 82.0 Mb and re-joining."""
        ms = MapState(self.pt).replay(self.st.edits)
        pos = np.empty(self.pt.npix, dtype=int)
        pos[ms.pix] = np.arange(self.pt.npix)
        s0 = self.pt.scaffold_start[self.pt.names.index("scaffold_2")]
        px = lambda mb: int(s0 + mb * 1e6 / self.pt.bp_per_pix)
        p83, p88, p81, p55 = (pos[px(m)] for m in (83, 88, 81, 55))
        if p83 > p88:                                     # the chain may be laid out either way round
            p83, p88, p81, p55 = -p83, -p88, -p81, -p55
        self.assertTrue(p83 < p88 < p81 < p55, (p83, p88, p81, p55))
        self.assertLess(abs(pos[px(88.6)] - pos[px(81.9)]), 20)      # 88.7 end now touches 82.0 end
        weak = [d for d in self.dec if "weak junction" in d.evidence]
        self.assertIn(sum(d.action == "cut" for d in weak), (1, 2, 3))
        self.assertGreaterEqual(sum(d.action == "flag_cut" for d in weak), 8)

    def test_borderline_haplotigs_are_highlighted_not_waypointed(self):
        vert = 1 << self.st.tag_id("vertPaint")
        n_low = sum("POSSIBLE HAPLOTIG" in d.evidence for d in self.dec)
        self.assertGreater(n_low, 20)
        self.assertEqual(sum(bool(m & vert) for _, m in self.st.meta_flags), n_low)
        # cuts, weak junctions, flagged joins, proposed insertions and overlaps only (2026-11-19: two
        # proposed overlaps at cuts in scaffold_2 and 3, without --fasta here, made it 91)
        self.assertLess(len(self.st.waypoints), 95)

    def test_contaminant_left_alone(self):
        d = {x.scaffold: x for x in self.dec}["scaffold_18"]
        self.assertEqual(d.action, "none")
        self.assertIn("no cis-level contact", d.evidence)

    def test_chromosome_scale_scaffolds_stay_put(self):
        self.assertTrue(all(d.length_px <= 150 for d in self.dec if d.action.startswith(("attach", "insert"))))

    def test_misjoin_cuts_match_curation(self):
        """Every applied cut lies within 12 px of a cut the curator made."""
        cur = SaveState.load(TYL + ".savestate_4")
        ms = MapState(self.pt).replay(cur.edits)
        where = np.empty(self.pt.npix, dtype=int)
        where[ms.pix] = np.arange(self.pt.npix)
        q = np.arange(self.pt.npix - 1)
        sc = self.pt.orig_scaffold
        curated = q[(sc[q] == sc[q + 1]) & (np.abs(where[q + 1] - where[q]) != 1)]
        applied = [d for d in self.dec if d.action == "cut"]
        self.assertGreaterEqual(len(applied), 12)
        for d in applied:
            px = self.pt.scaffold_start[self.pt.names.index(d.scaffold)] + d.target_bp / self.pt.bp_per_pix
            self.assertLessEqual(np.abs(curated - px).min(), 12, d)

    def test_review_fixes_2026_11_19(self):
        """Curator review of xgTylRafi: scaffold_13 stays one chromosome after its haplotig contig
        is cut out; scaffold_30 loses its haplotig end contig at the gap, and its first part goes
        into scaffold_3 at 22.97 Mb after the reversed 46.65-47.17 Mb piece, as curated."""
        ms = MapState(self.pt).replay(self.st.edits)
        paint = dict(self.st.scaffs)
        chrom = np.zeros(self.pt.npix, dtype=int)
        pos = np.empty(self.pt.npix, dtype=int)
        pos[ms.pix] = np.arange(self.pt.npix)
        for i, c in enumerate(ms.contigs):
            chrom[ms.pix[c.start:c.start + c.length]] = paint.get(i, 0)
        px = lambda s, mb: int(self.pt.scaffold_start[self.pt.names.index(s)] + mb * 1e6 / self.pt.bp_per_pix)
        self.assertEqual(chrom[px("scaffold_13", 5)], chrom[px("scaffold_13", 40)])
        self.assertGreater(chrom[px("scaffold_13", 5)], 0)
        cut = [d for d in self.dec if d.scaffold == "scaffold_30" and d.action == "cut"]
        self.assertEqual(len(cut), 1)
        self.assertTrue(any(d.scaffold == "scaffold_30" and d.action in ("attach_haplotig", "tag_haplotig")
                            for d in self.dec))
        a, b, c = pos[px("scaffold_3", 46.8)], pos[px("scaffold_30", 0.2)], pos[px("scaffold_3", 23.5)]
        self.assertTrue(a < b < c and c - a < 60)

    def test_only_clear_insertions_by_default(self):
        ins = [d for d in self.dec if d.action == "insert"]
        self.assertTrue(6 <= len(ins) <= 10)
        got = {(d.scaffold, d.target, d.strand) for d in ins}
        # insertions the curator also made (reviewed savestate), all at reduced-contact gaps
        for want in (("scaffold_5", "scaffold_16", "+"), ("scaffold_2", "scaffold_5", "+"),
                     ("scaffold_6", "scaffold_2", "+"), ("scaffold_21", "scaffold_10", "-"),
                     ("scaffold_6", "scaffold_6", "+")):
            self.assertIn(want, got)
        self.assertTrue(all(d.length_px >= 15 and d.strand in "+-" for d in ins))
        # the 2 Mb piece of scaffold_3 goes into scaffold_1 near 166.66 Mb, reversed, as curated
        d = [x for x in ins if x.scaffold == "scaffold_3" and x.target == "scaffold_1"]
        self.assertEqual(len(d), 1)
        self.assertEqual(d[0].strand, "-")
        self.assertLess(abs(d[0].target_bp - 166663868), 300000)

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            cli.main(["run", "--map", TYL, "--out-dir", d])
            st = SaveState.load(os.path.join(d, os.path.basename(TYL) + ".savestate_autocurate"))
            self.assertEqual(st.header_hash, header_hash(TYL))
            self.assertEqual(st.payload(), self.st.payload())
            for f in ("autocurate.decisions.tsv", "autocurate.agp"):
                self.assertGreater(os.path.getsize(os.path.join(d, f)), 1000)


if __name__ == "__main__":
    unittest.main()
