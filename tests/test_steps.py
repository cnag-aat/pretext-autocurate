"""Tests per step of the allhaps pass (`allhaps.driver`) on synthetic toy genomes: no
fixtures, so they run everywhere (GitHub Actions too).

Each toy genome has chromosomes A (100 px), B (80 px) and C (60 px) at 20 kb per pixel in
two haplotypes (`synth.toy_map`), with alignment blocks between the haplotypes taken from
the same layout (`synth.toy_blocks`, what minimap2 would give). `curate_allhaps(...,
stop_after=step)` returns the run's state after that step, so a step's own result can be
checked; the final layout is read back through the toy's truth (chromosome, haplotype,
position per pixel)."""
import os
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, MapState                   # noqa: E402
from pretext_autocurate.allhaps import curate_allhaps              # noqa: E402
from pretext_autocurate.synth import toy_blocks, toy_map           # noqa: E402

BPP = 20000
CHROMS = {"A": 100, "B": 80, "C": 60, "R": 6}
REST = [("H1_scaffold_2", [("B", 1, 0, 80, False)]), ("H1_scaffold_3", [("C", 1, 0, 60, False)]),
        ("H2_scaffold_2", [("B", 2, 0, 80, False)]), ("H2_scaffold_3", [("C", 2, 0, 60, False)])]
H2_A = [("H2_scaffold_1", [("A", 2, 0, 100, False)])]


class Toy:
    """A toy genome on disk, its map, truth and alignment blocks."""

    def __init__(self, scaffolds, blocks=True, drop_blocks=(), npix=512):
        self.d = tempfile.TemporaryDirectory()
        path = os.path.join(self.d.name, "toy.pretext")
        self.truth = toy_map(path, scaffolds, CHROMS, npix=npix, bp_per_px=BPP)
        self.pt = Pretext(path)
        bl = toy_blocks(scaffolds, BPP) if blocks else None
        self.blocks = [b for b in bl if b.q not in drop_blocks and b.t not in drop_blocks] if bl else None

    def run(self, **kw):
        return curate_allhaps(self.pt, phasing="hic", blocks=self.blocks, **kw)

    def layout(self, st):
        """The final map as runs of (chromosome + haplotype, first position, last position);
        a run continues while positions step by at most 2 px in either direction."""
        ms = MapState(self.pt).replay(st.edits)
        self.is_permutation = np.array_equal(np.sort(ms.pix), np.arange(self.pt.npix))
        seq = []
        for px in ms.pix:
            c = self.truth[int(px)]
            if c and seq and seq[-1][0] == f"{c[0]}{c[1]}" and abs(c[2] - seq[-1][2]) <= 2:
                seq[-1][2] = c[2]
            elif c:
                seq.append([f"{c[0]}{c[1]}", c[2], c[2]])
        return [tuple(x) for x in seq]

    def close(self):
        self.pt.close()
        self.d.cleanup()


def whole(layout, copy, n):
    """`copy` (e.g. 'A1') laid out in one piece, forward or reverse, covering 0..n-1."""
    runs = [x for x in layout if x[0] == copy]
    return len(runs) == 1 and sorted(runs[0][1:]) == [0, n - 1]


def rows_of(rows, action):
    return [r for r in rows if r["action"] == action]


class TestSites(unittest.TestCase):
    def test_misjoin_cut(self):
        """Two chromosomes in one scaffold: cut at the gap between them."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 100, False), ("B", 1, 0, 80, False)]),
                 ("H1_scaffold_3", [("C", 1, 0, 60, False)])] + H2_A + REST[2:])
        try:
            run_ = t.run(stop_after="_setup")
            self.assertEqual([t.pt.names[t.pt.orig_scaffold[q]] for q, _, _ in run_.cuts], ["H1_scaffold_1"])
        finally:
            t.close()

    def test_gap_inside_a_chromosome_kept(self):
        """A gap between two adjacent pieces of one chromosome is not cut."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 50, False), ("A", 1, 50, 100, False)])] + H2_A + REST)
        try:
            run_ = t.run(stop_after="_setup")
            self.assertEqual(run_.cuts, [])
        finally:
            t.close()


class TestPairing(unittest.TestCase):
    def test_pairs(self):
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 100, False)])] + H2_A + REST)
        try:
            run_ = t.run(stop_after="_setup")
            names = {frozenset((t.pt.names[p.a], t.pt.names[p.b])) for p in run_.pairs}
            self.assertEqual(names, {frozenset((f"H1_scaffold_{i}", f"H2_scaffold_{i}")) for i in (1, 2, 3)})
        finally:
            t.close()

    def test_split_in_both_haplotypes_merged(self):
        """A chromosome in two scaffolds in both haplotypes: the two pairs are one chromosome;
        each copy ends up whole. (1024 px: merge_pairs looks only at copies of 60 px or more.)"""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 60, False)]), ("H1_scaffold_4", [("A", 1, 60, 100, False)]),
                 ("H2_scaffold_1", [("A", 2, 0, 60, False)]), ("H2_scaffold_4", [("A", 2, 60, 100, False)])] + REST,
                npix=1024)
        try:
            run_ = t.run(stop_after="_merge")
            self.assertEqual({(t.pt.names[v], t.pt.names[u]) for v, u in run_.host_of.items()},
                             {("H1_scaffold_4", "H1_scaffold_1"), ("H2_scaffold_4", "H2_scaffold_1")})
            rows, st, _ = t.run()
            self.assertEqual([(r["scaffold"], r["partner"]) for r in rows_of(rows, "merge")],
                             [("H1_scaffold_4", "H2_scaffold_4")])
            lay = t.layout(st)
            self.assertTrue(whole(lay, "A1", 100) and whole(lay, "A2", 100), lay)
        finally:
            t.close()


class TestSmall(unittest.TestCase):
    def test_unique_signal_unloc(self):
        """A small piece (no alignment) whose MQ0 signal is unique to one pair: Unloc of that
        chromosome, laid out with it."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 94, False)]), ("H1_scaffold_9", [("A", 1, 94, 100, False)])]
                + H2_A + REST, drop_blocks=("H1_scaffold_9",))
        try:
            run_ = t.run(stop_after="_small")
            calls = {t.pt.names[c["u"]]: c for c in run_.small_calls}
            self.assertEqual(calls["H1_scaffold_9"]["action"], "unloc")
            self.assertEqual(t.pt.names[calls["H1_scaffold_9"]["copy"]], "H1_scaffold_1")
            rows, st, _ = t.run()
            self.assertEqual([r["tag"] for r in rows if r["scaffold"] == "H1_scaffold_9"], ["Unloc"])
        finally:
            t.close()

    def test_unrelated_piece_not_placed(self):
        """A small piece with no contact to any chromosome gets no call."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 100, False)]), ("H1_scaffold_9", [("R", 1, 0, 6, False)])]
                + H2_A + REST)
        try:
            run_ = t.run(stop_after="_small")
            self.assertNotIn("H1_scaffold_9", [t.pt.names[c["u"]] for c in run_.small_calls])
        finally:
            t.close()

    def test_small_piece_in_gap(self):
        """A small piece that belongs in a gap goes back into it (the gap is cut; the piece
        aligns to the other haplotype between the two sides)."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 47, False), ("A", 1, 53, 100, False)]),
                 ("H1_scaffold_9", [("A", 1, 47, 53, False)])] + H2_A + REST)
        try:
            rows, st, _ = t.run()
            lay = t.layout(st)
            self.assertTrue(t.is_permutation)
            self.assertTrue(whole(lay, "A1", 100), lay)
        finally:
            t.close()


class TestArrange(unittest.TestCase):
    def test_inverted_piece_turned_back(self):
        """A piece inverted in the assembly (at a gap) is cut and turned back."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 60, False), ("A", 1, 60, 100, True)])] + H2_A + REST)
        try:
            rows, st, _ = t.run()
            self.assertEqual(len(rows_of(rows, "cut")), 1)
            lay = t.layout(st)
            self.assertTrue(whole(lay, "A1", 100), lay)
        finally:
            t.close()

    def test_large_piece_into_its_gap(self):
        """A large piece missing from a gap of its chromosome goes back into the gap."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 40, False), ("A", 1, 70, 100, False)]),
                 ("H1_scaffold_5", [("A", 1, 40, 70, False)])] + H2_A + REST)
        try:
            rows, st, _ = t.run()
            self.assertIn("H1_scaffold_5", [r["scaffold"] for r in rows_of(rows, "join")])
            self.assertTrue(whole(t.layout(st), "A1", 100))
        finally:
            t.close()


class TestFollow(unittest.TestCase):
    def test_other_copy_follows_reference(self):
        """The other haplotype's copy, assembled in two scaffolds in the wrong order and one
        reversed, takes the reference copy's order and orientation."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 100, False)]), ("H2_scaffold_1", [("A", 2, 50, 100, False)]),
                 ("H2_scaffold_5", [("A", 2, 0, 50, True)])] + REST)
        try:
            rows, st, _ = t.run()
            self.assertEqual(len(rows_of(rows, "rearrange")), 1)
            lay = t.layout(st)
            self.assertTrue(whole(lay, "A2", 100), lay)
            run1 = [x for x in lay if x[0] == "A1"][0]
            run2 = [x for x in lay if x[0] == "A2"][0]
            self.assertEqual(run1[1] < run1[2], run2[1] < run2[2])        # same orientation
        finally:
            t.close()


class TestLayout(unittest.TestCase):
    def test_savestate(self):
        """The edits are a permutation; the copies of each chromosome sit side by side; every
        piece is painted and tagged HAP1 / HAP2."""
        t = Toy([("H1_scaffold_1", [("A", 1, 0, 100, False), ("B", 1, 0, 80, False)]),
                 ("H1_scaffold_3", [("C", 1, 0, 60, False)])] + H2_A + REST[2:])
        try:
            rows, st, _ = t.run()
            lay = t.layout(st)
            self.assertTrue(t.is_permutation)
            copies = [x[0] for x in lay]
            for ch in "ABC":
                i = copies.index(ch + "1")
                self.assertIn(ch + "2", copies[max(0, i - 1):i + 2], copies)
            ms = MapState(t.pt).replay(st.edits)
            flags, paint = dict(st.meta_flags), dict(st.scaffs)
            for i in range(len(ms.contigs)):
                self.assertTrue(set(st.tag_names(flags.get(i, 0))) & {"HAP1", "HAP2"})
                self.assertIn(i, paint)
        finally:
            t.close()


if __name__ == "__main__":
    unittest.main()
