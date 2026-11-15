"""Sex chromosomes from coverage and alignment (sexchrom): hemizygous stretches, PARs, sex
chromosomes joined to another chromosome, and what curate_allhaps does with them."""
import collections
import os
import random
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, SaveState, MapState
from pretext_autocurate import sexchrom
from pretext_autocurate.allhaps import curate_allhaps
from pretext_autocurate.homology import read_paf, find_tool, align_within

FIX = os.path.join(ROOT, "seed", "fixtures")
TRO = os.path.join(FIX, "qqTroRaco.hfsm.asm.bp.yhs_mq10.diploid_mq0.extensions.HR.pretext")
TRO_PAF = os.path.join(FIX, "qqTroRaco.H2_vs_H1.asm5.paf")
OWL = os.path.join(FIX, "bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid_mq0.extensions.HR.pretext")
OWL_PAF = os.path.join(FIX, "bAegFun_asm2_run1.H2_vs_H1.asm5.paf")


def _mb(pt, px):
    return pt.rel_coord[px] * pt.bp_per_pix / 1e6


def _stretches(pt, res):
    """{scaffold: (start Mb, end Mb, name of the stretch it aligns into or '')}"""
    return {pt.names[d["s"]]: (round(_mb(pt, d["a"]), 1), round(_mb(pt, d["b"]) + pt.bp_per_pix / 1e6, 1),
                               pt.names[d["via"][0]] if d["via"] else "") for d in res["stretches"]}


def _tags_by_bp(pt, st, wanted):
    """{tag: {scaffold: bp}} for the tags in `wanted`, from a savestate on map `pt`."""
    ms = MapState(pt).replay(st.edits)
    flags = dict(st.meta_flags)
    out = collections.defaultdict(collections.Counter)
    for i, c in enumerate(ms.contigs):
        for t in st.tag_names(flags.get(i, 0)):
            if t in wanted:
                out[t][pt.names[c.scaffold]] += c.length * pt.bp_per_pix
    return out


class TestLabels(unittest.TestCase):
    """sexchrom.sex_labels: the species' sex system decides the tags (curator, 2026-11-15)."""
    size = {"a": 90e6, "b": 30e6, "c": 50e6, "d": 10e6}
    L = staticmethod(lambda *p: {frozenset(x): 1.0 for x in p})

    def test_two_chromosomes(self):
        tags, note = sexchrom.sex_labels("XY", ["a", "b"], self.size, self.L(("a", "b")))
        self.assertEqual(tags, {"a": "X", "b": "Y"})
        tags, _ = sexchrom.sex_labels("ZW", ["b", "a"], self.size, {})        # unrelated: by size
        self.assertEqual(tags, {"a": "Z", "b": "W"})
        tags, note = sexchrom.sex_labels("XY", ["c"], self.size, {})
        self.assertEqual(tags, {"c": "X"})
        self.assertIn("Y not found", note)
        self.assertEqual(sexchrom.sex_labels(None, ["a"], self.size, {}), ({}, ""))

    def test_single_chromosome(self):
        tags, note = sexchrom.sex_labels("X0", ["a", "c"], self.size, {})       # an X in two pieces
        self.assertEqual(tags, {"a": "X", "c": "X"})
        _, note = sexchrom.sex_labels("Z0", ["a", "b"], self.size, self.L(("a", "b")))
        self.assertIn("check the sex system", note)

    def test_three_chromosomes(self):
        # X1X2Y: the Y (b) carries parts of both X's
        tags, _ = sexchrom.sex_labels("X1X2Y", ["a", "b", "c"], self.size, self.L(("a", "b"), ("b", "c")))
        self.assertEqual(tags, {"b": "Y", "a": "X1", "c": "X2"})
        tags, _ = sexchrom.sex_labels("ZW1W2", ["a", "b", "c", "d"], self.size, self.L(("a", "b"), ("a", "d")))
        self.assertEqual(tags, {"a": "Z", "b": "W1", "d": "W2"})
        tags, note = sexchrom.sex_labels("X1X2Y", ["a", "b"], self.size, self.L(("a", "b")))
        self.assertEqual(tags, {})
        self.assertIn("not tagged", note)


class TestRuns(unittest.TestCase):
    def test_short_runs_are_absorbed(self):
        st = np.array([1] * 10 + [0] * 2 + [1] * 10 + [0] * 20 + [3] * 1 + [0] * 5)
        runs = sexchrom._runs(st, 3)
        self.assertEqual([r[0] for r in runs], [1, 0])
        self.assertEqual((runs[0][1], runs[0][2], runs[1][1], runs[1][2]), (0, 21, 22, 47))


@unittest.skipUnless(os.path.exists(TRO) and os.path.exists(TRO_PAF), "missing fixture")
class TestTroglobisium(unittest.TestCase):
    """The curator's reading (savestate_final): X = H1_scaffold_1 25-119.5 Mb with
    H2_scaffold_46; Y = H2_scaffold_4 0-33.3 Mb; an autosome joined to both (H1_scaffold_1
    0-25 Mb, H2_scaffold_4 33-58 Mb)."""

    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(TRO)
        cls.blocks = read_paf(TRO_PAF, min_mapq=0)
        cls.res = sexchrom.analyse(cls.pt, cls.blocks)

    def test_stretches(self):
        s = _stretches(self.pt, self.res)
        self.assertEqual(set(s), {"H1_scaffold_1", "H2_scaffold_4", "H2_scaffold_46"})
        self.assertEqual(s["H1_scaffold_1"][:2], (25.0, 119.5))           # X: hemizygous and gametologous parts
        self.assertEqual(s["H2_scaffold_4"], (0.0, 33.4, "H1_scaffold_1"))  # Y: aligns into the X, both half coverage
        self.assertEqual(s["H2_scaffold_46"][:2], (0.0, 23.1))
        # the divergent heterozygous stretches of autosomes (half coverage, aligned) are not sex-linked
        self.assertNotIn("H1_scaffold_4", s)

    def test_autosome_joined_not_par(self):
        """The normal-coverage ends align to each other like PARs would, but share no
        chromosome-level contact with the sex chromosome: another chromosome joined to it."""
        ends = {self.pt.names[e["s"]]: e for e in self.res["ends"]}
        self.assertEqual(set(ends), {"H1_scaffold_1", "H2_scaffold_4"})
        for e in ends.values():
            self.assertEqual(e["kind"], "joined")
            self.assertLess(e["contact"], 0.03)
        self.assertFalse((self.res["par_mask"] > 0).any())

    def test_levels_and_map_only(self):
        for lv in (0.0, 1.0):
            self.assertEqual(set(_stretches(self.pt, sexchrom.analyse(self.pt, self.blocks, level=lv))),
                             {"H1_scaffold_1", "H2_scaffold_4", "H2_scaffold_46"})
        # without alignments: the same stretches (the Y now from coverage alone), no gametolog link
        s = _stretches(self.pt, sexchrom.analyse(self.pt))
        self.assertEqual(set(s), {"H1_scaffold_1", "H2_scaffold_4", "H2_scaffold_46"})
        self.assertFalse(any(v[2] for v in s.values()))
        # both sex-linked scaffolds of H2 are named for an alignment within the haplotype
        self.assertEqual(sexchrom.same_hap_candidates(self.pt, self.blocks), {"H2": ["H2_scaffold_4", "H2_scaffold_46"]})

    def test_curation(self):
        rows, st, _ = curate_allhaps(self.pt, phasing="hic", blocks=self.blocks, sex_system="XY")
        by = collections.defaultdict(list)
        for r in rows:
            by[r["action"]].append(r)
        # H2_scaffold_46 (in no pair) goes into the X's chain by contact
        self.assertEqual([(r["scaffold"], r["partner"].split(":")[0]) for r in by["sex_chromosome_join"]],
                         [("H2_scaffold_46", "H1_scaffold_1")])
        self.assertFalse(by["sex_chromosome_candidate"])
        # the cuts at the coverage steps say what they are
        steps = [r for r in by["cut"] if "a sex chromosome joined to another chromosome" in r["evidence"]]
        self.assertEqual({(r["scaffold"], round(r["length_mb"])) for r in steps}, {("H1_scaffold_1", 25), ("H2_scaffold_4", 33)})
        self.assertFalse(by["par"] or by["sex_boundary"])
        self.assertEqual(len(by["sex_linked_region"]), 3)
        self.assertIn("X the larger and Y the smaller", by["sex_summary"][0]["evidence"])
        # X and Y tags where the curator put them (bp of each scaffold carrying the tag)
        tags = _tags_by_bp(self.pt, st, {"X", "Y"})
        cur = _tags_by_bp(self.pt, SaveState.load(TRO + ".savestate_final"), {"X", "Y"})
        for t in "XY":
            self.assertEqual(set(tags[t]), set(cur[t]), t)
            for sc in cur[t]:
                self.assertAlmostEqual(tags[t][sc] / cur[t][sc], 1.0, delta=0.03)
        # one X: H1_scaffold_1's sex-linked pieces and H2_scaffold_46 painted alike
        ms = MapState(self.pt).replay(st.edits)
        paint = dict(st.scaffs)
        x_ids = {paint.get(i) for i, c in enumerate(ms.contigs)
                 if self.pt.names[c.scaffold] == "H2_scaffold_46"
                 or (self.pt.names[c.scaffold] == "H1_scaffold_1" and _mb(self.pt, ms.pix[c.start]) > 30)}
        self.assertEqual(len(x_ids), 1)
        self.assertNotIn(None, x_ids)


@unittest.skipUnless(os.path.exists(OWL) and os.path.exists(OWL_PAF), "missing fixture")
class TestOwl(unittest.TestCase):
    """Boreal owl, a female: Z = H1_scaffold_3, W = H2_scaffold_8, a PAR at the start of the Z
    and the end of the W (normal coverage, aligned to each other)."""

    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(OWL)
        cls.blocks = read_paf(OWL_PAF, min_mapq=0)
        cls.res = sexchrom.analyse(cls.pt, cls.blocks)

    def test_stretches_and_pars(self):
        s = _stretches(self.pt, self.res)
        self.assertEqual(set(s), {"H1_scaffold_3", "H2_scaffold_8"})
        pars = {self.pt.names[e["s"]]: e for e in self.res["ends"] if e["kind"] == "par"}
        self.assertEqual({k: v["partner"] for k, v in pars.items()}, {"H1_scaffold_3": "H2_scaffold_8", "H2_scaffold_8": "H1_scaffold_3"})
        z0, z1 = pars["H1_scaffold_3"]["par"]
        w0, w1 = pars["H2_scaffold_8"]["par"]
        self.assertLess(z0, 0.1e6)
        self.assertAlmostEqual(z1 / 1e6, 6.0, delta=0.3)
        self.assertAlmostEqual(w0 / 1e6, 28.5, delta=0.3)
        self.assertAlmostEqual(w1 / 1e6, 34.5, delta=0.1)

    def test_nothing_hemizygous(self):
        """Homogametic sample or homomorphic sex chromosomes: flat coverage, nothing found."""
        cov = self.pt.graphs["coverage"]
        try:
            self.pt.graphs["coverage"] = np.full_like(cov, 50)
            res = sexchrom.analyse(self.pt, self.blocks)
        finally:
            self.pt.graphs["coverage"] = cov
        self.assertEqual((res["stretches"], res["ends"]), ([], []))

    def test_par_from_hic_alone(self):
        """Without alignments the two PARs are found by their MQ0 contact (reads map to both)."""
        res = sexchrom.analyse(self.pt)
        pars = {self.pt.names[e["s"]]: e for e in res["ends"] if e["kind"] == "par"}
        self.assertEqual(set(pars), {"H1_scaffold_3", "H2_scaffold_8"})
        self.assertTrue(all(e["by"] == "MQ0 Hi-C" for e in pars.values()))

    def test_tracks(self):
        with tempfile.TemporaryDirectory() as d:
            got = sexchrom.write_tracks(self.pt, self.res, d)
            self.assertEqual(set(got), {"sex_chromosome", "par"})
            def lines(p):
                with open(p, encoding="utf-8") as f:
                    return [l.split() for l in f]
            par = lines(got["par"])
            self.assertEqual({x[0] for x in par}, {"H1_scaffold_3", "H2_scaffold_8"})
            self.assertTrue(all(x[3] == "2" for x in par))
            sx = lines(got["sex_chromosome"])
            self.assertEqual({x[0] for x in sx}, {"H1_scaffold_3", "H2_scaffold_8"})

    def test_curation(self):
        rows, st, _ = curate_allhaps(self.pt, phasing="hic", blocks=self.blocks, sex_system="ZW")
        pars = [r for r in rows if r["action"] == "par"]
        self.assertEqual({r["scaffold"] for r in pars}, {"H1_scaffold_3", "H2_scaffold_8"})
        self.assertFalse([r for r in rows if r["action"] in ("sex_boundary", "sex_chromosome_join")])
        tags = _tags_by_bp(self.pt, st, {"Z", "W"})
        self.assertEqual(set(tags["Z"]), {"H1_scaffold_3"})
        # the W chain: H2_scaffold_8, and H2_scaffold_34 that end contact puts in front of it
        self.assertLessEqual(set(tags["W"]), {"H2_scaffold_8", "H2_scaffold_34"})
        self.assertGreater(tags["W"]["H2_scaffold_8"], 30e6)
        summary = [r for r in rows if r["action"] == "sex_summary"]
        self.assertEqual(len(summary), 1)
        self.assertIn("2 PARs", summary[0]["evidence"])
        # without --sex-system: no Z/W tags
        _, st0, _ = curate_allhaps(self.pt, phasing="hic", blocks=self.blocks, do_place=False, do_arrange=False)
        self.assertFalse(_tags_by_bp(self.pt, st0, {"X", "Y", "Z", "W"}))


@unittest.skipUnless(find_tool("minimap2"), "minimap2 not found")
class TestAlignWithin(unittest.TestCase):
    def test_shared_region_found(self):
        rnd = random.Random(1)
        seq = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))
        par, x, y, z = seq(60000), seq(100000), seq(80000), seq(100000)
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "a.fa")
            with open(fa, "w") as f:
                f.write(f">H1_x\n{par}{x}\n>H1_y\n{y}{par}\n>H1_z\n{z}\n")
            paf = align_within(fa, ["H1_x", "H1_y"], os.path.join(d, "w"), threads=1)
            b = [b for b in read_paf(paf, min_len=10000, min_mapq=0) if b.q != b.t]
            self.assertTrue(b)
            self.assertEqual({frozenset((x.q, x.t)) for x in b}, {frozenset(("H1_x", "H1_y"))})
            self.assertEqual(paf, align_within(fa, ["H1_y", "H1_x"], os.path.join(d, "w")))      # cached


if __name__ == "__main__":
    unittest.main()
