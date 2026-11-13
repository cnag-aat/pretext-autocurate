"""allhaps step 1 on the bird (bAegFun): pairing, orientation, tags."""
import collections
import glob
import os
import sys
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import Pretext, SaveState, MapState
from pretext_autocurate.allhaps import curate_allhaps

M0 = os.path.join(ROOT, "seed", "fixtures", "bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid_mq0.extensions.HR.pretext")


@unittest.skipUnless(os.path.exists(M0), "missing fixture")
class TestAllhaps(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(M0)
        cls.hi = Pretext(M0.replace("diploid_mq0", "diploid_mq10"), load_graphs=False)
        cls.rows, cls.st, cls.info = curate_allhaps(cls.pt, phasing="hic", hi=cls.hi)
        cur = SaveState.load(M0 + ".savestate_6")
        ms = MapState(cls.pt).replay(cur.edits)
        paint = dict(cur.scaffs)
        groups = collections.defaultdict(list)
        for i, c in enumerate(ms.contigs):
            if paint.get(i):
                groups[paint[i]].append((c.length, cls.pt.names[c.scaffold], c.inverted))
        cls.truth = {}
        for k in groups:
            if k % 2 == 1 and k + 1 in groups:
                a, b = max(groups[k]), max(groups[k + 1])
                cls.truth[frozenset((a[1], b[1]))] = a[2] != b[2]
        cls.split = {n for n, c in collections.Counter(x[1] for g in groups.values() for x in g).items() if c > 1}

    @staticmethod
    def base(name):
        return name.split(":")[0]            # pieces of a cut scaffold are named parent:k

    def test_pairs_are_curator_pairs(self):
        tool = {frozenset((self.base(r["scaffold"]), self.base(r["partner"])))
                for r in self.rows if r["action"].startswith("pair")}
        self.assertGreaterEqual(len(tool), 36)
        # H1_scaffold_31 / H2_scaffold_37 is not in the original curation but the curator paired
        # them in savestate_curator_1 (2026-11-13)
        self.assertEqual(tool - set(self.truth) - {frozenset(("H1_scaffold_31", "H2_scaffold_37"))}, set())

    def test_telomere_sites(self):
        """Telomere-at-a-gap sites: misjoins are cut (also at the homologous point of the other
        haplotype), scaffold_2 stays joined with the H1 end piece inverted (curator, 2026-11-13)."""
        act = {(r["scaffold"], round(r["length_mb"])): r["action"] for r in self.rows
               if r["action"] in ("cut", "invert_tip", "keep")}
        for site in (("H1_scaffold_1", 135), ("H2_scaffold_1", 135), ("H1_scaffold_24", 9),
                     ("H2_scaffold_30", 3), ("H2_scaffold_31", 3)):
            self.assertEqual(act[site], "cut", site)
        self.assertEqual(act[("H1_scaffold_2", 166)], "invert_tip")
        self.assertEqual(act[("H2_scaffold_2", 5)], "keep")
        why = {r["scaffold"]: r["evidence"] for r in self.rows if r["action"] == "cut"}
        self.assertIn("homologous", why["H2_scaffold_1"])
        self.assertIn("different homologs", why["H2_scaffold_30"])
        self.assertEqual(len(self.st.waypoints), len(act))
        # the inverted end piece of H1_scaffold_2 is still attached to its scaffold, reversed
        ms = MapState(self.pt).replay(self.st.edits)
        pos = np.empty(self.pt.npix, dtype=int)
        pos[ms.pix] = np.arange(self.pt.npix)
        s0 = self.pt.scaffold_start[self.pt.names.index("H1_scaffold_2")]
        px = lambda mb: int(s0 + mb * 1e6 / self.pt.bp_per_pix)
        body = pos[px(160)] - pos[px(150)]
        tip = pos[px(169.5)] - pos[px(167.5)]
        self.assertEqual(np.sign(body), -np.sign(tip))
        self.assertLess(abs(pos[px(170.4)] - pos[px(166.2)]), 10)      # far end of the tip now touches the body

    def test_orientation_is_parallel_as_curated(self):
        n = 0
        for r in self.rows:
            key = frozenset((self.base(r["scaffold"]), self.base(r["partner"])))
            if (r["action"] == "pair" and key in self.truth and not key & self.split
                    and ":" not in r["scaffold"] + r["partner"] and "H1_scaffold_2" not in key):
                self.assertEqual(r["strand"] == "-", self.truth[key], r)
                n += 1
        self.assertGreaterEqual(n, 22)

    def test_sex_chromosomes_flagged_and_both_hap1(self):
        r = [r for r in self.rows if r["action"] == "pair_unequal"]
        self.assertEqual({(x["scaffold"], x["partner"]) for x in r}, {("H1_scaffold_3", "H2_scaffold_8")})   # Z, W
        self.assertEqual(r[0]["tag"], "HAP1/HAP1")
        self.assertEqual(r[0]["strand"], "+")

    def test_no_prefix_mixing_and_layout(self):
        ms = MapState(self.pt).replay(self.st.edits)
        self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(self.pt.npix)))
        n_sites = sum(r["action"] in ("cut", "invert_tip") for r in self.rows)
        self.assertEqual(len(ms.contigs), int((self.pt.scaffold_npix > 0).sum()) + n_sites)
        paint, flags = dict(self.st.scaffs), dict(self.st.meta_flags)
        ids = collections.defaultdict(set)
        for i, c in enumerate(ms.contigs):
            if paint.get(i):
                ids[paint[i]].add(self.pt.names[c.scaffold].split("_")[0])
        self.assertTrue(all(len(v) == 1 for v in ids.values()))                        # one prefix per superscaffold
        # the two members of a pair are adjacent in the map
        pos = {self.pt.names[c.scaffold]: i for i, c in enumerate(ms.contigs)}
        for r in self.rows:
            if r["action"] == "pair" and ":" not in r["scaffold"] + r["partner"]:
                self.assertLessEqual(abs(pos[r["scaffold"]] - pos[r["partner"]]), 2, r)   # 2: an inverted end piece
        self.assertEqual(len(flags), len(ms.contigs))                                  # every piece has a HAP tag

    def test_trio_keeps_prefix_tags(self):
        rows, st, _ = curate_allhaps(self.pt, phasing="trio", hi=self.hi)
        for r in rows:
            if r["action"] == "pair":
                self.assertTrue(r["scaffold"].startswith("H1_") and r["tag"] == "HAP1/HAP2", r)


PAF = os.path.join(ROOT, "seed", "fixtures", "bAegFun_asm2_run1.H2_vs_H1.asm5.paf")


@unittest.skipUnless(os.path.exists(M0) and os.path.exists(PAF), "missing fixture")
class TestAllhapsSequence(unittest.TestCase):
    """Same pass with homology from a minimap2 alignment of H2 against H1 (--paf / --fasta)."""

    @classmethod
    def setUpClass(cls):
        from pretext_autocurate.homology import read_paf
        cls.pt = Pretext(M0)
        hi = Pretext(M0.replace("diploid_mq0", "diploid_mq10"), load_graphs=False)
        cls.rows, cls.st, cls.info = curate_allhaps(cls.pt, phasing="hic", hi=hi, blocks=read_paf(PAF, min_mapq=0))

    def test_sites(self):
        act = {(r["scaffold"], round(r["length_mb"])): r["action"] for r in self.rows
               if r["action"] in ("cut", "invert_tip", "keep")}
        cuts = {k for k, v in act.items() if v == "cut"}
        self.assertEqual(cuts, {("H1_scaffold_1", 135), ("H2_scaffold_1", 135), ("H1_scaffold_24", 9),
                                ("H2_scaffold_30", 3), ("H2_scaffold_31", 3)})          # all curator cuts
        self.assertEqual(act[("H1_scaffold_2", 166)], "invert_tip")
        self.assertEqual(act[("H2_scaffold_2", 5)], "keep")

    def test_pairs(self):
        base = lambda n: n.split(":")[0]
        got = {frozenset((base(r["scaffold"]), base(r["partner"]))) for r in self.rows
               if r["action"] in ("pair", "pair_unequal", "place_next_to_homolog")}
        for a, b in (("H1_scaffold_1", "H2_scaffold_1"), ("H1_scaffold_42", "H2_scaffold_36"),
                     ("H1_scaffold_37", "H2_scaffold_31"), ("H1_scaffold_34", "H2_scaffold_35"),
                     ("H1_scaffold_36", "H2_scaffold_40"), ("H1_scaffold_5", "H2_scaffold_39"),
                     ("H1_scaffold_3", "H2_scaffold_8")):
            self.assertIn(frozenset((a, b)), got)
        self.assertGreaterEqual(len(got), 38)
        self.assertEqual(self.info["source"], "sequence alignment")
        # H2_scaffold_37 is too repetitive to assign (curator, 2026-11-13): unpaired, with the reason
        self.assertFalse(any("H2_scaffold_37" in k for k in got))
        why = [r["evidence"] for r in self.rows if r["action"] == "unpaired" and r["scaffold"] == "H2_scaffold_37"]
        self.assertTrue(why and "too repetitive" in why[0])

    def test_busco(self):
        """BUSCO counts per haplotype; in this assembly the copies of a pair hardly differ, so the
        HAP1 choice still rests on gaps, telomeres and length for most pairs."""
        from pretext_autocurate import busco
        from pretext_autocurate.homology import read_paf
        mp = os.path.join(ROOT, "seed", "fixtures", "bAegFun_asm2_run1.aves_odb12.ancestral.miniprot.paf")
        if not os.path.exists(mp):
            self.skipTest("missing miniprot fixture")
        hits = busco.read_hits(mp)
        self.assertGreater(len(hits), 6000)
        rows, st, _ = curate_allhaps(self.pt, phasing="hic", blocks=read_paf(PAF, min_mapq=0), busco_hits=hits, do_cuts=False)
        summ = {r["scaffold"]: r["evidence"] for r in rows if r["action"] == "busco_summary"}
        self.assertEqual(set(summ), {"H1", "H2"})
        pct = lambda t: float(t.split("(")[1].split("%")[0])
        self.assertGreater(pct(summ["H1"]), 97)
        self.assertGreater(pct(summ["H1"]), pct(summ["H2"]))          # H2 lacks the Z chromosome's genes
        by_busco = sum("BUSCO genes is HAP1" in r["evidence"] for r in rows if r["action"] == "pair")
        self.assertLess(by_busco, 10)

    def test_busco_full_table(self):
        import tempfile
        from pretext_autocurate import busco
        text = ("# BUSCO version is: 6.0.0\n# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
                "1at8782\tDuplicated\tH1_scaffold_2\t200\t100\t-\t204.7\t87\n"
                "1at8782\tDuplicated\tH2_scaffold_2\t300\t400\t+\t201.0\t87\n"
                "2at8782\tComplete\tH1_scaffold_3\t10\t90\t+\t374.0\t165\n"
                "3at8782\tMissing\n4at8782\tFragmented\tH1_scaffold_2\t5\t9\t+\t50.0\t20\n")
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as f:
            f.write(text)
        try:
            hits = busco.read_full_table(f.name)
        finally:
            os.remove(f.name)
        self.assertEqual(set(hits), {"1at8782", "2at8782"})
        self.assertEqual(hits["1at8782"][0][:4], ("H1_scaffold_2", 100, 200, "-"))
        haps, per, _ = busco.summarise(hits, lambda n: n.split("_")[0], n_genes=4)
        self.assertEqual((haps["H1"]["single"], haps["H1"]["duplicated"], haps["H1"]["missing"]), (2, 0, 2))
        self.assertEqual(haps["H2"]["single"], 1)                     # one copy per haplotype is not a duplication

    def test_busco_per_haplotype_tables_and_false_duplicates(self):
        """Tables made on each haplotype alone get the prefix added; a small scaffold whose genes
        all have a second copy on one region of another scaffold of the same haplotype is reported."""
        import tempfile
        from pretext_autocurate import busco
        head = "# Busco id\tStatus\tSequence\tGene Start\tGene End\tStrand\tScore\tLength\n"
        h1 = head + "".join(f"{g}at8782\tDuplicated\tscaffold_14\t{13700000 + 1000 * g}\t{13700500 + 1000 * g}\t-\t400\t185\n"
                            f"{g}at8782\tDuplicated\tscaffold_192\t{1000 * g}\t{1000 * g + 500}\t-\t400\t185\n" for g in range(1, 6))
        h1 += "".join(f"{g}at8782\tComplete\tscaffold_14\t{1000 * g}\t{1000 * g + 500}\t+\t300\t100\n" for g in range(10, 30))
        h2 = head + "1at8782\tComplete\tscaffold_14\t5\t500\t+\t400\t185\n"
        d = tempfile.mkdtemp()
        try:
            for name, text in (("a.tsv", h1), ("b.tsv", h2)):
                with open(os.path.join(d, name), "w") as f:
                    f.write(text)
            hits = busco.read_full_tables([f"H1={d}/a.tsv", f"H2={d}/b.tsv"])
        finally:
            import shutil
            shutil.rmtree(d)
        self.assertEqual({l[0] for l in hits["1at8782"]}, {"H1_scaffold_14", "H1_scaffold_192", "H2_scaffold_14"})
        fd = busco.false_duplicates(hits, lambda n: n.split("_")[0])
        self.assertEqual([(r[0], r[1], r[4]) for r in fd], [("H1_scaffold_192", "H1_scaffold_14", 5)])
        # and on the real bird (miniprot fixture): the three the curator's BUSCO tables show
        mp = os.path.join(ROOT, "seed", "fixtures", "bAegFun_asm2_run1.aves_odb12.ancestral.miniprot.paf")
        if os.path.exists(mp):
            real = {(r[0], r[1]) for r in busco.false_duplicates(busco.read_hits(mp), lambda n: n.split("_")[0])}
            self.assertEqual(real, {("H1_scaffold_192", "H1_scaffold_14"), ("H2_scaffold_43", "H2_scaffold_3"),
                                    ("H1_scaffold_209", "H1_scaffold_12")})

    def test_curator_busco_tables(self):
        """The curator's own per-haplotype BUSCO full tables (names not prefixed; hap in the file name)."""
        import glob as _glob
        from pretext_autocurate import busco
        t1 = _glob.glob(os.path.join(ROOT, "seed", "fixtures", "*hap1*full_table.tsv"))
        t2 = _glob.glob(os.path.join(ROOT, "seed", "fixtures", "*hap2*full_table.tsv"))
        if not (t1 and t2):
            self.skipTest("missing BUSCO tables")
        hits = busco.read_full_tables([f"H1={t1[0]}", f"H2={t2[0]}"])
        names = {l[0] for loci in hits.values() for l in loci}
        self.assertTrue(names <= set(self.pt.names))                    # every prefixed name is in the map
        haps, _, _ = busco.summarise(hits, lambda n: n.split("_")[0], n_genes=6251)
        self.assertGreater(haps["H1"]["complete"], 6200)
        self.assertGreater(haps["H1"]["complete"], haps["H2"]["complete"])
        fd = {(r[0], r[1]) for r in busco.false_duplicates(hits, lambda n: n.split("_")[0])}
        self.assertEqual(fd, {("H1_scaffold_192", "H1_scaffold_14"), ("H2_scaffold_43", "H2_scaffold_3"),
                              ("H1_scaffold_209", "H1_scaffold_12")})

    def test_tools_are_found_or_reported(self):
        from pretext_autocurate.homology import find_tool
        self.assertIsNone(find_tool("no-such-aligner-xyz"))


if __name__ == "__main__":
    unittest.main()
