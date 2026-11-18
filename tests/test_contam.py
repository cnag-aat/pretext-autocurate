"""Non-target sequence (contam.py): lists, Tiara output, the Hi-C check, tags, organelle split."""
import os
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from pretext_autocurate import MapState, Pretext, contam       # noqa: E402
from pretext_autocurate.solver import solve                    # noqa: E402
from pretext_autocurate.homology import Block                  # noqa: E402
from pretext_autocurate.synth import toy_blocks, toy_map       # noqa: E402
from test_synthetic import CHROMS                              # noqa: E402

# F: a foreign genome (one copy) in three pieces, shared out over both haplotypes as hifiasm does;
# H1_scaffold_6 is a nuclear piece (the end of C) that Tiara would wrongly call prokaryotic
FOREIGN = [("H1_scaffold_1", [("A", 1, 0, 100, False)]), ("H1_scaffold_2", [("B", 1, 0, 80, False)]),
           ("H1_scaffold_3", [("C", 1, 0, 40, False)]),
           ("H1_scaffold_4", [("F", 1, 0, 12, False)]), ("H1_scaffold_6", [("C", 1, 40, 50, False)]),
           ("H1_scaffold_5", [("F", 1, 24, 36, False)]),
           ("H2_scaffold_1", [("A", 2, 0, 100, False)]), ("H2_scaffold_2", [("B", 2, 0, 80, False)]),
           ("H2_scaffold_3", [("C", 2, 0, 50, False)]), ("H2_scaffold_4", [("F", 1, 12, 24, False)])]
TIARA = {"H1_scaffold_4": "bacteria", "H1_scaffold_5": "bacteria", "H2_scaffold_4": "bacteria",
         "H1_scaffold_6": "archaea"}


class TestContam(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.dir = tempfile.TemporaryDirectory()
        p = os.path.join(cls.dir.name, "foreign.pretext")
        toy_map(p, FOREIGN, dict(CHROMS, F=36))
        cls.pt = Pretext(p)

    @classmethod
    def tearDownClass(cls):
        cls.dir.cleanup()

    def test_read_list(self):
        p = os.path.join(self.dir.name, "list.txt")
        with open(p, "w") as f:
            f.write("# non-target\nctg1\nctg2,bacteria\nctg3 mitochondrion\n\nctg4\tplastid\textra\n")
        self.assertEqual(contam.read_list(p), {"ctg1": "listed", "ctg2": "bacteria", "ctg3": "mitochondrion",
                                               "ctg4": "plastid"})
        # this tool's own table: the class column is named in the header
        with open(p, "w") as f:
            f.write("scaffold\tlength_mb\tclass\tsource\nctg5\t0.022\tmitochondrion\tTiara\n")
        self.assertEqual(contam.read_list(p), {"ctg5": "mitochondrion"})

    def test_read_tiara(self):
        p = os.path.join(self.dir.name, "tiara.txt")
        with open(p, "w") as f:
            f.write("sequence_id\tclass_fst_stage\tclass_snd_stage\n"
                    "s1\teukarya\tn/a\ns2\tbacteria\tn/a\ns3 desc\torganelle\tmitochondrion\n"
                    "s4\torganelle\tplastid\ns5\tunknown\tn/a\ns6\tarchaea\tn/a\ns7\torganelle\tunknown\n")
        self.assertEqual(contam.read_tiara(p), {"s2": "bacteria", "s3": "mitochondrion", "s4": "plastid",
                                                "s6": "archaea", "s7": "organelle"})

    def test_confirm(self):
        # the foreign pieces contact each other, not the chromosomes; the nuclear piece its chromosome
        kept, unconf = contam.confirm(self.pt, TIARA, {n: "Tiara" for n in TIARA}, big_bp=300_000)
        self.assertEqual(set(kept), {"H1_scaffold_4", "H1_scaffold_5", "H2_scaffold_4"})
        self.assertEqual([u[0] for u in unconf], ["H1_scaffold_6"])
        self.assertIn("H1_scaffold_3", unconf[0][4])
        # a listed name is kept as given
        kept, unconf = contam.confirm(self.pt, TIARA, {n: "list.txt" for n in TIARA}, big_bp=300_000)
        self.assertEqual((len(kept), unconf), (4, []))
        # an organelle of unknown kind is never enough; a mitochondrion needs organelle copy number
        # (these are at haploid coverage); a call that aligns to a chromosome is nuclear
        calls = {"H1_scaffold_4": "organelle", "H1_scaffold_5": "mitochondrion", "H2_scaffold_4": "bacteria"}
        on_c = Block("H2_scaffold_4", 0, 100_000, "H1_scaffold_1", 0, 100_000, True, 100_000, 60)
        kept, unconf = contam.confirm(self.pt, calls, dict.fromkeys(calls, "Tiara"), big_bp=300_000, blocks=[on_c])
        why = {u[0]: u[4] for u in unconf}
        self.assertEqual(kept, {})
        self.assertIn("unknown kind", why["H1_scaffold_4"])
        self.assertIn("coverage", why["H1_scaffold_5"])
        self.assertIn("aligns to H1_scaffold_1", why["H2_scaffold_4"])

    def test_solver_keeps_them_apart_and_tag(self):
        pt = self.pt
        off = {"H1_scaffold_4", "H1_scaffold_5", "H2_scaffold_4"}
        rows, st, info = solve(pt, blocks=toy_blocks(FOREIGN), exclude=off)
        for r in rows:
            if r["action"] in ("join", "pair"):
                self.assertFalse({r["scaffold"].split(":")[0], r["partner"].split(":")[0]} & off, r)
        ms = MapState(pt).replay(st.edits)
        first = lambda c: pt.names[int(sum(ms.pix[c.start] >= pt.scaffold_start) - 1)]
        # unplaced shrapnel carries its prefix's tag (H1_ HAP1, H2_ HAP2), non-target nothing yet
        # (a pair's copies are HAP1 / HAP2 by which is better, whatever their prefix)
        flags, painted = dict(st.meta_flags), dict(st.scaffs)
        for i, c in enumerate(ms.contigs):
            if painted.get(i):
                continue
            want = set() if first(c) in off else {"HAP" + first(c)[1]}
            self.assertTrue(want <= set(st.tag_names(flags.get(i, 0))), (first(c), st.tag_names(flags.get(i, 0))))
        # as the rules engine would: a haplotype tag on the non-target pieces too, which tag() replaces
        st.meta_flags = [(i, flags.get(i, 0) | 1 << st.tag_id("HAP" + first(c)[1])) for i, c in enumerate(ms.contigs)]
        mixed = contam.tag(pt, st, dict.fromkeys(off, "bacteria"))
        self.assertEqual(mixed, [])
        flags = dict(st.meta_flags)
        bit = 1 << st.tag_id("Contaminant")
        tagged = {pt.names[int(sum(ms.pix[c.start] >= pt.scaffold_start) - 1)]
                  for i, c in enumerate(ms.contigs) if flags.get(i, 0) & bit}
        self.assertEqual(tagged, off)
        # Contaminant and nothing else: no HAP1/HAP2, no paint
        painted = dict(st.scaffs)
        for i in range(len(ms.contigs)):
            if flags.get(i, 0) & bit:
                self.assertEqual(st.tag_names(flags[i]), ["Contaminant"])
                self.assertFalse(painted.get(i))

    def test_split_organelles(self):
        d = os.path.join(self.dir.name, "asm")
        os.makedirs(d, exist_ok=True)
        fa = os.path.join(d, "xxTest1.1.contaminant.fa")
        with open(fa, "w") as f:
            f.write(">scaffold_1\nACGT\n>scaffold_2\nGGCC\n>scaffold_3\nTTAA\n")
        with open(fa + ".agp", "w") as f:
            f.write("##agp-version 2.1\n"
                    "scaffold_1\t1\t4\t1\tW\tH1_scaffold_9\t1\t4\t+\n"
                    "scaffold_2\t1\t4\t1\tW\tH2_scaffold_7\t1\t4\t+\n"
                    "scaffold_3\t1\t4\t1\tW\tH1_scaffold_8\t1\t4\t+\n")
        out = contam.split_organelles(d, {"H1_scaffold_9": "bacteria", "H2_scaffold_7": "mitochondrion",
                                          "H1_scaffold_8": "plastid"})
        mito = os.path.join(d, "xxTest1.1.contaminant.mito.fa")
        self.assertEqual(out, {mito: 1, os.path.join(d, "xxTest1.1.contaminant.plastid.fa"): 1})
        self.assertEqual(open(mito).read(), ">scaffold_2\nGGCC\n")
        self.assertEqual(open(fa).read(), ">scaffold_1\nACGT\n")
        self.assertIn("H2_scaffold_7", open(mito + ".agp").read())
        self.assertNotIn("H2_scaffold_7", open(fa + ".agp").read())


if __name__ == "__main__":
    unittest.main()
