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


TRO = os.path.join(ROOT, "seed", "fixtures", "qqTroRaco.hfsm.asm.bp.yhs_mq10.diploid_mq0.extensions.HR.pretext")


@unittest.skipUnless(os.path.exists(TRO), "missing fixture")
class TestAllhapsPseudoscorpion(unittest.TestCase):
    """Troglobisium racovitzai: an allhaps map that needed heavy curation (38 cuts, joins, 47
    unlocs), MQ0 only, no sequence. Scored against the curator's final savestate."""

    @classmethod
    def setUpClass(cls):
        cls.pt = Pretext(TRO)
        cls.rows, cls.st, _ = curate_allhaps(cls.pt, phasing="hic")
        cur = SaveState.load(TRO + ".savestate_final")
        ms = MapState(cls.pt).replay(cur.edits)
        paint, flags = dict(cur.scaffs), dict(cur.meta_flags)
        cls.groups = collections.defaultdict(collections.Counter)      # scaffold -> {painted id: px}
        size, tag = collections.Counter(), {}
        for i, c in enumerate(ms.contigs):
            k = paint.get(i)
            if k:
                cls.groups[cls.pt.names[c.scaffold]][k] += c.length
                tags = cur.tag_names(flags.get(i, 0))
                if "Unloc" not in tags:
                    size[k] += c.length
                tag.setdefault(k, "HAP1" if "HAP1" in tags else "HAP2")
        cls.partner, ks, i = {}, sorted(size), 0
        while i < len(ks) - 1:
            a, b = ks[i], ks[i + 1]
            if tag[a] != tag[b] and min(size[a], size[b]) > 0.6 * max(size[a], size[b]):
                cls.partner[a], cls.partner[b] = b, a
                i += 2
            else:
                i += 1

    def test_pairs(self):
        self.assertEqual(len(self.partner), 120)
        good = bad = 0
        for r in self.rows:
            if r["action"] in ("pair", "pair_unequal"):
                ga, gb = self.groups[r["scaffold"].split(":")[0]], self.groups[r["partner"].split(":")[0]]
                # also homologous in the curation, though not painted as a HAP1/HAP2 pair of equal size:
                # X (1) with Y (2); the two halves of H1_scaffold_1's autosome (3, 4) with H2_scaffold_4's
                # (5); H1_scaffold_33 (70) with the pieces of H2_scaffold_6 (72, 73)
                also = {(1, 2), (3, 5), (4, 5), (70, 72), (70, 73)}
                if any(self.partner.get(x) in gb for x in ga) or any((min(x, y), max(x, y)) in also for x in ga for y in gb):
                    good += 1
                else:
                    bad += 1
        self.assertGreaterEqual(good, 55)
        self.assertEqual(bad, 0)
        xy = [r for r in self.rows if r["action"] == "pair_unequal"]
        self.assertTrue(any(r["tag"] == "HAP1/HAP1" and "H2_scaffold_4" in r["scaffold"] + r["partner"] for r in xy))    # X with Y

    def test_unlocs_join_the_curators_chromosome(self):
        un = [r for r in self.rows if r["action"] in ("unloc", "join")]
        self.assertGreaterEqual(len(un), 9)
        for r in un:
            g = set(self.groups[r["partner"].split(":")[0]])
            g |= {self.partner.get(x) for x in g}                       # either copy of the chromosome
            self.assertTrue(set(self.groups[r["scaffold"].split(":")[0]]) & g, r)

    def test_rearrangements(self):
        """Pieces of one chromosome copy are ordered and oriented by end contact (the other
        haplotype breaks ties). These are arrangements the curator made."""
        re = {r["scaffold"]: r["evidence"] for r in self.rows if r["action"] == "rearrange"}
        self.assertGreaterEqual(len(re), 8)
        # X: 25-58 Mb reversed, then 58-119 Mb
        self.assertRegex(re["H1_scaffold_1:2"], r"new order: H1_scaffold_1 25\.0-58\.4 Mb reversed, H1_scaffold_1 58\.4-")
        # the autosome half of H2_scaffold_4: its two pieces swapped (the curator inverted each in place)
        self.assertRegex(re["H2_scaffold_4:2"], r"new order: H2_scaffold_4 45\.\d-58\.4 Mb, H2_scaffold_4 33\.\d-45\.\d Mb;")
        # joins of whole scaffolds in front of their chromosome
        self.assertTrue(re["H1_scaffold_32"].startswith("new order: H1_scaffold_65 0.0-9.7 Mb, H1_scaffold_32"))
        self.assertTrue(re["H2_scaffold_31"].startswith("new order: H2_scaffold_61 0.0-12.1 Mb, H2_scaffold_31"))
        joins = {r["scaffold"] for r in self.rows if r["action"] == "join"}
        self.assertTrue({"H1_scaffold_65", "H2_scaffold_61", "H2_scaffold_68"} <= joins)
        ms = MapState(self.pt).replay(self.st.edits)
        self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(self.pt.npix)))
        # switched off, nothing is rearranged
        rows, _, _ = curate_allhaps(self.pt, phasing="hic", do_arrange=False)
        self.assertFalse(any(r["action"] in ("rearrange", "join") for r in rows))

    def test_pairs_of_one_chromosome_are_merged(self):
        """Two homolog pairs whose ends join in both haplotypes are one chromosome. Each merge
        here is one chromosome in the curation (H1_scaffold_2 was cut by the tool and the curator)."""
        merged = {frozenset((r["scaffold"].split(":")[0], r["partner"].split(":")[0])): r for r in self.rows if r["action"] == "merge"}
        self.assertEqual(set(merged), {frozenset(("H1_scaffold_2", "H2_scaffold_59")), frozenset(("H2_scaffold_64", "H1_scaffold_66")),
                                       frozenset(("H2_scaffold_62", "H1_scaffold_10"))})
        for r in merged.values():
            host = r["evidence"].split("the pair ")[1].split(" are parts")[0].split(" / ")
            ga = set(self.groups[r["scaffold"].split(":")[0]]) | set(self.groups[r["partner"].split(":")[0]])
            gh = set(self.groups[host[0].split(":")[0]]) | set(self.groups[host[1].split(":")[0]])
            self.assertTrue(ga & gh, r)
        # merged units are painted like the copy they joined
        ms = MapState(self.pt).replay(self.st.edits)
        paint = dict(self.st.scaffs)
        ids = collections.defaultdict(set)
        for i, c in enumerate(ms.contigs):
            if c.length >= 60:
                ids[self.pt.names[c.scaffold]].add(paint.get(i))
        self.assertEqual(ids["H2_scaffold_59"], ids["H2_scaffold_8"])
        self.assertTrue(ids["H2_scaffold_59"] and None not in ids["H2_scaffold_59"])

    def test_alignment_guides_orientation(self):
        """With alignment blocks, a piece's direction in the other copy comes from the sequence."""
        from pretext_autocurate.allhaps import _aligned_dir, _block_index
        from pretext_autocurate.homology import Block
        pt, bpp = self.pt, self.pt.bp_per_pix
        ix = {n: i for i, n in enumerate(pt.names)}
        rng = lambda n: (int(pt.scaffold_start[ix[n]]), int(pt.scaffold_start[ix[n]] + pt.scaffold_npix[ix[n]] - 1))
        a, b = rng("H2_scaffold_67"), rng("H1_scaffold_61")
        L = int((a[1] - a[0] + 1) * bpp)
        fwd = _block_index([Block("H2_scaffold_67", 0, L, "H1_scaffold_61", 1000000, 1000000 + L, True, L, 60)])
        rev = _block_index([Block("H2_scaffold_67", 0, L, "H1_scaffold_61", 1000000, 1000000 + L, False, L, 60)])
        self.assertEqual(_aligned_dir(pt, fwd, a, [b]), 1)
        self.assertEqual(_aligned_dir(pt, rev, a, [b]), -1)
        self.assertEqual(_aligned_dir(pt, rev, a, [rng("H1_scaffold_7")]), 0)          # aligned, but not to this copy
        self.assertEqual(_aligned_dir(pt, {}, a, [b]), 0)
        paf = os.path.join(ROOT, "seed", "fixtures", "qqTroRaco.H2_vs_H1.asm5.paf")
        if os.path.exists(paf):
            from pretext_autocurate.homology import read_paf
            rows, _, _ = curate_allhaps(self.pt, phasing="hic", blocks=read_paf(paf, min_mapq=0))
            re = {r["scaffold"]: r["evidence"] for r in rows if r["action"] == "rearrange"}
            self.assertIn("sequence alignment", re["H2_scaffold_67"])
            # the alignment leads for order and orientation: several copies are laid out along their homolog
            self.assertGreaterEqual(sum("order and orientation follow the other haplotype" in v for v in re.values()), 3)
            # as the curator has it: H2_scaffold_70 forward, then H2_scaffold_67 reversed (here the whole chain is turned)
            self.assertRegex(re["H2_scaffold_67"], r"new order: H2_scaffold_67 0\.0-9\.1 Mb, H2_scaffold_70 0\.0-4\.4 Mb reversed;")

    def _sites(self, rows):
        """(action, pixel) of every break site, and the curator's cut pixels."""
        ix = {n: i for i, n in enumerate(self.pt.names)}
        px = lambda r: int(self.pt.scaffold_start[ix[r["scaffold"]]] + round(r["length_mb"] * 1e6 / self.pt.bp_per_pix))
        return [(r["action"], px(r)) for r in rows if r["action"] in ("cut", "keep", "invert_tip")]

    def test_cuts_without_telomeres(self):
        """Contact-based sites (the telomere track is unusable here): every site is at a place
        the curator cut, at the default level and at the conservative end."""
        cur = SaveState.load(TRO + ".savestate_final")
        ms = MapState(self.pt).replay(cur.edits)
        truth = np.array(sorted(c.start for c in ms.contigs if self.pt.rel_coord[c.start] > 0)
                         + sorted(c.start + c.length - 1 for c in ms.contigs))
        orig = np.array([ms.pix[c.start] if not c.inverted else ms.pix[c.start + c.length - 1] for c in ms.contigs])
        orig = orig[self.pt.rel_coord[orig] > 0]                       # piece starts that are not scaffold starts
        near = lambda q: np.abs(orig - q).min() <= 12
        sites = self._sites(self.rows)
        cuts = [q for a, q in sites if a == "cut"]
        self.assertGreaterEqual(len(cuts), 11)
        self.assertTrue(all(near(q) for q in cuts))                    # no false cut
        self.assertTrue(all(near(q) for _, q in sites))                # flagged sites are real too
        self.assertGreaterEqual(sum(a == "invert_tip" for a, _ in sites), 3)
        # the level is one axis: conservative cuts are a subset, aggressive ones a superset
        lo = [q for a, q in self._sites(curate_allhaps(self.pt, phasing="hic", level=0.0)[0]) if a == "cut"]
        hi = [q for a, q in self._sites(curate_allhaps(self.pt, phasing="hic", level=1.0)[0]) if a == "cut"]
        self.assertTrue(8 <= len(lo) < len(cuts) < len(hi))
        self.assertGreaterEqual(len(set(lo) & set(cuts)), len(lo) - 1)
        self.assertLessEqual(sum(not near(q) for q in hi), 2)
        self.assertTrue(all(near(q) for q in lo))


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

    def test_small_scaffolds(self):
        """With the small-against-all alignment: unlocs where sequence and Hi-C agree (or Hi-C is
        strong on both maps), candidates only flagged, and the same-prefix pair H1_33/H1_40."""
        from pretext_autocurate.homology import read_paf
        sp = os.path.join(ROOT, "seed", "fixtures", "bAegFun_asm2_run1.small_vs_all.asm5.paf")
        if not os.path.exists(sp):
            self.skipTest("missing small-scaffold alignment")
        hi = Pretext(M0.replace("diploid_mq0", "diploid_mq10"), load_graphs=False)
        rows, st, _ = curate_allhaps(self.pt, phasing="hic", hi=hi, blocks=read_paf(PAF, min_mapq=0),
                                     small_blocks=read_paf(sp, min_mapq=0))
        un = {r["scaffold"]: r["partner"] for r in rows if r["action"] in ("unloc", "join", "false_duplicate")}
        joins = {r["scaffold"] for r in rows if r["action"] == "join"}
        # sequence already present in both copies is a false duplicate (curator, 2026-11-14): at the end
        # of its copy, tagged FalseDuplicate instead of Unloc
        fd = {r["scaffold"] for r in rows if r["action"] == "false_duplicate"}
        self.assertTrue({"H2_scaffold_70", "H1_scaffold_209", "H2_scaffold_43"} <= fd, fd)
        # H1_scaffold_44 is not one (curator, 2026-11-14): a 1.7 Mb repeat array aligns loosely to both
        # copies, but its 262 kb of non-repetitive sequence aligns to the start of H2_scaffold_40 only,
        # at its own end. It continues H2_scaffold_40 and goes into the H2 superscaffold (more
        # parsimonious than a segmental duplication in H1 against a deletion in H2).
        self.assertNotIn("H1_scaffold_44", fd)
        self.assertIn("H1_scaffold_44", joins)
        # H2_scaffold_34 (3.9 Mb) is large enough to be placed by end contact: joined in front of W
        re = {r["scaffold"]: r["evidence"] for r in rows if r["action"] == "rearrange"}
        self.assertTrue(re["H2_scaffold_8"].startswith("new order: H2_scaffold_34"))
        # H1_scaffold_168 (0.5 Mb, in no pair) aligns to H2_scaffold_5 where H1_scaffold_6 has a gap with
        # no contact across and no sequence: it goes in there. Scaffolds that align where their own
        # copy already has sequence (H2_scaffold_43, H2_scaffold_70) are duplicates and stay out.
        self.assertEqual(set(re), {"H2_scaffold_8", "H1_scaffold_6", "H2_scaffold_40"})
        self.assertTrue(re["H2_scaffold_40"].startswith("new order: H1_scaffold_44 0.0-2.0 Mb, H2_scaffold_40 0.0-2.6 Mb;"))
        # H1_scaffold_31: its homologs are both H1_ (H1_scaffold_33 / H1_scaffold_40); it goes to the
        # shorter copy, which is tagged HAP2, to balance the haplotypes (curator, 2026-11-14)
        self.assertEqual(un.get("H1_scaffold_31"), "H1_scaffold_40")
        why44 = [r["evidence"] for r in rows if r["scaffold"] == "H1_scaffold_44" and r["action"] == "join"][0]
        self.assertIn("non-repetitive sequence (262 kb", why44)
        self.assertIn("H1_scaffold_168 0.0-0.5 Mb", re["H1_scaffold_6"])
        self.assertEqual(un.get("H1_scaffold_168"), "H1_scaffold_6")
        self.assertEqual(un.get("H2_scaffold_70"), "H2_scaffold_30:2")
        # the alignment disagrees with two arrangements; both are flagged, neither is changed:
        # the inverted end of H1_scaffold_2 (the curator's instruction) and H2_scaffold_39 inside H2_scaffold_4
        self.assertEqual({r["scaffold"] for r in rows if r["action"] == "arrangement_differs"},
                         {"H1_scaffold_2", "H2_scaffold_4"})
        # the curator's own placements (savestate_6) that have evidence, with the same copy; a small
        # piece may go to the other prefix's copy (curator, 2026-11-13)
        for s, copy in (("H2_scaffold_34", "H2_scaffold_8"), ("H2_scaffold_61", "H2_scaffold_5"),
                        ("H2_scaffold_70", "H2_scaffold_30:2"), ("H1_scaffold_44", "H2_scaffold_40"),
                        ("H1_scaffold_209", "H2_scaffold_13")):
            self.assertEqual(un.get(s), copy, s)
        self.assertLessEqual(len(un), 20)
        t_rows, _, _ = curate_allhaps(self.pt, phasing="trio", hi=hi, blocks=read_paf(PAF, min_mapq=0),
                                      small_blocks=read_paf(sp, min_mapq=0), do_cuts=False)
        self.assertTrue(all(r["scaffold"][:2] == r["partner"][:2] for r in t_rows
                            if r["action"] == "unloc" and "scaffold_8" not in r["partner"]))      # trio: no mixing
        cand = {r["scaffold"]: r["evidence"] for r in rows if r["action"] == "unloc_candidate"}
        self.assertIn("H2_scaffold_53", cand)                          # sequence only: flagged, not moved
        self.assertIn("collapsed repeats", cand["H1_scaffold_62"])     # repeat-like coverage: not moved
        pairs = {frozenset((r["scaffold"], r["partner"])): r["evidence"] for r in rows if r["action"] == "pair"}
        self.assertIn("same haplotype file", pairs[frozenset(("H1_scaffold_33", "H1_scaffold_40"))])
        # layout: every unloc is painted like its copy, tagged Unloc, and sits right after it
        ms = MapState(self.pt).replay(st.edits)
        paint, flags = dict(st.scaffs), dict(st.meta_flags)
        name = lambda c: self.pt.names[c.scaffold]
        idx = {}
        for i, c in enumerate(ms.contigs):
            idx.setdefault(name(c), []).append(i)
        for s, copy in un.items():
            i = idx[s][0]
            if s in fd:
                self.assertIn("FalseDuplicate", st.tag_names(flags[i]), s)
                self.assertNotIn("Unloc", st.tag_names(flags[i]), s)
            elif s not in joins:
                self.assertIn("Unloc", st.tag_names(flags[i]), s)
            js = idx[copy.split(":")[0]]
            self.assertTrue(any(paint.get(j) == paint[i] for j in js), s)
            self.assertLessEqual(min(abs(i - j) for j in js), len(un), s)
        self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(self.pt.npix)))

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
