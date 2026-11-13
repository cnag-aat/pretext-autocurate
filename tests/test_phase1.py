"""Phase 1 tests: reader, hash, savestate, simulator, planner, layout.

Run with:  python3 -m unittest discover -s tests -v
Map-dependent tests are skipped when the fixture is missing.
"""
import glob
import os
import struct
import sys
import tempfile
import unittest

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
FX = os.path.join(ROOT, "seed", "fixtures")

from pretext_autocurate import (Pretext, SaveState, MapState, Edit, Waypoint, plan_edits,
                                header_hash, check_same_layout, BREAK_DELTA)
from pretext_autocurate.pretextio import inflate
from pretext_autocurate.layout import build_layout, write_agp
from pretext_autocurate import cli

BIRD0 = os.path.join(FX, "bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid_mq0.extensions.HR.pretext")
BIRD10 = os.path.join(FX, "bAegFun_asm2_run1.hfsm.asm.hic.yhs_mq10.diploid_mq10.extensions.HR.pretext")
TYL = os.path.join(FX, "xgTylRafi_asm4_run1.hfsm.asm.hic.hap2.yhs_mq10_mq0.extensions.HR.pretext")
BIRD_SS, TYL_SS = BIRD0 + ".savestate_6", TYL + ".savestate_4"
_cache = {}


def pretext(path):
    if not os.path.exists(path):
        raise unittest.SkipTest(f"missing fixture {os.path.basename(path)}")
    if path not in _cache:
        _cache[path] = Pretext(path)
    return _cache[path]


def raw_payload(path):
    with open(path, "rb") as f:
        d = f.read()
    off = 14 if d[5] == 3 else 6
    nc, nf = struct.unpack_from("<II", d, off)
    return inflate(d[off + 8:off + 8 + nc], nf)


class FakeMap:
    """Small synthetic map: scaffolds of the given pixel lengths."""
    def __init__(self, lengths):
        self.n_scaffolds = len(lengths)
        self.npix = sum(lengths)
        self.orig_scaffold = np.repeat(np.arange(len(lengths)), lengths)
        self.rel_coord = np.concatenate([np.arange(l) for l in lengths])


class TestSavestate(unittest.TestCase):
    def test_payload_roundtrip_all_fixtures(self):
        files = glob.glob(os.path.join(FX, "*savestate*"))
        self.assertGreaterEqual(len(files), 5)
        for f in files:
            with self.subTest(f=os.path.basename(f)):
                self.assertEqual(SaveState.load(f).payload(), raw_payload(f))

    def test_save_reload(self):
        st = SaveState.load(TYL_SS)
        with tempfile.TemporaryDirectory() as d:
            out = os.path.join(d, "x.savestate")
            st.save(out)
            st2 = SaveState.load(out)
        self.assertEqual(st2.payload(), st.payload())
        self.assertEqual(st2.header_hash, st.header_hash)

    def test_edit_count_multiple_of_8_is_read(self):
        st = SaveState.load(BIRD_SS)      # 264 edits: last flag bit lies past the flag array
        self.assertEqual(len(st.edits), 264)
        self.assertFalse(st.edits[-1].invert)

    def _template(self, edits, n_wp=0):
        st = SaveState.load(TYL_SS)
        st.edits = edits
        st.waypoints = [Waypoint(0.0, 0.0, 30.0, i) for i in range(n_wp)]
        st.scaffs, st.meta_flags = [], []
        return st

    def test_quirk_last_edit_inverted(self):
        edits = [Edit(10 * i, 10 * i + 5, 0, True) for i in range(7)] + [Edit(100, 120, 50, True)]
        st = self._template(edits)
        safe = st.safe_edits()
        self.assertEqual(len(safe), 9)
        st2 = SaveState()
        st2._parse(st.payload())
        self.assertEqual(st2.edits, safe)
        a = MapState(FakeMap([300])).replay(edits).pix
        b = MapState(FakeMap([300])).replay(st2.edits).pix
        self.assertTrue(np.array_equal(a, b))

    def test_quirk_odd_waypoints(self):
        edits = [Edit(10 * i, 10 * i + 5, 20, False) for i in range(8)]
        st = self._template(edits, n_wp=1)
        st2 = SaveState()
        st2._parse(st.payload())
        self.assertEqual(len(st2.edits), 9)
        self.assertEqual(len(st2.waypoints), 1)
        a = MapState(FakeMap([300])).replay(edits).pix
        b = MapState(FakeMap([300])).replay(st2.edits).pix
        self.assertTrue(np.array_equal(a, b))

    def test_quirk_not_triggered(self):
        edits = [Edit(10 * i, 10 * i + 5, 20, False) for i in range(8)]
        self.assertEqual(self._template(edits).safe_edits(), edits)


class TestHash(unittest.TestCase):
    def test_known_answers(self):
        for path, want in ((BIRD0, "a4fd387b13578ea7"), (TYL, "06712b247229d374")):
            if not os.path.exists(path):
                self.skipTest("missing fixture")
            self.assertEqual(header_hash(path).hex(), want)

    def test_matches_curator_savestates(self):
        for m, s in ((BIRD0, BIRD_SS), (TYL, TYL_SS)):
            if not os.path.exists(m):
                self.skipTest("missing fixture")
            self.assertEqual(SaveState.load(s).header_hash, header_hash(m))

    def test_cae71_hashes_recorded(self):
        g = os.path.join(FX, "CAE71_hap_merged_assembly_hic_map_mapq%s.pretext.savestate_%s")
        self.assertEqual(SaveState.load(g % ("0", "1")).header_hash.hex(), "3ed8d2ab9a753cec")
        self.assertEqual(SaveState.load(g % ("30", "claude_scaffold1")).header_hash.hex(), "2f42a0604797c6e9")


class TestReader(unittest.TestCase):
    def test_header(self):
        pt = pretext(TYL)
        self.assertEqual((pt.total_len, pt.n_scaffolds, pt.npix), (1551929974, 360, 65536))
        self.assertEqual(set(pt.graphs), {"gap", "telomere", "coverage"})
        self.assertEqual(int((pt.scaffold_npix > 0).sum()), 328)
        self.assertEqual(int(pt.scaffold_npix.sum()), pt.npix)

    def test_same_layout(self):
        check_same_layout([pretext(BIRD0), pretext(BIRD10)])
        with self.assertRaises(ValueError):
            check_same_layout([pretext(BIRD0), pretext(TYL)])

    def test_matrix_symmetric_with_diagonal(self):
        pt = pretext(TYL)
        m = pt.matrix(pt.n_mips - 1)
        self.assertEqual(m.shape, (pt.npix >> (pt.n_mips - 1),) * 2)
        self.assertTrue(np.array_equal(m, m.T))
        self.assertGreater(np.diag(m).mean(), 4 * m.mean())


class TestMapState(unittest.TestCase):
    def test_replay_curated(self):
        for m, s, n0, n1 in ((BIRD0, BIRD_SS, 469, 478), (TYL, TYL_SS, 328, 400)):
            with self.subTest(m=os.path.basename(m)[:8]):
                pt, st = pretext(m), SaveState.load(s)
                ms = MapState(pt)
                self.assertEqual(len(ms.contigs), n0)
                ms.replay(st.edits)
                self.assertEqual(len(ms.contigs), n1)
                self.assertTrue(np.array_equal(np.sort(ms.pix), np.arange(pt.npix)))
                self.assertLess(max(c for c, _ in st.scaffs), n1)
                self.assertLess(max(c for c, _ in st.meta_flags), n1)
                self.assertEqual(sum(c.length for c in ms.contigs), pt.npix)

    def test_matches_seed_simulator(self):
        """The packaged simulator agrees with the session-1 reference, edit by edit."""
        sys.path.insert(0, os.path.join(ROOT, "seed", "scripts"))
        import pretext as seed
        fm = FakeMap([40, 3, 25, 1, 60, 17])
        rng = np.random.default_rng(1)
        ref, new = seed.MapState.__new__(seed.MapState), MapState(fm)
        ref.n, ref.n_orig = fm.npix, fm.n_scaffolds
        ref.pix, ref.orig, ref.rel = np.arange(fm.npix), fm.orig_scaffold.copy(), fm.rel_coord.copy()
        ref.update_contigs()
        for _ in range(300):
            k = rng.integers(3)
            a = int(rng.integers(0, fm.npix - 1))
            b = int(rng.integers(a, min(fm.npix, a + 30)))
            if k == 0:
                e = (a, 0, BREAK_DELTA, False)
            else:
                d = int(rng.integers(-a, fm.npix - b)) if k == 1 else 0
                e = (a, b, d, bool(rng.integers(2)) or k == 2)
            ref.apply(e)
            new.apply(Edit(*e))
            self.assertTrue(np.array_equal(ref.pix, new.pix))
            self.assertEqual([tuple(c) for c in new.contigs], [tuple(c) for c in ref.contigs])

    def test_break(self):
        ms = MapState(FakeMap([10, 10]))
        ms.apply(Edit(4, 0, BREAK_DELTA))
        self.assertEqual([(c.start, c.length, c.scaffold) for c in ms.contigs], [(0, 5, 0), (5, 5, 0), (10, 10, 1)])


class TestSegmentedMap(unittest.TestCase):
    def test_split(self):
        from pretext_autocurate.pretextio import SegmentedMap
        fm = FakeMap([10, 5])
        fm.names = ["a", "b"]
        fm.bp_per_pix = 100.0
        v = SegmentedMap(fm, [3])
        self.assertEqual(v.names, ["a:1", "a:2", "b"])
        self.assertEqual(v.scaffold_start.tolist(), [0, 4, 10])
        self.assertEqual(v.scaffold_npix.tolist(), [4, 6, 5])
        self.assertEqual(v.parent_of("a:2"), ("a", 4))
        self.assertEqual(v.rel_coord[4], 0)
        self.assertEqual(v.was_cut.tolist(), [True, True, False])


class TestPlanner(unittest.TestCase):
    def test_random_targets(self):
        rng = np.random.default_rng(7)
        for _ in range(20):
            ms = MapState(FakeMap([50, 30, 20]))
            for _ in range(15):
                a = int(rng.integers(0, 99)); b = int(rng.integers(a, min(100, a + 25)))
                ms.apply(Edit(a, b, int(rng.integers(-a, 100 - b)), bool(rng.integers(2))), update=False)
            edits = plan_edits(ms.pix)
            got = MapState(FakeMap([50, 30, 20])).replay(edits)
            self.assertTrue(np.array_equal(got.pix, ms.pix))

    def test_reproduces_curation(self):
        pt, st = pretext(TYL), SaveState.load(TYL_SS)
        want = MapState(pt).replay(st.edits)
        edits = plan_edits(want.pix)
        got = MapState(pt).replay(edits)
        self.assertTrue(np.array_equal(got.pix, want.pix))
        self.assertEqual(got.contigs, want.contigs)
        self.assertLessEqual(len(edits), len(want.contigs))

    def test_continue_from_state(self):
        ms = MapState(FakeMap([100])).replay([Edit(10, 20, 30, True)])
        target = np.arange(100)
        edits = plan_edits(target, current=ms.pix)
        self.assertTrue(np.array_equal(ms.replay(edits).pix, target))


class TestLayout(unittest.TestCase):
    def test_layout_covers_every_visible_scaffold_once(self):
        pt, st = pretext(TYL), SaveState.load(TYL_SS)
        pieces = build_layout(pt, st)
        self.assertEqual(len(pieces), 400)
        lens = dict(zip(pt.names, pt.scaffold_len_bp().tolist()))
        by = {}
        for p in pieces:
            by.setdefault(p.name, []).append((p.start_bp, p.end_bp))
        self.assertEqual(len(by), 328)
        for name, iv in by.items():
            iv.sort()
            self.assertEqual(iv[0][0], 1)
            self.assertEqual(iv[-1][1], lens[name])
            for (a, b), (c, d) in zip(iv, iv[1:]):
                self.assertEqual(c, b + 1)
        self.assertEqual(len({p.painted_id for p in pieces if p.painted_id}), 18)
        self.assertEqual(sum("Haplotig" in p.tags for p in pieces), 132)

    def test_agp_and_cli(self):
        pretext(TYL); pretext(BIRD10)
        with tempfile.TemporaryDirectory() as d:
            agp = os.path.join(d, "t.agp")
            cli.main(["convert", "--map", TYL, "--state", TYL_SS, "--to", "agp", "--out", agp])
            with open(agp) as f:
                rows = [l.rstrip("\n").split("\t") for l in f if not l.startswith("#")]
            w = [r for r in rows if r[4] == "W"]
            self.assertEqual(len(w), 400)
            for r in rows:
                self.assertEqual(int(r[2]) - int(r[1]) + 1, int(r[5]) if r[4] == "U" else int(r[7]) - int(r[6]) + 1)
            # same edits, bound to the MQ10 map of the same assembly
            out = os.path.join(d, os.path.basename(BIRD10) + ".savestate_x")
            cli.main(["convert", "--map", BIRD0, "--state", BIRD_SS, "--to", "state",
                      "--target-map", BIRD10, "--out", out])
            st = SaveState.load(out)
            self.assertEqual(st.header_hash, header_hash(BIRD10))
            self.assertEqual(st.payload(), SaveState.load(BIRD_SS).payload())
            with self.assertRaises(SystemExit):
                cli.main(["convert", "--map", BIRD0, "--state", BIRD_SS, "--to", "state",
                          "--target-map", TYL, "--out", out])


if __name__ == "__main__":
    unittest.main()
