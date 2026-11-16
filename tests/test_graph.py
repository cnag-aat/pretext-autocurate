"""hifiasm contig graph + YaHS AGPs (graph.py): links and bridges at gaps, overlaps checked in
the sequence, bridges put into a PretextView-style AGP, gaps closed in pretext-to-asm output."""
import os
import random
import sys
import tempfile
import unittest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, "src"))
from pretext_autocurate import graph as G          # noqa: E402

FIX = os.path.join(ROOT, "seed", "fixtures")
GDIR = os.path.join(FIX, "qqTroRaco_hifiasm_gfas_yahs_agps")
TRO_FA = os.path.join(FIX, "qqTroRaco.hfsm.asm.bp.yhs_mq10.diploid.fa")


def _rand(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


class TestSynthetic(unittest.TestCase):
    """Three contig joins of a toy assembly:
      S1 = A, gap, B   with A+ -> B- linked by a 50 bp overlap (B stored reversed by hifiasm)
      S3 = D, gap, E   no link, but D -> C -> E with 40 bp overlaps, and C alone in S2."""

    @classmethod
    def setUpClass(cls):
        rng = random.Random(1)
        cls.d = tempfile.mkdtemp()
        A = _rand(300, rng)
        B = A[-50:] + _rand(250, rng)                 # as laid out in S1
        D = _rand(300, rng)
        C = D[-40:] + _rand(100, rng)
        E = C[-40:] + _rand(260, rng)
        cls.A, cls.B, cls.C, cls.D, cls.E = A, B, C, D, E
        N = "N" * 200
        cls.S = {"H1_scaffold_1": A + N + B, "H1_scaffold_2": C, "H1_scaffold_3": D + N + E}
        cls.fa = os.path.join(cls.d, "asm.fa")
        with open(cls.fa, "w") as f:
            for n, s in cls.S.items():
                f.write(f">{n}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
        cls.agp = os.path.join(cls.d, "x.hap1.yhs.agp")
        with open(cls.agp, "w") as f:
            f.write("scaffold_1\t1\t300\t1\tW\tctgA\t1\t300\t+\n"
                    "scaffold_1\t301\t500\t2\tN\t200\tscaffold\tyes\tproximity_ligation\n"
                    "scaffold_1\t501\t800\t3\tW\tctgB\t1\t300\t-\n"
                    "scaffold_2\t1\t140\t1\tW\tctgC\t1\t140\t+\n"
                    "scaffold_3\t1\t300\t1\tW\tctgD\t1\t300\t+\n"
                    "scaffold_3\t301\t500\t2\tN\t200\tscaffold\tyes\tproximity_ligation\n"
                    "scaffold_3\t501\t800\t3\tW\tctgE\t1\t300\t+\n")
        cls.gfa = os.path.join(cls.d, "x.hap1.p_ctg.noseq.gfa")
        with open(cls.gfa, "w") as f:
            for n, ln in (("ctgA", 300), ("ctgB", 300), ("ctgC", 140), ("ctgD", 300), ("ctgE", 300)):
                f.write(f"S\t{n}\t*\tLN:i:{ln}\n")
            # hifiasm writes each link from both sides
            f.write("L\tctgA\t+\tctgB\t-\t50M\nL\tctgB\t+\tctgA\t-\t50M\n"
                    "L\tctgD\t+\tctgC\t+\t40M\nL\tctgC\t-\tctgD\t-\t40M\n"
                    "L\tctgC\t+\tctgE\t+\t40M\nL\tctgE\t-\tctgC\t-\t40M\n")

    def test_junctions_and_overlaps(self):
        js, lay, gr = G.load([self.gfa], [self.agp], self.fa, names=list(self.S))
        by = {j.scaffold: (j, v) for j, v in js}
        self.assertEqual(set(by), {"H1_scaffold_1", "H1_scaffold_3"})
        j, v = by["H1_scaffold_1"]
        self.assertEqual((j.kind, j.left, j.left_o, j.right, j.right_o, j.overlap, v), ("link", "ctgA", "+", "ctgB", "-", 50, (50, 0)))
        self.assertEqual((j.gap_start, j.gap_end), (300, 500))
        j, v = by["H1_scaffold_3"]
        self.assertEqual((j.kind, j.via, j.via_scaffold, v), ("bridge", "ctgC", "H1_scaffold_2", (40, 40)))

    def test_bridge_insert_and_gap_closing(self):
        js, lay, gr = G.load([self.gfa], [self.agp], self.fa, names=list(self.S))
        # a PretextView-style AGP: S1 and S3 painted, S2 on its own
        pv = os.path.join(self.d, "pv.agp")
        with open(pv, "w") as f:
            f.write("# HiC MAP RESOLUTION: 100.0 bp/texel\n"
                    "Scaffold_1\t1\t800\t1\tW\tH1_scaffold_1\t1\t800\t+\tPainted\tHAP1\n"
                    "Scaffold_2\t1\t800\t1\tW\tH1_scaffold_3\t1\t800\t-\tPainted\tHAP1\n"
                    "Scaffold_3\t1\t140\t1\tW\tH1_scaffold_2\t1\t140\t+\tHAP1\n")
        out = os.path.join(self.d, "pv.bridged.agp")
        rep = G.insert_bridges(pv, out, [j for j, _ in js], lay)
        self.assertTrue(all(done for _, done, _ in rep))
        _, objs = G.read_agp(out)
        self.assertEqual(set(objs), {"Scaffold_1", "Scaffold_2"})
        # S3 reversed in the layout: E, C, D, all reversed, with the host's tags
        w = [(r.name, r.beg, r.end, r.orient, r.extra) for r in objs["Scaffold_2"] if r.kind == "W"]
        self.assertEqual(w, [("H1_scaffold_3", 501, 800, "-", ("Painted", "HAP1")),
                             ("H1_scaffold_2", 1, 140, "-", ("Painted", "HAP1")),
                             ("H1_scaffold_3", 1, 300, "-", ("Painted", "HAP1"))])
        # pretext-to-asm-like output: components split at gaps; one object reversed
        rc = G.revcomp
        o_agp, o_fa = os.path.join(self.d, "o.agp"), os.path.join(self.d, "o.fa")
        rows = {"SUPER_1": [("H1_scaffold_1", 1, 300, "+"), 200, ("H1_scaffold_1", 501, 800, "+")],
                "SUPER_2": [("H1_scaffold_1", 501, 800, "-"), 200, ("H1_scaffold_1", 1, 300, "-")],
                "SUPER_3": [("H1_scaffold_3", 1, 300, "+"), 200, ("H1_scaffold_2", 1, 140, "+"), 200,
                            ("H1_scaffold_3", 501, 800, "+")],
                "OTHER": [("H1_scaffold_2", 1, 140, "-")]}
        with open(o_agp, "w") as fa_, open(o_fa, "w") as ff:
            for obj, comps in rows.items():
                pos, seq = 1, []
                for k, c in enumerate(comps, 1):
                    if isinstance(c, int):
                        fa_.write(f"{obj}\t{pos}\t{pos + c - 1}\t{k}\tU\t{c}\tscaffold\tyes\tproximity_ligation\n")
                        seq.append("N" * c)
                        pos += c
                    else:
                        n, b, e, o = c
                        fa_.write(f"{obj}\t{pos}\t{pos + e - b}\t{k}\tW\t{n}\t{b}\t{e}\t{o}\n")
                        s = self.S[n][b - 1:e]
                        seq.append(s if o == "+" else rc(s))
                        pos += e - b + 1
                s = "".join(seq)
                ff.write(f">{obj}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
        report = []
        n = G.close_gaps(o_agp, o_fa, G.Faidx(self.fa), lay, gr, report)
        self.assertEqual(n, 4)
        got = {}
        with open(o_fa) as f:
            for line in f:
                if line.startswith(">"):
                    cur = line[1:].strip()
                    got[cur] = ""
                else:
                    got[cur] += line.strip()
        A, B, C, D, E = self.A, self.B, self.C, self.D, self.E
        self.assertEqual(got["SUPER_1"], A + B[50:])
        self.assertEqual(got["SUPER_2"], rc(B) + rc(A)[50:])
        self.assertEqual(got["SUPER_3"], D + C[40:] + E[40:])
        self.assertEqual(got["OTHER"], rc(C))                # untouched
        _, objs = G.read_agp(o_agp)
        self.assertEqual([r.kind for r in objs["SUPER_3"]], ["W", "W", "W"])
        self.assertEqual(int(objs["SUPER_1"][1].beg), 551)  # B lost its first 50 bp
        # a link whose overlap is not in the sequence is reported and left alone
        bad = os.path.join(self.d, "bad.gfa")
        with open(self.gfa) as f, open(bad, "w") as g:
            g.write(f.read().replace("50M", "45M"))
        js2, _, _ = G.load([bad], [self.agp], self.fa, names=list(self.S))
        link = [v for j, v in js2 if j.kind == "link"][0]
        self.assertEqual(link, (50, 0))                     # found within the slack: the real 50 bp


class TestUnitigFill(unittest.TestCase):
    """True sequence T = A + X + B; YaHS laid A, gap, B; the unitig graph has U1 (A's end and
    the start of X), U2 (inside X), U3 (end of X and B's start) linked with 100 bp overlaps.
    The contig ends reach the unitigs through shared reads (GFA A lines)."""

    def test_fill_from_unitig_path(self):
        rng = random.Random(2)
        d = tempfile.mkdtemp()
        A, X, B = _rand(3000, rng), _rand(500, rng), _rand(3000, rng)
        T = A + X + B
        U = {"U1": T[2000:3200], "U2": T[3100:3450], "U3": T[3350:5000]}
        S1 = A + "N" * 200 + B
        fa = os.path.join(d, "asm.fa")
        with open(fa, "w") as f:
            f.write(">H1_scaffold_1\n" + "".join(S1[i:i + 60] + "\n" for i in range(0, len(S1), 60)))
        agp = os.path.join(d, "x.hap1.agp")
        with open(agp, "w") as f:
            f.write("scaffold_1\t1\t3000\t1\tW\tctgA\t1\t3000\t+\n"
                    "scaffold_1\t3001\t3200\t2\tN\t200\tscaffold\tyes\tproximity_ligation\n"
                    "scaffold_1\t3201\t6200\t3\tW\tctgB\t1\t3000\t+\n")
        ctg = os.path.join(d, "x.hap1.p_ctg.noseq.gfa")
        with open(ctg, "w") as f:
            f.write("S\tctgA\t*\tLN:i:3000\nS\tctgB\t*\tLN:i:3000\n"
                    "A\tctgA\t0\t+\tr0\t0\t1000\nA\tctgA\t2000\t+\tr1\t0\t1000\n"
                    "A\tctgB\t0\t+\tr2\t0\t1000\nA\tctgB\t2000\t+\tr3\t0\t1000\n")
        utg = os.path.join(d, "x.p_utg.gfa")
        with open(utg, "w") as f:
            for n, s in U.items():          # U2 at twice the haploid depth: collapsed, so it may fill
                f.write(f"S\t{n}\t{s}\tLN:i:{len(s)}\trd:i:{52 if n == 'U2' else 26}\n")
            f.write("L\tU1\t+\tU2\t+\t100M\nL\tU2\t-\tU1\t-\t100M\nL\tU2\t+\tU3\t+\t100M\nL\tU3\t-\tU2\t-\t100M\n"
                    "A\tU1\t0\t+\tr1\t0\t1000\nA\tU2\t0\t+\tr9\t0\t300\nA\tU3\t150\t+\tr2\t0\t1000\n")
        js, lay, gr = G.load([ctg], [agp], fa)
        self.assertEqual(js, [])                          # no contig-graph link: only the unitigs join them
        ug = G.UnitigGraph(utg)
        ug.haploid = 26.0                   # too few long unitigs here to estimate it
        ups = G.unitig_paths(lay, gr, ug)
        self.assertEqual(len(ups), 1)
        up = ups[0]
        self.assertEqual(([u for u, _ in up.path], up.n_paths, up.bp, up.source), (["U1", "U2", "U3"], 1, 350, "no contig 100%"))
        self.assertEqual(G.path_fill(G.Faidx(fa), lay, ug, up), ("fill", X))
        # pretext-to-asm-like output, closed with and without --fill-gaps
        for fill in (False, True):
            o_agp, o_fa = os.path.join(d, f"o{fill}.agp"), os.path.join(d, f"o{fill}.fa")
            with open(o_agp, "w") as f:
                f.write("SUPER_1\t1\t3000\t1\tW\tH1_scaffold_1\t1\t3000\t+\n"
                        "SUPER_1\t3001\t3200\t2\tU\t200\tscaffold\tyes\tproximity_ligation\n"
                        "SUPER_1\t3201\t6200\t3\tW\tH1_scaffold_1\t3201\t6200\t+\n")
            with open(o_fa, "w") as f:
                f.write(">SUPER_1\n" + "".join(S1[i:i + 60] + "\n" for i in range(0, len(S1), 60)))
            rows, fills = [], {}
            n = G.close_gaps(o_agp, o_fa, G.Faidx(fa), lay, gr, rows, utg=ug, fill=fill, fills=fills)
            with open(o_fa) as f:
                got = "".join(l.strip() for l in f if not l.startswith(">"))
            if fill:
                self.assertEqual((n, got, fills), (1, T, {"gapfill_1": X}))
                self.assertTrue(rows[0]["action"].startswith("gap_filled: gapfill_1, 500 bp; collapsed"))
            else:
                self.assertEqual((n, got), (0, S1))
                self.assertTrue(rows[0]["action"].startswith("kept: 500 bp could fill it"))
        # the same path at haploid depth is one haplotype's sequence: never filled
        ug.depth["U2"] = 26
        rows = []
        o_agp, o_fa = os.path.join(d, "oTrue.agp"), os.path.join(d, "oTrue.fa")
        with open(o_agp, "w") as f:
            f.write("SUPER_1\t1\t3000\t1\tW\tH1_scaffold_1\t1\t3000\t+\n"
                    "SUPER_1\t3001\t3200\t2\tU\t200\tscaffold\tyes\tproximity_ligation\n"
                    "SUPER_1\t3201\t6200\t3\tW\tH1_scaffold_1\t3201\t6200\t+\n")
        with open(o_fa, "w") as f:
            f.write(">SUPER_1\n" + "".join(S1[i:i + 60] + "\n" for i in range(0, len(S1), 60)))
        n = G.close_gaps(o_agp, o_fa, G.Faidx(fa), lay, gr, rows, utg=ug, fill=True, fills={})
        self.assertEqual(n, 0)
        self.assertIn("haploid", rows[0]["action"])


class TestPlaceContig(unittest.TestCase):
    """T = A + X + B; YaHS laid A, gap, B in scaffold_1 and contig C (= X) alone in scaffold_2,
    same haplotype. The unitig path across the gap runs through C: C belongs in the gap
    (placed, not copied)."""

    def test_place(self):
        rng = random.Random(3)
        d = tempfile.mkdtemp()
        A, X, B = _rand(3000, rng), _rand(2500, rng), _rand(3000, rng)
        T = A + X + B
        U = {"U1": T[2000:3200], "U2": T[3100:5600], "U3": T[5500:7000]}
        S = {"H1_scaffold_1": A + "N" * 200 + B, "H1_scaffold_2": X}
        fa = os.path.join(d, "asm.fa")
        with open(fa, "w") as f:
            for n, s in S.items():
                f.write(f">{n}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
        agp = os.path.join(d, "x.hap1.agp")
        with open(agp, "w") as f:
            f.write("scaffold_1\t1\t3000\t1\tW\tctgA\t1\t3000\t+\n"
                    "scaffold_1\t3001\t3200\t2\tN\t200\tscaffold\tyes\tproximity_ligation\n"
                    "scaffold_1\t3201\t6200\t3\tW\tctgB\t1\t3000\t+\n"
                    "scaffold_2\t1\t2500\t1\tW\tctgC\t1\t2500\t+\n")
        ctg = os.path.join(d, "x.hap1.p_ctg.noseq.gfa")
        with open(ctg, "w") as f:
            f.write("S\tctgA\t*\tLN:i:3000\nS\tctgB\t*\tLN:i:3000\nS\tctgC\t*\tLN:i:2500\n"
                    "A\tctgA\t2000\t+\tr1\t0\t1000\nA\tctgB\t0\t+\tr2\t0\t1000\n"
                    "A\tctgC\t0\t+\tc1\t0\t1000\nA\tctgC\t1500\t+\tc2\t0\t1000\n")
        utg = os.path.join(d, "x.p_utg.gfa")
        with open(utg, "w") as f:
            for n, s in U.items():
                f.write(f"S\t{n}\t{s}\tLN:i:{len(s)}\trd:i:26\n")
            f.write("L\tU1\t+\tU2\t+\t100M\nL\tU2\t-\tU1\t-\t100M\nL\tU2\t+\tU3\t+\t100M\nL\tU3\t-\tU2\t-\t100M\n"
                    "A\tU1\t0\t+\tr1\t0\t1000\nA\tU2\t0\t+\tc1\t0\t1000\nA\tU2\t1500\t+\tc2\t0\t1000\n"
                    "A\tU3\t0\t+\tr2\t0\t1000\n")
        js, lay, gr = G.load([ctg], [agp], fa)
        ug = G.UnitigGraph(utg)
        ug.haploid = 26.0
        up, = G.unitig_paths(lay, gr, ug)
        kind, val, note = G.choose_fill(G.Faidx(fa), lay, gr, ug, up)
        self.assertEqual((kind, val), ("place", ("ctgC", "+")))
        self.assertIn("ctgC (H1_scaffold_2)", note)
        # a junction between two pieces of a PretextView-style AGP: C goes between them
        pv = os.path.join(d, "pv.agp")
        with open(pv, "w") as f:
            f.write("# HiC MAP RESOLUTION: 100.0 bp/texel\n"
                    "Scaffold_1\t1\t3000\t1\tW\tH1_scaffold_1\t1\t3000\t+\tPainted\tHAP1\n"
                    "Scaffold_1\t3001\t3100\t2\tU\t100\tscaffold\tyes\tproximity_ligation\n"
                    "Scaffold_1\t3101\t6100\t3\tW\tH1_scaffold_1\t3201\t6200\t+\tPainted\tHAP1\n"
                    "Scaffold_2\t1\t2500\t1\tW\tH1_scaffold_2\t1\t2500\t+\tHAP1\n")
        rep = G.place_between(pv, lay, gr, ug, G.Faidx(fa))
        self.assertEqual([r[3] for r in rep], ["ctgC"])
        _, objs = G.read_agp(pv)
        self.assertEqual(list(objs), ["Scaffold_1"])
        self.assertEqual([(r.name, r.beg, r.end, r.orient) for r in objs["Scaffold_1"] if r.kind == "W"],
                         [("H1_scaffold_1", 1, 3000, "+"), ("H1_scaffold_2", 1, 2500, "+"),
                          ("H1_scaffold_1", 3201, 6200, "+")])
        # C abuts both neighbours exactly along the path: closed with nothing added
        ups = G.unitig_paths(lay, gr, ug, pairs=[("x", "ctgA", "+", "ctgC", "+", 0, 0),
                                                 ("x", "ctgC", "+", "ctgB", "+", 0, 0)])
        self.assertEqual([G.choose_fill(G.Faidx(fa), lay, gr, ug, u)[:2] for u in ups], [("overlap", 0)] * 2)


@unittest.skipUnless(os.path.isdir(GDIR) and os.path.exists(TRO_FA), "missing fixture")
class TestTroglobisium(unittest.TestCase):
    def test_links_and_bridges(self):
        gfa = [os.path.join(GDIR, f"qqTroRaco.hfsm.asm.bp.hap{h}.p_ctg.noseq.gfa") for h in (1, 2)]
        agp = [os.path.join(GDIR, f"qqTroRaco.hfsm.asm.bp.hap{h}.yhs_mq10.agp") for h in (1, 2)]
        js, lay, _ = G.load(gfa, agp, TRO_FA)
        self.assertEqual(len(lay.scaffolds), 765)
        links = {(j.scaffold, j.overlap) for j, v in js if j.kind == "link" and v[0] == j.overlap}
        self.assertEqual(links, {("H1_scaffold_18", 17443), ("H1_scaffold_38", 13943), ("H1_scaffold_44", 13284),
                                 ("H2_scaffold_23", 17829)})
        bridges = {(j.scaffold, j.via_scaffold) for j, v in js if j.kind == "bridge" and all(v)}
        self.assertEqual(bridges, {("H1_scaffold_3", "H1_scaffold_373"), ("H1_scaffold_54", "H1_scaffold_229")})
        # unitig graph (structure only): about half the YaHS gaps are crossed by a unitig path,
        # including the end piece of H1_scaffold_16 the curator kept (the tool used to invert it)
        utg = os.path.join(GDIR, "qqTroRaco.hfsm.asm.bp.p_utg.noseq.gfa")
        if os.path.exists(utg):
            ups = G.unitig_paths(lay, G.Graph(gfa), G.UnitigGraph(utg))
            self.assertEqual(len(ups), 106)
            # one path only (alternatives searched up to 150 kb longer than the shortest)
            self.assertEqual(sum(u.n_paths == 1 for u in ups), 76)
            self.assertIn(("H1_scaffold_16", 450), {(u.scaffold, round(u.gap_start / 1e3, -1)) for u in ups})


if __name__ == "__main__":
    unittest.main()
