"""Do hifiasm's assembly graphs support the junctions between contigs that YaHS, the
curator and pretext-autocurate made?

Usage: python3 tools/graph_junctions.py MAP.pretext GFA_DIR CURATOR.savestate [TOOL.savestate ...]

GFA_DIR holds hifiasm's noseq GFAs (`*.hap1.p_ctg.noseq.gfa`, `*.hap2.p_ctg.noseq.gfa`,
`*.p_utg.noseq.gfa`, `*.r_utg.noseq.gfa`) and YaHS's AGPs (`*.hap1.*.agp`, `*.hap2.*.agp`;
scaffold_N becomes H1_scaffold_N / H2_scaffold_N). No sequence is needed: contigs are
placed on the unitig graphs through the reads they share (the GFA `A` lines).

A junction is two contig ends laid next to each other (contig A's exit, contig B's entry,
in the orientation the layout gives them). For each one:
  link     hifiasm's haplotype contig graph links the two ends directly (L line)
  comp     both ends are in one connected component of the unitig graph
  path     length (bp of unitigs walked, oriented) of the shortest path from A's exit to B's
           entry in the unitig graph, up to MAX_BP; 0 = the ends share a unitig
  tip      A's exit has no outgoing link (a dead end of the graph)
Junctions come from: YaHS (adjacent contigs within each scaffold of the AGP), the curator's
layout (consecutive pieces painted alike) and any tool savestate given. Exploratory, not
part of the package.
"""
import bisect
import collections
import glob
import heapq
import os
import sys

import numpy as np

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate import Pretext, SaveState, MapState          # noqa: E402

MAX_BP = 2_000_000
END_READS = 8


def gfa_parts(path, want_a=True):
    """S lengths, oriented L edges, and A lines {read: (segment, offset, strand)} or, per
    segment, [(offset, end, strand, read)]."""
    seg_len, links = {}, []
    by_read, by_seg = {}, collections.defaultdict(list)
    with open(path) as f:
        for line in f:
            t = line[0]
            if t == "S":
                c = line.split("\t")
                seg_len[c[1]] = int(c[3].split(":")[2])
            elif t == "L":
                c = line.split("\t")
                links.append((c[1], c[2], c[3], c[4]))
            elif t == "A" and want_a:
                c = line.split("\t", 7)
                off, rs, re_ = int(c[2]), int(c[5]), int(c[6])
                by_read[c[4]] = (c[1], off, c[3])
                by_seg[c[1]].append((off, off + re_ - rs, c[3], c[4]))
    return seg_len, links, by_read, by_seg


flip = {"+": "-", "-": "+"}


def main():
    pt = Pretext(sys.argv[1], load_graphs=False)
    d = sys.argv[2]
    states = sys.argv[3:]
    bpp = pt.bp_per_pix
    lens = dict(zip(pt.names, pt.scaffold_len_bp()))

    # --- YaHS AGPs: scaffold -> [(start, end, contig, orient)], contig -> scaffold placement
    agp, place = {}, {}
    for h in (1, 2):
        for path in glob.glob(os.path.join(d, f"*hap{h}*.agp")):
            for line in open(path):
                c = line.rstrip("\n").split("\t")
                if len(c) < 9 or c[4] != "W":
                    continue
                sc = f"H{h}_{c[0]}"
                agp.setdefault(sc, []).append((int(c[1]) - 1, int(c[2]), c[5], c[8]))
                place[c[5]] = (sc, int(c[1]) - 1, int(c[2]), c[8])
    bad = [s for s in agp if s in lens and abs(agp[s][-1][1] - lens[s]) > bpp]
    print(f"AGP: {sum(len(v) for v in agp.values())} contigs in {len(agp)} scaffolds; "
          f"{sum(s in lens for s in agp)} scaffolds in the map; length mismatches: {len(bad)}")
    starts = {s: [x[0] for x in v] for s, v in agp.items()}

    # --- haplotype contig graphs: ends and direct links
    ctg_reads, ctg_len, ctg_links = {}, {}, set()
    for h in (1, 2):
        path = glob.glob(os.path.join(d, f"*hap{h}.p_ctg.noseq.gfa"))[0]
        sl, links, _, by_seg = gfa_parts(path)
        ctg_len.update(sl)
        for s, v in by_seg.items():
            v.sort()
            ctg_reads[s] = v
        for a, oa, b, ob in links:
            ctg_links.add(((a, oa), (b, ob)))
            ctg_links.add(((b, flip[ob]), (a, flip[oa])))

    # --- unitig graph (p_utg: bubbles popped); oriented adjacency and components
    results = {}
    for gname in ("p_utg", "r_utg"):
        path = glob.glob(os.path.join(d, f"*.{gname}.noseq.gfa"))[0]
        ulen, ulinks, read_utg, _ = gfa_parts(path)
        adj = collections.defaultdict(list)
        parent = {u: u for u in ulen}

        def find(x, parent=parent):         # bound now: `parent` is rebound for the next graph
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for a, oa, b, ob in ulinks:
            adj[(a, oa)].append((b, ob))
            adj[(b, flip[ob])].append((a, flip[oa]))
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb
        comp_size = collections.Counter(find(u) for u in ulen)

        def end_node(ctg, side, read_utg=read_utg):
            """The unitig (and its orientation along the contig) at the contig's head or
            tail, from the outermost reads that are in a unitig."""
            v = ctg_reads.get(ctg, [])
            seq = v[:END_READS] if side == "h" else sorted(v, key=lambda x: -x[1])[:END_READS]
            for off, end, st, rd in seq:
                if rd in read_utg:
                    u, _, ust = read_utg[rd]
                    return (u, "+" if st == ust else "-")
            return None

        ends = {}
        for c in ctg_reads:
            ends[(c, "h")] = end_node(c, "h")
            ends[(c, "t")] = end_node(c, "t")
        results[gname] = (ulen, adj, find, comp_size, ends)

    def junction_info(A, oa, B, ob, gname):
        """A then B in the layout, with their orientations relative to hifiasm's contigs."""
        ulen, adj, find, comp_size, ends = results[gname]
        # exit of A: its tail if forward, its head (flipped) if reversed
        ex = ends.get((A, "t" if oa == "+" else "h"))
        en = ends.get((B, "h" if ob == "+" else "t"))
        if ex is None or en is None:
            return dict(comp=None, path=None, tip=None)
        ex = ex if oa == "+" else (ex[0], flip[ex[1]])
        en = en if ob == "+" else (en[0], flip[en[1]])
        comp = find(ex[0]) == find(en[0])
        tip = not adj.get(ex)
        # Dijkstra by bp entered
        dist = None
        if comp:
            if ex == en:
                dist = 0
            else:
                seen, pq = {ex: 0}, [(0, ex)]
                while pq:
                    dd, n = heapq.heappop(pq)
                    if n == en:
                        dist = dd
                        break
                    if dd > seen.get(n, 1e18) or dd > MAX_BP:
                        continue
                    for m in adj.get(n, ()):
                        nd = dd + (0 if m == en else ulen[m[0]])
                        if nd < seen.get(m, 1e18):
                            seen[m] = nd
                            heapq.heappush(pq, (nd, m))
        return dict(comp=comp, path=dist, tip=tip, csize=comp_size[find(ex[0])])

    # --- junctions of a layout: from a savestate (None = the original YaHS scaffolds)
    def contigs_of_piece(sc, lo_px, hi_px):
        """AGP contigs (with orientation) a map piece covers, snapping its ends to the
        nearest contig boundary within 1.5 px; ('cut inside', ...) when there is none."""
        v = agp.get(sc)
        if v is None:
            return None
        b0, b1 = lo_px * bpp, (hi_px + 1) * bpp
        bounds = [x[0] for x in v] + [v[-1][1]]
        snap = lambda b: min(bounds, key=lambda x: abs(x - b))
        s0, s1 = snap(b0), snap(b1)
        inside = (abs(s0 - b0) > 1.5 * bpp, abs(s1 - b1) > 1.5 * bpp)
        if inside[0]:
            s0 = v[bisect.bisect_right(starts[sc], b0) - 1][0]
        if inside[1]:
            s1 = v[bisect.bisect_right(starts[sc], b1 - 1) - 1][1]
        return [x for x in v if x[0] >= s0 and x[1] <= s1], inside

    def layout_junctions(state_path):
        out = []
        if state_path is None:
            for sc, v in agp.items():
                for x, y in zip(v, v[1:]):
                    out.append((x[2], x[3], y[2], y[3], sc, x[1] / 1e6, False))
            return out
        st = SaveState.load(state_path)
        ms = MapState(pt).replay(st.edits)
        paint = dict(st.scaffs)
        seq = []          # per piece: paint id, [(contig, orient in layout)], inside flags
        for i, c in enumerate(ms.contigs):
            sc = pt.names[c.scaffold]
            got = contigs_of_piece(sc, c.lo, c.hi)
            if got is None:
                seq.append((paint.get(i, 0), [], (False, False), sc, i))
                continue
            cs, inside = got
            lay = [(x[2], x[3]) for x in cs]
            if c.inverted:
                lay = [(n, flip[o]) for n, o in reversed(lay)]
                inside = inside[::-1]
            seq.append((paint.get(i, 0), lay, inside, sc, i))
        # inside a piece: its own contig junctions; between pieces: only if painted alike
        for p, lay, inside, sc, i in seq:
            for (a, oa), (b, ob) in zip(lay, lay[1:]):
                out.append((a, oa, b, ob, sc, place[a][2] / 1e6, False))
        for (p1, l1, in1, sc1, _), (p2, l2, in2, sc2, _) in zip(seq, seq[1:]):
            if p1 and p1 == p2 and l1 and l2:
                (a, oa), (b, ob) = l1[-1], l2[0]
                out.append((a, oa, b, ob, f"{sc1}|{sc2}", 0.0, in1[1] or in2[0]))
        return out

    def key(j):
        a, oa, b, ob = j[:4]
        x = ((a, "t" if oa == "+" else "h"), (b, "h" if ob == "+" else "t"))
        return tuple(sorted(x))

    yahs = layout_junctions(None)
    cur = layout_junctions(states[0])
    yk = {key(j) for j in yahs}
    ck = {key(j) for j in cur}

    def describe(j, gname="p_utg"):
        a, oa, b, ob = j[:4]
        info = junction_info(a, oa, b, ob, gname)
        info["link"] = ((a, oa), (b, ob)) in ctg_links
        return info

    def summary(label, js):
        js = [j for j in js if not j[6]]
        if not js:
            print(f"{label:46} n=0")
            return
        line = f"{label:46} n={len(js):4}"
        for gname in ("p_utg", "r_utg"):
            inf = [describe(j, gname) for j in js]
            link = sum(i["link"] for i in inf)
            comp = sum(bool(i["comp"]) for i in inf)
            path = sum(i["path"] is not None for i in inf)
            short = sum(i["path"] is not None and i["path"] <= 200_000 for i in inf)
            tip = sum(bool(i["tip"]) for i in inf)
            line += f" | {gname}: link {link} comp {comp} path {path} (<=200kb {short}) tip {tip}"
        print(line)

    print("\n== YaHS junctions (contig ends adjacent in the AGP), by the curator's verdict")
    summary("kept by the curator", [j for j in yahs if key(j) in ck])
    summary("broken by the curator", [j for j in yahs if key(j) not in ck])
    print("== curator junctions not made by YaHS")
    summary("new in the curation", [j for j in cur if key(j) not in yk])
    for sp in states[1:]:
        tl = layout_junctions(sp)
        tk = {key(j) for j in tl}
        print(f"== {os.path.basename(sp)}: junctions not made by YaHS")
        summary("in the curation", [j for j in tl if key(j) not in yk and key(j) in ck])
        summary("not in the curation", [j for j in tl if key(j) not in yk and key(j) not in ck])
        summary("YaHS junctions it broke, curator kept", [j for j in yahs if key(j) not in tk and key(j) in ck])
        summary("YaHS junctions it broke, curator broke", [j for j in yahs if key(j) not in tk and key(j) not in ck])

    # --- the graph as a proposer: from each contig end, the nearest contig end of the same
    # haplotype reachable in the unitig graph (entering it as a layout would), and whether
    # the curator has those two ends together
    ulen, adj, find, comp_size, ends = results["p_utg"]
    entry_of = collections.defaultdict(list)       # oriented node -> [(contig, its orientation)]
    for (c, side), n in ends.items():
        if n is None:
            continue
        if side == "h":
            entry_of[n].append((c, "+"))
        else:
            entry_of[(n[0], flip[n[1]])].append((c, "-"))
    hp = lambda c: c[:2]
    cur_ends = {}
    for k in ck:
        cur_ends.setdefault(k[0], set()).add(k[1])
        cur_ends.setdefault(k[1], set()).add(k[0])
    yahs_ends = {}
    for k in yk:
        yahs_ends.setdefault(k[0], set()).add(k[1])
        yahs_ends.setdefault(k[1], set()).add(k[0])
    small = lambda c: ctg_len.get(c, 0) < 1e6
    tally = collections.Counter()
    rows = []
    for (A, side), n in ends.items():
        if n is None:
            continue
        ex = n if side == "t" else (n[0], flip[n[1]])     # leaving A through this end
        seen, pq, found = {ex: 0}, [(0, ex)], None
        while pq and found is None:
            dd, x = heapq.heappop(pq)
            if dd > seen.get(x, 1e18) or dd > MAX_BP:
                continue
            for B, ob in entry_of.get(x, ()):
                if B != A and hp(B) == hp(A) and x != ex:
                    found = (dd, B, ob)
                    break
            for m in adj.get(x, ()):
                nd = dd + ulen[m[0]]
                if nd < seen.get(m, 1e18):
                    seen[m] = nd
                    heapq.heappush(pq, (nd, m))
        if found is None:
            tally["no contig end reached"] += 1
            continue
        dd, B, ob = found
        endB = (B, "h" if ob == "+" else "t")
        mine = (A, side)
        in_cur = endB in cur_ends.get(mine, ())
        in_yahs = endB in yahs_ends.get(mine, ())
        cls = ("YaHS and curator" if in_cur and in_yahs else "curator only" if in_cur
               else "YaHS only (curator broke)" if in_yahs else "neither")
        tally[cls] += 1
        tally[f"{cls}, <=200 kb"] += dd <= 200_000
        if min(ctg_len.get(A, 0), ctg_len.get(B, 0)) >= 1e6:
            tally[f"[both >= 1 Mb] {cls}"] += 1
        elif max(ctg_len.get(A, 0), ctg_len.get(B, 0)) >= 1e6:
            tally[f"[one >= 1 Mb] {cls}"] += 1
        if cls in ("curator only", "neither"):
            same_sc = place.get(A, ("",))[0] == place.get(B, ("",))[0]
            rows.append((cls, dd, A, side, B, endB[1], ctg_len.get(A, 0), ctg_len.get(B, 0), same_sc,
                         bool(cur_ends.get(mine)), bool(cur_ends.get(endB))))
    print("\n== graph as proposer: nearest same-haplotype contig end from each contig end (p_utg)")
    for k in sorted(tally):
        print(f"  {k:40} {tally[k]}")
    print("  'neither' and 'curator only' cases (kb; lengths Mb; same YaHS scaffold; each end used by the curator):")
    for r in sorted(rows, key=lambda r: (r[0], r[1])):
        print(f"    {r[0]:13} {r[1] / 1e3:7.0f}  {r[2]}:{r[3]} -> {r[4]}:{r[5]}  {r[6] / 1e6:6.2f} {r[7] / 1e6:6.2f}  "
              f"same_scaffold={r[8]} A_joined={r[9]} B_joined={r[10]}")

    # detail: every YaHS junction the curator broke and every new curator junction
    print("\n== detail (p_utg; path in kb)")
    extra = []
    for sp in states[1:]:
        tl = layout_junctions(sp)
        tk = {key(j) for j in tl}
        extra += [("tool_ok", [j for j in tl if key(j) not in yk and key(j) in ck and not j[6]]),
                  ("tool_bad", [j for j in tl if key(j) not in yk and key(j) not in ck and not j[6]]),
                  ("tool_cut_bad", [j for j in yahs if key(j) not in tk and key(j) in ck]),
                  ("tool_cut_ok", [j for j in yahs if key(j) not in tk and key(j) not in ck])]
    for label, js in [("broken", [j for j in yahs if key(j) not in ck]),
                      ("new", [j for j in cur if key(j) not in yk and not j[6]])] + extra:
        for j in js:
            i = describe(j)
            p = "-" if i["path"] is None else f"{i['path'] / 1e3:.0f}"
            print(f"{label:7} {j[4]:28} {j[5]:7.2f}  {j[0]}{j[1]} {j[2]}{j[3]}  link={int(i['link'])} "
                  f"comp={i['comp']} path={p} tip={i['tip']} csize={i.get('csize')}")


main()
