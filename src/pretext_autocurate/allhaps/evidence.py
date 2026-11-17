"""hifiasm graph evidence at break sites: contig-graph links and bridges
(`_graph_evidence`) and unitig paths (`_unitig_evidence`)."""


def _hold_mate(sites, x):
    """A site the hifiasm graph backs: its homologous site in the other haplotype (a kept one)
    is not used to rearrange either; where the two haplotypes then disagree, that is left as
    a difference between them for the curator (Troglobisium H1_scaffold_16 / H2_scaffold_17
    0.45 Mb: the end piece aligns in opposite orientations; the curator kept both)."""
    m = next((y for y in sites if y["q"] == x.get("mate")), None)
    if m is not None and m["action"] == "keep" and not m.get("graph") and not m.get("hold"):
        m["hold"] = True
        m["reason"] += ("; not used to rearrange: the homologous site's join is backed by hifiasm's graph, so a "
                        "disagreement here is a difference between the haplotypes")


def _unitig_evidence(raw, sites, unitig_paths):
    """hifiasm unitig-graph paths (`graph.unitig_paths`) between the two contig ends of YaHS
    gaps: weaker than a contig-graph link (on Troglobisium 107 of the 197 gaps the curator
    kept have one, and 1 of the 24 they broke: a haplotype switch where the alignment to the
    other haplotype breaks). So a path only stops a gap from being used to rearrange: a kept
    site there is no split point, and a tip inversion or a cut on contact alone becomes a
    keep; a cut backed by the alignment or a telomere stands. Returns rows: graph_keep where
    a decision changed, graph_path otherwise."""
    idx = {n: i for i, n in enumerate(raw.names)}
    bpp = raw.bp_per_pix
    rows = []
    for up in unitig_paths:
        s = idx.get(up.scaffold)
        if s is None:
            continue
        q0 = int(raw.scaffold_start[s] + up.gap_start // bpp)
        q1 = int(raw.scaffold_start[s] + max(up.gap_end - 1, up.gap_start) // bpp)
        undone, stands = [], []
        for x in sites:
            if x["s"] != s or not q0 - 2 <= x["q"] <= q1 + 2 or x.get("graph"):
                continue
            strong = "alignment" in x["reason"] or "telomere" in x["reason"]
            if x["action"] == "cut" and strong:
                stands.append("cut")
                continue
            if x["action"] != "keep":
                undone.append(x["action"].replace("_", " "))
            x["action"], x["graph"] = "keep", True
            _hold_mate(sites, x)
        what = (f"a unitig path of hifiasm's graph joins the two contig ends YaHS put here ({up.left}{up.left_o} .. "
                f"{up.right}{up.right_o}: {len(up.path)} unitigs, {up.bp / 1e3:.0f} kb between them, "
                + ("one path" if up.n_paths == 1 else "with alternatives (bubbles)")
                + (f"; reads of the inner unitigs: {up.source}" if up.source else "") + ")")
        if undone:
            what += f": the {', '.join(sorted(set(undone)))} proposed here is not made"
        elif stands:
            what += ": the cut stands (the alignment or a telomere says so; e.g. a haplotype switch)"
        else:
            what += ": not used as a point to rearrange at"
        rows.append(dict(action="graph_keep" if undone else "graph_path", scaffold=up.scaffold,
                         length_mb=round(up.gap_start / 1e6, 2), partner="", strand="", tag="", score=up.bp,
                         next_best=up.n_paths, evidence=what))
    return rows


def _graph_evidence(raw, sites, graph_junctions):
    """hifiasm graph evidence (`graph.load`) at YaHS gaps. A gap whose two contig ends the
    contig graph links, or bridges through one contig that is a scaffold of its own, is a
    right join: a site there is turned into a keep marked `graph` (never cut, tip-inverted or
    used as a split point). Returns one row per junction: graph_link (asm closes the gap)
    or graph_bridge (asm puts the bridging contig into the gap)."""
    idx = {n: i for i, n in enumerate(raw.names)}
    bpp = raw.bp_per_pix
    rows = []
    for j, ver in graph_junctions:
        s = idx.get(j.scaffold)
        if s is None:
            continue
        q0 = int(raw.scaffold_start[s] + j.gap_start // bpp)
        q1 = int(raw.scaffold_start[s] + max(j.gap_end - 1, j.gap_start) // bpp)
        undone = []
        for x in sites:
            if x["s"] == s and q0 - 2 <= x["q"] <= q1 + 2:
                if x["action"] != "keep":
                    undone.append(x["action"].replace("_", " "))
                x["action"], x["graph"] = "keep", True
                _hold_mate(sites, x)
        if ver is None:
            seq = "not checked in the sequence (no --fasta)"
        elif j.kind == "link":
            seq = "the same bases in the assembly" if ver[0] else "NOT found in the assembly sequence"
        else:
            seq = "both the same bases in the assembly" if ver[0] and ver[1] else "NOT both found in the assembly sequence"
        if j.kind == "link":
            what = (f"hifiasm's contig graph links {j.left}{j.left_o} to {j.right}{j.right_o}, the two contigs YaHS "
                    f"joined here (overlap {j.overlap:,} bp, {seq}): the join is right; asm closes the gap")
        else:
            what = (f"the only graph link out of {j.left}{j.left_o} and the only one into {j.right}{j.right_o} both go "
                    f"to {j.via}{j.via_o} ({j.via_scaffold}, a scaffold of its own; overlaps {j.overlap:,} and "
                    f"{j.overlap2:,} bp, {seq}): it "
                    "fills this gap; asm puts it there and closes both sides")
        if j.alternatives and j.kind == "link":
            what += f" ({j.alternatives} other link(s) at these ends, e.g. a bubble branch; YaHS chose the linked one)"
        elif j.alternatives:
            what += (f" (the bridging contig has {j.alternatives} other link(s), e.g. to a bubble branch; the gap's two "
                     "ends have this one only)")
        if undone:
            what += f"; the {', '.join(sorted(set(undone)))} proposed here is not made"
        rows.append(dict(action="graph_" + j.kind, scaffold=j.scaffold, length_mb=round(j.gap_start / 1e6, 2),
                         partner=j.via_scaffold, strand="", tag="", score=j.overlap, next_best=j.overlap2 or "",
                         evidence=what))
    return rows
