"""Final assembly files from a curated layout, with pretext-to-asm (sanger-tol/agp-tpf-utils).

pretext-to-asm turns the input assembly FASTA and PretextView's AGP export into the curated
assembly (FASTA + AGP per haplotype, chromosome lists, unlocs, haplotigs, contaminants ...),
naming chromosomes by the Tree of Life rules. It is called as an external program rather
than reimplemented: its naming rules are maintained there, and it streams the sequence
through a .fai index (this tool never holds the sequence in memory). The AGP it is given is
written by `layout.write_agp`, which mirrors PretextView's export exactly.

Install: pip install git+https://github.com/sanger-tol/agp-tpf-utils
(or uv tool install 'git+https://github.com/sanger-tol/agp-tpf-utils').
"""
import collections
import os
import re
import shlex
import subprocess
import sys

from .homology import find_tool

INSTALL = "pip install git+https://github.com/sanger-tol/agp-tpf-utils"


def default_name(fasta):
    """Assembly name for the output files: the FASTA's base name up to the first dot
    (qqTroRaco.hfsm.asm.bp.diploid.fa -> qqTroRaco)."""
    base = os.path.basename(fasta)
    return re.split(r"\.", base, maxsplit=1)[0] or "assembly"


def run_pretext_to_asm(fasta, agp, out_dir, name=None, version=1, exe=None, extra=""):
    """Run pretext-to-asm on `fasta` and the PretextView-style `agp`; outputs go to
    `out_dir` as <name>.<version>.* (pretext-to-asm names them per haplotype and kind).
    `extra`: further pretext-to-asm options as one string (e.g. "--autosome-prefix SUPER_").
    Returns the command run. Raises RuntimeError with the install hint if the program is
    missing, or with its log tail if it fails. Note: pretext-to-asm writes a .fai (and an
    .agp) next to the FASTA if they are missing or older."""
    path = find_tool("pretext-to-asm", exe)
    if path is None:
        raise RuntimeError(f"pretext-to-asm not found (on PATH or in the pretext-autocurate conda env); install it with: {INSTALL}")
    os.makedirs(out_dir, exist_ok=True)
    out = os.path.join(out_dir, f"{name or default_name(fasta)}.{version}.fa")
    cmd = [path, "--assembly", fasta, "--pretext", agp, "--output", out] + shlex.split(extra or "")
    print("running: " + " ".join(shlex.quote(c) for c in cmd), file=sys.stderr)
    r = subprocess.run(cmd, stdin=subprocess.DEVNULL, capture_output=True, text=True)
    if r.returncode != 0:
        tail = "\n".join((r.stderr or r.stdout or "").strip().splitlines()[-15:])
        raise RuntimeError(f"pretext-to-asm failed (exit {r.returncode}):\n{tail}")
    if r.stderr.strip():
        print(r.stderr.strip(), file=sys.stderr)       # it lists the files it wrote on STDERR
    return cmd


def build(fasta, agp, out_dir, name=None, exe=None, extra="", gfa=(), yahs_agp=(), names=None, utg_gfa=None,
          fill_gaps=False):
    """The assembly files from a PretextView-style AGP, with the hifiasm graph when given
    (`gfa`: contig GFAs, `yahs_agp`: YaHS AGPs, both lists): bridge contigs are put into
    their gaps before pretext-to-asm, and gaps between graph-linked contig ends are closed in
    its output after (overlap checked in the sequence). With `utg_gfa` (unitig graph, with
    sequence) contigs that overlap along a single unitig path are joined too, and with
    `fill_gaps` a gap with a single unitig path is filled with the path's sequence between
    the two contigs (OUT_DIR/<name>.gapfills.fa holds those sequences; the AGPs name them).
    Every graph edit is listed in OUT_DIR/<name>.graph_edits.tsv. Returns that list."""
    from . import graph as G
    name = name or default_name(fasta)
    os.makedirs(out_dir, exist_ok=True)
    rows = []
    if gfa and yahs_agp:
        fa = G.Faidx(fasta)
        lay = G.Layout(yahs_agp, names=set(names) if names is not None else None)
        bad = lay.check({n: fa.length(n) for n in lay.scaffolds if n in (names or lay.scaffolds)} if not fa.gz else {})
        if bad:
            raise RuntimeError("the YaHS AGPs do not match the assembly: " + "; ".join(bad[:5]))
        gr = G.Graph(gfa)
        juncs = G.junctions(lay, gr)
        utg = None
        if utg_gfa:
            print(f"reading the unitig graph {utg_gfa} ...", file=sys.stderr)
            utg = G.UnitigGraph(utg_gfa)
            if fill_gaps and not utg.has_seq:
                raise RuntimeError("--fill-gaps needs the unitig GFA with sequence")
            if fill_gaps and not utg.haploid:
                raise RuntimeError("--fill-gaps needs hifiasm's read depth (rd:i: tags) in the unitig GFA")
            print(f"unitig graph: haploid read depth {utg.haploid:.0f} (collapsed sequence at about twice that)"
                  if utg.haploid else "unitig graph: no read depth", file=sys.stderr)
        if utg is not None and utg.has_seq:
            # contigs of the gap's own haplotype that lie wholly on a unitig path across it:
            # put into the gap like a bridge (a contig named for two gaps is left alone)
            done = {(j.scaffold, j.gap_start) for j in juncs}
            place = []
            for up in G.unitig_paths(lay, gr, utg):
                if (up.scaffold, up.gap_start) in done:
                    continue
                kind, val, note = G.choose_fill(fa, lay, gr, utg, up)
                if kind == "place":
                    c, oc = val
                    place.append(G.Junction("bridge", up.scaffold, up.gap_start, up.gap_end, up.left, up.left_o,
                                            up.right, up.right_o, 0, c, oc, 0, 0, lay.place[c].scaffold))
            named = collections.Counter(j.via for j in place)
            juncs = juncs + [j for j in place if named[j.via] == 1]
            print(f"unitig graph: {sum(named[j.via] == 1 for j in place)} contig(s) of the gap's own haplotype "
                  "to put into gaps", file=sys.stderr)
        bridged = os.path.join(out_dir, name + ".pretext.bridged.agp")
        for j, done, note in G.insert_bridges(agp, bridged, juncs, lay):
            # overlap 0: placed along a unitig path (its junctions are closed or filled below)
            what = "bridge" if j.overlap else "contig_placed"
            rows.append(dict(file=os.path.basename(agp), object=j.scaffold, left=f"{j.left}{j.left_o}",
                             right=f"{j.right}{j.right_o}",
                             graph_overlap=f"{j.overlap},{j.overlap2}" if j.overlap else "",
                             overlap="", action=f"{what}_{'inserted' if done else 'skipped'}: {note}"
                                                + ("" if j.overlap else f" (unitig path through {j.via}{j.via_o})")))
        if utg is not None and utg.has_seq:
            # the same for joins between two pieces of the layout (made by the tool or curator)
            for o, left, right, c, note in G.place_between(bridged, lay, gr, utg, fa):
                rows.append(dict(file=os.path.basename(agp), object=o, left=left, right=right, graph_overlap="",
                                 overlap="", action=f"contig_placed_inserted: {lay.place[c].scaffold} between "
                                                    f"{left} and {right} (unitig path through {c})"))
        agp = bridged
    run_pretext_to_asm(fasta, agp, out_dir, name, exe=exe, extra=extra)
    if gfa and yahs_agp:
        fills = {}
        for path in sorted(os.listdir(out_dir)):
            if not path.endswith(".agp") or ".pretext." in path:
                continue
            stem = os.path.join(out_dir, path[:-4])
            fa_out = next((stem + e for e in (".fa", ".fa.gz") if os.path.exists(stem + e)), None)
            if fa_out:
                G.close_gaps(os.path.join(out_dir, path), fa_out, fa, lay, gr, rows, utg=utg, fill=fill_gaps,
                             fills=fills)
        if fills:
            with open(os.path.join(out_dir, name + ".gapfills.fa"), "w") as f:
                for k, s in fills.items():
                    f.write(f">{k}\n" + "".join(s[i:i + 60] + "\n" for i in range(0, len(s), 60)))
        rep = os.path.join(out_dir, name + ".graph_edits.tsv")
        cols = ["file", "object", "left", "right", "graph_overlap", "overlap", "action"]
        with open(rep, "w") as f:
            f.write("\t".join(cols) + "\n")
            for r in rows:
                f.write("\t".join(str(r[c]) for c in cols) + "\n")
        n = sum(r["action"].startswith("gap_closed") for r in rows)
        nf = sum(r["action"].startswith("gap_filled") for r in rows)
        nb = sum(r["action"].startswith(("bridge_inserted", "contig_placed_inserted")) for r in rows)
        print(f"hifiasm graph: {nb} bridge contig(s) put into gaps, {n} gap(s) closed, {nf} filled; see {rep}",
              file=sys.stderr)
    return rows
