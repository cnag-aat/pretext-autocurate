"""Command line: info, dump, convert. (run and simulate come in later phases.)"""
import argparse
import collections
import sys

from . import __version__
from .pretextio import Pretext, check_same_layout
from .savestate import SaveState
from .layout import build_layout, invisible_scaffolds, write_agp, write_tsv


def cmd_info(a):
    pt = Pretext(a.map)
    print(f"file\t{pt.path}")
    print(f"hash\t{pt.hash.hex()}")
    print(f"total_bp\t{pt.total_len}")
    print(f"pixels\t{pt.npix}\t({pt.n_tex} x {pt.tex_res}, {pt.n_mips} mip levels)")
    print(f"bp_per_pixel\t{pt.bp_per_pix:.1f}")
    inv = invisible_scaffolds(pt)
    print(f"scaffolds\t{pt.n_scaffolds}\t({len(inv)} shorter than a pixel, {sum(l for _, l in inv)} bp)")
    for name, g in pt.graphs.items():
        print(f"track\t{name}\tmin {int(g.min())}\tmax {int(g.max())}\tnonzero px {int((g != 0).sum())}")
    if a.scaffolds:
        lens = pt.scaffold_len_bp()
        print("#scaffold\tlength_bp\tstart_px\tn_px")
        for i, name in enumerate(pt.names):
            print(f"{name}\t{lens[i]}\t{pt.scaffold_start[i]}\t{pt.scaffold_npix[i]}")


def cmd_dump(a):
    st = SaveState.load(a.state)
    pt = Pretext(a.map, load_graphs=False)
    ok = st.header_hash == pt.hash
    print(f"hash\t{st.header_hash.hex() if st.header_hash else None}\t{'matches map' if ok else 'DOES NOT MATCH MAP'}")
    kinds = collections.Counter(
        "break" if e.is_break else "invert" if e.delta == 0 else "move+invert" if e.invert else "move"
        for e in st.edits)
    print(f"edits\t{len(st.edits)}\t" + ", ".join(f"{k} {v}" for k, v in sorted(kinds.items())))
    print(f"waypoints\t{len(st.waypoints)}")
    print(f"painted_contigs\t{len(st.scaffs)}\tscaffold ids {len({s for _, s in st.scaffs})}")
    used = collections.Counter(t for _, m in st.meta_flags for t in st.tag_names(m))
    print("tags_used\t" + ", ".join(f"{k} {v}" for k, v in used.items()))
    if a.edits:
        print("#edit\tx\ty\tdelta\tinvert")
        for i, e in enumerate(st.edits):
            print(f"{i}\t{e.x}\t{e.y}\t{'break' if e.is_break else e.delta}\t{int(e.invert)}")
    for w in st.waypoints:
        print(f"waypoint\t{w.id + 1}\tpx {round((w.x + 0.5) * pt.npix)}")
    print("#contig\tout_scaffold\tpainted\tscaffold\tstart_bp\tend_bp\tstrand\tmap_px\tlen_px\ttags")
    for p in build_layout(pt, st):
        print("\t".join(map(str, p[:-1])) + "\t" + ",".join(p.tags))


def cmd_convert(a):
    pt = Pretext(a.map, load_graphs=False)
    st = SaveState.load(a.state)
    if st.header_hash != pt.hash:
        sys.exit(f"error: {a.state} was not saved from {a.map} (hash mismatch)")
    if a.to == "state":
        if not a.target_map:
            sys.exit("error: --to state needs --target-map (another map of the same assembly)")
        other = Pretext(a.target_map, load_graphs=False)
        try:
            check_same_layout([pt, other])
        except ValueError as e:
            sys.exit(f"error: {e}")
        st.save(a.out, other.path)
    else:
        pieces = build_layout(pt, st)
        (write_agp if a.to == "agp" else write_tsv)(pieces, a.out)
    print(f"wrote {a.out}", file=sys.stderr)


def cmd_haplotigs(a):
    import os
    from .haplotig import detect_haplotigs, write_calls, Params
    pt = Pretext(a.map)
    if "coverage" not in pt.graphs:
        print("warning: no coverage track in the map; calls rest on Hi-C alone and none will be tagged",
              file=sys.stderr)
    calls = detect_haplotigs(pt, Params(cov_max=a.cov_max, line_excess=a.line_excess))
    os.makedirs(a.out_dir, exist_ok=True)
    base = os.path.join(a.out_dir, os.path.basename(a.map))
    write_calls(calls, base + ".haplotigs.tsv")
    tiers = {"high": ("haplotig_high",), "medium": ("haplotig_high", "haplotig_medium")}[a.tag]
    st = SaveState.blank()
    bit = 1 << st.tag_id("Haplotig")
    index = {int(s): i for i, s in enumerate(sorted(pt.scaffold_start[pt.scaffold_npix > 0]))}
    st.meta_flags = sorted((index[c.start_px], bit) for c in calls if c.call in tiers)
    st.save(base + ".savestate_autocurate", pt.path)
    n = collections.Counter(c.call for c in calls)
    print("\t".join(f"{k} {v}" for k, v in sorted(n.items())))
    print(f"wrote {base}.haplotigs.tsv and {base}.savestate_autocurate ({len(st.meta_flags)} Haplotig tags)",
          file=sys.stderr)


def cmd_run(a):
    import os
    from .curate import curate_single, write_decisions
    maps = [Pretext(m, load_graphs=(i == 0)) for i, m in enumerate(a.map)]
    try:
        check_same_layout(maps)
    except ValueError as e:
        sys.exit(f"error: {e}")
    pt = maps[0]
    if a.mode == "allhaps":
        from .allhaps import curate_allhaps, write_rows
        blocks = None
        if a.paf or a.fasta:
            from .homology import align_haplotypes, read_paf, fasta_lengths
            paf = a.paf
            if not paf:
                missing = set(pt.names) - set(fasta_lengths(a.fasta))
                if missing:
                    sys.exit(f"error: {len(missing)} scaffolds of the map are not in {a.fasta} (e.g. {sorted(missing)[0]}); "
                             "give the assembly the map was made from")
                try:
                    paf = align_haplotypes(a.fasta, a.work_dir or os.path.join(a.out_dir, "work"), threads=a.threads,
                                           minimap2=a.minimap2)
                except RuntimeError as e:
                    sys.exit(f"error: {e}")
            blocks = read_paf(paf, min_mapq=0)      # pairing uses mapq >= 30; the rest shows repeats
            print(f"sequence homology: {len(blocks)} alignment blocks from {paf}", file=sys.stderr)
        hits = None
        if a.busco_table:
            from . import busco
            hits = busco.read_full_tables(a.busco_table)
            known = set(pt.names)
            if hits and not any(l[0] in known for loci in hits.values() for l in loci):
                sys.exit(f"error: no sequence name in {a.busco_table} matches the map; BUSCO must be run on the same assembly")
            print(f"BUSCO: {len(hits)} genes with complete copies from {a.busco_table}", file=sys.stderr)
        elif a.miniprot_paf or a.busco_lineage:
            from . import busco
            mp = a.miniprot_paf
            if not mp:
                if not a.fasta:
                    sys.exit("error: --busco-lineage needs --fasta (or give --miniprot-paf)")
                try:
                    lin = busco.resolve_lineage(a.busco_lineage, download=a.busco_download)
                    mp = busco.run_miniprot(a.fasta, busco.lineage_proteins(lin, a.busco_variants),
                                            os.path.join(a.work_dir or os.path.join(a.out_dir, "work"),
                                                         os.path.basename(lin.rstrip("/")) + ".miniprot.paf"),
                                            threads=a.threads, miniprot=a.miniprot)
                except RuntimeError as e:
                    sys.exit(f"error: {e}")
            hits = busco.read_hits(mp)
            print(f"BUSCO: complete hits for {len(hits)} genes from {mp}", file=sys.stderr)
        rows, st, info = curate_allhaps(pt, phasing=a.phasing, hi=maps[1] if len(maps) > 1 else None, blocks=blocks,
                                        busco_hits=hits)
        os.makedirs(a.out_dir, exist_ok=True)
        prefix = os.path.join(a.out_dir, a.prefix)
        write_rows(rows, prefix + ".decisions.tsv")
        write_agp(build_layout(pt, st), prefix + ".agp")
        for m in maps:
            st.save(os.path.join(a.out_dir, os.path.basename(m.path) + ".savestate_autocurate"), m.path)
        n = collections.Counter(r["action"] for r in rows)
        print("\t".join(f"{k} {v}" for k, v in sorted(n.items())))
        print(f"allhaps so far: telomere-supported cuts, pairing, orientation, tags (no joins, unlocs or haplotigs yet). "
              f"edits {len(st.edits)}, painted {len(st.scaffs)}, tags {len(st.meta_flags)}; outputs in {a.out_dir}",
              file=sys.stderr)
        return
    decisions, st = curate_single(pt, apply_inserts=a.apply_inserts)
    os.makedirs(a.out_dir, exist_ok=True)
    prefix = os.path.join(a.out_dir, a.prefix)
    write_decisions(decisions, prefix + ".decisions.tsv")
    pieces = build_layout(pt, st)
    write_agp(pieces, prefix + ".agp")
    for m in maps:
        st.save(os.path.join(a.out_dir, os.path.basename(m.path) + ".savestate_autocurate"), m.path)
    n = collections.Counter(d.action for d in decisions)
    print("\t".join(f"{k} {v}" for k, v in sorted(n.items())))
    print(f"edits {len(st.edits)}, painted contigs {len(st.scaffs)}, tags {len(st.meta_flags)}, "
          f"waypoints {len(st.waypoints)}; outputs in {a.out_dir}", file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(prog="pretext-autocurate")
    ap.add_argument("--version", action="version", version=__version__)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("info", help="header, tracks and savestate hash of a .pretext file")
    p.add_argument("map")
    p.add_argument("--scaffolds", action="store_true", help="also list every scaffold")
    p.set_defaults(fn=cmd_info)
    p = sub.add_parser("dump", help="human-readable savestate and the layout it produces")
    p.add_argument("state")
    p.add_argument("--map", required=True)
    p.add_argument("--edits", action="store_true", help="also list every edit")
    p.set_defaults(fn=cmd_dump)
    p = sub.add_parser("convert", help="savestate -> AGP, TSV, or a savestate for another map of the same assembly")
    p.add_argument("--map", required=True)
    p.add_argument("--state", required=True)
    p.add_argument("--to", required=True, choices=["agp", "tsv", "state"])
    p.add_argument("--target-map", help="for --to state: the map the new savestate is for")
    p.add_argument("--out", required=True)
    p.set_defaults(fn=cmd_convert)
    p = sub.add_parser("haplotigs", help="single mode: call haplotigs among small scaffolds; TSV + tagged savestate")
    p.add_argument("map")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--tag", choices=["high", "medium"], default="medium", help="lowest tier that gets a Haplotig tag")
    p.add_argument("--cov-max", type=float, default=0.8, help="coverage ratio below this counts as reduced")
    p.add_argument("--line-excess", type=float, default=0.05)
    p.set_defaults(fn=cmd_haplotigs)
    p = sub.add_parser("run", help="first curation pass: haplotigs + placement -> savestate, AGP, decisions")
    p.add_argument("--map", action="append", required=True, help="first map (MQ0) is analysed; repeat for other maps of the same assembly")
    p.add_argument("--mode", choices=["single", "allhaps"], default="single")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--prefix", default="autocurate")
    p.add_argument("--phasing", choices=["none", "hic", "trio"], default="hic",
                   help="allhaps only: how freely sequence may move between haplotypes. "
                        "trio: never; hic: whole chromosomes only, no H1_/H2_ mixing within one; none: freely")
    p.add_argument("--fasta", help="allhaps: the assembly the maps were made from (plain or .gz); "
                                   "haplotypes are aligned with minimap2 and homology is read from that")
    p.add_argument("--paf", help="allhaps: an existing alignment between haplotypes, used in place of --fasta")
    p.add_argument("--threads", type=int, default=8)
    p.add_argument("--work-dir", help="where alignments are kept and reused (default: OUT_DIR/work)")
    p.add_argument("--minimap2", help="path to minimap2 (default: PATH, then conda env pretext-autocurate)")
    p.add_argument("--busco-table", action="append",
                   help="allhaps: full_table.tsv of an existing BUSCO run (nothing is aligned; BUSCO evidence is optional). "
                        "Either one table made on the diploid assembly, or one per haplotype as PREFIX=path, repeated: "
                        "--busco-table H1=hap1.full_table.tsv --busco-table H2=hap2.full_table.tsv (the prefix is added "
                        "to the sequence names)")
    p.add_argument("--busco-lineage", help="allhaps: BUSCO lineage, as the path to an unpacked lineage directory "
                                           "or a name such as aves_odb12 (looked up in ~/.cache/pretext-autocurate/busco)")
    p.add_argument("--busco-download", action="store_true",
                   help="allow downloading a lineage given by name from busco-data.ezlab.org (file and size are printed)")
    p.add_argument("--busco-variants", action="store_true", help="align all protein variants (several times slower)")
    p.add_argument("--miniprot-paf", help="allhaps: an existing miniprot alignment of BUSCO proteins, used in place of running it")
    p.add_argument("--miniprot", help="path to miniprot (default: PATH, then conda env pretext-autocurate)")
    p.add_argument("--apply-inserts", action="store_true", help="insert localised pieces instead of only proposing them")
    p.set_defaults(fn=cmd_run)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
