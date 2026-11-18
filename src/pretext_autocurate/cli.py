"""Command line: run (the pre-curation pass), info, dump, convert, asm, haplotigs."""
import argparse
import collections
import os
import sys

from . import __version__, asm
from .pretextio import Pretext, check_same_layout
from .savestate import SaveState
from .layout import build_layout, invisible_scaffolds, unplaced_hap_tags, write_agp, write_tsv


def cmd_info(a):
    """`info`: header, savestate hash, tracks (and with --scaffolds each scaffold's bp and px)."""
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
    """`dump`: a savestate in readable form (edits, tags, waypoints) and the layout it gives."""
    st = SaveState.load(a.state)
    pt = Pretext(a.map, load_graphs=False)
    ok = st.header_hash == pt.hash              # PretextView loads a savestate only on its own map
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
    for w in st.waypoints:                      # x is stored in screen units, -0.5..0.5 over the map
        print(f"waypoint\t{w.id + 1}\tpx {round((w.x + 0.5) * pt.npix)}")
    print("#contig\tout_scaffold\tpainted\tscaffold\tstart_bp\tend_bp\tstrand\tmap_px\tlen_px\ttags")
    for p in build_layout(pt, st):
        print("\t".join(map(str, p[:-1])) + "\t" + ",".join(p.tags))


def cmd_convert(a):
    """`convert`: a savestate to AGP or TSV, or re-hashed for another map of the same assembly
    (same scaffolds and pixels, e.g. the high-MQ map), so PretextView loads it there."""
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
        if a.to == "agp":
            write_agp(pt, st, a.out, unplaced=unplaced_hap_tags(pt, st))
        else:
            write_tsv(build_layout(pt, st), a.out)
    print(f"wrote {a.out}", file=sys.stderr)


def cmd_haplotigs(a):
    """`haplotigs`: single-mode haplotig calls among small scaffolds; writes a TSV and a
    savestate with Haplotig tags for calls of tier `--tag` and above."""
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
    # Tags are keyed by contig index in map order, counting only scaffolds with at least a pixel.
    index = {int(s): i for i, s in enumerate(sorted(pt.scaffold_start[pt.scaffold_npix > 0]))}
    st.meta_flags = sorted((index[c.start_px], bit) for c in calls if c.call in tiers)
    st.save(base + ".savestate_autocurate", pt.path)
    n = collections.Counter(c.call for c in calls)
    print("\t".join(f"{k} {v}" for k, v in sorted(n.items())))
    print(f"wrote {base}.haplotigs.tsv and {base}.savestate_autocurate ({len(st.meta_flags)} Haplotig tags)",
          file=sys.stderr)


# Decisions a curator should look at in allhaps mode, in review order (dict order is the sort
# key), with why each is listed: everything else was applied with clear evidence.
REVIEW_ALLHAPS = {"rephase": "a copy of a same-haplotype scaffold, moved to the other haplotype (all-vs-all PAF)",
                  "keep": "site flagged, nothing changed","arrangement_differs": "the two haplotypes disagree",
                  "graph_bridge": "small contig that fills a gap (hifiasm graph): put there by asm",
                  "graph_keep": "unitig path across the gap: rearrangement there dropped",
                  "graph_join_candidate": "scaffold ends only the hifiasm graph joins: not made",
                  "merge_candidate": "chromosome-level contact says one chromosome: not merged at this level",
                  "self_duplication": "two copies of one region in a scaffold (cut only if haplotypic, at a gap)",
                  "sex_boundary": "sex chromosome meets another chromosome, not cut",
                  "unloc_candidate": "chromosome suggested, not moved",
                  "insert_candidate": "missing piece of a chromosome, belongs inside it: not moved",
                  "false_duplicate_candidate": "BUSCO duplicate, not moved",
                  "sex_chromosome_candidate": "hemizygous chromosome", "pair_unequal": "possible sex chromosome pair",
                  "sex_chromosome_join": "hemizygous piece placed by contact", "par_candidate": "possible PAR, not confirmed",
                  "par": "pseudoautosomal region: check its extent", "sex_linked_region": "coverage and alignment of a sex chromosome",
                  "sex_summary": "sex chromosomes found (or why none)",
                  "place_next_to_homolog": "position along the homolog not resolved", "merge": "two pairs laid out as one chromosome",
                  "false_duplicate": "tagged FalseDuplicate", "unpaired": "no homolog found"}


def _write_review(rows, path, mode):
    """The decisions that need a curator's eye, most important first, one line each.

    `rows` are allhaps row dicts or single-mode Decision objects, by `mode`. Returns the number
    of lines written."""
    if mode == "allhaps":
        # Unpaired scaffolds under 1 Mb (or with no length) are too many and too small to list.
        sel = [(list(REVIEW_ALLHAPS).index(r["action"]), r) for r in rows if r["action"] in REVIEW_ALLHAPS
               and not (r["action"] == "unpaired" and (r["length_mb"] == "" or float(r["length_mb"]) < 1.0))]
        sel.sort(key=lambda x: x[0])
        with open(path, "w") as f:
            f.write("action\twhy_listed\tscaffold\tposition_or_length_mb\tpartner\tevidence\n")
            for _, r in sel:
                f.write(f"{r['action']}\t{REVIEW_ALLHAPS[r['action']]}\t{r['scaffold']}\t{r['length_mb']}\t{r['partner']}\t"
                        f"{r['evidence']}\n")
    else:
        sel = [d for d in rows if d.confidence in ("low", "medium") or d.action.startswith("flag")]
        sel.sort(key=lambda d: (d.confidence != "low", d.scaffold))     # low confidence first
        with open(path, "w") as f:
            f.write("action\tconfidence\tscaffold\tlength_bp\ttarget\tevidence\n")
            for d in sel:
                f.write(f"{d.action}\t{d.confidence}\t{d.scaffold}\t{d.length_bp}\t{d.target}\t{d.evidence}\n")
    return len(sel)


def _write_run_info(a, path, mode, extra=""):
    """Version, date, command line, mode and every option, so a run can be repeated; `extra`
    is appended as is (tab-separated lines, e.g. the level of a sweep)."""
    import datetime
    import platform
    with open(path, "w") as f:
        f.write(f"pretext-autocurate {__version__}\n")
        f.write(f"date\t{datetime.datetime.now().isoformat(timespec='seconds')}\n")
        f.write(f"python\t{platform.python_version()}\n")
        f.write(f"command\t{' '.join(sys.argv)}\n")
        f.write(f"mode\t{mode}\n")
        for k, v in sorted(vars(a).items()):
            if k != "fn":                       # the subcommand's function, not an option
                f.write(f"option\t{k}\t{v}\n")
        f.write(extra)


def _detect_mode(pt, hap_regex):
    """allhaps if at least two haplotype prefixes each hold a fifth of the sequence.

    Returns (mode, sorted prefixes found). The fifth keeps a few prefixed unplaced scaffolds
    (e.g. a lone H2_ contig in a primary assembly) from switching the mode."""
    import re
    lens = pt.scaffold_len_bp()
    by = collections.Counter()
    for n, l in zip(pt.names, lens):
        m = re.match(hap_regex, n)
        if m:
            by[m.group(1)] += int(l)
    big = [h for h, l in by.items() if l >= 0.2 * sum(lens)]
    return ("allhaps" if len(big) >= 2 else "single"), sorted(by)


def _load_tracks(a, pt):
    """Telomere strand and repeat tracks: from the map, from files, or computed from --fasta.

    Puts them in pt.graphs (see tracks.load) and keeps the files in `a.track_files` for
    --add-tracks. Exits on a tool error."""
    import os
    from . import tracks
    files = {"telomere_5p": a.telomere_5p, "telomere_3p": a.telomere_3p, "repeat": a.repeats}
    # Given files win; the tools run only for what is missing (and their outputs are cached).
    work = a.work_dir or os.path.join(a.out_dir, "work")
    try:
        if a.telomere_motif and not (a.telomere_5p and a.telomere_3p):
            if not a.fasta:
                sys.exit("error: --telomere-motif needs --fasta (or give --telomere-5p/--telomere-3p)")
            print(f"telomeres: running telo_scan.py for {a.telomere_motif} (cached in {work})", file=sys.stderr)
            files.update(tracks.run_telo_scan(a.fasta, work, a.telomere_motif, a.telo_scan, a.threads))
        if a.red and not a.repeats:
            if not a.fasta:
                sys.exit("error: --red needs --fasta (or give --repeats BED)")
            print(f"repeats: running Red on each haplotype (cached in {work})", file=sys.stderr)
            files["repeat"] = tracks.run_red(a.fasta, work, a.hap_regex)
    except (RuntimeError, OSError) as e:
        sys.exit(f"error: {e}")
    except Exception as e:                                  # an external tool failed: say which
        sys.exit(f"error: {e}")
    a.track_files = {k: v for k, v in files.items() if v}
    got = tracks.load(pt, files)
    for k, src in sorted(got.items()):
        print(f"track {k}: from {src}", file=sys.stderr)
    for k, path in files.items():
        if path and k not in got:
            print(f"warning: no sequence name in {path} matches the map; track {k} not used", file=sys.stderr)
    if any(src != "map" for src in got.values()) and not a.add_tracks:
        # Tracks from files are used by the rules but not visible in PretextView: say how to add them.
        names ={"telomere_5p": "5p_telomere", "telomere_3p": "3p_telomere", "repeat": "repeat_density"}
        print("to see these tracks in PretextView, add them to a copy of the map with PretextGraph, e.g.\n  PretextGraph -i "
              + a.map[0] + " -o with_tracks.pretext "
              + " ".join(f"-t{i + 1} {p} -n{i + 1} {names[k]}" for i, (k, p) in enumerate(
                  (k, p) for k, p in files.items() if p and k in got)), file=sys.stderr)


def _with_tracks(a, hits):
    """--add-tracks: a copy of every map with the evidence tracks, for PretextView. Returns the
    paths written; the caller saves a savestate for each of them too (the savestate hash covers
    the file name, so the original map's savestate does not load on the copy).

    `hits` are BUSCO gene copies (None: no BUSCO tracks). A failure on one map is a warning."""
    import os
    import re
    from . import tracks
    if not a.add_tracks:
        return []
    # Graph names as telo_scan/PretextGraph users know them; tracks.from_map reads them back.
    names = {"telomere_5p": "5p_telomere", "telomere_3p": "3p_telomere", "repeat": "repeat_density"}
    files = {names[k]: v for k, v in getattr(a, "track_files", {}).items() if k in names}
    if hits:
        hap_of = lambda n: (re.match(a.hap_regex, n).group(1) if re.match(a.hap_regex, n) else "")
        files.update(tracks.busco_tracks(hits, hap_of, os.path.join(a.out_dir, "tracks")))
    files.update(getattr(a, "sex_track_files", {}))          # sex_chromosome, par (sexchrom.write_tracks)
    if not files:
        print("warning: --add-tracks: no tracks to add (give --repeats/--red, --telomere-*, or BUSCO)", file=sys.stderr)
        return []
    made = []
    for m in a.map:                      # every map: on a high-MQ map a repeat track tells repeat from contaminant
        base = os.path.basename(m)
        out = os.path.join(a.out_dir, (base[:-len(".pretext")] if base.endswith(".pretext") else base) + ".tracks.pretext")
        try:
            tracks.add_tracks(m, files, out)
        except (RuntimeError, OSError) as e:
            print(f"warning: --add-tracks failed for {m}: {e}", file=sys.stderr)
            continue
        except Exception as e:
            print(f"warning: --add-tracks failed for {m}: {e} (see {out}.pretextgraph.log)", file=sys.stderr)
            continue
        print(f"map with tracks ({', '.join(files)}): {out}", file=sys.stderr)
        made.append(out)
    return made


def cmd_run(a):
    """`run`: the pre-curation pass. Analyses the first (MQ0) map in single or allhaps mode;
    writes per output directory (one per level with --sweep) decisions.tsv, review.tsv, an AGP,
    run.txt, and a savestate for every map given (and every --add-tracks copy)."""
    from .curate import curate_single, write_decisions
    if bool(a.gfa) != bool(a.yahs_agp):
        sys.exit("error: --gfa and --yahs-agp go together (the AGPs say where the graph's contigs are)")
    if a.pretext_to_asm:
        if not a.fasta:
            sys.exit("error: --pretext-to-asm needs --fasta (the assembly the maps were made from)")
        if asm.find_tool("pretext-to-asm", a.pretext_to_asm_path) is None:     # fail before the long run
            sys.exit(f"error: pretext-to-asm not found; install it with: {asm.INSTALL}")
    for m in a.map:
        if not os.path.exists(m):
            sys.exit(f"error: {m} not found")
    maps = [Pretext(m, load_graphs=(i == 0)) for i, m in enumerate(a.map)]   # tracks: first map only
    # One savestate serves every map only if they share scaffolds and pixels.
    try:
        check_same_layout(maps)
    except ValueError as e:
        sys.exit(f"error: {e}")
    pt = maps[0]
    for g in ("gap", "coverage"):
        if g not in pt.graphs:
            print(f"warning: {a.map[0]} has no '{g}' track; calls that need it are weaker or skipped "
                  "(add it with PretextGraph)", file=sys.stderr)
    mode, prefixes = _detect_mode(pt, a.hap_regex)
    if a.mode != "auto":
        if a.mode == "allhaps" and mode != "allhaps":
            sys.exit(f"error: --mode allhaps needs scaffold names with two haplotype prefixes matching {a.hap_regex!r} "
                     f"(found: {', '.join(prefixes) or 'none'}); set --hap-regex or use --mode single")
        mode = a.mode
    print(f"mode: {mode}" + (f" (haplotype prefixes {', '.join(prefixes)})" if mode == "allhaps" else ""), file=sys.stderr)
    a.mode = mode
    _load_tracks(a, pt)
    if a.mode == "allhaps":
        from .allhaps import curate_allhaps, write_rows
        blocks = small = None
        if a.small_paf or (a.fasta and not a.no_small_align):
            from .homology import align_small, read_paf
            sp = a.small_paf
            if not sp:
                print("aligning small scaffolds against the whole assembly (slow the first time; --no-small-align skips it)",
                      file=sys.stderr)
                try:
                    sp = align_small(a.fasta, a.work_dir or os.path.join(a.out_dir, "work"), threads=a.threads,
                                     minimap2=a.minimap2)
                except RuntimeError as e:
                    sys.exit(f"error: {e}")
            small = read_paf(sp, min_mapq=0)        # mapq 0 kept: copies within a haplotype score 0
            print(f"small scaffolds: {len(small)} alignment blocks from {sp}", file=sys.stderr)
        if a.paf or a.fasta:
            from .homology import align_haplotypes, read_paf, fasta_lengths
            paf = a.paf
            if not paf:
                # Catch a wrong --fasta before a long alignment, not after.
                missing = set(pt.names) - set(fasta_lengths(a.fasta))
                if missing:
                    sys.exit(f"error: {len(missing)} scaffolds of the map are not in {a.fasta} (e.g. {sorted(missing)[0]}); "
                             "give the assembly the map was made from")
                try:
                    paf = align_haplotypes(a.fasta, a.work_dir or os.path.join(a.out_dir, "work"), threads=a.threads,
                                           minimap2=a.minimap2, hap_regex=a.hap_regex)
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
        # Sex chromosomes (sexchrom): two sex-linked scaffolds in one haplotype are not in the
        # alignment between haplotypes, so their PARs are looked for by aligning them to each other.
        same = None
        if "coverage" in pt.graphs:
            import re
            from . import sexchrom
            groups = sexchrom.same_hap_candidates(pt, blocks, a.hap_regex, level=a.level) if a.fasta and blocks else {}
            if groups:
                from .homology import align_within, read_paf
                names = sorted(n for v in groups.values() for n in v)
                print(f"sex chromosomes: aligning {', '.join(names)} to each other within their haplotype (PARs)",
                      file=sys.stderr)
                try:
                    wp = align_within(a.fasta, names, a.work_dir or os.path.join(a.out_dir, "work"), threads=a.threads,
                                      minimap2=a.minimap2)
                except RuntimeError as e:
                    sys.exit(f"error: {e}")
                pre = lambda n: (re.match(a.hap_regex, n).group(1) if re.match(a.hap_regex, n) else "")
                same = [b for b in read_paf(wp, min_mapq=0) if pre(b.q) == pre(b.t)]
            # bedgraphs of the sex-linked stretches and PARs (added to the map with --add-tracks)
            res = sexchrom.analyse(pt, blocks, same, level=a.level)
            a.sex_track_files = sexchrom.write_tracks(pt, res, os.path.join(a.out_dir, "tracks")) if res["stretches"] else {}
        tracked = _with_tracks(a, hits)             # once: the tracks do not depend on the level
        gj, ups, ejoins = None, None, None
        if a.gfa:
            from . import graph
            try:
                gj, lay, gr = graph.load(a.gfa, a.yahs_agp, a.fasta, names=pt.names)
            except ValueError as e:
                sys.exit(f"error: {e}")
            nl = sum(j.kind == "link" for j, _ in gj)
            print(f"hifiasm graph: {nl} gap(s) between linked contig ends, {len(gj) - nl} filled by a bridging contig",
                  file=sys.stderr)
            utg = None
            if a.utg_gfa:
                print(f"reading the unitig graph {a.utg_gfa} ...", file=sys.stderr)
                utg = graph.UnitigGraph(a.utg_gfa)
                done = {(j.scaffold, j.gap_start) for j, _ in gj}
                ups = [u for u in graph.unitig_paths(lay, gr, utg) if (u.scaffold, u.gap_start) not in done]
                print(f"hifiasm unitig graph: {len(ups)} more gap(s) crossed by a unitig path", file=sys.stderr)
            # joins only the graph sees, between scaffold ends: review rows per level (below)
            ejoins = graph.end_joins(lay, gr, utg)
            print(f"hifiasm graph: {len(ejoins)} join(s) between scaffold ends that only the graph sees", file=sys.stderr)
        levels =[float(x) for x in a.sweep.split(",")] if a.sweep else [a.level]
        summary = []
        for lv in levels:
            try:
                if a.engine == "adjacency":
                    from .solver import solve
                    rows, st, info = solve(pt, blocks=blocks, level=lv, hap_regex=a.hap_regex, graph_junctions=gj,
                                           no_yahs=a.no_yahs, fine_px=a.fine_px,
                                           end_joins=graph.scaffold_end_joins(lay, ejoins) if ejoins else None)
                    print("adjacency solver: " + ", ".join(f"{k} {v}" for k, v in info.items()), file=sys.stderr)
                else:
                    rows, st, info = curate_allhaps(
                        pt, phasing=a.phasing, hi=maps[1] if len(maps) > 1 else None, blocks=blocks, busco_hits=hits,
                        small_blocks=small, level=lv, hap_regex=a.hap_regex, same_hap_blocks=same,
                        sex_system=a.sex_system, graph_junctions=gj, unitig_paths=ups,
                        end_joins=graph.scaffold_end_joins(lay, ejoins) if ejoins else None)
            except ValueError as e:
                sys.exit(f"error: {e}")
            out_dir = os.path.join(a.out_dir, f"level_{lv:.2f}") if a.sweep else a.out_dir
            os.makedirs(out_dir, exist_ok=True)
            prefix = os.path.join(out_dir, a.prefix)
            write_rows(rows, prefix + ".decisions.tsv")
            write_agp(pt, st, prefix + ".agp", unplaced=unplaced_hap_tags(pt, st))
            _maybe_asm(a, prefix + ".agp", out_dir, pt.names)
            # One savestate per map, each with that map's hash (it covers the file name).
            for path in [m.path for m in maps] + tracked:
                name = os.path.basename(path) + ".savestate_autocurate" + (f"_level{lv:.2f}" if a.sweep else "")
                st.save(os.path.join(out_dir, name), path)
            n = collections.Counter(r["action"] for r in rows)
            summary.append((lv, n, len(st.edits)))
            nrev = _write_review(rows, prefix + ".review.tsv", "allhaps")
            _write_run_info(a, prefix + ".run.txt", "allhaps", f"level\t{lv}\n")
            print(f"level {lv:.2f}\t" + "\t".join(f"{k} {v}" for k, v in sorted(n.items())))
            print(f"edits {len(st.edits)}, painted {len(st.scaffs)}, tags {len(st.meta_flags)}, waypoints {len(st.waypoints)}; "
                  f"{nrev} items to review in {prefix}.review.tsv; outputs in {out_dir}", file=sys.stderr)
        if a.sweep:
            acts = sorted({k for _, n, _ in summary for k in n})
            with open(os.path.join(a.out_dir, a.prefix + ".sweep.tsv"), "w") as f:
                f.write("level\tedits\t" + "\t".join(acts) + "\n")
                for lv, n, e in summary:
                    f.write(f"{lv:.2f}\t{e}\t" + "\t".join(str(n.get(k, 0)) for k in acts) + "\n")
            print(f"sweep summary: {os.path.join(a.out_dir, a.prefix + '.sweep.tsv')}", file=sys.stderr)
        return
    decisions, st = curate_single(pt, apply_inserts=a.apply_inserts)
    os.makedirs(a.out_dir, exist_ok=True)
    prefix = os.path.join(a.out_dir, a.prefix)
    write_decisions(decisions, prefix + ".decisions.tsv")
    write_agp(pt, st, prefix + ".agp", unplaced=unplaced_hap_tags(pt, st))
    _maybe_asm(a, prefix + ".agp", a.out_dir, pt.names)
    tracked = _with_tracks(a, None)
    for path in [m.path for m in maps] + tracked:
        st.save(os.path.join(a.out_dir, os.path.basename(path) + ".savestate_autocurate"), path)
    n = collections.Counter(d.action for d in decisions)
    nrev = _write_review(decisions, prefix + ".review.tsv", "single")
    _write_run_info(a, prefix + ".run.txt", "single")
    print("\t".join(f"{k} {v}" for k, v in sorted(n.items())))
    print(f"edits {len(st.edits)}, painted contigs {len(st.scaffs)}, tags {len(st.meta_flags)}, "
          f"waypoints {len(st.waypoints)}; {nrev} items to review in {prefix}.review.tsv; outputs in {a.out_dir}",
          file=sys.stderr)


def _asm_args(p):
    """Options shared by `asm` and `run --pretext-to-asm`."""
    p.add_argument("--asm-name", help="name of the assembly files (default: the FASTA name up to its first dot)")
    p.add_argument("--asm-args", default="", help="further pretext-to-asm options, quoted as one string, e.g. "
                                                  "'--autosome-prefix SUPER_ --min-contig-length 1000'")
    p.add_argument("--pretext-to-asm-path", help="path to pretext-to-asm (default: PATH, then conda env pretext-autocurate)")


def _graph_args(p):
    """hifiasm graph and YaHS AGPs (`run` and `asm`)."""
    p.add_argument("--gfa", action="append", default=[],
                   help="hifiasm contig graph of one haplotype (*.hap1.p_ctg.gfa, with or without sequence); repeat "
                        "for each haplotype. With --yahs-agp: gaps between contig ends the graph links are kept "
                        "(and closed by asm), and small contigs that fill a gap unambiguously are put there by asm")
    p.add_argument("--yahs-agp", action="append", default=[],
                   help="YaHS AGP of one haplotype (how its scaffolds were built from the GFA's contigs), as "
                        "H1=path or a path whose name contains hap1/hap2; repeat for each haplotype")
    p.add_argument("--utg-gfa", help="hifiasm unitig graph (*.p_utg.gfa; uncompressed, with sequence for --fill-gaps), "
                                     "with --gfa and --yahs-agp: a gap whose contig ends a unitig path connects is not "
                                     "rearranged (run), and contigs that overlap along a single path are joined (asm)")
    p.add_argument("--fill-gaps", action="store_true",
                   help="asm: fill a gap that a single unitig path crosses with the path's sequence between the two "
                        "contigs (often sequence the other haplotype also has; each fill is listed)")


def _maybe_asm(a, agp, out_dir, names):
    """`run --pretext-to-asm`: the assembly from this run's AGP, in OUT_DIR/assembly."""
    if not a.pretext_to_asm:
        return
    try:
        asm.build(a.fasta, agp, os.path.join(out_dir, "assembly"), a.asm_name, exe=a.pretext_to_asm_path,
                  extra=a.asm_args, gfa=a.gfa, yahs_agp=a.yahs_agp, names=names, utg_gfa=a.utg_gfa,
                  fill_gaps=a.fill_gaps)
    except RuntimeError as e:
        sys.exit(f"error: {e}")


def cmd_asm(a):
    """`asm`: a savestate on its map -> the AGP PretextView would export -> pretext-to-asm."""
    pt = Pretext(a.map, load_graphs=False)
    st = SaveState.load(a.state)
    if st.header_hash != pt.hash:
        sys.exit(f"error: {a.state} was not saved from {a.map} (hash mismatch)")
    os.makedirs(a.out_dir, exist_ok=True)
    name = a.asm_name or asm.default_name(a.fasta)
    agp = os.path.join(a.out_dir, name + ".pretext.agp")
    write_agp(pt, st, agp, unplaced=unplaced_hap_tags(pt, st))
    if bool(a.gfa) != bool(a.yahs_agp) or ((a.utg_gfa or a.fill_gaps) and not a.gfa):
        sys.exit("error: --gfa and --yahs-agp go together (the AGPs say where the graph's contigs are), and "
                 "--utg-gfa / --fill-gaps need both")
    if a.fill_gaps and not a.utg_gfa:
        sys.exit("error: --fill-gaps needs --utg-gfa (the unitig graph with sequence)")
    try:
        asm.build(a.fasta, agp, a.out_dir, name, exe=a.pretext_to_asm_path, extra=a.asm_args,
                  gfa=a.gfa, yahs_agp=a.yahs_agp, names=pt.names, utg_gfa=a.utg_gfa, fill_gaps=a.fill_gaps)
    except RuntimeError as e:
        sys.exit(f"error: {e}")


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
    p = sub.add_parser("run", help="pre-curation pass -> savestate per map, decisions, review list, AGP")
    p.add_argument("--map", action="append", required=True,
                   help="the MQ0 map first (it is analysed); repeat for other maps of the same assembly (in allhaps mode "
                        "the second is used as the high-MQ map); each gets its own savestate")
    p.add_argument("--mode", choices=["auto", "single", "allhaps"], default="auto",
                   help="auto (default): allhaps when scaffold names carry two haplotype prefixes (see --hap-regex)")
    p.add_argument("--hap-regex", default=r"^(H\d+)_",
                   help="regular expression whose first group is the haplotype prefix of a scaffold name [^(H\\d+)_]")
    p.add_argument("--out-dir", default=".")
    p.add_argument("--prefix", default="autocurate")
    p.add_argument("--phasing", choices=["none", "hic", "trio"], default="hic",
                   help="allhaps only: how freely sequence may move between haplotypes. "
                        "trio: never; hic: whole chromosomes only, no H1_/H2_ mixing within one; none: freely")
    p.add_argument("--level", type=float, default=0.5,
                   help="allhaps: how readily to act, 0 (conservative: only the clearest evidence, flag the rest) to 1 "
                        "(aggressive); moves every act-or-flag threshold together [0.5]")
    p.add_argument("--engine", choices=["rules", "adjacency"], default="rules",
                   help="allhaps: rules = the step-by-step pass; adjacency = the experimental adjacency solver "
                        "(contig ends scored and joined best first; solver.py), for comparison [rules]")
    p.add_argument("--no-yahs", action="store_true",
                   help="adjacency engine: do not use the YaHS layout as evidence (scaffold from the contigs; its "
                        "joins are still reported for comparison)")
    p.add_argument("--fine-px", type=int, default=0,
                   help="adjacency engine: score contigs too small for the binned contact windows at full map "
                        "resolution, down to this many pixels, instead of merging them into a neighbour (0: off) [0]")
    p.add_argument("--sweep", help="allhaps: comma-separated levels, e.g. 0,0.25,0.5,0.75,1: one output directory and "
                                   "savestate per level (OUT_DIR/level_X.XX) and a summary table, to pick a starting point")
    p.add_argument("--telomere-motif", help="run telo_scan.py (CNAG, bundled) on --fasta for this motif (e.g. TTAGGG for "
                                            "vertebrates, TTAGG for most arthropods) and use telomeres by strand")
    p.add_argument("--telo-scan", help="use this telo_scan.py instead of the bundled copy")
    p.add_argument("--telomere-5p", help="existing telo_scan.py 5p_telomere bedgraph (with --telomere-3p)")
    p.add_argument("--telomere-3p", help="existing telo_scan.py 3p_telomere bedgraph (with --telomere-5p)")
    p.add_argument("--repeats", help="repeat intervals (BED: name, start, end), e.g. from Red")
    p.add_argument("--red", action="store_true", help="run Red on --fasta (each haplotype) for a repeat track")
    p.add_argument("--add-tracks", action="store_true",
                   help="also write a copy of each map with the evidence tracks (repeats, telomeres by strand, "
                        "BUSCO duplicate pairs and genes missing from one haplotype) for PretextView, with its own "
                        "savestate; needs PretextGraph")
    p.add_argument("--fasta", help="allhaps: the assembly the maps were made from (plain or .gz); "
                                   "haplotypes are aligned with minimap2 and homology is read from that")
    p.add_argument("--paf", help="allhaps: an existing alignment between haplotypes, used in place of --fasta")
    p.add_argument("--small-paf", help="allhaps: an existing alignment of the small scaffolds against the whole assembly "
                                       "(minimap2 -D, secondary alignments kept), used in place of running it from --fasta")
    p.add_argument("--no-small-align", action="store_true",
                   help="allhaps with --fasta: do not align the small scaffolds against the whole assembly")
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
    p.add_argument("--sex-system", choices=["XY", "ZW", "X0", "Z0", "X1X2Y", "XY1Y2", "Z1Z2W", "ZW1W2"],
                   help="allhaps: the species' sex chromosome system, if known; the sex chromosome copies found from "
                        "coverage are then tagged: XY/ZW the larger X (Z); X0/Z0 all X (Z); X1X2Y etc. the copy homologous "
                        "to both others Y (W), the others X1, X2 by size. Without it they are tagged HAP1 only")
    p.add_argument("--apply-inserts", action="store_true", help="insert localised pieces instead of only proposing them")
    _asm_args(p)
    _graph_args(p)
    p.add_argument("--pretext-to-asm", action="store_true",
                   help="also build the assembly from this pre-curation (needs --fasta and pretext-to-asm): "
                        "OUT_DIR/assembly/NAME.1.* (with --sweep, one per level directory)")
    p.set_defaults(fn=cmd_run)
    p = sub.add_parser("asm", help="savestate -> PretextView-style AGP -> assembly FASTA/AGP files with pretext-to-asm "
                                   "(sanger-tol/agp-tpf-utils)")
    p.add_argument("--map", required=True)
    p.add_argument("--state", required=True, help="the savestate, e.g. the curator's after review")
    p.add_argument("--fasta", required=True, help="the assembly the map was made from (plain or .gz)")
    p.add_argument("--out-dir", default=".")
    _asm_args(p)
    _graph_args(p)
    p.set_defaults(fn=cmd_asm)
    a = ap.parse_args(argv)
    a.fn(a)


if __name__ == "__main__":
    main()
