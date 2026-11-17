# Changelog

## Unreleased

- Alignment blocks are shared out over the cut pieces they cover (`common.split_block`)
  instead of being credited to the pieces at their midpoints: one 56 Mb Troglobisium block
  across three cuts had left H1_scaffold_1's first 25 Mb and its homolog H2_scaffold_4
  33-58 Mb unpaired (now one chromosome, as the curator laid it out).
- Chromosome-level joins (`allhaps/chromjoin.py`, curator: chromosome arms apart across
  centromeric repeats touch only weakly, over their whole length; common in beetles, hard to
  see): two pairs are merged when, with repeat bins masked (coverage, telomere signal,
  promiscuous bins), their copies are each other's best partner by `chrom_margin`, with
  contact spread over `chrom_spread` of the block (a hot spot from a shared repeat does not
  count), confirmed in the other haplotype, never across a cut the pass made, and each side
  at least `chrom_min_mb` (small chromosomes cluster). Level 0: review rows
  `merge_candidate` only. Troglobisium: the H2_scaffold_18 / 66 chromosome; owl unchanged.

- Joins only hifiasm's graph sees are made (curator's review of the Troglobisium rows: all
  should be placed): `allhaps/graphjoin.py` puts the smaller scaffold at the larger one's
  end as the graph orients it, ahead of what the tool had put there; rows `graph_join`, or
  `graph_join_candidate` saying why not (under one pixel, the larger side an Unloc or
  unplaced). Troglobisium: 10 of 15 made.
- Pairs that are parts of one chromosome (`merge_pairs`) are also merged with weaker end
  contact (`merge_link_weak`) when their whole-block contact is `merge_block_min` times the
  trans background in both haplotypes, or when hifiasm's graph joins their ends. The
  curator found these joins obvious on the map: end contact 0.25-0.5, whole blocks 7-36x,
  against 3.5x for a false candidate. Troglobisium at level 0.5 against the curation:
  junction recall 0.57 -> 0.64, chromosome grouping 89% -> 94%. `tools/missed_joins.py`
  lists a reference curation's joins the tool does not make.

- Two copies of one region in a scaffold (`allhaps/selfdup.py`): from the alignment between
  the haplotypes, review rows `self_duplication` classed by read depth (coverage track) as
  repeat, duplication or haplotypic. Curator's rule: only a haplotypic copy (it and the other
  haplotype's copy at about half depth) with a gap between the copies is cut out, at the
  gaps around it; the scaffold's sides are joined again as assembled and the copy is tagged
  FalseDuplicate. Troglobisium: 12 rows, none haplotypic; owl: none. `tools/self_dups.py`
  lists such copies from a PAF. `synth.toy_map` takes coverage per region; `toy_blocks`
  reads a scaffold's haplotype from its name.

- `allhaps.py` split into the package `allhaps/` (common, pairing, sites, evidence, small,
  arrange, driver, state); `curate_allhaps` is now nine steps over a shared `Run`
  (`driver.STEPS`, `stop_after=` for tests), and every level-dependent threshold is named
  (`common.THRESHOLDS`, `Thresholds`). Outputs byte-identical on all golden runs
  (`tools/golden.py`). New `tests/test_steps.py` (per-step tests on toy genomes) and
  `synth.toy_blocks` (alignment blocks for a toy genome). Fixed on the way: `place_small`
  crashed on a map with a single homolog pair.

- Joins only hifiasm's graph sees, between the outer ends of two scaffolds (a contig link
  that is the only one at both ends, or mutually nearest ends along the unitig graph),
  where one scaffold is chromosome-size and the layout does not make the join: review rows
  `graph_join_candidate` with the Hi-C end contact against the best other end. Not made
  automatically. Troglobisium: 18, nearly all a small piece at a chromosome end.

- Tests that run without the private fixtures: `pretext_autocurate/synth.py` writes small
  synthetic `.pretext` maps of a toy genome, `tests/test_synthetic.py` runs the whole pass
  and the command line on them; GitHub Actions runs these on every push. Fixed on the way:
  the decay curve crashed on maps with textures under 200 pixels.

- `tools/score_all.py` (with `tools/score_cases.json` and `pretext_autocurate/score.py`):
  runs the tool on every test genome and scores it against the curators' savestates
  (junctions, breaks, chromosome grouping, tags by bp, edits, review items) in one table.
  Baseline for v1.2.0 in `docs/PROGRESS.md`.
- Small scaffolds (under 20 px, no sequence evidence) are tagged Unloc when their high-MQ
  signal is unique to one homolog pair (curator's rule; MQ0 uniqueness strengthens it, and
  stands alone without a high-MQ map), with a floor on the signal so that a few reads do not
  count. Before, Hi-C alone was never tried below 20 px.
- Large attached pieces are placed in their chromosome by end contact (when the place and the
  contact are clear) instead of being left as Unlocs at its end, also when the chromosome's
  own pieces are not rearranged. Troglobisium: 1 more junction agrees with the curation at
  each level, none lost.

## 1.2.0 (2026-11-16)

Sex chromosomes, hifiasm's graphs, and the final assembly files.

Sex chromosomes (`sexchrom.py`)
- Hemizygous stretches (half coverage, no alignment to the other haplotype) mark the
  sex-specific part of a sex chromosome; a half-coverage stretch aligned into one is the
  other sex chromosome's homologous part. Divergent heterozygous stretches of autosomes,
  also at half coverage but aligned to each other, are no longer mistaken for them.
- Normal-coverage sequence next to a sex-linked stretch is a PAR when it has
  chromosome-level contact with it and aligns to (or, without alignments, shares MQ0 reads
  with) the other sex chromosome's end; with no contact it is another chromosome joined
  there, and the cut row says so. Two sex-linked scaffolds in one haplotype are aligned to
  each other with `--fasta` (`homology.align_within`). Hemizygous scaffolds in no pair join
  the chain of the sex chromosome they contact.
- `--sex-system` (XY, ZW, X0, Z0, X1X2Y, XY1Y2, Z1Z2W, ZW1W2) tags the copies; without it
  they are not labelled (that needs species knowledge). A `sex_summary` row reports what was
  found, and when nothing was, that a homogametic sample or homomorphic sex chromosomes look
  like autosomes. New rows `sex_linked_region`, `par`, `par_candidate`, `sex_boundary`,
  `sex_chromosome_join`; new tracks `sex_chromosome` and `par`.

hifiasm's graphs (`graph.py`; `--gfa`, `--yahs-agp`, `--utg-gfa`, `--fill-gaps`)
- The YaHS AGPs place every contig in its scaffold, and the GFAs' read lists place contigs
  on the unitig graph; no graph alignment is needed, and GFAs with or without sequence work.
- A gap whose two contig ends the contig graph links (`graph_link`), or that one small
  contig bridges (`graph_bridge`), is never cut. A gap a unitig path crosses is not used to
  rearrange (`graph_path`; a tip inversion or a contact-only cut there is dropped,
  `graph_keep`; cuts backed by the alignment or a telomere stand), and neither is its
  homologous kept site.
- When the assembly files are built (`asm`): bridges and contigs of the gap's own haplotype
  that lie wholly on a unitig path are put into the gap; gaps between linked contig ends,
  or contigs that overlap or abut along a unitig path, are closed by trimming the overlap
  (checked in the sequence); with `--fill-gaps`, a gap is filled only with collapsed
  homozygous sequence (unitig read depth about twice the haploid, none of it in this
  haplotype's other contigs, not next to a circular contig), never with the other
  haplotype's own allele. `<name>.graph_edits.tsv` lists every gap and why;
  `<name>.gapfills.fa` holds the fills.
- Troglobisium: 2 bridges and 2 contigs placed (all pieces the curator had left
  unplaced), 16 gaps closed, 15 filled (1.4 Mb; each fill of 1 kb or more is, by
  alignment, only in the other haplotype); H1_scaffold_16's end piece no longer inverted.
  Every edited chromosome checked against pretext-to-asm's own sequence.

Assembly files (`asm.py`)
- The AGP is written exactly as PretextView exports it after loading the savestate
  (`# HiC MAP RESOLUTION` header, pixel-derived coordinates, painted scaffolds first, "-"
  only inside painted scaffolds of several pieces), without the paint-display tags
  (`vertPaint` used to reach it), and with the scaffolds too short for a pixel appended,
  tagged HAP1/HAP2 by prefix (else pretext-to-asm makes separate h1/h2 assemblies).
- New command `asm` and `run --pretext-to-asm`: the assembly files from a savestate with
  pretext-to-asm (sanger-tol/agp-tpf-utils, installed separately); `--asm-name`,
  `--asm-args`, `--pretext-to-asm-path`.

Layout
- Contigs outrank scaffold joins: the copy of a chromosome that leads the other's order and
  orientation is the one whose arrangement the other copy's contigs contradict less; when
  each contradicts the other, the pair is flagged `arrangement_differs`.
- Extra pieces of a chromosome copy are split at their kept gaps before order and
  orientation are searched, so one piece can go into another's gap. Troglobisium: the X is
  laid out as the curator did; 24 of 32 tool junctions in the curation (21 of 30 before).

Measured and not built (PROGRESS): a shared Hi-C backbone of both copies
(`tools/backbone_gaps.py`: the copies give nearly the same contact, no gain) and merges
checked by alignment (no sequence spans a merge join).

Other
- Every module commented inline; `docs/HANDOFF.md` lists what is left.
- `tools/graph_junctions.py`: graph evidence at the junctions of any layout.

## 1.1.0 (2026-11-14)

Evidence tracks: telomeres by strand, repeats, BUSCO pairs.

- Telomeres by strand from telo_scan.py (CNAG), now bundled in `pretext_autocurate/vendor/`:
  `--telomere-motif MOTIF` with `--fasta` runs it (a few minutes for a 4 Gb assembly);
  `--telomere-5p` / `--telomere-3p` take existing telo_scan bedgraphs; maps that already
  carry its `5p_telomere` / `3p_telomere` tracks are read directly. At a gap: telomeres
  facing it from both sides, or from the long side, mean two chromosome ends joined (cut);
  from the end piece only, an inverted end piece (inverted, with the evidence in the
  decisions table). The check includes the gap's own pixel, since telomeres are shorter
  than a bin.
- Repeats: `--repeats BED`, `--red` (runs Red on each haplotype of `--fasta`), or a
  `repeat_density` track in the map. Alignment breaks inside repeats are ignored; a
  scaffold placed by Hi-C alone is not moved when it is 90% repeat. With sequence evidence
  the decision rests on unique alignments and repeat content does not count.
- `--add-tracks`: a copy of each map with the evidence tracks for PretextView
  (`5p_telomere`, `3p_telomere`, `repeat_density`, `busco_dup_pairs`, `busco_missing_hap`)
  and a savestate for each copy (needs PretextGraph). On a high-MQ map a repeat track tells
  repeat from contaminant. BUSCO duplicate pairs share one colour per pair of scaffolds, so
  both places of a possible haplotig light up together; genes missing from one haplotype
  are shown where the other has them.
- Results on the test genomes: on the boreal owl the strands confirm the scaffold_1 cut and
  the H1_scaffold_2 tip, and invert the end pieces of H1_scaffold_5 (as in the original
  curation) and H2_scaffold_2; on Troglobisium (TTAGG) four flagged sites become cuts, all
  at places the curator broke.
- Tools and notes: `tools/bench_aligners.py` (minimap2 against FastGA; minimap2 stays the
  aligner), `docs/RELATED.md` (hapscotch), `CLAUDE.md` (run commands under bash).

## 1.0.1 (2026-11-14)

- License: PolyForm Noncommercial License 1.0.0 (`LICENSE.md`). Free for noncommercial use,
  including by educational institutions and public research organizations; commercial use
  requires a license from CNAG (`NOTICE`).
- Repository moved to https://github.com/cnag-aat/pretext-autocurate.
- No change to the code or its results.

## 1.0.0 (2026-11-14)

First release for testing by curators ("pre-curation").

Two-haplotype mode (`allhaps`)
- Homolog pairing from sequence alignment (minimap2, `--fasta` / `--paf`) or MQ0 Hi-C;
  pairs laid out side by side and parallel, better copy tagged HAP1 (`--phasing`).
- Cuts at gaps with reduced contact, at alignment breaks and at interstitial telomeres,
  carried to the matching gap of the other haplotype; sharp breaks away from gaps.
  Telomere signal is supporting evidence only.
- Order and orientation of pieces within each chromosome copy, along the other haplotype
  when aligned, by end contact otherwise; merging of pairs that are parts of one
  chromosome; scaffolds joined where they align.
- Small scaffolds: unlocs where sequence and Hi-C agree, joins when a piece continues a
  copy, `FalseDuplicate` for second copies, decisions on non-repetitive sequence only.
- Sex chromosome candidates from half coverage.
- BUSCO tables as an optional tie-breaker and for false-duplication candidates.
- `--level` (conservative to aggressive) and `--sweep` (one savestate per level).

Both modes
- `--mode auto` (default) from the scaffold names; `--hap-regex`.
- New outputs: `autocurate.review.tsv` (what to look at) and `autocurate.run.txt`
  (version, command and options).
- Clear errors for missing maps, maps of different assemblies and missing prefixes.

Known limitations
- Hemizygous chromosomes built from several scaffolds are not chained.
- Two haplotypes of one chromosome inside one scaffold are not separated.
- No repeat track yet (Red / PretextGraph planned); repeats are inferred from alignments.
- The AGP has not been checked against PretextView's own export.
