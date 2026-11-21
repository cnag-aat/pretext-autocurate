# pretext-autocurate

Pre-curation of Hi-C contact maps, straight from PretextMap `.pretext` files. It does the
routine first pass of manual curation and writes a PretextView savestate that a curator
opens, reviews and finishes, together with a table explaining every decision and a short
list of the things that need a curator's eye. It does not replace the curator: marginal
calls are flagged, never hidden.

**Version 1.2.** Developed and tested on three genomes: a gastropod (single haplotype),
the boreal owl *Aegolius funereus* (two haplotypes, MQ0 + MQ10) and the pseudoscorpion
*Troglobisium racovitzai* (two haplotypes, heavily rearranged, MQ0 only). Thresholds are
relative to each map's own signal. Always review the output before accepting it.

## Install

    git clone https://github.com/cnag-aat/pretext-autocurate.git && cd pretext-autocurate
    conda env create -f environment.yml
    conda activate pretext-autocurate
    pretext-autocurate --version

Or, with Python >= 3.10 and numpy available: `pip install -e .`. minimap2 (and miniprot for
`--busco-lineage`) are only needed when the tool aligns sequence itself; they are looked
for on `PATH`, then in a conda environment named `pretext-autocurate`. To build the final
assembly files (`asm`, `run --pretext-to-asm`), also install pretext-to-asm:
`pip install git+https://github.com/sanger-tol/agp-tpf-utils`.

Tiara (`--tiara`) pins pytorch 1.7, so it goes in an environment of its own, found by name
(`tiara`, or `tiara-x86`), or give `--tiara-path`:

    conda create -n tiara -c conda-forge -c bioconda tiara
    # Apple silicon (no arm64 pytorch 1.7): an Intel environment, run under Rosetta;
    # libcxx 14 because Tiara's llvmlite needs a symbol newer libc++ dropped
    CONDA_SUBDIR=osx-64 conda create -n tiara-x86 -c pytorch -c conda-forge -c bioconda tiara "libcxx=14"

The biocontainer (`quay.io/biocontainers/tiara`) works too: put a `tiara` wrapper on `PATH`.

## Quick start

Two haplotypes (scaffolds named `H1_...` / `H2_...`), the recommended way: give the maps
and the assembly they were made from.

    pretext-autocurate run \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext \
        --fasta asm.diploid.fa.gz --threads 16 --out-dir precur

With hifiasm's graphs and the YaHS AGPs as well (recommended when you have them; see
"Optional evidence"), and the assembly files built at the end:

    pretext-autocurate run \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext \
        --fasta asm.diploid.fa --threads 16 \
        --gfa asm.hap1.p_ctg.gfa --gfa asm.hap2.p_ctg.gfa \
        --yahs-agp asm.hap1.yahs.agp --yahs-agp asm.hap2.yahs.agp \
        --utg-gfa asm.p_utg.gfa --fill-gaps --pretext-to-asm --out-dir precur

After your review in PretextView, the same files from your savestate:

    pretext-autocurate asm --map asm.diploid_mq0.pretext --state my.savestate \
        --fasta asm.diploid.fa --gfa ... --yahs-agp ... --utg-gfa asm.p_utg.gfa --fill-gaps \
        --out-dir final

Single haplotype (or a collapsed assembly):

    pretext-autocurate run --map asm.mq0.pretext --out-dir precur

The mode is detected from the scaffold names (`--mode single|allhaps` to force it). The
first `--map` is analysed and should be the MQ0 map; it should carry the `gap` and
`coverage` tracks (and ideally `telomere`). Every further map of the same assembly gets a
savestate with the same edits; in two-haplotype mode the second map is read as the high-MQ
map.

Then, in PretextView: load the map, load `precur/<map name>.savestate_autocurate`, and go
through `precur/autocurate.review.tsv`.

## Outputs

| File | Content |
|---|---|
| `<map name>.savestate_autocurate` | load in PretextView after loading that map (one per `--map`) |
| `autocurate.review.tsv` | what to look at, most important first, with the reason |
| `autocurate.decisions.tsv` | every action with its evidence (cut, join, pairing, unloc, flag ...) |
| `autocurate.agp` | the proposed assembly, written as PretextView would export it after loading the savestate (input for pretext-to-asm) |
| `autocurate.run.txt` | version, command line and every option, for reproducibility |
| `autocurate.breaks.tsv` | with `--fasta` (both modes): cuts placed to the base by the sequence and gapless joins, for `asm --breaks` |
| `work/` | cached alignments (reused on the next run) |
| `tracks/` | bedgraphs of the evidence tracks (sex chromosomes, PARs, telomeres by strand, ...) |
| `assembly/` | with `--pretext-to-asm`: pretext-to-asm's files (FASTA and AGP per haplotype, chromosome lists, ...), plus `<name>.graph_edits.tsv` (every gap closed, filled or left, and contigs placed, with the reason) and `<name>.gapfills.fa` (gap fills) when the hifiasm graphs are given; `asm` writes the same into its `--out-dir` |

In the savestate:

- **Paint**: one colour per chromosome copy; unlocs and false duplicates are painted with
  their chromosome and sit at its end.
- **Tags**: `HAP1` / `HAP2` on everything in two-haplotype mode; `Unloc`, `Haplotig`,
  `FalseDuplicate`; sex chromosome candidates are all `HAP1` (and `X`, `Y`, `Z`, `W`,
  `X1` ... with `--sex-system`).
- **Waypoints**: positions inside contigs to look at (cuts, flagged sites, proposed
  insertion points).
- **vertPaint** (vertical band): a standalone piece that needs a look (an unloc candidate,
  a BUSCO-only duplicate, a borderline haplotig).

## How conservative: `--level` and `--sweep`

One setting moves every act-or-flag threshold together: `--level 0` acts only on the
clearest evidence and flags the rest, `--level 1` acts on weaker evidence too (default
0.5). To compare, write one savestate per level and pick the one to start from:

    pretext-autocurate run ... --sweep 0,0.25,0.5,0.75,1 --out-dir precur

This writes `precur/level_0.00/` ... `precur/level_1.00/` and `precur/autocurate.sweep.tsv`
(counts of each action per level). After the first level, each further one takes a
fraction of the time. `--sweep` applies to two-haplotype mode (and to `--engine adjacency`). In
single-haplotype mode the rules pass reads `--level` for its insertion and inversion thresholds
(`curate.SINGLE_THRESHOLDS`):

| threshold | level 0 | 0.5 | 1 |
|---|---|---|---|
| orientation ratio to propose an insertion (waypoint at the gap) | 1.3 | 1.2 | 1.15 |
| orientation ratio to apply an insertion | 1.4 | 1.2 | 1.15 |
| smallest piece (px) whose insertion is applied, not only proposed | 20 | 15 | 8 |
| inversion: contact across the reversed ends needed to apply (1 = adjacent) | 1.2 | 1.0 | 0.8 |
| inversion: contact across the assembled ends allowed (0 = as far apart as the segment is long) | 0.2 | 0.3 | 0.4 |
| inversion at a scaffold end (half a bow-tie): applied (two haplotypes: only if the homolog alignment agrees) | no | no | at 1.0 |

An inversion is a segment assembled backwards inside a scaffold. Its two ends touch the
wrong neighbours, which shows as a bow-tie of off-diagonal contact at the segment's corners.
The tool calls one only when the reversed ends touch like adjacent sequence **and** the
assembled ends lose contact. When both stay high (both orientations in the reads, as with a
heterozygous inversion or a repeat), nothing is called. Weaker calls are `flag_invert` rows
with waypoints at both ends.

## What it does

Two haplotypes (`allhaps`):

- **Cuts** at gaps where contact across is reduced, where the alignment to the other
  haplotype stops being colinear, or where a telomere sits inside a scaffold; and, away
  from gaps, only where contact drops to background. A cut found in one haplotype is
  carried to the matching gap of the other.
- **Pairs** homologous scaffolds and lays each pair side by side, parallel. The better copy
  (fewer gaps, telomeres, length; BUSCO as a tie-breaker) is `HAP1`.
- **Orders and orients** the pieces of each chromosome copy: along the other haplotype when
  the alignment is available (the copy whose contigs say so leads: contigs outrank scaffold
  joins), by end-to-end contact otherwise; pairs that are parts of one
  chromosome are merged; scaffolds that align into a copy are joined where they align.
- **Small scaffolds**: attached as `Unloc` where sequence and Hi-C agree, or, for pieces too
  small to weigh contact ratios, where their high-MQ signal is unique to one homolog pair
  (MQ0 without a high-MQ map); joined when they continue a copy, and large ones placed in
  their chromosome where their ends touch it clearly best; a second copy of sequence its haplotype already has is tagged
  `FalseDuplicate`; repeat-rich scaffolds are judged on their non-repetitive part.
- **Sex chromosomes** (from the `coverage` track, and the alignment when given): a
  stretch at about half coverage that does not align to the other haplotype is
  hemizygous, i.e. the sex-specific part of a sex chromosome; a half-coverage stretch that
  aligns into one is the other sex chromosome's homologous part. Next to it, normal-coverage
  sequence that shares chromosome-level contact with it and aligns to the other sex
  chromosome's end is a PAR (recombines like an autosome; never a haplotig); with no
  contact it is another chromosome joined there (the cut row says so). Hemizygous scaffolds
  in no pair join the sex chromosome they contact. Sex chromosomes are all `HAP1`. Which is
  which needs knowledge of the species, so they are only labelled with `--sex-system`:
  `XY` / `ZW` (the larger of the pair `X` / `Z`), `X0` / `Z0` (all `X` / `Z`), or `X1X2Y`,
  `XY1Y2`, `Z1Z2W`, `ZW1W2` (the copy homologous to both others is the single one).
  Only what is hemizygous in the sequenced individual shows: a sample of the homogametic
  sex (XX, ZZ) or homomorphic sex chromosomes look like autosomes and are paired as such;
  the `sex_summary` row says when nothing was found.
- **hifiasm's graphs** (`--gfa`, `--yahs-agp`, `--utg-gfa`): a gap whose two contig ends
  the contig graph links, or that one small contig bridges, is a right join and is never
  cut; a gap a unitig path crosses is not used to rearrange. When the assembly files are
  built (`asm`), linked contigs are joined by trimming their overlap, bridging contigs and
  contigs of the gap's own haplotype that lie on a unitig path are put into the gap, and,
  with `--fill-gaps`, collapsed homozygous sequence fills the gap (never sequence the
  other haplotype already holds as its own allele). PretextView cannot express these
  changes, so they are made on the assembly files, each listed with its evidence. Joins
  only the graph sees between scaffold ends (usually a small piece at a chromosome end)
  are listed for review (`graph_join_candidate`), not made.
- **With `--fasta`, as in single-haplotype mode** (below): an end of one scaffold that copies an
  end of another scaffold of the same haplotype is cut off at the base the alignment gives and
  tagged `FalseDuplicate` (`overlap`, a review row: re-phase it instead if the other haplotype
  lacks that sequence); gaps with a copy on both sides are closed by `asm` (`close_gap`); both
  go to `autocurate.breaks.tsv`.
- **Inversions** inside a scaffold are turned in place (`inversion`) when the bow-tie passes the
  level's thresholds and the segment aligns to the other haplotype in the other direction from
  its flanks; at a scaffold end that alignment is required. The same pattern at the same place
  in both haplotypes is folding of the chromosome, not a misassembly (`flag_inversion`).
- **Shrapnel** specific to one homolog pair (MQ0 contact with it at least 3x that with any other
  pair) goes with that pair as an `Unloc` of its same-haplotype copy, even when nothing says
  where in it: pieces in no pair (including a piece cut off a scaffold), further pieces
  homologous to a scaffold, and small pairs of their own (both copies under 10 Mb) whose two
  copies both belong to a larger pair. A curator who wants to localise them finds them next
  to their chromosome.
- **Which way the second copy runs** is read from the two chains as arranged: their alignment,
  else a clear Hi-C homology line (one vote per piece, weighted by rows; a piece whose own line
  runs the other way from the majority is turned in place), else the pairing's strand. Sex
  chromosomes are laid out with their PAR first.
- Flags (does not change): sites where both sides belong together, places where the two
  haplotypes disagree, unloc candidates with one line of evidence, scaffolds with no homolog.

`--phasing` sets how freely sequence may change haplotype: `trio` never; `hic` (default)
whole chromosomes may swap so the better copy is HAP1, and small pieces may go to the other
prefix's copy when the alignment says so; `none` currently behaves like `hic`.

Single haplotype: misjoin cuts (preferably at gaps), joins and orientation of large pieces,
rearrangements inside a scaffold (tentative cuts kept only if the pieces fit better
elsewhere), inversions (the bow-tie, see the level table), placement of small scaffolds
(inserted at a gap when clear, `Unloc` otherwise), haplotigs (off-diagonal line plus reduced
coverage; borderline ones highlighted), **displaced stretches**, and, with `--fasta`, the
items after them.

- **Displaced stretches** (`move`): a stretch inside a scaffold with almost no contact with
  its own two sides and a line off the diagonal that points at one place elsewhere. One pass
  over the map (about 100 s on a 65,536 px map, once per `--sweep`) finds them anywhere. Each
  is cut out at both ends and goes in at the junction (gap, cut or weak point) its line points
  at, its sides joined back together. With `--fasta` the stretch is aligned there: the
  alignment gives the strand, a copy at one end says it goes right next to that copy (cut off
  as a `Haplotig` when 5 px or more, joined without a gap by `asm`), and a stretch that is
  mostly a copy is a haplotig of that place. Without the sequence and with no clear
  orientation it is attached as an `Unloc`, with a waypoint at the place. The bar
  (`move_line`, `move_ratio`, `move_spec`) moves with `--level`. In allhaps mode the scan
  looks within each haplotype, cuts the stretch out (`displaced` in the review list), and the
  pairing and placement steps put it where it belongs.

- **Overlapping scaffold ends** (`overlap` / `flag_overlap`): one scaffold's end is the other
  haplotype of another scaffold's end, and the copy runs to that scaffold's end (hifiasm left
  both haplotypes at two contig ends). On the map this is a line of contact as bright as the
  diagonal pairing the two ends inwards; the sequence confirms it (the ends align in order at
  haplotype-level divergence) and says where the copy stops. The copy on the longer scaffold
  is cut off at that base and tagged `Haplotig`, and the other scaffold's end is joined to the
  cut. Without `--fasta`, or when the sequence and the map disagree on the copy, it is proposed.
- **Gaps with a copy on both sides** (`close_gap`): the sequence after a gap aligns in order
  with the sequence before it, up to the gap. Nothing changes on the map; `asm` closes the gap.
- **Candidate haplotigs are always aligned** to where their Hi-C copy line points. The
  alignment says what each is: a copy (a borderline haplotig call is then confirmed), nothing
  (said in the row), or a **gap fill** (`fill_gap`): its two ends are copies of the two sides
  of a gap there and its own sequence between them fills it. The piece goes into the gap in
  the strand the alignment gives, copies of 5 px or more are cut off at the base and tagged
  `Haplotig`, and `asm` joins both junctions without a gap. A fill needs at least 5 px of its
  own sequence (a piece that is nearly all copy is a haplotig) and copies of at least 2 px at
  both ends; otherwise it is proposed (`flag_fill_gap`). A collapsed tandem repeat can look
  like a copy; telling it from a real haplotype needs read phasing, which the tool does not
  take yet.

Cuts the sequence placed, and gapless joins, are written to `autocurate.breaks.tsv`. `asm
--breaks` (and `run --pretext-to-asm`) uses them: a cut moves from its pixel (tens of kb) to
the base, and a join at an overlap or a haplotypic gap is made without a gap, switching from
one copy to the other inside a stretch of at least 1 kb where the two are identical; the other
path through the bubble (the other copy of each side) is written as a haplotig, so no sequence
is lost. They apply only where the savestate still has that cut or join: whatever you change
in PretextView keeps its pixel position. `<name>.exact_edits.tsv` lists what was applied.

## Optional evidence

| Option | What it adds |
|---|---|
| `--fasta ASM.fa[.gz]` | aligns the haplotypes (minimap2) and the small scaffolds against everything; homology, order, orientation, duplicates. The first run takes long (both alignments for the 2.7 Gb owl assembly about 15 min; the haplotype alignment alone for a 4.4 Gb assembly about 1.5 h, 10 GB, 6 threads); alignments are cached in `--work-dir`. The small-scaffold alignment indexes the whole assembly: on a machine with little memory, skip it with `--no-small-align` |
| `--paf`, `--small-paf` | use existing alignments instead (H2 against H1, and small scaffolds against the whole assembly) |
| `--busco-table H1=hap1.full_table.tsv --busco-table H2=hap2.full_table.tsv` | existing BUSCO results (one per haplotype, or one on the diploid assembly): tie-breaker for HAP1, false-duplication candidates. Never decisive alone |
| `--busco-lineage`, `--miniprot-paf` | let the tool find BUSCO genes itself (slow) |
| `--repeats BED`, `--red` | a repeat track (e.g. from Red); `--red` runs Red on each haplotype of `--fasta` (about 25 min for a 2.7 Gb assembly) |
| `--telomere-motif MOTIF` | the telomere motif for the bundled telo_scan.py (CNAG), e.g. TTAGGG vertebrates, TTAGG most arthropods. **With `--fasta`, telo_scan runs and decides by default** (telomeres by strand; a few minutes for a 4 Gb assembly, cached); without this option the motif is detected from the ends of the largest scaffolds (the unit forming the most tandem arrays at both starts and ends; the run prints the evidence and records the motif in `run.txt`). The map's own telomere track is kept in the map and used only with `--telomere-from-map`, with no `--fasta`, or when no motif is found. `--telomere-5p/3p` take existing telo_scan bedgraphs |
| `--add-tracks` | write a copy of each map with these tracks, the BUSCO duplicate pairs / missing genes, and the sex chromosome tracks (`sex_chromosome`: 2 hemizygous, 1 the rest of a sex-linked stretch; `par`: 2 PAR, 1 candidate), for viewing in PretextView (needs PretextGraph). The bedgraphs are in `OUT_DIR/tracks/` in any case |
| `--sex-system SYSTEM` | the species' sex chromosome system, if known (`XY`, `ZW`, `X0`, `Z0`, `X1X2Y`, `XY1Y2`, `Z1Z2W`, `ZW1W2`): tags the sex chromosome copies `X`/`Y`, `Z`/`W`, `X1`/`X2` ... |
| `--gfa GFA --yahs-agp AGP` (each once per haplotype) | hifiasm's contig graph (`*.hap1.p_ctg.gfa`; with or without sequence) and the YaHS AGP that built the scaffolds from those contigs (`H1=path`, or a name containing `hap1`). Where the graph links the two contig ends YaHS put next to each other, the join is right: the gap is never cut (row `graph_link`), and `asm` closes it, dropping the overlap the graph reports after checking it in the sequence. Where the only link out of the left contig and the only one into the right go to one small contig that is a scaffold of its own, that contig fills the gap (row `graph_bridge`); `asm` puts it there and closes both sides. `OUT_DIR/<name>.graph_edits.tsv` lists every change |
| `--utg-gfa P_UTG.gfa` (with `--gfa`, `--yahs-agp`) | hifiasm's unitig graph (uncompressed; the `noseq` file is enough for `run` and loads fast). A gap whose two contig ends a unitig path joins is not used to rearrange (row `graph_path`; a tip inversion or a cut on contact alone there is dropped, `graph_keep`; a cut the alignment or a telomere backs stands). In `asm`, contigs that overlap along a single path are joined |
| `--fill-gaps` (`asm`, `run --pretext-to-asm`; needs the unitig GFA with sequence and hifiasm's read depth) | fill gaps from the unitig graph without duplicating what the assembly already has for that haplotype. Every path across the gap is examined: if the contigs overlap or abut along one, the gap is closed; if a path runs through another contig of this haplotype, that contig belongs there: it is put into the gap when it lies wholly on the path and is a scaffold of its own (also at joins the tool or you made; this needs the unitig GFA with sequence but not `--fill-gaps`), and is named in the report otherwise (nothing is copied); otherwise the gap is filled only from a path whose unitigs are at about twice the haploid read depth (collapsed homozygous sequence: two copies in the genome, one in the assembly, in the other haplotype). Haploid paths (the other haplotype's allele), repeats and gaps next to a circular contig (a tandem repeat unit) are not filled. Every decision is in `<name>.graph_edits.tsv`; the fills are in `OUT_DIR/<name>.gapfills.fa` (named `gapfill_N` in the output AGPs) |
| `--contaminants FILE` | non-target sequences (contaminants, cobionts, organelles): a text file with a sequence name in the first column (comma, space or tab delimited; an optional second column is the class, e.g. `bacteria`, `mitochondrion`, `plastid`). They are tagged `Contaminant` and nothing else, no `HAP1`/`HAP2` (also when under one pixel, in the tool's AGP), never joined, cut, paired, painted or placed as unlocs. `OUT.contaminants.tsv` lists them with their class |
| `--tiara` (with `--fasta`) | classifies the scaffolds up to `--tiara-max-mb` (default 10) with [Tiara](https://github.com/ibe-uw/tiara) and treats non-eukaryotic ones like `--contaminants`. Tiara alone also calls nuclear sequence non-eukaryotic (repeat-rich pieces "prokaryote", a high-copy AT-rich repeat "organelle", NUMT-like pieces "mitochondrion"), so each call needs other evidence: a prokaryote (or plastid) call counts only if Hi-C sets it apart (its second-best contact with another called scaffold at least twice its peak contact anywhere on a chromosome, best 1 Mb window: the pieces of one foreign genome contact each other, a nuclear piece one place on its chromosome; contact among the calls only supports them); a mitochondrion call only at organelle copy number (coverage ≥ 5× haploid); an organelle of unknown kind never; and no call that aligns to a chromosome along ≥ 20% of its length. Calls without that evidence are rows `flag_contaminant` for review, not tagged. About 3–4 min for the small scaffolds of a 3–4 Gb assembly (cached in `--work-dir`) |
| `--hap-regex` | how to read the haplotype from a scaffold name (default `^(H\d+)_`) |

Run `pretext-autocurate run --help` for every option.

## Reporting problems

Please send: `autocurate.run.txt`, the decisions and review tables, and, if you corrected
the result, your savestate (saved from PretextView on top of the tool's). Corrected
savestates are the best way to improve the rules.

## Other commands

    pretext-autocurate info MAP.pretext                 # header, tracks, scaffolds
    pretext-autocurate dump STATE --map MAP.pretext     # a savestate in readable form
    pretext-autocurate convert --map MAP --state STATE --to agp|tsv|state --out FILE
    pretext-autocurate asm --map MAP --state STATE --fasta ASM.fa --out-dir final [--breaks precur/autocurate.breaks.tsv]

`asm` builds the assembly files from a savestate (typically yours, after review) without
exporting from PretextView: it writes the AGP PretextView would export and runs
pretext-to-asm from [sanger-tol/agp-tpf-utils](https://github.com/sanger-tol/agp-tpf-utils)
on it (FASTA and AGP per haplotype, chromosome lists, unlocs, haplotigs; Tree of Life
naming). `run --pretext-to-asm` does the same for the pre-curation itself, in
`OUT_DIR/assembly/`. Options: `--asm-name` (default: the FASTA name up to its first dot),
`--asm-args "..."` (passed on, e.g. `--autosome-prefix SUPER_`). pretext-to-asm is not
bundled: `pip install git+https://github.com/sanger-tol/agp-tpf-utils`; it writes a `.fai`
next to the FASTA if there is none. Scaffolds too short to have a pixel in the map are not
in a PretextView export; the tool's AGP lists them at the end, tagged `HAP1`/`HAP2` by
their `H1_`/`H2_` prefix, so they join their haplotype instead of forming separate `h1`
and `h2` assemblies.

## License

Free for noncommercial use, including by universities and public research institutes:
[PolyForm Noncommercial License 1.0.0](LICENSE.md). Commercial use requires a license from
CNAG (https://www.cnag.eu); see [NOTICE](NOTICE).

## Development

    python3 -m unittest discover -s tests

Most tests need the development maps and curated savestates in `seed/fixtures/` (not in
the repository: large and unpublished) and skip themselves without them. The rest run
anywhere in about a second, among them `tests/test_synthetic.py`, which writes small
synthetic `.pretext` maps of a toy genome (`pretext_autocurate/synth.py`: deliberate
misjoins, two haplotypes) and runs the whole pass and the command line on them, and
`tests/test_steps.py`, which tests each step of the allhaps pass on its own (sites, pairing,
merges, small scaffolds, arrangement, the other copy following the reference, layout) with
synthetic alignments (`synth.toy_blocks`) and `curate_allhaps(..., stop_after=step)`. GitHub
Actions runs those on every push (`.github/workflows/tests.yml`).

The allhaps pass is the package `pretext_autocurate/allhaps/`: `driver.py` runs the steps
listed in `STEPS` over a shared state (`state.Run`), and every level-dependent threshold is
named in `common.THRESHOLDS`. A refactoring must leave the outputs unchanged:
`python3 tools/golden.py check` compares them byte for byte with the recorded runs
(`record` first; needs the fixtures).

    python3 tools/score_all.py                   # run every test genome and score it

scores the tool against the curators' savestates on every genome listed in
`tools/score_cases.json` (junctions, breaks, chromosome grouping, tags) and prints one
table; add a genome by adding an entry there.

- `CHANGELOG.md`: what changed per version
- `docs/HANDOFF.md`: state of the code, conventions, open problems
- `docs/PROGRESS.md`: dated log of methods tried and results on the test genomes
- `docs/FORMATS.md`: `.pretext` and savestate formats (PretextView commit in
  `docs/pretextview_commit.txt`)
