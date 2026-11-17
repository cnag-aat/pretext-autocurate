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
fraction of the time. `--level` applies to two-haplotype mode.

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
- Flags (does not change): sites where both sides belong together, places where the two
  haplotypes disagree, unloc candidates with one line of evidence, scaffolds with no homolog.

`--phasing` sets how freely sequence may change haplotype: `trio` never; `hic` (default)
whole chromosomes may swap so the better copy is HAP1, and small pieces may go to the other
prefix's copy when the alignment says so; `none` currently behaves like `hic`.

Single haplotype: misjoin cuts (preferably at gaps), joins and orientation of large pieces,
rearrangements inside a scaffold (tentative cuts kept only if the pieces fit better
elsewhere), placement of small scaffolds (inserted at a gap when clear, `Unloc` otherwise),
haplotigs (off-diagonal line plus reduced coverage; borderline ones highlighted).

## Optional evidence

| Option | What it adds |
|---|---|
| `--fasta ASM.fa[.gz]` | aligns the haplotypes (minimap2) and the small scaffolds against everything; homology, order, orientation, duplicates. The first run takes long (both alignments for the 2.7 Gb owl assembly about 15 min; the haplotype alignment alone for a 4.4 Gb assembly about 1.5 h, 10 GB, 6 threads); alignments are cached in `--work-dir`. The small-scaffold alignment indexes the whole assembly: on a machine with little memory, skip it with `--no-small-align` |
| `--paf`, `--small-paf` | use existing alignments instead (H2 against H1, and small scaffolds against the whole assembly) |
| `--busco-table H1=hap1.full_table.tsv --busco-table H2=hap2.full_table.tsv` | existing BUSCO results (one per haplotype, or one on the diploid assembly): tie-breaker for HAP1, false-duplication candidates. Never decisive alone |
| `--busco-lineage`, `--miniprot-paf` | let the tool find BUSCO genes itself (slow) |
| `--repeats BED`, `--red` | a repeat track (e.g. from Red); `--red` runs Red on each haplotype of `--fasta` (about 25 min for a 2.7 Gb assembly) |
| `--telomere-motif MOTIF` | runs the bundled telo_scan.py (CNAG) on `--fasta` for telomeres by strand (e.g. TTAGGG vertebrates, TTAGG most arthropods); tells two chromosome ends joined from an inverted end piece. A few minutes for a 4 Gb assembly. `--telomere-5p/3p` take existing telo_scan bedgraphs; maps that already carry its 5p/3p tracks are read directly |
| `--add-tracks` | write a copy of each map with these tracks, the BUSCO duplicate pairs / missing genes, and the sex chromosome tracks (`sex_chromosome`: 2 hemizygous, 1 the rest of a sex-linked stretch; `par`: 2 PAR, 1 candidate), for viewing in PretextView (needs PretextGraph). The bedgraphs are in `OUT_DIR/tracks/` in any case |
| `--sex-system SYSTEM` | the species' sex chromosome system, if known (`XY`, `ZW`, `X0`, `Z0`, `X1X2Y`, `XY1Y2`, `Z1Z2W`, `ZW1W2`): tags the sex chromosome copies `X`/`Y`, `Z`/`W`, `X1`/`X2` ... |
| `--gfa GFA --yahs-agp AGP` (each once per haplotype) | hifiasm's contig graph (`*.hap1.p_ctg.gfa`; with or without sequence) and the YaHS AGP that built the scaffolds from those contigs (`H1=path`, or a name containing `hap1`). Where the graph links the two contig ends YaHS put next to each other, the join is right: the gap is never cut (row `graph_link`), and `asm` closes it, dropping the overlap the graph reports after checking it in the sequence. Where the only link out of the left contig and the only one into the right go to one small contig that is a scaffold of its own, that contig fills the gap (row `graph_bridge`); `asm` puts it there and closes both sides. `OUT_DIR/<name>.graph_edits.tsv` lists every change |
| `--utg-gfa P_UTG.gfa` (with `--gfa`, `--yahs-agp`) | hifiasm's unitig graph (uncompressed; the `noseq` file is enough for `run` and loads fast). A gap whose two contig ends a unitig path joins is not used to rearrange (row `graph_path`; a tip inversion or a cut on contact alone there is dropped, `graph_keep`; a cut the alignment or a telomere backs stands). In `asm`, contigs that overlap along a single path are joined |
| `--fill-gaps` (`asm`, `run --pretext-to-asm`; needs the unitig GFA with sequence and hifiasm's read depth) | fill gaps from the unitig graph without duplicating what the assembly already has for that haplotype. Every path across the gap is examined: if the contigs overlap or abut along one, the gap is closed; if a path runs through another contig of this haplotype, that contig belongs there: it is put into the gap when it lies wholly on the path and is a scaffold of its own (also at joins the tool or you made; this needs the unitig GFA with sequence but not `--fill-gaps`), and is named in the report otherwise (nothing is copied); otherwise the gap is filled only from a path whose unitigs are at about twice the haploid read depth (collapsed homozygous sequence: two copies in the genome, one in the assembly, in the other haplotype). Haploid paths (the other haplotype's allele), repeats and gaps next to a circular contig (a tandem repeat unit) are not filled. Every decision is in `<name>.graph_edits.tsv`; the fills are in `OUT_DIR/<name>.gapfills.fa` (named `gapfill_N` in the output AGPs) |
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
    pretext-autocurate asm --map MAP --state STATE --fasta ASM.fa --out-dir final

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
misjoins, two haplotypes) and runs the whole pass and the command line on them. GitHub
Actions runs those on every push (`.github/workflows/tests.yml`).

    python3 tools/score_all.py                   # run every test genome and score it

scores the tool against the curators' savestates on every genome listed in
`tools/score_cases.json` (junctions, breaks, chromosome grouping, tags) and prints one
table; add a genome by adding an entry there.

- `CHANGELOG.md`: what changed per version
- `docs/HANDOFF.md`: state of the code, conventions, open problems
- `docs/PROGRESS.md`: dated log of methods tried and results on the test genomes
- `docs/FORMATS.md`: `.pretext` and savestate formats (PretextView commit in
  `docs/pretextview_commit.txt`)
