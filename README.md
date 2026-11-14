# pretext-autocurate

Pre-curation of Hi-C contact maps, straight from PretextMap `.pretext` files. It does the
routine first pass of manual curation and writes a PretextView savestate that a curator
opens, reviews and finishes, together with a table explaining every decision and a short
list of the things that need a curator's eye. It does not replace the curator: marginal
calls are flagged, never hidden.

**Version 1.0.** Developed and tested on three genomes: a gastropod (single haplotype),
the boreal owl *Aegolius funereus* (two haplotypes, MQ0 + MQ10) and the pseudoscorpion
*Troglobisium racovitzai* (two haplotypes, heavily rearranged, MQ0 only). Thresholds are
relative to each map's own signal. Always review the output before accepting it.

## Install

    git clone <this repository> && cd pretext-autocurate
    conda env create -f environment.yml
    conda activate pretext-autocurate
    pretext-autocurate --version

Or, with Python >= 3.10 and numpy available: `pip install -e .`. minimap2 (and miniprot for
`--busco-lineage`) are only needed when the tool aligns sequence itself; they are looked
for on `PATH`, then in a conda environment named `pretext-autocurate`.

## Quick start

Two haplotypes (scaffolds named `H1_...` / `H2_...`), the recommended way: give the maps
and the assembly they were made from.

    pretext-autocurate run \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext \
        --fasta asm.diploid.fa.gz --threads 16 --out-dir precur

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
| `autocurate.agp` | the proposed assembly (format not yet checked against PretextView's own export) |
| `autocurate.run.txt` | version, command line and every option, for reproducibility |
| `work/` | cached alignments (reused on the next run) |

In the savestate:

- **Paint**: one colour per chromosome copy; unlocs and false duplicates are painted with
  their chromosome and sit at its end.
- **Tags**: `HAP1` / `HAP2` on everything in two-haplotype mode; `Unloc`, `Haplotig`,
  `FalseDuplicate`; sex chromosome candidates are both `HAP1`.
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
  the alignment is available, by end-to-end contact otherwise; pairs that are parts of one
  chromosome are merged; scaffolds that align into a copy are joined where they align.
- **Small scaffolds**: attached as `Unloc` where sequence and Hi-C agree, or joined when
  they continue a copy; a second copy of sequence its haplotype already has is tagged
  `FalseDuplicate`; repeat-rich scaffolds are judged on their non-repetitive part.
- **Sex chromosomes**: a pair or scaffold at half coverage is reported as a candidate, both
  copies `HAP1`.
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

## Development

    python3 -m unittest discover -s tests

The tests need the development maps and curated savestates in `seed/fixtures/` (not in the
repository: large and unpublished); without them most tests are skipped.

- `CHANGELOG.md`: what changed per version
- `docs/HANDOFF.md`: state of the code, conventions, open problems
- `docs/PROGRESS.md`: dated log of methods tried and results on the test genomes
- `docs/FORMATS.md`: `.pretext` and savestate formats (PretextView commit in
  `docs/pretextview_commit.txt`)
