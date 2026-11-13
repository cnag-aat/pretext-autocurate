# pretext-autocurate

Automatic first-pass curation of Hi-C contact maps, straight from PretextMap `.pretext`
files. It writes a PretextView savestate (plus an AGP and a table explaining every decision)
that a curator opens, reviews and finishes. The aim is to take over the routine part of
manual curation, not to replace the curator: marginal calls are flagged, never hidden.

**Status: research prototype.** Developed on two genomes (a gastropod, single haplotype; a
bird, two haplotypes). Thresholds are relative to each map's own signal but have not been
tested widely. Check the output before trusting it.

## What it does

Single-haplotype maps (`--mode single`):

- cuts misjoins, preferably at gaps, and re-joins and orients the large pieces
- handles rearrangements inside a scaffold by cutting tentatively at weak junctions and
  keeping a cut only if the pieces fit better another way
- assigns small scaffolds to a chromosome; inserts them at a gap when position and
  orientation are clear, otherwise attaches them as `Unloc` at the chromosome end
- calls haplotigs (off-diagonal line of contact plus reduced coverage); borderline ones are
  highlighted, not tagged
- never joins at an end that carries a telomere signal

Two-haplotype maps (`--mode allhaps`, scaffolds named `H1_...` / `H2_...`):

- pairs homologous scaffolds and lays each pair side by side in the same orientation,
  tagged `HAP1` / `HAP2`
- decides telomere-bearing gaps inside a scaffold as a misjoin (cut, also at the matching
  gap of the other haplotype), an inverted end piece (inverted in place), or neither
- flags a hemizygous, lopsided pair as possible sex chromosomes (both `HAP1` by convention)
- leaves scaffolds that look repetitive (high coverage) unplaced, and says why
- optional evidence: sequence homology from the assembly (minimap2) and BUSCO genes

Not done yet in allhaps mode: joins, unlocs, haplotigs, same-haplotype pairs.

## Install

    conda env create -f environment.yml
    conda activate pretext-autocurate

or, with Python >= 3.10 and numpy already available: `pip install -e .`
minimap2 and miniprot are only needed for `--fasta` and `--busco-lineage`. They are looked
for on `PATH`, then in a conda environment named `pretext-autocurate`.

## Usage

    pretext-autocurate run --map MAP.pretext [--map HIGH_MQ_MAP.pretext] --mode single|allhaps --out-dir DIR [options]
    pretext-autocurate info MAP.pretext
    pretext-autocurate dump STATE --map MAP.pretext
    pretext-autocurate convert --map MAP.pretext --state STATE --to agp|tsv|state --out FILE
    pretext-autocurate haplotigs MAP.pretext --out-dir DIR

The first `--map` is analysed and should be the MQ0 map (it must carry `gap`, `coverage`
and, ideally, `telomere` tracks). Further maps of the same assembly get their own savestate
with the same edits; in allhaps mode the second map is used as the high-MQ map.

Outputs in `--out-dir`:

| File | Content |
|---|---|
| `<map name>.savestate_autocurate` | load in PretextView after loading the map (one per `--map`) |
| `autocurate.decisions.tsv` | every cut, join, insertion, pairing and flag, with its evidence |
| `autocurate.agp` | the proposed assembly (format not yet checked against PretextView's own export) |

In the savestate: waypoints mark positions inside contigs to look at (cuts, weak junctions,
flagged joins, proposed insertion points); standalone pieces that need a look carry the
`vertPaint` tag, which PretextView draws as a vertical band.

### Examples

Single haplotype:

    pretext-autocurate run --map asm.mq0.pretext --mode single --out-dir out

Two haplotypes, Hi-C evidence only:

    pretext-autocurate run --mode allhaps --phasing hic \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext --out-dir out

Two haplotypes, with the assembly (haplotypes are aligned with minimap2; cached in `out/work`):

    pretext-autocurate run --mode allhaps --phasing hic \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext \
        --fasta asm.diploid.fa.gz --threads 8 --out-dir out

With an existing alignment and existing BUSCO results (one table per haplotype; the prefix
is added to the sequence names):

    pretext-autocurate run --mode allhaps --phasing hic \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext \
        --paf H2_vs_H1.paf \
        --busco-table H1=hap1.full_table.tsv --busco-table H2=hap2.full_table.tsv --out-dir out

Letting the tool run the BUSCO protein search (slow; lineage directory, or a name with
`--busco-download` to fetch it from busco-data.ezlab.org):

    pretext-autocurate run --mode allhaps ... --fasta asm.diploid.fa --busco-lineage /path/to/aves_odb12

### Options worth knowing

| Option | Meaning |
|---|---|
| `--phasing trio\|hic\|none` | allhaps: how freely sequence may change haplotype. `trio`: never. `hic` (default): whole chromosomes may swap so the better copy is HAP1, but `H1_`/`H2_` are never mixed within one. `none`: as `hic` for now |
| `--fasta`, `--paf` | allhaps: sequence homology between haplotypes (run minimap2, or use an alignment) |
| `--level`, `--sweep` | allhaps: one setting from 0 (conservative) to 1 (aggressive) that moves every act-or-flag threshold; `--sweep 0,0.25,0.5,0.75,1` writes a savestate per level to choose from |
| `--small-paf`, `--no-small-align` | allhaps: alignment of the small scaffolds against the whole assembly (run from `--fasta` by default, slow the first time); used to attach small scaffolds as unlocs and to find homologs within one haplotype file |
| `--busco-table`, `--busco-lineage`, `--miniprot-paf` | allhaps: BUSCO genes as a guide for the HAP1 choice and for false-duplication candidates. Entirely optional; nothing acts on a single gene |
| `--apply-inserts` | single: apply every proposed insertion, not only the clear ones |
| `--threads`, `--work-dir` | for the external aligners |

Run `pretext-autocurate run --help` for the full list.

## Tests

    python3 -m unittest discover -s tests

The tests need the development maps and curated savestates in `seed/fixtures/`, which are
not in this repository yet (too large, and unpublished). Without them most tests are
skipped and the savestate round-trip test fails. A small public test set is planned.

## Documentation

- `docs/HANDOFF.md`: state of the code, module by module, conventions, known gaps
- `docs/PROGRESS.md`: dated log of methods tried and results on the two genomes
- `docs/FORMATS.md`: `.pretext` and savestate binary formats
- `PLAN.md`: section 0 holds the curation conventions the tool follows; the rest is the
  original plan and partly out of date

Requires maps made with PretextMap and savestates for PretextView
(sanger-tol/PretextView; the commit the formats were read from is in
`docs/pretextview_commit.txt`).
