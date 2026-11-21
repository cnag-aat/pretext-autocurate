# Testing pretext-autocurate on your genomes

Thank you for trying the tool. It has been developed against curations of four genomes so
far, so what we need most is to see how its proposals compare with **your** curation of
genomes it has never seen: other taxa, other assemblers, other map settings. This page says
how to run it, what to look at, and what to send back.

## 1. Install

Follow "Install" in the [README](../README.md) (a conda environment; minimap2 comes with it).
Check the version you test, and give it in every report:

    pretext-autocurate --version        # 1.3.0 or later

## 2. Choose a genome

Best: one you have **already curated by hand** (or are about to), so that your savestate can
be compared with the tool's. Any of these is useful:

- two haplotypes in one map (scaffolds named `H1_...` / `H2_...`), or a single haplotype /
  primary assembly;
- taxa we have not tested: plants, fungi, most invertebrates, fish, reptiles, anything with
  holocentric chromosomes, many microchromosomes, or a hard karyotype;
- maps from any pipeline, as long as they are PretextMap `.pretext` files.

What the tool needs:

| Input | Needed? |
|---|---|
| MQ0 `.pretext` map (first `--map`) with `gap` and `coverage` tracks | yes (older names GAPs / ONTcov are read too; with `--fasta` a missing gap track is made from the runs of N) |
| high-MQ map of the same assembly (second `--map`), two haplotypes | recommended |
| the assembly FASTA (`--fasta`) | strongly recommended: alignments between haplotypes, telomeres, sequence-placed cuts and joins |
| hifiasm GFAs (`--gfa`, `--utg-gfa`) and YaHS AGPs (`--yahs-agp`) | if you have them |
| the telomere motif (`--telomere-motif`) | optional: detected from the scaffold ends if not given |

## 3. Run it at three levels

`--level` sets how much the tool does on its own, from 0 (conservative: only clear edits,
more left for review) to 1 (aggressive). `--sweep` runs several levels in one go and reuses
the alignments:

    # two haplotypes
    pretext-autocurate run \
        --map asm.diploid_mq0.pretext --map asm.diploid_mq10.pretext \
        --fasta asm.diploid.fa --threads 16 --sweep 0,0.5,1 --out-dir precur

    # single haplotype / primary assembly
    pretext-autocurate run --map asm.mq0.pretext --fasta asm.fa --threads 16 --out-dir precur

Add the graph options from the README's "Quick start" when you have hifiasm's files. Note
the wall time and peak memory (for example with `/usr/bin/time -v` on Linux or
`/usr/bin/time -l` on macOS). Alignments are cached in `precur/work/`, so a second run is
much faster.

## 4. Review it in PretextView

Load the map, then the savestate `precur/<map name>.savestate_autocurate` (with `--sweep`,
one per level under `precur/level_*/`). Go through `autocurate.review.tsv` first: it lists
what the tool was unsure of, most important first. Every edit, applied or not, is explained
in `autocurate.decisions.tsv`.

If PretextView's file picker cuts the long names off, copy or link the map under a short
name and re-save the savestate for it (`tools/mark_disagreements.py` below does this for
you).

## 5. Compare with your curation

With your own savestate for the same map:

    python3 tools/mark_disagreements.py asm.diploid_mq0.pretext my.savestate cmp \
        precur/level_0.50/*.savestate_autocurate_level0.50 --label auto

It writes a short-named copy of the map with a track marking each disagreement (extra break,
wrong join, missed break, missed join), the savestates bound to that copy, and
`auto.disagreements.tsv` listing them. Open the copy with the tool's savestate and step
through the waypoints.

## 6. What to send back

Open an issue with the **"Tester feedback"** template on GitHub
(https://github.com/cnag-aat/pretext-autocurate/issues/new/choose), one per genome. Most
useful, in order:

1. **Your curated savestate** for the same map, and whether you curated from scratch or
   corrected the tool's output (the second is biased towards the tool). With it we can score
   the tool the same way as on our own genomes.
2. **The tool's mistakes you would not accept**: scaffold, position (Mb), what the tool did,
   what is right, and why (what you see on the map). A screenshot helps.
3. **What it got right that you would have had to do by hand**, and roughly how much time it
   saved or cost you.
4. Version, `autocurate.run.txt` (command line and options), level(s), map size, taxon,
   one or two haplotypes, wall time and peak memory, and any error message in full.

**Unpublished data:** do not attach assemblies, maps or savestates of unpublished genomes to a
public GitHub issue. Describe the case there and share the files with the CNAG team privately.

## Known limits (1.3.0)

- Thresholds were tuned on four genomes (a gastropod, an insect, a bird, a pseudoscorpion).
- Telomere tracks are unreliable in some taxa; outside vertebrates and butterflies, expect
  more to review.
- A collapsed tandem repeat can look like a haplotypic copy; telling them apart needs read
  phasing, which the tool does not take yet.
- Displaced stretches (a piece with a clear line to another place) are moved in single
  haplotype mode; in two-haplotype mode they are only cut out and placed by the usual steps,
  and this has not yet been seen on real data.
- `--engine adjacency` is experimental; the default rules pass is the one to test.
