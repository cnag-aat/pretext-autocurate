# Handoff: state of pretext-autocurate (written 2026-11-13)

Read this first, then `PLAN.md` section 0 (the curator's priorities and conventions; the
rest of PLAN.md is the original plan and partly out of date), then the tail of
`docs/PROGRESS.md` (dated log with every result and number).

## What exists (Python, `src/pretext_autocurate/`)

| Module | Purpose |
|---|---|
| `pretextio.py` | `.pretext` reader, mip-level matrix, header hash, `SegmentedMap` (view of the map split at cuts; detectors run on it) |
| `savestate.py` | PretextView savestate read/write, blank template (`template.bin`), nEdits%8 quirk |
| `mapstate.py` | pixel simulator of PretextView edits, `plan_edits` (target permutation -> edits) |
| `layout.py` | savestate + map -> pieces, AGP, TSV |
| `signal.py` | decay curve, full-resolution `strip`/`block` readers, `insulation` (contact across each boundary / expectation), coverage ratio |
| `misjoin.py` | `detect_misjoins` (cuts), `weak_junctions` (candidate breakpoints) |
| `haplotig.py` | `detect_haplotigs` (whole small scaffolds, tiers high/medium/low), `detect_contig_haplotigs` (gap-to-gap contigs) |
| `place.py` | `place_scaffolds`: anchor (chromosome) for each small piece; localise only at reduced-contact junctions |
| `join.py` | `find_joins`: mutual-best end contact between pieces > 150 px, telomeric ends never joined, chains with orientation |
| `curate.py` | `curate_single`: the whole single-haplotype pass -> decisions + savestate |
| `allhaps.py` | `Level` (the one conservative-to-aggressive setting), `arrange_pieces` (order/orientation within a copy), and the allhaps pass: telomere sites (cut / invert_tip / keep), homolog pairing (Hi-C or sequence), same-prefix pairs, small scaffolds as unlocs (`place_small`), parallel layout, HAP tags |
| `homology.py` | minimap2 between haplotypes (`align_haplotypes`) and of small scaffolds against everything (`align_small`), PAF reading, external tool lookup |
| `busco.py` | BUSCO lineage (path or name, opt-in download), miniprot, gene counts per haplotype and scaffold |
| `cli.py` | `info`, `dump`, `convert`, `haplotigs`, `run` |

Order inside `curate_single`: misjoin cuts -> contig haplotigs (cut at gaps) -> weak
junctions cut tentatively and kept only if pieces fit better elsewhere -> `SegmentedMap`
-> joins/chains -> haplotig calls -> placement -> layout per chain (pieces, insertions,
then unlocs, then haplotigs at the chromosome end) -> `plan_edits` -> paint/tags/waypoints.

## How to run and check

    PYTHONPATH=src python3 -m pretext_autocurate.cli run --map MAP.pretext --out-dir out
    python3 -m unittest discover -s tests            # about 10 minutes; tests.test_allhaps alone about 5
    ... run --mode allhaps --map MQ0 --map MQ10 --paf H2_vs_H1.paf --small-paf small_vs_all.paf --busco-table H1=... --busco-table H2=...
    PYTHONPATH=src python3 -m pretext_autocurate.cli run --mode allhaps --phasing hic --map MQ0 --map MQ10 --out-dir out_bird
    python3 tools/compare_layouts.py MAP TOOL.savestate CURATOR.savestate

Outputs: `<map>.savestate_autocurate`, `autocurate.decisions.tsv` (every action with its
evidence), `autocurate.agp`. pytest is not installed; tests are unittest.

## Fixtures (`seed/fixtures/`)

- xgTylRafi (gastropod, single haplotype, MQ0 only, 65536 px): the development genome.
  `.savestate_4` = curator's original curation; `.savestate_3_curator_on_autocurate` =
  curator's 30 edits on top of an earlier tool output; `.savestate_4_curator_reviewed_tags`
  = the same plus reviewed Haplotig/Unloc tags. **The last is the reference for scoring.**
  Its paint ids split one chromosome (scaffold_3) into two ids; ignore that.
- bAegFun (boreal owl, *Aegolius funereus*; allhaps H1_/H2_, MQ0 + MQ10, 65536 px) with
  curator `.savestate_6`. "The bird" in the notes is this genome.
- qqTroRaco (Troglobisium racovitzai, a pseudoscorpion; allhaps, MQ0 only) with curator
  `.savestate_final` and `.savestate_curator_1_on_autocurate_level1` (the curator's
  corrections of the tool's map-only level 1.0 output, 2026-11-14), the FASTA, the
  H2-against-H1 PAF and per-haplotype BUSCO tables. Added 2026-11-13 as the harder allhaps test.
- CAE71 savestates only; the CAE71 maps are not in the repo.
- Wanted (curator, 2026-11-13): an allhaps test that needs extensive curation, as the
  gastropod does. The bird is too well behaved for that. The curator will check whether
  the other haplotype of xgTylRafi can be recovered (S3 restore, read download and
  remapping: slow, so not expected soon) and will look for another genome among the lab's
  recent assemblies. Until then allhaps cuts/joins/placement are tested on the bird only.

## Conventions agreed with the curator (also in PLAN.md section 0)

- The curator's savestates are a guide, not exact truth. Flag when marginal, never hide an edit.
- Cuts preferably at gaps; a cut inside a contig needs a stated reason.
- Waypoints only for positions inside contigs (cuts, weak junctions, flagged joins,
  proposed insertion loci). Standalone pieces needing a look get the `vertPaint` tag.
- Unlocs and haplotigs go at the end of their chromosome, painted with it.
- Telomeric/subtelomeric end contacts are common at MQ0 and are not join evidence.
- Telomere tracks are not reliable in every species: telomere signal is supporting
  evidence, never a requirement. Vertebrates and butterflies are easy; other arthropods,
  other invertebrates and other phyla are the hard cases to design for.
- One setting, `--level` (0 conservative .. 1 aggressive), moves every act-or-flag
  threshold; new thresholds must be read through `allhaps.Level`. `--sweep` writes one
  savestate per level for the curator to choose from, so keep a level fast (shared caches).
- Equal score inverted or not = probably heterozygous inversion: flag, do not edit.
- Under `--phasing hic` large scaffolds never mix prefixes within a superscaffold, but a
  small piece (unloc) may go to the other prefix's copy. `trio`: no mixing at all.
- Sex chromosomes: PAR is not a haplotig; in allhaps both X/Y (Z/W) go to HAP1.
- allhaps with two haplotypes: pairs side by side, better one (fewer gaps, more complete,
  telomeres) labelled HAP1; more than two: HAP1, HAP2, then Haplotig/FalseDuplicate.

## allhaps spec from the curator (details in PLAN.md section 0)

- Step 1: put the two haplotypes of each chromosome side by side, same orientation
  (parallel, not an X), tag HAP1/HAP2, identify sex chromosomes if possible.
- `--phasing none|hic|trio` controls haplotype reassignment: `trio` none at all; `hic` no
  mixing of `H1_`/`H2_` within a superscaffold but whole superscaffolds may swap (the
  bird); `none` free mixing with HAP1 best and balanced. The flag exists in `run` but
  **does nothing yet**.
- A collapsed chromosome (low heterozygosity) has no partner: report it, do not force a pair.
- One pass preferred. Pair homologs on MQ0; cuts/joins/placement within each haplotype on
  the high-MQ map at an adaptive resolution. A second pass (haploid maps, single mode) is
  optional for hard cases and needs no new code.

## Verified in PretextView by the curator

Savestates from the blank template, break edits, paint, tags and waypoints all load.
Not verified: the AGP format (written from memory of PretextView's style), the CAE71
`claude_scaffold1` fixtures, and tool savestates generated after 2026-11-12.

## Current numbers on xgTylRafi (against the reviewed curation)

28 cuts (13 flagged), 13 joins (4 flagged), 8 insertions, 128 unlocs, 86 haplotigs
attached. Confident haplotig calls: 83 of 84 right; 74 borderline pieces highlighted.
Tool junctions in the curation: 23 of 29. Curator junctions reproduced: 33 of 45.

## To think about next (curator, end of session 2026-11-13)

- **Coverage for sex chromosomes.** The coverage track often separates sex chromosomes
  when they are not homomorphic. The tool uses it in two places only (a pair with both
  copies at about half coverage; unpaired scaffolds of 5 Mb or more at half coverage).
  Use it more: for hemizygous chromosomes built from several scaffolds (Troglobisium X is
  five pieces from two scaffolds and cannot be chained by two-haplotype evidence), for
  the boundary between a sex chromosome and an autosome joined in one scaffold
  (H1_scaffold_1, H2_scaffold_4), and to label X/Y or Z/W.
- **Sex chromosome identification, references from the curator** (not read yet):
  - homomorphic sex chromosomes: SCINKD, https://github.com/DrPintoThe2nd/SCINKD, paper
    https://academic.oup.com/mbe/article/43/4/msag067/8664413
  - heteromorphic sex chromosomes: HalfDeep, https://github.com/makovalab-psu/HalfDeep
  See what signal each uses and whether it can be computed from what this tool already
  has (coverage track, haplotype alignments, the FASTA).
  If k-mer counting is needed for any of this, the curator notes that FASTK is faster
  than meryl: prefer FASTK.
- **FASTGA** (https://github.com/thegenemyers/FASTGA) for the pairwise alignment between
  haplotypes: benchmark against minimap2 asm5. On Troglobisium minimap2 took 94 minutes
  and 10 GB for H2 against H1, and the small-against-all alignment does not fit in this
  machine's memory. Compare run time, memory, aligned fraction and what the tool does
  with each alignment (pairs, orientation guide).
- **hapscotch** by Chenxi Zhou (https://github.com/c-zhou/hapscotch), who uses FASTGA:
  possibly related to what this tool attempts. It starts from contigs, not scaffolds, but
  some ideas look similar. Read the code and compare approaches (homology between
  haplotypes as scaffolding evidence, how joins and switches are decided) before building
  more of the alignment-based steps here. Neither repository has been looked at yet.

## Known gaps, in rough priority

0. **allhaps with sequence.** With `--paf` / `--fasta` the alignment now leads: breaks of
   colinearity at gaps are cut/split candidates (`alignment_breaks`), each copy's pieces are
   ordered and oriented along the other copy (`_project`), and loose scaffolds go into the
   chain where they align. End contact is the check and Hi-C alone is the fallback
   (`merge_pairs`, `arrange_pieces`). See the last PROGRESS entries for scores against the
   two Troglobisium curator states. Still to do: which copy is the reference when the two
   disagree (a contig should outrank a scaffold join), merges checked by alignment,
   hemizygous chromosomes built from several scaffolds, two haplotypes joined in one
   scaffold, pieces under 4 px.
1. **allhaps mode is partial** (`allhaps.py`: telomere-at-gap sites decided as cut /
   invert_tip / keep with dual-scaffold logic, pairing, parallel layout, HAP tags, sex
   chromosome candidates). Sequence homology is implemented (`homology.py`; `run --fasta` / `--paf`; minimap2 and
   miniprot live in conda env `pretext-autocurate`) and is the better source for pairing.
   BUSCO (`--busco-table`, `--busco-lineage`) is wired into `run` and is a tie-breaker for
   HAP1. Small scaffolds: `--small-paf` / `align_small`, unlocs where sequence and Hi-C
   agree, same-prefix pairs (see the last PROGRESS entry). Next: position of a piece along
   its homolog (filling a hole in one copy, as the curator did with H1_scaffold_44 and
   H2_scaffold_39), joins, Haplotig/FalseDuplicate tags for third copies, the three kept
   telomere sites the curator cut (H1_scaffold_5, H2_scaffold_4, H2_scaffold_36).
   Decided by the curator (2026-11-13): a small piece may go to the other prefix's copy
   under `--phasing hic` (never under `trio`); `place_small` chooses the copy.
   File mip levels saturate: use `Signal.binned` for anything quantitative at coarse
   resolution.
2. Chromosome 2 of xgTylRafi is two groups (a join flagged because of an interstitial
   telomere signal and a tie). Grouping pieces of one chromosome whose ends do not touch.
3. Small local rearrangements with weak orientation signal; nested inversions (iterate).
4. Second piece at an already-used junction; pieces under 15 px.
5. `Contaminant` tag (scaffold_18 is correctly left unplaced but not tagged); circular
   bacterial pattern. p-to-q orientation of chromosomes.
6. Haplotigs with normal coverage (4 missed outright); contig-level haplotigs (2 of 14).
7. Thresholds are tuned on one genome. They are map-relative, but untested elsewhere.
8. Speed: `curate_single` reads full-resolution tiles repeatedly (about 60 s per run).

## allhaps: what is known and what to do first

- Everything in single mode assumes one copy per locus. On an allhaps MQ0 map homologous
  scaffolds show an off-diagonal line (like a haplotig), and placement/joins see contact
  with both haplotypes. `detect_haplotigs` must not be run as is.
- The bird has an MQ10 map: at MQ10 homolog lines vanish and only real contacts remain.
  Plan: do cuts, joins and placement on MQ10 (or MQ0 restricted to the same haplotype by
  the `H1_`/`H2_` name prefix), and use MQ0 minus MQ10 to pair homologs.
- Curator's bird result: 82 painted units = 41 pairs incl. Z and W (both HAP1); 9 cuts,
  all at gaps (H1_scaffold_1 and H2_scaffold_1 are each two chromosomes, ~134 + ~130 Mb);
  6 unlocs; HAP1/HAP2 assigned per pair by quality, not by name prefix.
- First measurements (last entry of PROGRESS.md): the single-mode cut detector fails on
  both bird maps. MQ0 hides misjoins (contact across them is 0.24-0.67); MQ10 is too sparse
  at full resolution (1133 false cuts) but has the right signal. Make the working
  resolution adaptive (mip level chosen from the decay curve) before anything else.
