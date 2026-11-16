# Handoff: state of pretext-autocurate (v1.1.0, updated 2026-11-14)

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
| `tracks.py` | evidence tracks: bedgraph/BED or map tracks -> per-pixel arrays (`to_pixels`, `from_map`, `load`); runs telo_scan.py (`run_telo_scan`) and Red (`run_red`, Intel build under Rosetta); BUSCO duplicate-pair / missing-gene tracks (`busco_tracks`); map copies with tracks via PretextGraph (`add_tracks`) |
| `sexchrom.py` | sex chromosomes from the coverage track and the alignments (`analyse`): sex-linked stretches (half coverage, not aligned, or aligned into one: the other sex chromosome), PARs (normal-coverage end with chromosome-level contact, aligned or MQ0-linked to the other sex chromosome's end), another chromosome joined to a sex chromosome; `same_hap_candidates` for `homology.align_within`; `sex_chromosome` / `par` bedgraphs (`write_tracks`) |
| `vendor/telo_scan.py` | the curator's telomere scanner (cnag-aat/scripts 666994c), bundled unchanged apart from a provenance header |
| `cli.py` | `info`, `dump`, `convert`, `haplotigs`, `run` (mode detection, `--level` / `--sweep`, review table, run info, track options) |

Order inside `curate_single`: misjoin cuts -> contig haplotigs (cut at gaps) -> weak
junctions cut tentatively and kept only if pieces fit better elsewhere -> `SegmentedMap`
-> joins/chains -> haplotig calls -> placement -> layout per chain (pieces, insertions,
then unlocs, then haplotigs at the chromosome end) -> `plan_edits` -> paint/tags/waypoints.

## How to run and check

    PYTHONPATH=src python3 -m pretext_autocurate.cli run --map MAP.pretext --out-dir out
    python3 -m unittest discover -s tests            # about 12 minutes; tests.test_allhaps alone about 6
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
- A second copy of sequence that its own haplotype already has is tagged `FalseDuplicate`
  and kept out of the chain, whatever produced it (curator, 2026-11-14).
- Telomeric/subtelomeric end contacts are common at MQ0 and are not join evidence.
- Telomere tracks are not reliable in every species: telomere signal is supporting
  evidence, never a requirement. Vertebrates and butterflies are easy; other arthropods,
  other invertebrates and other phyla are the hard cases to design for.
- One setting, `--level` (0 conservative .. 1 aggressive), moves every act-or-flag
  threshold; new thresholds must be read through `allhaps.Level`. `--sweep` writes one
  savestate per level for the curator to choose from, so keep a level fast (shared caches).
- Equal score inverted or not = probably heterozygous inversion: flag, do not edit.
- When the evidence is clear (e.g. telomere strands), act and explain the evidence in the
  decisions table; the curator can reverse it with everything in front of them (curator,
  2026-11-14, on H2_scaffold_2).
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
  bird), and small pieces may go to the other prefix's copy; `none` free mixing with HAP1
  best and balanced. `trio` and `hic` are implemented; `none` currently behaves like `hic`.
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
- **Contigs outrank scaffold joins** (curator, 2026-11-14). The scaffolds are from YaHS,
  which is not perfect; contigs are built from unitigs and are not perfect either. When
  the two copies of a chromosome disagree, the one that is contiguous sequence at that
  point is the reference (not yet implemented: the reference is the copy in fewer pieces).
  For a future version the curator can provide the hifiasm unitig GFAs of hap1 and hap2;
  sequence names differ from the assembly, so matching needs a graph alignment.
- **Repeats** (curator, 2026-11-14). Repetitive sequence is already "duplicated"
  sequence: variable amounts of an array in the two haplotypes should not weigh in a
  duplicate or placement decision; judge on the non-repetitive part. Repeats can be
  inferred from contacts, from the hap1-to-hap2 alignment (done for small scaffolds: loose
  alignments against a dense one), and from a repeat track: Red
  (https://github.com/BioinformaticsToolsmith/Red) is fast and its output can be added to
  the map with PretextGraph (https://github.com/sanger-tol/PretextGraph). To do: read a
  repeat track when the map has one (maps already carry gap, telomere, coverage), offer to
  run Red on `--fasta`, and use it wherever sequence is weighed (pairing, duplicates,
  unlocs, alignment breaks).
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

## Next, in the curator's order (2026-11-14)

1. ~~Red repeat mask as a map track (PretextGraph) and as evidence wherever sequence is
   weighed~~ done in 1.1.0: `--red`, `--repeats`, `--add-tracks`. On the owl Red changes no
   decision (sequence evidence is already judged on unique alignments); it matters for
   Hi-C-only placement and for display.
2. ~~telo_scan.py (cnag-aat/scripts) for telomeres from `--fasta`, both strands~~ done in
   1.1.0 and bundled (`--telomere-motif`).
3. ~~Benchmark minimap2 against FASTGA~~ done 2026-11-14 (PROGRESS): FastGA needs over
   50 GB of free disk for its index on a 4.4 Gb assembly and gives no mapping quality.
   Curator's decision: stick with minimap2 (the tool runs on a MacBook).
4. ~~Read hapscotch (Chenxi Zhou) and compare approaches~~ done: `docs/RELATED.md`. The curator trusts the work of
   Richard Durbin's group and of anyone who has worked with him: treat their methods as
   references to learn from.

## Suggestions for further development (2026-11-14, after v1.1.0)

Proposals, not curator decisions; roughly in order of expected value for curators. Items
the curator has already asked for are marked (curator).

Curation quality, two-haplotype mode
1. ~~**Contig outranks scaffold join** (curator)~~ done 2026-11-15 (`_contig_conflicts`,
   per pair, see PROGRESS). Still open: a pair where each copy is right in a different
   place gets one reference only (flagged `arrangement_differs`); a merged layout taking
   each disagreement from the contiguous copy would settle those.
2. ~~**Hemizygous and sex chromosomes from coverage** (curator)~~ done 2026-11-15
   (`sexchrom.py`, see PROGRESS): sex-linked stretches, PARs, coverage steps where another
   chromosome is joined, hemizygous pieces chained by contact, `sex_chromosome` / `par`
   tracks, tags with `--sex-system` (XY, ZW, X0, Z0, X1X2Y, XY1Y2, Z1Z2W, ZW1W2); chained
   pieces split at their kept gaps so they can interleave (the Troglobisium X as the
   curator laid it out). Still open: only gaps with a `keep` site are split points, so a
   piece belonging in a gap the cut detection did not look at stays whole; phase 1 searches
   at most 9 pieces, and extras beyond that go in whole. Which copy is which without `--sex-system` is left to the
   curator (needs species knowledge). Homomorphic sex chromosomes, and any sample of the
   homogametic sex, show no hemizygous stretch and are paired as autosomes (the
   `sex_summary` row says so); SCINKD-style k-mers would need reads of both sexes.
3. ~~**Shared Hi-C backbone** (from hapscotch)~~ measured 2026-11-16, not built
   (`tools/backbone_gaps.py`, see PROGRESS): the two copies give nearly the same contact
   across each gap (correlation 0.97-0.99 on MQ0 maps), so pooling does not separate the
   curator's cut and kept gaps any better (AUC 0.962 own copy, 0.969 pooled, Troglobisium).
   The weak gaps are those with no mapped homolog (AUC 0.69 / 0.60), where there is
   nothing to pool; that is where gap decisions need another signal.
4. **Merges checked by alignment.** `merge_pairs` decides on end contact; when `--paf`
   is given, require the two pairs to be colinear pieces of one homolog, as `_project`
   already does for order.
5. **Two haplotypes in one scaffold**: a scaffold that aligns to itself (an off-diagonal
   line inside its own block at MQ0, plus a self-alignment) should be split and the halves
   paired.
6. **`--phasing none` as an optimisation** (from hapscotch): assign pieces to haplotypes
   to minimise inter-haplotype contact within a chromosome; a principled way to decide
   cases like owl H1_scaffold_31 / H1_scaffold_44, which now rest on fixed rules.
7. **Small pieces** (under 4 px in allhaps, 15 px in single mode): read full-resolution
   strips for them instead of the working mip level.
8. **Contaminant tag**: on the high-MQ map with the repeat track, a scaffold with no
   contact to anything, no homolog and no repeats is a contaminant candidate; tag it
   rather than leave it unplaced and silent.

Evidence
9. **Unitig GFAs** (curator, hifiasm hap1/hap2): independent evidence for where a contig
   is trustworthy. Needs a graph alignment (names differ from the assembly); minigraph or
   GraphAligner are candidates.
10. **FastGA on a cluster**: keep minimap2 as the default (curator), but accept a FastGA
    PAF via `--paf` once the tool has its own uniqueness filter (best chain per query
    stretch), since FastGA gives no mapping quality.
11. **telo_scan threshold**: optionally drop windows below telo_scan's own score (100)
    before PretextGraph; removes weak interstitial runs (measured on the owl, optional).

Validation, so the rules keep improving
12. **One scoring script over all fixtures** (`tools/score_all.py`): for every genome
    with a curator savestate, run the tool at each level and report junction precision /
    recall, cuts, pairs and tags against the curation, as a table to keep in PROGRESS.
    Today the numbers come from separate tests and ad hoc runs.
13. **More genomes from the team, especially hard taxa** (not vertebrates, not
    butterflies): every corrected savestate (tool output + curator's edits) is a new
    test. Ask curators to send `autocurate.run.txt` with it (README already says so).
14. **AGP** checked against PretextView's own export (still unverified).

Engineering
15. **Split `allhaps.py`** (about 1,700 lines; `curate_allhaps` alone about 600): telomere
    sites, small scaffolds, merging/arrangement and the driver into separate modules, each
    step returning its rows. Easier to test step by step.
16. **Tests that run without the fixtures**: a few small synthetic `.pretext` maps (or
    pickled `SegmentedMap`s) so GitHub Actions can run on every push; keep the 12-minute
    fixture suite for releases.
17. **`--level` in single mode**: its thresholds are fixed; reading them through `Level`
    would let a sweep cover both modes.
18. **Packaging**: a bioconda recipe (with minimap2 and miniprot as dependencies), after
    checking that bioconda accepts a noncommercial license; the Intel-only Red / PretextGraph stay optional.

## Known gaps, in rough priority

0. **allhaps with sequence.** With `--paf` / `--fasta` the alignment now leads: breaks of
   colinearity at gaps are cut/split candidates (`alignment_breaks`), each copy's pieces are
   ordered and oriented along the other copy (`_project`), and loose scaffolds go into the
   chain where they align. End contact is the check and Hi-C alone is the fallback
   (`merge_pairs`, `arrange_pieces`). See the last PROGRESS entries for scores against the
   two Troglobisium curator states. The reference copy is the one whose arrangement the other
   copy's contigs contradict less (2026-11-15). Hemizygous pieces in no pair join the sex
   chromosome they contact (`sexchrom`, 2026-11-15). Still to do: merges checked by
   alignment, two haplotypes joined in one scaffold, pieces under 4 px.
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
