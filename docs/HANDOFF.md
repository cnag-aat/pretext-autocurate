# Handoff: state of pretext-autocurate (v1.2.0, updated 2026-11-16)

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
| `allhaps/` | the allhaps pass, a package (split 2026-11-17, `docs/ALLHAPS_SPLIT.md`): |
| `allhaps/common.py` | `Level` (the one conservative-to-aggressive setting), `THRESHOLDS` + `Thresholds` (every act-or-flag threshold by name, at levels 0 / 0.5 / 1), `Pair`, `write_rows` |
| `allhaps/pairing.py` | homolog pairing (`pair_homologs`, Hi-C or sequence), `same_prefix_pairs`, `merge_pairs` (pairs that are parts of one chromosome) |
| `allhaps/sites.py` | `telomere_sites` (cut / invert_tip / keep at gaps), `alignment_breaks` |
| `allhaps/evidence.py` | hifiasm graph evidence at sites (`_graph_evidence`, `_unitig_evidence`) |
| `allhaps/chromjoin.py` | chromosome-level joins: pairs merged on specific, spread-out whole-block contact with repeat bins masked (arms apart across centromeric repeats); guards against cuts the pass made and small clustering chromosomes |
| `allhaps/graphjoin.py` | joins only hifiasm's graph sees, between scaffold ends: `links` (evidence for `merge_pairs`), `place` (the smaller scaffold put at the larger one's end) |
| `allhaps/selfdup.py` | two copies of one region in a scaffold: `self_duplications` (repeat / duplication / haplotypic by read depth), `apply_cuts` (a haplotypic copy cut out at the gaps around it), `rows` |
| `allhaps/small.py` | small scaffolds (`place_small`): unloc, candidate, join, false duplicate |
| `allhaps/arrange.py` | `arrange_pieces` (order/orientation within a copy), `_project` (the other copy follows the reference), `_contig_conflicts` |
| `allhaps/driver.py` | `curate_allhaps`: the steps in `STEPS` (`_setup`, `_small`, `_merge`, `_sex_chain`, `_arrange`, `_reference`, `_follow`, `_rows`, `_layout`), each a function of a shared `Run` (`allhaps/state.py`); `stop_after=` returns the state after a step |
| `homology.py` | minimap2 between haplotypes (`align_haplotypes`) and of small scaffolds against everything (`align_small`), PAF reading, external tool lookup |
| `busco.py` | BUSCO lineage (path or name, opt-in download), miniprot, gene counts per haplotype and scaffold |
| `tracks.py` | evidence tracks: bedgraph/BED or map tracks -> per-pixel arrays (`to_pixels`, `from_map`, `load`); runs telo_scan.py (`run_telo_scan`) and Red (`run_red`, Intel build under Rosetta); BUSCO duplicate-pair / missing-gene tracks (`busco_tracks`); map copies with tracks via PretextGraph (`add_tracks`) |
| `sexchrom.py` | sex chromosomes from the coverage track and the alignments (`analyse`): sex-linked stretches (half coverage, not aligned, or aligned into one: the other sex chromosome), PARs (normal-coverage end with chromosome-level contact, aligned or MQ0-linked to the other sex chromosome's end), another chromosome joined to a sex chromosome; `same_hap_candidates` for `homology.align_within`; `sex_chromosome` / `par` bedgraphs (`write_tracks`) |
| `graph.py` | hifiasm graphs + YaHS AGPs: `Layout` (contigs in scaffolds), `Graph` (contig links, contig-end reads), `UnitigGraph` (unitig links, reads, read depth, sequence by file offset), `junctions` (links, bridges), `unitig_paths`, `choose_fill` (close / place / fill / leave, without duplicating the other haplotype), `insert_bridges` and `place_between` (on the PretextView-style AGP), `close_gaps` (on pretext-to-asm's output), `Faidx` |
| `asm.py` | `build`: PretextView-style AGP -> graph edits -> pretext-to-asm -> gaps closed / filled; `run_pretext_to_asm` |
| `vendor/telo_scan.py` | the curator's telomere scanner (cnag-aat/scripts 666994c), bundled unchanged apart from a provenance header |
| `cli.py` | `info`, `dump`, `convert`, `asm`, `haplotigs`, `run` (mode detection, `--level` / `--sweep`, review table, run info, track and graph options, `--pretext-to-asm`) |

Order inside `curate_single`: misjoin cuts -> contig haplotigs (cut at gaps) -> weak
junctions cut tentatively and kept only if pieces fit better elsewhere -> `SegmentedMap`
-> joins/chains -> haplotig calls -> placement -> layout per chain (pieces, insertions,
then unlocs, then haplotigs at the chromosome end) -> `plan_edits` -> paint/tags/waypoints.

## How to run and check

    PYTHONPATH=src python3 -m pretext_autocurate.cli run --map MAP.pretext --out-dir out
    python3 -m unittest discover -s tests            # 80 tests, about 25-45 minutes on the MacBook
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
4. ~~**Merges checked by alignment**~~ measured 2026-11-16, not built (see PROGRESS): no
   sequence spans a merge join, and the two false Troglobisium merges align end to end
   (their ends are real chromosome ends in both haplotypes) while the real ones have the
   joined ends 1.5-4 Mb inside the homolog, so a colinearity check would keep the false
   and reject the real. Open idea, too few cases yet: ends that coincide in both
   haplotypes as a guard against merging chromosome ends without a telomere signal.
5. ~~**Two haplotypes in one scaffold**~~ done 2026-11-17 (`allhaps/selfdup.py`, curator's
   rule: split only at a gap, and only with evidence that the copies are haplotypic, not
   paralogous; read-depth classes after Flagger): two stretches of one scaffold that align to
   the same stretch of the other haplotype are a `self_duplication` review row, classed by
   the coverage track as `repeat` (about 2x or more), `duplication` (all copies normal) or
   `haplotypic` (that copy and the other haplotype's at about half depth, the other copy
   normal). Only `haplotypic` with a gap between the copies acts: the half-depth copy is cut
   out at the gaps around it (then tagged FalseDuplicate by the small-scaffold step) and the
   scaffold's two sides joined again as assembled; a gap the hifiasm graph backs is not cut.
   Troglobisium: 12 rows, none haplotypic (fold-backs at gaps are repeats, H2_scaffold_17's
   end piece is a duplication); owl: none (`tools/self_dups.py` lists them from a PAF).
   Still open: the H2-against-H1 PAF (primary only) finds a double copy in H2 only; H1 needs
   the reverse alignment. The hifiasm unitig graph (two copies on the two sides of one
   bubble) would be a second, independent test of "haplotypic".
6. **`--phasing none` as an optimisation** (from hapscotch): assign pieces to haplotypes
   to minimise inter-haplotype contact within a chromosome; a principled way to decide
   cases like owl H1_scaffold_31 / H1_scaffold_44, which now rest on fixed rules.
7. **Small pieces** (under 4 px in allhaps, 15 px in single mode): read full-resolution
   strips for them instead of the working mip level.
8. **Contaminant tag**: on the high-MQ map with the repeat track, a scaffold with no
   contact to anything, no homolog and no repeats is a contaminant candidate; tag it
   rather than leave it unplaced and silent.

Evidence
9. **hifiasm graphs + YaHS AGPs** (curator): explored 2026-11-16 (`tools/graph_junctions.py`,
   see PROGRESS). No graph alignment is needed: the AGPs place each contig in its
   scaffold, and the noseq GFAs' read lists place contigs on the unitig graph. An oriented
   path of at most 2 Mb between the two contig ends of a YaHS gap: 107 of 197 gaps the
   curator kept, 1 of 24 they broke (a haplotype switch at a bubble); it would have saved
   the tool's one wrong cut with a path (H1_scaffold_16 at 0.45 Mb) and touched none of its
   right ones. As a proposer of joins it is weaker: for contigs >= 1 Mb the nearest graph
   neighbour is the curator's in 36 of 46 ends (7 of them joins YaHS missed), but paths
   also exist for 8 of the tool's 29 wrong joins, and for small contigs most proposals are
   pieces the curator left unplaced.
   Done in 1.2.0 (`graph.py`, `--gfa`, `--yahs-agp`, `--utg-gfa`, `--fill-gaps`; PROGRESS
   2026-11-16): contig links and bridges keep and close gaps; unitig paths stop
   rearrangement at a gap; `asm` closes overlaps, places bridges and own-haplotype contigs
   that lie wholly on a path, and fills only collapsed homozygous sequence (curator: never
   duplicate the other haplotype's allele; collapsed sequence must be duplicated).
   Still open:
   (a) joins only the graph makes: MADE since 2026-11-17 (curator's review of the 15
       Troglobisium rows: all should be placed). `allhaps/graphjoin.py`: the smaller scaffold
       goes at the larger one's end, before whatever the tool had put there, oriented as the
       graph says, taken from where the tool had it; a small homolog pair goes when both its
       scaffolds are graph-joined into the two copies of one chromosome; graph links between
       copy ends (direct, or through one small scaffold) also let `merge_pairs` merge two
       pairs. Rows `graph_join` (made) or `graph_join_candidate` with why not: under one
       pixel, the larger side is an Unloc or unplaced (PretextView cannot join Unlocs; the
       curator would like the tool to), or a cut piece. Troglobisium level 0.5: 10 of the 15
       made. Earlier note: done 2026-11-17 as review rows (`graph.end_joins`,
       `join_rows`; row `graph_join_candidate`, in the review list): a contig link that is
       the only one at both ends, or mutually nearest scaffold ends along the unitig graph
       with no close alternative, one side >= 1 Mb, not made by the layout; each row gives
       the Hi-C end contact against the best other large scaffold end. Not acted on: on
       Troglobisium 18 such joins, nearly all a small piece at a chromosome end; Hi-C prefers
       the graph's end for only 3 (small pieces touch many ends alike), 6 are under 2 px, and
       the curator made 1. Next: the curator's verdict on some rows decides whether the
       tool should make them (e.g. graph unique + Hi-C best, or graph alone under 2 px);
   (b) single mode does not read the graph;
   (c) own-haplotype contigs that extend beyond the path (14 of 18 on Troglobisium) are
       only named: placing them would mean cutting the contig where it leaves the path;
   (d) H2_scaffold_17 0.45 Mb: with H1's join backed by the graph, the tool now inverts H2's
       end piece on telomere + alignment evidence, which the curator kept; one case;
   (e) only tested on Troglobisium (the only genome with GFAs); `run --utg-gfa` reads the
       noseq file in about 10 s, `asm` the sequence file in 2-4 min.
10. **FastGA on a cluster**: keep minimap2 as the default (curator), but accept a FastGA
    PAF via `--paf` once the tool has its own uniqueness filter (best chain per query
    stretch), since FastGA gives no mapping quality.
11. **telo_scan threshold**: optionally drop windows below telo_scan's own score (100)
    before PretextGraph; removes weak interstitial runs (measured on the owl, optional).

Validation, so the rules keep improving
12. ~~**One scoring script over all fixtures**~~ done 2026-11-16: `tools/score_all.py`
    with `tools/score_cases.json` and `pretext_autocurate/score.py` (junctions, breaks,
    chromosome grouping, tags by bp, edits, review items); baseline table in PROGRESS. To add
    a genome (next: icLarUrsu, diploid with GFAs): put the map(s), inputs and the curated
    savestate in `seed/fixtures/`, add an entry to the cases file, run
    `python3 tools/score_all.py --only NAME`. What the baseline points at:
    (a) Unloc tags of the tool and the curator almost never coincide on the two-haplotype
        genomes (0% on Troglobisium, 0-11% on the owl, with similar totals). Looked at
        2026-11-16 (PROGRESS): the tool's few Unlocs on Troglobisium were large pieces the
        curator joins (now placed by end contact when clear), and the curator's are small
        pieces the tool never weighed (under 20 px without sequence evidence). Curator's rule
        now implemented: a small piece is an Unloc when its high-MQ signal is unique to one
        homolog pair (MQ0 uniqueness stronger; alone without a high-MQ map). Troglobisium
        has no high-MQ map in the fixtures; the curator's small Unlocs there are not unique
        at MQ0, so the scores will only move with that map. The curations predate the rule
        (curator: they are not ground truth);
    (b) HAP1 agrees with the owl curator's own curation (savestate_6) on 53% of the bp: not
        a disagreement. That curation kept each haplotype as assembled; the per-chromosome
        choice of the better copy as HAP1 is the curator's rule implemented in the tool and
        not yet applied in their own curations. Against the curator's corrections of the
        tool (which follow the rule): 91%. HAP1 agreement is only meaningful against
        references made with the rule;
    (c) Troglobisium at level 0: sex tags 41% (H2_scaffold_46 not chained at that level).
13. **More genomes from the team, especially hard taxa** (not vertebrates, not
    butterflies): every corrected savestate (tool output + curator's edits) is a new
    test. Ask curators to send `autocurate.run.txt` with it (README already says so).
14. **AGP**: since 2026-11-16 `layout.write_agp` mirrors PretextView's `GenerateAGP` line
   for line (resolution header, pixel-derived coordinates, painted scaffolds first, "-" only
   inside multi-piece painted scaffolds, paint-display tags dropped; before, `vertPaint`
   reached the AGP, which pretext-to-asm would read as a haplotype). One deliberate
   departure: scaffolds too short for a pixel are appended, unpainted, tagged HAP1/HAP2
   by their H1_/H2_ prefix; otherwise pretext-to-asm puts them in extra "h1"/"h2"
   assemblies (curator: a known nuisance with PretextView exports too). pretext-to-asm
   1.4.1 runs end to end on the tool's Troglobisium and owl results (PROGRESS). Still to
   do: compare with a real PretextView export of one savestate (base-pair rounding).
   On the curator's Troglobisium savestate_final pretext-to-asm stops at a 1-pixel
   unpainted piece inside a contig (H2_scaffold_6 25.83-25.90 Mb): "Expecting 2 abutting
   sub fragments but got 1"; a PretextView export would contain the same piece. Either
   remove the piece in PretextView or report it to sanger-tol/agp-tpf-utils.

Engineering
15. ~~**Split `allhaps.py`**~~ done 2026-11-17 (`docs/ALLHAPS_SPLIT.md`): the package
    `allhaps/` (module table above), named thresholds (`Thresholds`), `curate_allhaps` as
    nine steps over a shared `Run`; outputs byte-identical on all five golden runs
    (`tools/golden.py check`). `tests/test_steps.py` tests the steps on toy genomes with
    synthetic alignments (`synth.toy_blocks`). Still to do when touching it: `_follow` (about
    200 lines) and `_layout` (about 220) are the longest steps; the steps still unpack many
    names from `Run`, which could shrink as the steps are given clearer inputs.
16. ~~**Tests that run without the fixtures**~~ done 2026-11-16: `synth.py` writes real
    `.pretext` files (BC4 textures, every mip, tracks) from a toy genome with misjoins and
    the MQ0 homology line; `tests/test_synthetic.py` runs the format, both passes, the
    command line (with a sweep) and scoring on them; fixture-dependent tests skip
    themselves. GitHub Actions (`.github/workflows/tests.yml`, Python 3.10 and 3.12) runs
    25 tests in about a second on every push; the fixture suite stays for releases. Next:
    toy genomes for more rules (sex chromosomes from the coverage track, small pieces,
    the hifiasm graph end to end with `asm` minus pretext-to-asm).
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
