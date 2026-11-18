# Progress log

## 2026-11-12: Phase 1 (I/O and state), Python

Decisions
- Python first (curator's choice); Rust port deferred until the methods settle.
- Priorities reordered: haplotigs, placement/unlocs, obvious cuts and joins, contaminants,
  p-to-q orientation. See PLAN.md section 0.
- New fixtures: bAegFun (bird, allhaps, MQ0 + MQ10) stands in for the butterfly;
  xgTylRafi (gastropod, hap2 only, MQ0) is the hard single-haplotype case.

Done
- `src/pretext_autocurate/`: `pretextio` (reader, mip-level matrix, header hash, layout
  check), `savestate` (read/write), `mapstate` (simulator, `plan_edits`), `layout`
  (AGP, TSV), `cli` (`info`, `dump`, `convert`).
- `tests/test_phase1.py`, 20 tests, all passing: `python3 -m unittest discover -s tests`.
  (pytest is not installed in this environment; the tests are plain unittest.)

Findings
- Both new maps are 65536 px (HR). A full matrix is 4.3 GB; mip level 3 (8192 px) decodes
  in about 3 s and is 67 MB.
- The header hash uses the first 64 bytes of the file name; verified on the 72/73-character
  names of the new fixtures.
- The seed savestate reader failed when nEdits is a multiple of 8 (bird: 264). Fixed in the
  package; the writer avoids the ambiguous case by adding one result-neutral edit.
- Curators never use explicit breaks: moving part of a scaffold splits it implicitly.
- `plan_edits` reproduces the curated layouts exactly with fewer edits than the curator
  used (bird 88 vs 264, gastropod 281 vs 357).
- Curator cuts: bird 9/9 at gaps; gastropod 55/72 within 3 px of a gap.
- In the bird, HAP1/HAP2 tags are assigned per pair by the curator and do not follow the
  `H1_`/`H2_` name prefix.

Not done / not verified
- AGP output follows PretextView's style from memory (100 bp `U` gaps, `Painted` and tag
  columns). It has not been compared with a real PretextView export or run through
  pretext-to-asm. A sample AGP exported from PretextView for one fixture would settle it.
- No AGP input (`convert --agp`), and no writing of a savestate from scratch (a template
  savestate is needed for colours and display settings; today an existing one is reused).
- Phase 0 is still open: nobody has confirmed that the CAE71 `claude_scaffold1` savestates
  load in PretextView, and the CAE71 `.pretext` maps are not in `seed/fixtures/`.
- Sub-pixel scaffolds (20 in the bird, 32 in the gastropod) are absent from every output.

Next
- Phase 2/3 start with priority 1: signal model (distance decay, background) at mip level,
  then haplotig detection scored against the 132 curator-tagged haplotigs in xgTylRafi.

## 2026-11-12 (later): haplotig detection v0, single mode

Guidance recorded in PLAN.md section 0: p-to-q orientation, biology as a prior,
sex-chromosome and PAR rules, and that the curators' savestates are a guide, not exact truth.

Done
- `signal.py`: distance-decay curve, full-resolution row strips read on demand (no 4.3 GB
  matrix), coverage normalised to the median of large scaffolds.
- `haplotig.py`: for every scaffold of 3 to 400 px, finds the best off-diagonal line (both
  strands, whole length and best half-length run) and combines three pieces of evidence:
  line brighter than adjacency would give, low coverage of the piece, low coverage of the
  partner locus. Tiers: `haplotig_high` (all three), `haplotig_medium` (line + low
  coverage), `review` (one of them), `none`, `too_small`.
- `pretext-autocurate haplotigs MAP --out-dir DIR`: TSV of every call plus a savestate
  with `Haplotig` tags (no edits). About 15 s on xgTylRafi.
- `tests/test_haplotig.py` (4 tests; 24 in total, all passing).

Result on xgTylRafi (117 of the curator's 132 haplotigs are whole scaffolds)
- high + medium: 86 calls, 78 tagged Haplotig by the curator, 5 tagged Unloc, 3 untagged.
- review: 110 scaffolds (29 curator haplotigs among them). Too many to be useful yet.
- none: 37, of which 10 are curator haplotigs.
- The contaminant (scaffold_18) is not called.

Findings
- Near-diagonal values are saturated (about 250 of 255 out to 5 px, 24 kb/px), so line
  brightness cannot separate "copy of locus R" from "belongs next to locus R" for pieces
  under about 8 px. Coverage carries those calls. A high-MQ map would separate them
  directly (a duplication line vanishes at high MQ, adjacency does not); xgTylRafi has
  only MQ0.
- Curator haplotigs have median coverage 0.68 of the norm, not 0.5; their partner loci
  0.70. Unlocs tend to have high coverage (median 1.2, up to 6), i.e. repeats.

Not done
- The 15 curated haplotigs that are parts of larger scaffolds (need cuts at gaps).
- Placing called haplotigs next to their partner (the partner locus is in the TSV).
- allhaps mode: homologous haplotypes give the same line, so this detector must not be
  run on allhaps maps as is. Sex-chromosome guard is only "never call on coverage alone".
- The tagged savestate has not been opened in PretextView.

## 2026-11-12 (later still): placement and unlocs v0, single mode; `run` command

Curator's notes on haplotigs and allhaps (cut haplotigs off contig ends; >2 haplotypes;
HAP1 = better assembled) are recorded in PLAN.md section 0 for when haplotigs are revisited.

Done
- `place.py`: for each scaffold shorter than 0.2 x the longest, scans the mean contact of
  its rows against all "anchor" scaffolds (at least 5 x longer, at least 100 px). The
  strongest 20 px window gives the anchor; it must reach half of what an adjacent piece
  would show (from the map's decay curve) and beat every other anchor by 1.3 x. A finer
  5 px window, the nearest gap or anchor end, and an end-versus-flank contact ratio give a
  proposed insertion point and strand.
- `curate.py`: combines haplotig calls and placement into one layout. Haplotigs and unlocs
  are moved to the end of their anchor, painted with it, and tagged. Edits come from
  `plan_edits` and are verified by replay.
- `pretext-autocurate run --map MAP --out-dir DIR`: savestate (one per map given),
  `autocurate.agp`, `autocurate.decisions.tsv`. About 25 s on xgTylRafi.
- `tests/test_run.py` (6 tests; 30 in total, all passing).

Result on xgTylRafi (whole scaffolds only; curator's label in columns)

    tool \ curator    placed  Unloc  Haplotig  unpainted  Contaminant  cut
    attach_unloc          27     44        37         14            0    5
    attach_haplotig        1      3        75          2            0    0
    tag_haplotig           0      2         3          0            0    0
    none                   5      8         2         82            1    0

- Chromosome: the anchor chosen lies on the curator's chromosome for 190 of 192 moved
  scaffolds that the curator painted. The largest scaffold moved is 69 px (1.6 Mb).
- The contaminant (scaffold_18) is left alone with "no cis-level contact".
- 196 edits, 225 painted contigs, 213 tags, 11 waypoints.

Findings
- Insertion (localising a piece at a gap) is not reliable on this map: of 11 pieces that
  passed every localisation test, 3 were placed there by the curator, 5 are curator
  haplotigs, 2 unlocs, 1 a piece the curator cut. The haplotig cases have their duplication
  line right next to the proposed gap. With saturated MQ0-only signal, "belongs at this
  gap" and "copy of the sequence beside this gap" look the same for small pieces. So
  insertions are only proposed by default (attached as unloc, waypoint at the locus);
  `--apply-inserts` applies them.
- Peak sharpness, gap proximity and coverage do not separate the curator's "placed" from
  "Unloc"; the curator's choice seems to rest on judgement the map statistics do not carry.

Not done / known limits
- Unlocs and haplotigs are attached at the end of the anchor *scaffold*. Where a scaffold
  holds two chromosomes (scaffold_2) or a chromosome is split over several scaffolds, that
  is the wrong end until cuts and joins exist (priority 3).
- 37 curator haplotigs are attached as Unloc (right chromosome, wrong tag); 27 scaffolds
  the curator localised are attached as Unloc.
- Pieces whose best contact is another small scaffold (chains of debris) are assigned via
  larger anchors only.
- Only scaffolds of 3 to 150 px are candidates. (A first version had no upper limit and
  attached the 32 Mb and 27 Mb chromosomes scaffold_16 and scaffold_17 to scaffold_2 as
  unlocs, because the adjacency expectation is meaningless beyond the 200 px decay range.
  Fixed, with a test.) Larger pieces need the cut/join logic.
- Scaffolds under 3 px are not moved. Parts of scaffolds are never cut out.
- Bird (allhaps) not run: placement there must sum contact over both haplotypes.
- AGP format and the written savestates are still unverified in PretextView.

Next
- Cuts and joins (priority 3), so that chromosome units are right before attaching.
- Then back to haplotigs: contig-end haplotigs (cut at gaps), allhaps HAP1/HAP2 logic.

## 2026-11-12: verified in PretextView

The curator loaded `out/xgTylRafi...savestate_autocurate` (written from the built-in blank
template, 196 edits, paint, tags, 11 waypoints) in PretextView: it loads and looks right.
This confirms the savestate writer, the header hash and the blank template on a real map.
Still unverified: the CAE71 `claude_scaffold1` fixtures and the AGP output.

## 2026-11-12: misjoin cuts v0, single mode

Done
- `signal.py`: `block()` reader and `insulation(w, skip)`: contact across every boundary
  relative to the decay expectation (integral images over full-resolution tiles, ~14 s).
- `misjoin.py`: a cut is called where contact across a boundary is below 0.2 of expectation
  on both a 10 px window and a 30 px window set back 10 px, while each side is internally
  connected (so compartments and no-signal regions do not trigger). Cuts snap to a gap
  within 10 px. Without a gap the cut is applied only below 0.1; otherwise flagged. Small
  foreign pieces (two sharp dips close together, far window normal) are cut only when both
  ends are at gaps, otherwise flagged.
- `SegmentedMap`: view of the map split at the cuts, so haplotig calls and placement run on
  the corrected pieces. `curate_single` writes cuts as break edits, paints the pieces of a
  cut scaffold separately, and puts waypoints at applied and flagged cuts.
- 32 tests, all passing.

Result on xgTylRafi
- The curator's 72 cuts split into: 25 between different chromosomes, 10 rearrangements
  within a chromosome, 25 made only to insert another piece, 12 to remove a haplotig/unloc.
- Tool: 25 calls, 24 within 12 px of a curator cut; the other is inside the contaminant
  (flagged only). 14 applied, all matching curator cuts (13 "different chromosome",
  1 rearrangement). 11 flagged.
- Cut-out pieces then find the home the curator gave them, e.g. scaffold_3 74.9-76.9 Mb
  goes to scaffold_1 at 166.66 Mb and scaffold_6 20.3-23.6 Mb to scaffold_2 at 23.1 Mb.
  They are attached as unlocs with a waypoint at the locus (insertions still not applied).
- Chromosome agreement for moved scaffolds: 196 of 198.

Not found (12 curator cuts between chromosomes or rearrangements): pieces under ~10 px,
cuts where contact is reduced but not at background (scaffold_4 16.6 Mb: 0.59), and
rearrangements within a chromosome, which need an orientation/order test, not insulation.

Unverified: break edits (delta INT32_MIN) have not been loaded in PretextView; the curators'
own savestates contain none.

Verified by the curator (2026-11-12): the savestate with 14 break edits (delta INT32_MIN,
x = y = pixel before the cut) loads in PretextView.

## 2026-11-12: clear insertions applied

A localised piece is now inserted (not just proposed) when it is at least 15 px and its
orientation ratio is at least 1.4, and it is not a haplotig call or low coverage. On
xgTylRafi that gives 3 insertions:
- scaffold_3 74.9-76.9 Mb -> scaffold_1 at 166.66 Mb, reversed: same place and strand as
  the curator.
- scaffold_21 -> scaffold_10 at 31.6 Mb, reversed: same place and strand as the curator.
- scaffold_6 20.3-23.6 Mb -> scaffold_2 at 23.1 Mb: same chromosome as the curator, and
  the curator cut scaffold_2 at 23.2 Mb for an insertion; exact neighbours not confirmed.
Smaller or less clearly oriented proposals stay as unlocs with a waypoint. 32 tests pass.
Totals: 14 cuts, 11 flagged cuts, 3 insertions, 130 unlocs, 83 haplotigs attached,
3 haplotigs tagged in place; 236 edits, 39 waypoints.

## 2026-11-12: curator's corrections on top of the tool's output (xgTylRafi)

Fixture: `seed/fixtures/xgTylRafi*.savestate_3_curator_on_autocurate` (the tool's 236 edits
followed by 30 curator edits, repainted into 18 chromosomes; 9.8 Mb left unpainted).

Kept unchanged by the curator: all 14 cuts, all 3 insertions, and every tag on whole
scaffolds except two Unloc -> placed.

What the curator had to do, i.e. what the tool should learn:
1. Join and orient the large pieces left by cuts. scaffold_2 became two chromosomes:
   [0-53.9] + [82.0-88.7] + [53.9-82.0 reversed], and [88.7-109.7 reversed] +
   scaffold_15[19.8-24.4 reversed] + [111.1-177.1]. Needs order/orientation of large
   segments, and candidates above the 150 px limit (scaffold_15 piece: 195 px).
2. Execute the flagged "small foreign piece" cuts. All of them were real, also where one
   end is not at a gap: scaffold_2 40.5-41.4 -> into scaffold_5 at 49.3; scaffold_9
   14.8-15.3 -> into scaffold_10 at 56.7; scaffold_3 0-0.43 -> into scaffold_1 at 3.5;
   scaffold_3 46.7-47.2 -> reversed to scaffold_3 at 23.0; scaffold_10 56.45-56.69 -> Unloc.
3. Insert large unlocs at a gap: scaffold_5 101.1-103.5 (104 px, attached as Unloc to
   scaffold_16) was inserted into scaffold_16 at 22.4 Mb.
4. Attach unlocs/haplotigs at the end of the *chromosome*. Groups the tool attached to the
   end of a cut piece ended up mid-chromosome and were moved (5 block moves).
5. Haplotigs inside scaffolds, at contig ends: scaffold_2 41.46-41.94, scaffold_3
   118.38-119.80, scaffold_9 15.30-15.61, scaffold_10 41.54-42.03, scaffold_13 23.63-23.80.
6. A local rearrangement in scaffold_13 (0-0.99 Mb moved, reversed, to 23.6-23.8 Mb).
7. Tag the contaminant (scaffold_18) `Contaminant`; it stays unpainted.
8. Paint each chromosome as one scaffold id including its unlocs and haplotigs; leave
   unplaced debris unpainted.

## 2026-11-12: haplotigs inside scaffolds (contig level), v0

- `detect_contig_haplotigs`: for each gap-to-gap contig of 8-150 px inside a
  chromosome-scale scaffold with coverage < 0.8, look for a line parallel to the diagonal
  within the same scaffold. Called only when the line is at least 0.75 of the diagonal
  value, at least 2x what distance explains, and the partner locus has coverage < 0.85.
  The contig is cut out at its two gaps and tagged Haplotig.
- xgTylRafi: 2 calls, both curated haplotigs (scaffold_3 119.82-121.3 Mb, scaffold_13
  23.63-23.85 Mb), 0 false calls. 12 other curated contig-level haplotigs are missed:
  looser thresholds give 6 of 14 at the cost of ~68 false calls, because the best line over
  many offsets is noisy (compartments) and small contigs are saturated.
- The scaffold_3 case is a pair of mutual copies either side of a gap. The curator removed
  the right copy in the first curation and the end of the left contig in the second; the
  tool removes the gap-delimited right copy (cuts only at gaps).
- The scaffold_3 haplotig is tagged in place rather than moved: after the cuts its two
  flanks are separate anchors with equal contact. Fix belongs with chromosome-level joins.
- Whole-scaffold tags: in the second curation the curator left all the tool's Unloc tags
  on whole scaffolds as they were, including 37 pieces tagged Haplotig in the first
  curation. Whether that is agreement or not-yet-reviewed is not known; ask.
- 33 tests pass. Totals now: 18 cuts, 11 flagged, 3 insertions, 130 unlocs, 84 haplotigs
  attached, 4 tagged in place.

## 2026-11-13: curator's review of the tags; borderline haplotig tier

Fixture: `seed/fixtures/xgTylRafi*.savestate_4_curator_reviewed_tags` (savestate_3 plus 60
edits and reviewed tags; annotated list in `out/review_unloc_or_haplotig_reviewed.xlsx`).

Review outcome
- Of 37 pieces the tool tagged Unloc that the first curation called Haplotig: 31 haplotig,
  4 unloc, 2 inserted into the chromosome.
- Of 8 tool-only haplotig calls: 6 confirmed (three are haplotigs of unloc sequence),
  2 rejected (scaffold_196, scaffold_237: repetitive unplaced sequence).
- So the first curation was right about most of the 37; the tool's coverage gate was the
  miss. Many true haplotigs have coverage 0.85-1.1 here.

Changes
- `haplotig_low` tier: not confident, but either coverage < 0.8, or the line is at least
  0.75 of the diagonal with coverage < 1.2. These keep their Unloc tag and get a waypoint
  ("POSSIBLE HAPLOTIG" in the decisions table). Curator's request: flag, do not tag.
- Coverage below 0.3 on the piece or its partner is treated as unmappable repeat and is
  not a haplotig signal (fixes scaffold_196 and scaffold_237).

Scores against the review (whole scaffolds): confident calls 84, of which 83 are reviewed
haplotigs. Of the 26 reviewed haplotigs not called confidently, 22 are in the borderline
tier (74 pieces in total). 4 are missed outright (coverage 0.92-1.3, weak line).
34 tests pass. Waypoints: 117 of 255.

Open: the curator asked for vertical-only waypoints. The savestate has no per-waypoint
line type (x, y, z, id); it is probably a display setting bit. Bytes 0-1 of the payload
differ between savestates (bits 0 and 5 of byte 1); which bit is which is not known.

## 2026-11-13: highlight convention (curator)

- Waypoints show both lines in PretextView; they are now used only for positions inside
  contigs: applied cuts, flagged cuts, proposed insertion loci.
- Standalone pieces that need a look (borderline haplotigs, attached at chromosome ends)
  get the `vertPaint` tag in addition to their Unloc tag. To the curator it means only
  "pay attention to this"; PretextView draws it as a vertical highlight band.

## 2026-11-13: joins, orientation, chromosome-end attachment

Done
- `join.py`: for pieces over 150 px, contact between terminal windows (100 px, outermost
  10 px skipped) relative to the adjacent-sequence expectation. Join when the two ends are
  mutual best, score >= 0.6, and the next best partner of either end is below 0.8 x the
  score. Ends with a telomere signal are never joined. Chains are oriented with their
  longest piece forward.
- Flagged cuts confirmed by the curator are now applied: small foreign pieces even when one
  end is not at a gap, and sharp breaks without a gap (both scales < 0.2, sides normal).
  Still flag-only: low contact over a wide region.
- Layout is built per chain: pieces in order and orientation, insertions inside them, then
  the chain's unlocs and haplotigs at the chromosome end. One paint id per chain.
- Decisions table has `join` and `flag_join` rows; flagged joins get waypoints.
- 36 tests pass.

Result on xgTylRafi (against the reviewed curation)
- 26 cuts, 3 flagged; 12 joins, 3 flagged; 4 insertions.
- Joins: scaffold_2 88.7-109.7 (reversed) + scaffold_15 19.5-24.4 (reversed) + scaffold_2
  111.1-177.1 reproduces the curator's chromosome 3 with the same orientations. The other
  joins re-connect flanks after a foreign piece or haplotig contig was cut out.
- Flagged, not joined: scaffold_11 end with scaffold_16 end (score 0.91, but scaffold_16's
  end is telomeric; the curator did not join them), scaffold_3 end with scaffold_8 start
  (both telomeric), scaffold_13 (two ends fit equally: the local rearrangement).
- Tool junctions between untagged pieces >= 15 px: 20, of which 15 are in the curation.
  The other 5 all differ by a small haplotig or unloc piece at the junction (a different
  copy removed, or a piece the curator moved later); none joins different chromosomes or
  has the wrong orientation.
- Chromosome groups: each curator chromosome maps to one tool group, except chromosome 2
  (scaffold_2 53.9-88.7 stays a separate group: its 53.9-82.0 part is inverted inside the
  piece) and the scaffold_13 chromosome (unresolved, see above).
- Curator junctions reproduced: 24 of 45. Missing ones are small-piece insertions (the tool
  attaches those as unlocs), rearrangements inside a chromosome, and the 2.4 Mb scaffold_5
  piece into scaffold_16.
- Tagged pieces: 171 of 184 curator unloc/haplotig pieces are painted into the tool group
  that holds most of their chromosome.

Not done
- Rearrangements and inversions inside a piece (scaffold_2 53.9-82.0, scaffold_6
  73.9-97.7, scaffold_13 0-1.0, chromosome 1 in the reviewed savestate).
- Grouping pieces that share a chromosome but whose ends do not touch.
- p-to-q orientation of whole chromosomes; contaminant tag; allhaps mode.

## 2026-11-13: rearrangements inside a piece (tentative cuts at weak junctions)

Curator's notes on MQ0, repeats and telomeric contacts are recorded in PLAN.md section 0.

Done
- `weak_junctions`: gaps in chromosome-scale scaffolds where contact across is reduced
  (near window < 0.5, far window < 0.6 of expectation) but not at background.
- `curate_single` cuts these tentatively, re-runs joins (and placement for pieces of
  15-150 px next to them), and keeps a cut only if a large piece joins a different end, or
  a small piece is clearly localised elsewhere or reversed. Otherwise the cut is dropped
  and the junction is reported as `flag_cut` ("possible insertion site") with a waypoint.
- Break edits are written only where both sides of a cut stay adjacent in the final layout.
- 37 tests pass (suite now takes about 5 minutes).

Result on xgTylRafi
- 13 weak junctions. 1 kept: scaffold_2 at 82.0 Mb. The piece 82.0-88.7 then joins by its
  88.7 end to the 82.0 end of 53.9-82.0, which is the curator's arrangement
  (82.0->88.7 followed by 82.0->53.9).
- 12 dropped and flagged. Five of them are sites where the curator inserted a piece
  (scaffold_1 166.7, scaffold_2 23.2, scaffold_5 49.3, scaffold_10 31.7, scaffold_16 22.4).
- Totals: 27 cuts, 14 flagged; 13 joins, 4 flagged; 4 insertions; 66 waypoints.
- Tool junctions in the curation: 16 of 21. Curator junctions reproduced: 25 of 45.

Not resolved
- The join of scaffold_2 41.5-53.9 to the 82.0 Mb end is flagged, not made: score 0.62 with
  another end at 0.61, and a telomere signal at one end (interstitial, since the curator
  joined there). So the curator's chromosome 2 is still two tool groups.
- Small local rearrangements: scaffold_6 96.2-97.7 (orientation ratio 1.11), scaffold_13
  0-1.0 (1.04), scaffold_11 61.8-62.5 (localised at 58.0 as in the first curation, ratio
  1.23). All below the 1.4 needed to act; they stay flagged.
- A weak junction is a useful hint that something belongs there; it is not yet used to
  support insertions.

## 2026-11-13: insertions at reduced-contact junctions

Curator's notes on inversions (equal score either way = probably heterozygous, do not
edit; nested inversions may need iteration) are recorded in PLAN.md section 0.

Finding: every insertion the curator made sits at a junction where contact across the gap
is reduced (near window <= 0.55 of expectation). The wrong proposals, mostly haplotigs,
sit at junctions with normal contact (0.7-1.3): there is no room there, so a piece that
seems to belong is a copy of the sequence beside it. Also, the tool had been testing the
gap nearest the contact peak, which was often a few pixels from the right one; at the
right gap the orientation ratio is clear (1.85 for the scaffold_5 piece, 2.0 for the
scaffold_2 piece).

Changes
- `place_scaffolds(junction_contact=...)`: a piece is localised only at a junction with
  contact below 0.6 (or a piece end), chosen among such junctions near the contact peak by
  the clearest orientation. Otherwise: unloc, "no junction with reduced contact".
- Insertion is applied for pieces >= 15 px with orientation ratio >= 1.3 (was 1.4), one
  piece per junction (the clearest), not at a chromosome end (end extensions stay unlocs).
- `Signal.insulation` is cached.

Result on xgTylRafi: 8 insertions, all consistent with the reviewed curation:
scaffold_5 101.1-103.5 -> scaffold_16 at 22.4 (+); scaffold_2 40.5-41.4 -> scaffold_5 at
49.3 (+); scaffold_6 20.3-23.6 -> scaffold_2 at 23.2 (+); scaffold_3 74.9-76.9 ->
scaffold_1 at 166.7 (-); scaffold_21 -> scaffold_10 at 31.7 (-); scaffold_3 46.7-47.2 ->
scaffold_3 at 23.0 (-); scaffold_6 96.3-97.7 -> scaffold_6 at 74.0 (+, a weak-junction cut
now validated); scaffold_20 -> scaffold_7 at 39.7 (-).
Tool junctions in the curation: 23 of 29; the
rest differ by a small haplotig/unloc piece at the junction. Curator junctions reproduced:
33 of 45 (was 25). 37 tests pass.

Still missing: scaffold_3 0-0.43 into scaffold_1 at 3.5; the curator's later
rearrangement of chromosome 1; scaffold_71 at scaffold_2 53.9 (15 px, orientation
unclear); scaffold_3 84.4-85.4 next to scaffold_21 (second piece at the same junction);
scaffold_9 14.8-15.3 into scaffold_10 at 56.7; scaffold_13 0-1.0 to 23.6.

## 2026-11-13: handoff written; first measurements on the bird (bAegFun, allhaps)

`docs/HANDOFF.md` describes the code, conventions, fixtures, numbers and gaps for a new
session. `tools/compare_layouts.py` is the scoring script used throughout.

Bird: the single-haplotype cut detector does not transfer, on either map.
- Curator's 9 cuts: H1_scaffold_1 134.5 Mb, H1_scaffold_5 5.3, H1_scaffold_24 9.0,
  H2_scaffold_1 134.6, H2_scaffold_4 83.2, H2_scaffold_29 4.9, H2_scaffold_30 2.8,
  H2_scaffold_31 2.5, H2_scaffold_36 0.5.
- MQ0 map (decay 254/247/191/113/57 at 0/1/5/20/100 px): 48 calls, 12 applied, none of
  them a curator cut. Contact across the curator's cuts is 0.24-0.67 of expectation, not
  background (0.45 at the H1_scaffold_1 misjoin): at MQ0 a false join between two
  chromosomes still shows contact, so the 0.2 threshold misses all of them, while 10
  other places in large scaffolds are cut wrongly.
- MQ10 map (decay 221/100/19/6/2: the signal is gone beyond a few pixels at 41 kb/px, and
  the far background is relatively high): 1332 calls, 1133 applied. The map is too sparse
  at full resolution for 10-30 px windows. But 6 of the 9 curator cuts have contact 0.0-0.17
  there, and 5 are among the applied calls, so the signal is right and the scale is wrong.
- Homolog signal at mip level 4: H1_scaffold_2 x H2_scaffold_2 mean 19.5 on MQ0 against
  0.01 for an unrelated pair; on MQ10 the same pair is 39 against 32 for an unrelated pair
  (MQ10 values are normalised differently; compare ratios within a map only).

What this means for allhaps mode (next session)
1. Work resolution must adapt to signal density: on a sparse map compute insulation and
   end contact on a mip level where the decay curve still has signal at the window size
   (choose the level from the decay curve, e.g. where g(window) > a few x background).
2. Cuts and joins from the high-MQ map (real contacts only); MQ0 for pairing homologs.
3. Thresholds assume a dense map with saturated diagonal; express them against each
   map's own background, not only against the decay expectation.
4. Nothing from this session has been run to completion on the bird; no bird savestate
   has been written.

## 2026-11-13: allhaps step 1 on the bird (pairing, orientation, tags)

Curator's specs recorded in PLAN.md section 0: allhaps first step (pairs side by side,
parallel not X, HAP tags, sex chromosomes if possible, one pass preferred) and
`--phasing none|hic|trio`.

Done
- `allhaps.py`: `pair_homologs` reads homology from the MQ0 map at mip level 3 (8192 px).
  For every pair of scaffolds from different haplotypes (name prefix, `^(H\d+)_`) it takes
  the mean along the diagonal and the anti-diagonal of their block. Partners must be
  mutual best by total line signal (line x length, so a few noisy bins of a tiny scaffold
  cannot outrank a full-length homolog), at least 5x the scaffold's own median line, and
  2x the next best. Diagonal vs anti-diagonal gives the relative orientation.
- `curate_allhaps`: pairs laid out largest first, the two copies adjacent, the second
  reversed when needed so the homolog line is parallel. Painted 2k-1 / 2k. Tags:
  `trio` = by prefix; `hic` and `none` = better copy (fewer gap pixels, telomeric ends,
  length) is HAP1, whole scaffolds only so prefixes are never mixed within a superscaffold.
  Every other scaffold gets HAP1/HAP2 by prefix.
- Unequal pairs (length ratio < 0.5) are reported as possible sex chromosomes, both tagged
  HAP1, not flipped (a shared region is too little to orient by). Unpaired scaffolds are
  reported with coverage relative to paired copies.
- `run --mode allhaps --phasing ...` writes savestates for every map given, AGP, decisions.
- `tests/test_allhaps.py`: 5 tests, pass in 5 s.

Result on bAegFun (`out_bird/`)
- 32 pairs, all of them pairs in the curator's savestate, none false. Includes
  H1_scaffold_1/H2_scaffold_1 (each still two chromosomes; no cuts in step 1).
- Relative orientation agrees with the curator for 27 of 27 pairs of unsplit scaffolds.
- Z and W found as the one unequal pair (H1_scaffold_3 96 Mb, H2_scaffold_8 34.5 Mb), both
  HAP1 as the curator has them. Their coverage is 0.49 and 0.51 of a paired copy, which is
  a second, independent sex-chromosome signal in this map.
- HAP1 choice: same copy as the curator in 21 of 32 pairs. The tool picked the H2_ copy
  as HAP1 in 14 of 31 pairs. What decided the curator's choice in the other 11 is not known.
- 8 curator pairs not found: all under 5 Mb or involving scaffolds the curator cut
  (H2_scaffold_30, _31, _36 hold two microchromosomes each). One extra piece
  (H2_scaffold_39) is placed next to its homolog, unpainted.
- 39 edits, 64 painted scaffolds, 469 tags. Not yet opened in PretextView.

Not done for allhaps: cuts (incl. the scaffold_1 misjoin in both haplotypes), joins, unlocs,
haplotigs/false duplications, microchromosomes under ~1.3 Mb (4 bins at level 3; needs full
resolution), `none` phasing beyond whole-scaffold swaps, collapsed-chromosome detection
(coverage ~2x is reported for unpaired scaffolds but none occurs here, so it is untested).

Verified by the curator (2026-11-13): the bird allhaps savestates in out_bird load in PretextView and look right.

## 2026-11-13: allhaps cuts on the bird (telomere at a gap, carried to the homolog)

Curator's note recorded in PLAN.md section 0: use telomere signal and the other
haplotype's contiguity (dual scaffolding) for misjoins and orientation.

Findings
- The file's mip levels are re-normalised and saturate (decay nearly flat at level 4), so
  they cannot be used for contact ratios. `Signal.binned(f)` averages full-resolution
  values into bins instead; `Signal.insulation_binned` works on that.
- Even binned, the bird's MQ10 map is too noisy to scan for misjoins by contact alone
  (218 regions below 0.25 at 329 kb bins; the curator's cuts range 0.0-1.8).
- The bird has only 27 gap pixels. 11 have a telomere peak within 3 px. The curator cut at
  8 of those; a 9th curator cut (H2_scaffold_1 134.7 Mb) has no telomere signal but is the
  homologous point of H1_scaffold_1's cut.
- Contact across (full resolution, 10 px window) separates the rest: the Z chromosome's
  telomere-at-gap site (H1_scaffold_3 78.9 Mb) has normal contact on MQ10 (1.82), i.e. an
  interstitial telomere repeat, and the curator did not cut it.

Done
- `telomere_cuts` in `allhaps.py`: cut at a gap inside a scaffold when a telomere peak is
  within 3 px and contact across is below 0.5 on the high-MQ map (MQ0 if there is none).
  Then carry each cut to the homologous point of the paired scaffold if that is a gap with
  reduced contact. `curate_allhaps(hi=...)` cuts first, then pairs the pieces.
- Cuts get waypoints and `cut` rows; break edits only where the sides stay adjacent.
- 6 bird tests pass (25 s). The full suite was not rerun after these changes.

Result on bAegFun
- 10 cuts. 8 are curator cuts (H1_scaffold_1, H1_scaffold_5, H1_scaffold_24, H2_scaffold_1
  via homology, H2_scaffold_4, H2_scaffold_29, H2_scaffold_30, H2_scaffold_31).
- 2 are not in the curation: H1_scaffold_2 at 166.3 Mb and H2_scaffold_2 at 5.1 Mb. They
  are homologous points (the two scaffolds are in opposite orientation), each a gap with a
  telomere peak (262, 130) and contact across of 0.13/0.15 on MQ0 and 0.0 on MQ10, lower
  than at the scaffold_1 misjoin. The 4.2 and 5.1 Mb end pieces pair with each other. To
  ask the curator: a missed misjoin, or a known interstitial telomere?
- Missed: H2_scaffold_36 at 0.5 Mb (telomere 610, but contact across is not reduced; the
  piece is 12 px, at the limit of the window).
- After cuts: 36 pairs + Z/W, all curator pairs, including the pieces
  (H1_scaffold_24:2 with H2_scaffold_30:1, H1_scaffold_35 with H2_scaffold_30:2,
  H1_scaffold_38 with H2_scaffold_31:1).
- Output in `out_bird/`: 67 edits, 74 painted, 10 waypoints. Not yet opened in PretextView.

## 2026-11-13: curator's bird corrections; telomere sites decided three ways; more pairs

Fixture: `seed/fixtures/bAegFun*.savestate_curator_1_on_autocurate` (13 curator edits on an
earlier tool output). The curator did not finish the microchromosomes in it; the original
`.savestate_6` is the fuller reference for those. Notes in PLAN.md section 0: scaffold_2,
dot chromosomes, sequence as input (minimap2 / miniprot).

What the curator said and did
- scaffold_2 is not a misjoin at its telomere-bearing gap: keep joined, invert the H1 end
  piece. Swapping the end pieces between haplotypes looks better but breaks `--phasing hic`.
- Paired small chromosomes by hand: H1_34/H2_35, H1_37/H2_31 (2.5-4.9), H1_36/H2_40,
  H1_42/H2_36, H1_31/H2_37, and H1_5 (0-5.3) with H2_39 + H2_4 (83.3-87.4).
- H2_29: the 2 Mb end piece moved to the other end, i.e. inverted in place.

Changes (`allhaps.py`)
- `telomere_sites` replaces `telomere_cuts`. Each telomere-at-gap site (and its homologous
  gap in the other haplotype) is decided as:
  - `cut` if the two sides have different homologs after a tentative cut, or share little
    chromosome-level contact (MQ0, 329 kb bins, up to 5 Mb each side, background removed,
    below 0.6), or it is the homologous site of a cut (dual scaffolding);
  - `invert_tip` if the sides belong together and the end piece (up to 40% of the scaffold)
    fits at least 1.3x better the other way round on the high-MQ map;
  - `keep` otherwise, with a waypoint.
- Pairing: margin over rivals lowered from 2.0 to 1.3. A scaffold whose own best partner is
  x counts as a piece of x's homolog, not a rival, provided the pieces fit in x's length
  (1.3x). Further pieces need 2x margin and a line comparable to pairs of similar size.

Result on bAegFun (`out_bird/`, 6 bird tests pass; full suite not rerun)
- 5 cuts, all curator cuts: H1_scaffold_1, H2_scaffold_1 (via homology), H1_scaffold_24,
  H2_scaffold_30, H2_scaffold_31.
- 2 tip inversions: H1_scaffold_2 (as the curator asked) and H2_scaffold_29 (as the
  curator arranged it). 3 kept and flagged: H2_scaffold_2 (correct), H1_scaffold_5 and
  H2_scaffold_4 (the curator cut these; in the original curation both pieces of each stay
  in one chromosome).
- 39 pairs + Z/W. Against the original curation: 38 of its 40 main pairs found; not found:
  H1_42/H2_36, and H1_33/H1_40 (both `H1_`: the two copies of one microchromosome ended up
  in the same haplotype file, which the pairing cannot see since it only pairs across
  prefixes). One pair is not in the original curation but is in the curator's new
  savestate (H1_31/H2_37).
- H2_scaffold_39 (homolog of the H1_scaffold_5 end piece together with the H2_scaffold_4
  end piece) is no longer placed: the Hi-C line of a 2.6 Mb scaffold against a 90 Mb one is
  diluted. Sliding the short scaffold along the long one fixed that but broke the pairing
  of cut pieces of scaffold_1, so it was reverted.

Why sequence homology is the right next step
- Margins here are thin: microchromosomes contact each other so much that the best homolog
  line is only 1.3-2x the next. minimap2 on the assembly would give homology directly, with
  coordinates and strand, including same-prefix pairs and small pieces.
- minimap2, miniprot and busco are not installed on this machine (checked with `which`).

## 2026-11-13: sequence homology (minimap2) for allhaps; BUSCO module

Curator's proposals recorded in PLAN.md section 0: assembly FASTA as input (minimap2
between haplotypes, BUSCO genes via miniprot), `telo_scan.py` as a future telomere option.

Environment
- conda env `pretext-autocurate` created with minimap2 2.31 and miniprot 0.18 (curator's
  request; base env untouched). `homology.find_tool` looks on PATH, then in that env.
- `aves_odb12` (308 MB archive, 1.9 GB unpacked) downloaded with the curator's approval to
  `~/.cache/pretext-autocurate/busco/aves_odb12`.
- Disk on this machine is nearly full (about 13 GB free). The split haplotype FASTAs
  (2.7 GB) are deleted after aligning; the PAF is kept.

Done
- `homology.py`: split the FASTA by haplotype prefix, `minimap2 -x asm5 --secondary=no`
  H2 against H1 (about 3 minutes for the bird with 6 threads), cached PAF, `read_paf`
  (blocks >= 10 kb, mapq >= 30: lower-quality repeat alignments, e.g. from W, compete with
  true homologs of small chromosomes).
- `pair_homologs(blocks=...)`: score = percent of the shorter scaffold aligned, ranking by
  aligned bases, strand by majority. Scaffolds >= 500 kb. Hi-C estimate remains the fallback.
- `telomere_sites(blocks=...)`: the matching gap in the partner is found through the
  alignments flanking it: a gap whose two sides are homologous to the two sides of the
  telomere site. (The crude proportional estimate failed for scaffold_1, where the two
  haplotypes joined the same two chromosomes in different arrangements.)
- A lopsided pair is called a sex-chromosome candidate only if both members have about half
  coverage; otherwise "one copy fragmented".
- CLI: `--fasta`, `--paf`, `--threads`, `--work-dir`, `--minimap2`. The FASTA's names are
  checked against the map.
- `busco.py`: lineage as a directory or a name (cache, opt-in download that prints file and
  size first), `run_miniprot`, `read_hits` (complete copies per gene), `summarise` (single /
  duplicated / missing per haplotype, per scaffold). Not yet wired into `run` or scored.
- Fixture `seed/fixtures/bAegFun_asm2_run1.H2_vs_H1.asm5.paf` (3.8 MB) and
  `TestAllhapsSequence`. 9 bird tests pass; the full suite was not rerun.

Result on bAegFun with `--paf` (`out_bird_seq/`)
- 5 cuts, all curator cuts; 2 tip inversions (H1_scaffold_2, H2_scaffold_29); 3 kept.
- 40 pairs + Z/W + 2 further pieces (H2_scaffold_39 with H1_scaffold_5, H2_scaffold_43 with
  H1_scaffold_4). 39 of the original curation's 40 main pairs found; the missing one is
  H1_scaffold_33/H1_scaffold_40, both `H1_` (the alignment is H2 against H1 only).
- Relative orientation agrees with the curator in 28 of 28 comparable pairs.
- One pair differs from the curator: H2_scaffold_37 goes with H1_scaffold_33 (52.7% of it
  aligns) where the curator's new savestate put it with H1_scaffold_31 (which shows little
  alignment to any H2 scaffold: best 11%). To ask.
- 139 scaffolds unpaired (the many small H1 scaffolds); not yet placed.

Next
- Score BUSCO completeness per haplotype copy and use it for the HAP1 choice.
- All-against-all alignment for same-prefix pairs and false duplications.
- Place small scaffolds along their homolog using block coordinates.

## 2026-11-13: BUSCO on the bird

- miniprot with all 50,767 `ancestral_variants` of aves_odb12 was heading for about 3.5 h
  on this machine (6 threads); stopped. With the 6,251 `ancestral` proteins (one per gene)
  it finished in about 25 minutes. `lineage_proteins` now defaults to `ancestral`;
  `--busco-variants` selects the full set.
- Fixture: `seed/fixtures/bAegFun_asm2_run1.aves_odb12.ancestral.miniprot.paf` (7.7 MB).
- Wired into `run --mode allhaps`: `--busco-lineage PATH_OR_NAME`, `--busco-download`,
  `--busco-variants`, `--miniprot-paf`, `--miniprot`. Adds `busco_summary` rows and the
  gene counts of each pair to the decisions table.
- Result: complete genes H1 6167 (99.3% of the 6213 found anywhere), H2 5948 (95.7%);
  duplicated 21 and 11. H2's deficit is mostly the Z chromosome, which is in H1.
- For the HAP1 choice BUSCO hardly helps here: the two copies of a pair differ by 0-3
  genes. Taking "more genes" as better agrees with the curator in 11 pairs, disagrees in 8,
  ties in 13; the gap/telomere/length rule agrees in 23 of 32. So BUSCO decides only when
  the copies differ by more than max(2 genes, 2%), which happens for 3 pairs.
- Microchromosomes carry few BUSCO genes in this set (0-20 per scaffold), so BUSCO adds
  little homology evidence for them; the sequence alignment does that job.
- 10 bird tests pass. The full suite was not rerun today after the allhaps work.

## 2026-11-13: repetitive scaffolds; BUSCO full_table input

- Curator: H2_scaffold_37 is too repetitive to assign uniquely (the original curation
  leaves it out of chromosomes). Alignment spread does not show that (82% of it aligns
  uniquely, 53% to H1_scaffold_33), but long-read coverage does: H2_scaffold_37 1.54,
  H1_scaffold_33 1.49, H1_scaffold_31 1.41 of a normal copy, against 0.5-1.05 for the
  small chromosomes the curator did pair. Rule: a scaffold with median coverage >= 1.3 is
  not paired; reported as "too repetitive to assign". (An alignment-based rule was tried
  first and wrongly rejected H1_scaffold_37 with H2_scaffold_31; dropped.)
- Bird with sequence homology now: 5 cuts, 2 tip inversions, 3 kept, 39 pairs + Z/W,
  2 further pieces, H1_scaffold_33 and H2_scaffold_37 left unplaced as repetitive.
  With Hi-C only: 38 pairs + Z/W.
- `--busco-table full_table.tsv` (`busco.read_full_table`): uses an existing BUSCO run on
  the same assembly; copies are counted per haplotype from the rows, so "Duplicated" in a
  diploid assembly is handled. BUSCO stays optional throughout. Curator's notes in PLAN.md.
- 11 bird tests pass.

## 2026-11-13: per-haplotype BUSCO tables; false-duplication candidates

- `--busco-table` can be repeated and takes `PREFIX=path` for a table made on one haplotype
  alone (`--busco-table H1=hap1.full_table.tsv --busco-table H2=hap2.full_table.tsv`); the
  prefix plus underscore is put in front of the sequence names. The Complete/Duplicated
  label is not used beyond keeping those rows; copies are counted from the rows.
- `busco.false_duplicates`: a scaffold with at least 3 BUSCO genes, at least 80% of which
  have a second copy on one other scaffold of the same haplotype, is reported as
  `false_duplicate_candidate` with the region on the other scaffold, and highlighted with
  `vertPaint` (not tagged: the curator had two of them as Unloc).
- On the bird (from the tool's own miniprot run on the diploid assembly) this gives exactly
  the three cases visible in the curator's per-haplotype BUSCO tables: H1_scaffold_192 on
  H1_scaffold_14 13.72-14.02 Mb (9 genes), H2_scaffold_43 on H2_scaffold_3 8.49-9.20 Mb
  (4 genes), H1_scaffold_209 on H1_scaffold_12 24.43-24.62 Mb (4 genes).
- Not handled: duplicated genes close together on the same scaffold (tandem duplicates or
  local false duplications, e.g. H1_scaffold_24 3.18-3.31 Mb). 12 bird tests pass.

## 2026-11-13: curator's per-haplotype BUSCO tables

Fixtures: `seed/fixtures/bAegFun_asm2_run1.hfsm.asm.hic.hap{1,2}.yhs_mq10.aves_odb12.full_table.tsv`
(BUSCO 6.0.0, aves_odb12, run on each haplotype; sequence names are not prefixed).
- Read with `--busco-table H1=...hap1... --busco-table H2=...hap2...`: all 68 prefixed
  sequence names are in the map and no gene lies beyond its scaffold's end.
- Complete genes: H1 6222 of 6251 (99.5%; 32 duplicated), H2 6036 (96.6%; 25 duplicated).
  The tool's own miniprot run gives 6167 and 5948: miniprot with one ancestral protein per
  gene is a little less sensitive than BUSCO, same ranking.
- False-duplication candidates are the same three from both sources: H1_scaffold_192 on
  H1_scaffold_14 (9 genes), H2_scaffold_43 on H2_scaffold_3 (5), H1_scaffold_209 on
  H1_scaffold_12 (4).
- 13 bird tests pass.

## 2026-11-13: small scaffolds in allhaps (alignment against everything, unlocs, same-prefix pair)

Note from the curator: an allhaps test that needs extensive curation (like the gastropod)
would be more robust than the bird. The other haplotype of xgTylRafi may be recoverable
(S3, reads, remapping: slow); the curator will also look for another recent genome.

Baseline before changes: 13 allhaps tests pass.

Done
- `homology.align_small`: every sequence under 10 Mb against the whole diploid assembly,
  `minimap2 -x asm5 -D -N 5 -p 0.5` (-D drops a sequence's alignment to itself; secondary
  alignments kept, marked in `Block.primary`). About 11 minutes for the bird with 6
  threads, 13 GB index. Fixture `seed/fixtures/bAegFun_asm2_run1.small_vs_all.asm5.paf`.
  CLI: run from `--fasta` unless `--no-small-align`; `--small-paf` takes an existing one.
- `allhaps.same_prefix_pairs`: two unpaired scaffolds of one prefix, each other's best
  match, at least half of the shorter aligned at mapq >= 30. Gives H1_scaffold_33 /
  H1_scaffold_40 (the pair the H2-against-H1 alignment could not see) and nothing else.
- `allhaps.place_small`: chromosome of each remaining small scaffold from two lines of
  evidence. Sequence: share of the scaffold aligned to a chromosome's copies (>= 50%, 2x
  the next chromosome). Hi-C: mean MQ0 contact per chromosome (best copy), best / next.
  - `unloc` (moved to the end of its copy, painted with it, tagged Unloc + HAP): sequence
    and Hi-C (>= 1.5x) agree; or, for pieces >= 20 px, Hi-C alone at >= 2.5x on MQ0 and
    >= 2x on the high-MQ map for the same chromosome.
  - `unloc_candidate` (not moved, `vertPaint`): sequence only, or coverage >= 1.5x a
    normal copy (repeat-like).
  - A piece aligned to both copies of the chromosome is noted as a third copy (haplotig or
    false duplication) and highlighted.
  - The copy is the one with the piece's own prefix (no mixing under any phasing for now);
    for a sex-chromosome pair it is the copy with the evidence.

What the measurements showed (392 unpaired small scaffolds)
- Most are repeats: at MQ0 they contact the dot chromosome H2_scaffold_30:1 /
  H1_scaffold_24:2 at 1.2-1.6x the next chromosome and show nothing at MQ10. A ratio
  under 2.5 at MQ0 without sequence support is therefore not evidence.
- Many align mostly to themselves (internal repeats); those alignments are ignored.

Result on bAegFun (`out_bird_seq/`, with both PAFs and the BUSCO tables)
- 41 pairs + Z/W (new: H1_scaffold_33/H1_scaffold_40), 5 cuts, 2 tip inversions, 3 kept.
- 8 unlocs. 5 are scaffolds the curator placed in the same chromosome: H2_scaffold_34
  (W; Hi-C alone, 3.6x and 6.9x), H2_scaffold_61, H2_scaffold_70, H1_scaffold_44,
  H1_scaffold_209. 3 are not in the curation: H1_scaffold_192 (also a BUSCO false
  duplication candidate), H1_scaffold_310 (1 px), H2_scaffold_86 (third copy).
- Differences in which copy: the curator put H1_scaffold_44 and H1_scaffold_209 with the
  H2_ copy (H1_44 inside the chromosome, not as Unloc); the tool keeps each with its own
  prefix. To ask: is mixing wanted for such pieces under `--phasing hic`?
- 15 candidates flagged, among them H2_scaffold_53 (the curator's Unloc; no Hi-C support).
- Curator placements with no evidence in either line: H1_scaffold_55 (38% self-aligned
  repeats, 3% to its chromosome), H1_scaffold_141, H1_scaffold_155.
- 14 allhaps tests pass (`tests.test_allhaps`, new `test_small_scaffolds`). The single-mode tests were not
  rerun; no single-mode code changed. Savestates in `out_bird_seq/` not yet opened in PretextView.

## 2026-11-13: small pieces may cross prefixes (curator's decision)

Curator: under `--phasing hic` a small piece may go to the other prefix's copy of its
chromosome. Never under `trio`. Large scaffolds still do not mix.
- `place_small` chooses the copy: sequence found in one copy only -> the copy that lacks
  it; found in both (third copy) -> the copy it contacts at least 1.5x more on the high-MQ
  map, else the shorter copy; Hi-C alone or `trio` -> its own prefix. The evidence text
  says how the copy was chosen and when prefixes are mixed.
- Bird: the five curator placements with evidence now sit on the curator's copy too
  (H1_scaffold_44 with H2_scaffold_40, H1_scaffold_209 with H2_scaffold_13, H2_scaffold_70
  with H2_scaffold_30:2 by high-MQ contact, H2_scaffold_61 with H2_scaffold_5,
  H2_scaffold_34 with W). The other three: H1_scaffold_192 with H2_scaffold_14,
  H2_scaffold_86 with H1_scaffold_21, H1_scaffold_310 (1 px) with H2_scaffold_5.
- Still Unloc at the end of the copy: the position inside the chromosome (H1_scaffold_44
  in the curation) is not resolved.

## 2026-11-13: first run on Troglobisium racovitzai (qqTroRaco, allhaps, MQ0 only)

Fixture from the curator: `seed/fixtures/qqTroRaco*.pretext` (4.36 Gb diploid, 66.6 kb/px,
765 scaffolds, textures of 2048 px) and `.savestate_final`: 129 painted units = 60 pairs
+ X, Y and a few single groups, 38 cuts in 22 scaffolds (27 at gaps), many joins, 47
unlocs, 8 superscaffolds with mixed prefixes. H1_scaffold_1 holds X plus both haplotypes
of another chromosome. No sequence, no high-MQ map.

Run as is (`run --mode allhaps --map ...`, about 1 minute, `out_tro/`)
- 61 pairs + 1 unequal. 60 are consistent with the curator's pairs; 56 of the curator's
  60 pairs are found. Wrong: H1_scaffold_1 with H2_scaffold_4 called an unequal pair
  (really X plus an autosome pair, against Y plus an autosome), and H1_scaffold_33 with
  H2_scaffold_6 (cut 7 times by the curator).
- Unlocs by Hi-C alone (pieces >= 20 px, 2.5x): 10, every one in the curator's chromosome
  and copy. The curator placed most of them inside the chromosome (joins), not as Unloc.
- New guard, found here: H1_scaffold_30 (32 Mb) and H1_scaffold_44 (25 Mb) were first
  attached as unlocs to the wrong chromosome. Each is the homolog of one half of a
  misjoined scaffold of the other haplotype (H2_scaffold_1, H2_scaffold_2). Such a piece
  touches only the other haplotype's copy: contact with its own haplotype's copy is
  0.05-0.06 of that, against 0.36-0.75 for the ten real pieces. Rule: below 0.2 the piece
  is not moved and is reported as pointing at a misjoin in the other copy.
- 0 cuts. The telomere-at-gap rule finds nothing: the telomere track is noisy here (99th
  percentile 37, 62880 non-zero pixels, threshold 185; only 2 curator cuts reach it), so
  the motif is probably not this species' telomere repeat. But contact across is clearly
  reduced at most curator cuts: 8 cuts not at gaps have 0.01-0.09 (30 px window), and most
  gap cuts 0.2-0.6. Unlike the bird's, this MQ0 map shows misjoins by contact.

Next for allhaps cuts (the main gap this genome exposes)
1. Cut on contact without telomere support: the single-mode detector on the allhaps MQ0
   map (windows 10 and 30 px), gaps first.
2. Confirm with the homolog: the two sides of a misjoin have different homologs, and an
   unpaired chromosome-sized piece that touches only the other copy marks the place.
3. Then joins of the cut pieces (most of the curator's 216 moves).
- Tests: `TestAllhapsPseudoscorpion` (2 tests, 1 minute).

## 2026-11-13: contact-based cuts for allhaps; one conservative-to-aggressive setting; sweep

Curator's guidance
- Telomere tracks are not always reliable. Vertebrates and butterflies are now easy;
  other arthropods, other invertebrates and various other phyla are hard. More examples
  will follow, to make the tool robust.
- Wanted: one global setting on the conservative-to-aggressive axis, and, if runs are fast,
  a sweep over it so the curator can pick the savestate to start from.
- More Troglobisium inputs arrived: per-haplotype BUSCO tables (arachnida_odb10) and the
  diploid FASTA. Not used yet (next: alignments, as for the bird).

Measurements on Troglobisium (217 gap pixels, 36 distinct curator cut sites, 26 at gaps)
- Contact across the gap (10 px window, minimum within 2 px) separates them: curator cuts
  at gaps are mostly 0.19-0.60; other gaps are 0.64 or more with four exceptions.
- 8 curator cuts away from gaps have contact at background on both scales (0.01-0.09).
- The Hi-C homolog of a 20 px window (best diagonal line elsewhere in the map) is too
  ambiguous at MQ0 to use alone: the best line is often only 1.0-1.3x the next, and it
  "jumps" at many gaps the curator left alone. Not used.
- The single-mode detector as it is finds 4-12 of them (it needs both windows low).

Done
- `allhaps.Level` and `curate_allhaps(level=...)`; CLI `--level` (0 conservative, 0.5
  default, 1 aggressive). It moves together: contact below which a gap is a candidate
  (0.45 / 0.6 / 0.75), sharp breaks away from gaps (0.05 / 0.1 / 0.2), cuts of short end
  pieces (0.35 / 0.5 / 0.65), pairing margin (1.5 / 1.3 / 1.15), unloc thresholds.
- `telomere_sites` now takes every gap with reduced contact as a candidate site; a
  telomere signal is supporting evidence only. Gap pixels within 5 px are one site. Sharp
  breaks away from gaps (both windows at background, normal contact each side, dense map
  only) are cut outright. The decision for gap sites is as before (different homologs or
  no chromosome-level contact: cut; end piece fits better inverted: invert_tip; else keep
  and flag), plus: when one side is too short to measure chromosome-level contact, cut on
  the contact across the gap. (Before, a NaN there fell through to "cut".)
- A pair whose two copies both have about half coverage is tagged as sex chromosomes
  (both HAP1) even when the copies are of similar size (X and Y here).
- `--sweep 0,0.25,0.5,0.75,1`: one directory and savestate per level plus `sweep.tsv`.
  Signal objects and the mip matrix are kept on the map object, so after the first level
  (about 90 s on Troglobisium) each further level takes about 20 s.

Result on Troglobisium, map only (36 curator cut sites)
| level | cuts | false | tip inversions | kept (flagged) |
|---|---|---|---|---|
| 0 | 9 | 0 | 2 | 2 |
| 0.25 | 9 | 0 | 4 | 8 |
| 0.5 | 12 | 0 | 4 | 9 |
| 0.75 | 14 | 1 | 3 | 11 |
| 1 | 15 | 1 | 4 | 14 (4 not at curator cuts) |
- At the default every one of the 25 sites is a place the curator cut. Three tip
  inversions are exactly what the curator did (H1_scaffold_22, H1_scaffold_26,
  H2_scaffold_26). The false cut at 0.75 and above is H1_scaffold_2 at 21.8 Mb.
- With the cuts, H2_scaffold_1 and H2_scaffold_2 are split and H1_scaffold_30 /
  H1_scaffold_44 pair with the right halves; X (H1_scaffold_1 25-58 Mb) pairs with Y
  (H2_scaffold_4 0-33 Mb). 68 pairs, all consistent with the curation.
- The 9 kept sites are all real: the curator cut there and rearranged pieces within the
  chromosome (H2_scaffold_6, H2_scaffold_46, H2_scaffold_47, H2_scaffold_4 45.3 Mb,
  H1_scaffold_1 109.8 Mb), or the two sides are the two haplotypes of one chromosome
  joined in one scaffold (H1_scaffold_1 at 12.4 Mb).
- Not found (11 sites): contact across is normal or one side is short (H1_scaffold_10
  34.9 Mb, H2_scaffold_8 3.5 Mb, H2_scaffold_10, H2_scaffold_44 11.4 Mb, H2_scaffold_64,
  H2_scaffold_6 46.7 Mb and the three small pieces at 25.8-26.1 Mb, H1_scaffold_118).

Next
1. Order and orientation of the pieces of one chromosome (the kept sites; most of the
   curator's 216 moves): cut every candidate, then chain the pieces by end contact within
   the pair, using the other haplotype as the guide.
2. Two haplotypes joined in one scaffold: a homolog line between the two sides of a site.
3. Troglobisium with sequence: minimap2 between haplotypes and small-against-all, BUSCO.

## 2026-11-13: order and orientation of pieces within a chromosome copy (allhaps)

Done
- `allhaps.arrange_pieces`: exact search (all orders and orientations, up to 9 pieces) for
  the arrangement with the most contact between neighbouring ends (`_adjacency`: end
  windows up to 100 px, outermost pixels skipped, relative to the decay expectation).
- In `curate_allhaps` each chromosome copy is split at its kept / tip-inverted sites;
  attached unlocs and "further pieces" of at least 40 px (60 / 40 / 20 by level) join the
  search. The result replaces the assembled arrangement only if the summed end contact
  rises by 0.15 (0.3 / 0.15 / 0.05). The largest piece keeps its orientation.
- Orientation: where flipping a piece changes the sum by less than 0.12 (0.2 / 0.12 /
  0.05) end contact cannot orient it. Then the other haplotype decides: `_homolog_dir`
  finds the piece's homology line in the other copy (sharpest diagonal or anti-diagonal,
  at least 2x the other direction) and the piece is turned to run the same way as the
  nearest clearly oriented piece. Without a clear line it keeps its orientation and the
  row says so.
- Rows: `rearrange` (new order, link scores, score as assembled) and `join` for an
  attached scaffold that went into the chain (no Unloc tag). `do_arrange=False` skips it.
- Fixed on the way: the contact-only candidate rule on the bird's sparse MQ10 map (a gap
  without a telomere must have reduced contact on the MQ0 map too); a site that is both a
  candidate on its own and the homologous point of a misjoin is cut again (H2_scaffold_1).

Why end contact alone is not enough for orientation: for pieces of 4-9 Mb the two ends
often touch the neighbour equally (H2_scaffold_70 with H2_scaffold_67: 1.24 and 1.23).
Which end of the neighbour they belong to is clear (1.2 against 0.5).

Result on Troglobisium, map only, scored with `tools/compare_layouts.py` (new junctions
between untagged contigs of 15 px or more; the curation has 49)
| level | cuts | rearranged copies | joins | tool junctions in the curation | curator junctions reproduced |
|---|---|---|---|---|---|
| 0 | 9 | 5 | 3 | 5 of 7 | 15 |
| 0.5 | 12 | 10 | 5 | 10 of 15 | 17 |
| 1 | 15 | 19 | 14 | 14 of 28 | 20 |
- Right at the default: X (25-58 Mb reversed, then 58-119 Mb), H2_scaffold_4 autosome
  half, the order of H2_scaffold_6 26-57 Mb, H1_scaffold_65 + H1_scaffold_32,
  H2_scaffold_61 + H2_scaffold_31, H2_scaffold_68 + H2_scaffold_11, H2_scaffold_71 +
  H2_scaffold_55, and the four tip inversions.
- Wrong or unresolved at the default: orientation of H2_scaffold_6 31-38 Mb and 0-8 Mb,
  of the 6 Mb piece of H2_scaffold_44, and of H2_scaffold_67 against H2_scaffold_70.
- Not attempted: pieces of one chromosome spread over several scaffolds that are not yet
  attached to the copy (H2_scaffold_46 into X between H1_scaffold_1 pieces, H2_scaffold_59
  / _69 / _8, H1_scaffold_2 with H1_scaffold_72), pieces under 40 px, the two haplotypes
  of one chromosome inside H1_scaffold_1 (0-12.4 and 12.4-25 Mb).
- Bird (with both PAFs): one rearrangement, H2_scaffold_34 joined in front of W
  (H2_scaffold_8), where the curator has it; cuts and tip inversions unchanged.

Troglobisium alignment (H2 against H1, minimap2 asm5, 6 threads): still running after
45 minutes when this was written; PAF will be `work/qqTroRaco/H2_vs_H1.asm5.paf`.

## 2026-11-13: Troglobisium with sequence (first look)

- 18 allhaps tests pass on the code with rearrangements (bird and Troglobisium; the run
  took 31 minutes because it shared the machine with minimap2). Single-mode tests not rerun.
- `tools/align_tro.py`: H2 against H1 with `minimap2 -x asm5`, 6 threads, finished after
  about 90 minutes. Fixture `seed/fixtures/qqTroRaco.H2_vs_H1.asm5.paf` (1.6 MB).
- The alignment is complete, not thin: 10596 alignments, 2161 Mb of H2 aligned, of which
  2147 Mb in the 1329 blocks of 10 kb or more (about 1.6 Mb per block). asm5 is adequate
  for this genome. Run time 94 minutes, peak memory 10 GB.
- Run with `--paf` and the two arachnida_odb10 BUSCO tables (`out_tro_seq/`): 69 pairs,
  12 cuts, 10 rearrangements, 15 pieces placed next to their homolog, 29 unpaired. The
  unpaired count is not comparable with the map-only run (3): with sequence every scaffold
  of 500 kb or more is considered, with Hi-C only those of 4 bins (about 2 Mb) or more.
  Layout score 9 of 15 tool junctions in the curation and 15 curator junctions
  reproduced, against 10 of 15 and 17 from the map alone. Why it is not better has not
  been looked at yet; cuts and rearrangements do not use the alignment at all so far.
- BUSCO: H1 2705 complete of 2726 found (211 duplicated), H2 2698 (200 duplicated). About
  200 duplicated genes per haplotype is far more than the bird's 25-32: worth a look for
  real duplications against false ones before using them as evidence here.

To do with sequence (the alignment is there and good; it is barely used yet)
1. Use alignment blocks as the orientation and order guide in `arrange_pieces` (now Hi-C
   only): this should settle the 4-9 Mb pieces whose ends tie.
2. Use them for cut decisions (the two sides of a site aligning to different places) and
   for chains across scaffolds of one chromosome.
3. Find out why the layout score with `--paf` is slightly below the map-only score.
4. Small-against-all alignment for Troglobisium (not run; the 4.4 Gb index needs more than
   this machine's 24 GB in one piece).

## 2026-11-13: alignment blocks as the orientation guide

Committed the work so far as da2a901 on branch `allhaps-cuts-joins-level`. The bird
(bAegFun) is the boreal owl, *Aegolius funereus* (curator).

- `_block_index`, `_aligned_dir`: direction of a piece relative to the other copy from the
  alignment blocks (mapq >= 30) that fall in the other copy's ranges: +1 / -1 when at
  least 70% of the aligned bases agree, and at least 100 kb or 10% of the piece aligns.
  Long blocks are clipped to the piece in proportion.
- In `arrange_pieces` the guide for orientation ties is now the alignment when `--paf` /
  `--fasta` is given, and the Hi-C homology line otherwise. End contact still decides
  when it is clear. The `rearrange` row names the pieces oriented this way and the source.
- Troglobisium with `--paf` at the default level: 11 of 15 tool junctions are in the
  curation (9 before the guide; 10 from the map alone), 17 curator junctions reproduced.
  Oriented by sequence: H2_scaffold_67, H2_scaffold_71, H2_scaffold_68, H2_scaffold_61,
  H2_scaffold_6 31.5-37.8 Mb.
- Still not in the curation: H2_scaffold_69 joined to H2_scaffold_8 (the curator has
  H2_scaffold_59, _77, _88, _81, _69 and the 3.5 Mb end of H2_scaffold_8 in a different
  order), H2_scaffold_6 0-8 Mb, the 6 Mb piece of H2_scaffold_44, and H2_scaffold_6
  26-31 Mb next to 38-57 Mb: the tool keeps these in one chromosome (end contact 1.13),
  the curator painted 26-38 Mb and 38-57 Mb as two units.
- Not done: alignments for the order of pieces, for cut decisions, and for chains across
  scaffolds of one chromosome.

## 2026-11-13: pairs of one chromosome merged (chains across scaffolds)

Why: after a correct cut (H1_scaffold_2 at 18.2 Mb) the two halves pair separately (with
H2_scaffold_59 and H2_scaffold_8) and were laid out as two chromosomes. The curator has
one. The same happens whenever a chromosome is in two scaffolds in both haplotypes.

Done
- `allhaps.merge_pairs`: for every two homolog pairs, the contact between an end of one
  copy and an end of the other, and the same for the corresponding ends of their two
  homologs (taken from each pair's strand); the score is the weaker of the two, best of
  the four end combinations. Computed on the map binned by 8 px (in memory, fast).
  A merge needs: score >= 0.5 (0.7 / 0.5 / 0.3 by level); no end with a telomere signal;
  score >= 3x (4 / 3 / 2) what either pair shows with its third-best partner; 1.2x any
  rival for the same end; each end used once.
- In `curate_allhaps` the smaller pair's units join the chains of the larger pair's copies
  (same prefix on each side), with their attached pieces; `arrange_pieces` then orders and
  orients everything (always applied for a merged chromosome). Row `merge`.
- At most 9 pieces per search: the smallest attached pieces are left as unlocs beyond that.

What separates real merges from clustering (measured, both genomes)
- Troglobisium, real: 0.92, 0.59, 0.58, 0.50; best false one 0.42. Third-best partner of
  the pairs involved: 0.03-0.13.
- Boreal owl microchromosomes (all false, MQ10 map): 0.96-2.60, but every one has a
  telomeric end, and the third-best partner is at 0.8-1.65: they touch many others alike.
- End specificity (best end combination against the second) points the wrong way: real
  merges of 9 Mb pieces have 0.59 against 0.54, the false telomere-to-telomere one 0.42
  against 0.03. Not used.

Result on Troglobisium, map only (sweep in `out_tro_sweep/`, with `compare.txt` per level)
| level | merges | tool junctions in the curation | curator junctions reproduced | curator chromosomes >= 90% in one tool group |
|---|---|---|---|---|
| 0 | 1 | 6 of 9 | 15 | 107 of 128 |
| 0.5 | 3 | 12 of 21 | 18 | 109 of 128 |
| 1 | 5 | 15 of 38 | 21 | 109 of 128 |
- The three merges at the default are one chromosome each in the curation: H1_scaffold_2
  halves with H2_scaffold_59 + H2_scaffold_8; H1_scaffold_10 halves with H2_scaffold_27 +
  H2_scaffold_62; H2_scaffold_60/H1_scaffold_63 with H2_scaffold_64/H1_scaffold_66.
- The strict junction score undercounts them: the curator has small scaffolds between the
  merged pieces (H2_scaffold_77, _88, _81, _69 between H2_scaffold_59 and H2_scaffold_8),
  and in the third case the curator's copies mix prefixes (H2_scaffold_64 with
  H1_scaffold_63) where the tool keeps one prefix per copy.
- With `--paf`: same three merges, 13 of 21 and 18.
- Boreal owl: no merge; output as before.

Not done: alignment positions for the order of pieces (order still comes from end contact
only), alignments for cut decisions, pieces under 40 px inside chains.

## 2026-11-13: end of session

Committed on branch `allhaps-cuts-joins-level` (da2a901, f1d1a6f, fb99581, and this note).
20 allhaps tests pass; the single-mode tests were not rerun this session. Curator's notes
for the next session (coverage for sex chromosomes, FASTGA, hapscotch) are in
`docs/HANDOFF.md` under "To think about next".

## 2026-11-14: curator's corrections of the level 1.0 sweep; alignment leads order and orientation

Fixture: `seed/fixtures/qqTroRaco*.savestate_curator_1_on_autocurate_level1` (the
curator's corrections on top of `out_tro_sweep/level_1.00`, which was a map-only run: the
alignment had not been given to the sweep). Against it the level 1.0 output had 21 of 38
new junctions right; the curator changed orientations of small and medium pieces and
undid wrong merges and joins (H1_scaffold_16 with H1_scaffold_54, H1_scaffold_24 with
H1_scaffold_48, H1_scaffold_70 between H1_scaffold_68 and H1_scaffold_5, ...).
Curator's point: the alignment blocks should set order and orientation; that would fix
most errors and guide further cuts and joins at gaps. Until now they only broke
orientation ties.

Done
- `allhaps._project`: position and direction of a piece along the other copy as laid out,
  from the alignment blocks (weighted by aligned bases; 70% strand majority; at least
  100 kb or 10% of the piece).
- Second phase in `curate_allhaps` when blocks are given: of the two copies of a
  chromosome the one in fewer pieces is the reference (arranged by end contact as
  before); the other copy's pieces are put in the order and orientation in which they
  align along it. Applied unless end contact prefers the contact arrangement by more than
  0.5 summed (0.3 / 0.5 / 0.8 by level); then the contact arrangement stays and the row
  says the haplotypes may differ there (`arrangement_differs` when nothing was moved).
  If the alignment says the assembled arrangement is the colinear one, a contact-based
  rearrangement is withdrawn.
- A tip inversion decided from the high-MQ map and telomere position is not undone by
  the alignment, only flagged (owl H1_scaffold_2, the curator's instruction).
- Tried and rejected: letting end contact overrule the alignment at a difference of 0.15.
  Troglobisium fell from 18 to 15 junctions in the curation.

Result on Troglobisium with `--paf`, default level
| | tool junctions in it | its junctions reproduced |
|---|---|---|
| curator's final curation (49 new junctions) | 18 of 23 | 21 |
| curator's corrections of level 1.0 (40) | 18 of 23 | 23 |
Before (alignment for ties only): 13 of 21 and 18 against the final curation.
- Owl: one rearrangement as before (H2_scaffold_34 in front of W); two flags, nothing
  moved: the inverted end of H1_scaffold_2, and H2_scaffold_39, which aligns between
  H2_scaffold_4 0-83.3 Mb and its 4 Mb end piece (where the curator put it) while end
  contact on the sparse MQ10 map prefers the assembled order.

Open case from the curator (2026-11-14): the half-coverage sex chromosome of Troglobisium
is hard to read. A second piece looks like a very large PAR, but there is no visible
non-recombining part and coverage stays at half throughout. Not understood; the tool
should flag these pieces, not label them. To do: a table of coverage and of the
alignments between the X/Y pieces for the curator.

Next
1. Alignments to guide cuts and joins at gaps (the curator's second point): a gap where
   the two sides align to different places or orientations in the other copy is a
   candidate even when contact is normal (11 curator cuts are still missed); pieces that
   align next to each other in the other copy are join candidates.
2. Reference choice when both copies are in several pieces (now simply the one in fewer).
3. Merges checked by alignment (the curator undid two at level 1.0).

## 2026-11-14: the alignment proposes cuts and joins

Measured first (`H2 against H1`, blocks of 100 kb or more, mapq >= 30)
- The PAF blocks are long chains (median 51 kb, 10% over 1 Mb, up to 56 Mb), and they end
  exactly at the curator's cuts: H1_scaffold_2 at 18.12 Mb (H2_scaffold_59 | H2_scaffold_69),
  H2_scaffold_1 at 43.36 Mb (H1_scaffold_7 reversed | forward), H2_scaffold_47 at 12.2 Mb
  (the two halves of H1_scaffold_50 in swapped order).
- Breakpoints between consecutive blocks of a scaffold, by where they fall:
  - at a gap: owl 6, all curator cuts; Troglobisium 32, of which 18 at curator cuts and
    several more within a few pixels of one;
  - inside a contig: owl 38, Troglobisium 125, almost never a curator cut. These mark
    where the other haplotype is in several scaffolds (along H1_scaffold_2: H2_scaffold_77,
    _88, _69, _59, the chain the curator built), so they say where those scaffolds belong.
- A first attempt at pixel resolution found nothing: the breaks sit within one pixel of
  the gap and the test window straddled them. Breakpoints must be taken in base pairs.

Done
- `allhaps.alignment_breaks`: breakpoints between consecutive long blocks (other side
  changes scaffold, strand, or jumps). A change of scaffold where both alignments run to
  the ends of their scaffolds (within 300 kb) is fragmentation of the other haplotype and
  is not a break.
- `telomere_sites`: a breakpoint within 2 px of a gap makes that gap a candidate site even
  with normal contact; the usual decision follows (different homologs or no
  chromosome-level contact: cut; otherwise the site splits the copy and the arrangement
  step orders the pieces along the other haplotype). The evidence text quotes both sides.
- Loose scaffolds (in no pair, 4 px or more; 8 / 4 / 2 by level) that align over half
  their length to one copy of a chromosome go into the chain of the other copy at the
  place where they align (`along` in `curate_allhaps`; rows `join` and `rearrange`). Also
  for the reference copy, provided its own pieces stay as they are.
- A loose scaffold that aligns where its own copy already has sequence (half or more of
  the stretch) is a second copy, not a missing piece: it stays out of the chain (owl
  H2_scaffold_43 and H2_scaffold_70, both Unloc in the curation).
- Tip-inversion protection now applies only when the inversion was decided on a high-MQ map.

Result on Troglobisium with `--paf`, default level
| | sites at a curator break | tool junctions in it | its junctions reproduced |
|---|---|---|---|
| final curation (49 new junctions) | 28 of 37 | 22 of 32 | 27 |
| curator's corrections of level 1.0 (40) | 28 of 37 | 22 of 32 | 27 |
Before this step: 25 sites (all at breaks), 18 of 23 junctions, 21 and 23 reproduced.
19 scaffolds joined by alignment. The two curator states differ from each other in a few
places (H1_scaffold_2, H2_scaffold_59 / _77 / _8), so neither can be matched fully.
- Sites not at a curator break (9): one cut (H1_scaffold_2 21.8 Mb, or H2_scaffold_6
  46.7 Mb depending on the state compared) and 8 kept sites that only split a copy
  (H1_scaffold_16 0.5 Mb, H1_scaffold_43, H1_scaffold_46, H1_scaffold_63, H2_scaffold_6
  6.6 Mb, H2_scaffold_17, H2_scaffold_32, H2_scaffold_60).
- Owl: cuts unchanged (every alignment break at a gap was already a site). One new join:
  H1_scaffold_168 (0.5 Mb) into H1_scaffold_6 at the gap at 33 Mb, where H2_scaffold_5 has
  sequence that H1_scaffold_6 lacks and contact across the gap is zero. Not in the
  curation; for the curator to judge.

Next: breakpoints inside contigs as evidence for the other haplotype (which of two
conflicting copies is the reference: a contig is trusted over a scaffold join); merges
checked by alignment; the same with FASTGA once benchmarked.

## 2026-11-14: false duplicates tagged; curator's comments

Curator
- A loose scaffold that is a second copy of sequence its own haplotype already has can be
  labelled `FalseDuplicate` (sequencing error, a tandem repeat or a segmental duplication
  in one copy; whatever the cause, if it is not in the contig it is kept out).
- Contigs have priority over scaffolds. The scaffolds were built with YaHS, which is not
  perfect. Contigs are not perfect either (built from unitigs), but that is harder to use.
  For a future version the curator could provide the hifiasm unitig GFAs for hap1 and
  hap2; names differ from the assembly, so matching would need a graph alignment.
- The curator agrees with the join of H1_scaffold_168 into H1_scaffold_6 at 33 Mb (owl).

Done
- Second copies are tagged `FalseDuplicate` (not Unloc, no vertPaint), laid out at the end
  of their copy and painted with it; row `false_duplicate`. Two sources: a small scaffold
  aligned to both copies of its chromosome (`place_small`, needs `--small-paf`), and a
  loose scaffold that aligns where its own copy already has sequence (`along`).
  BUSCO-only candidates stay `false_duplicate_candidate` with vertPaint (a guide only).
- Owl with both PAFs and the BUSCO tables: 6 false duplicates (H1_scaffold_44, _192, _209,
  H2_scaffold_43, _70, _86), 2 unlocs, 2 joins (H2_scaffold_34, H1_scaffold_168).
  To ask: H1_scaffold_44 (2 Mb) aligns 70% to both H1_scaffold_36 and H2_scaffold_40 and is
  now a false duplicate; the original curation has it inside the chromosome next to
  H2_scaffold_40.

## 2026-11-14: duplicates and placement judged on non-repetitive sequence (owl H1_scaffold_44)

Curator: H1_scaffold_44 belongs inside the chromosome (flip it, next to H2_scaffold_40).
A large repeat array was assembled in different amounts in the two haplotypes; the
non-repetitive sequence after it is haplotypic. Repetitive sequence is already
"duplicated", so it should not weigh in the decision. Repeats can also be inferred from
contacts and from a repeat track (Red, added with PretextGraph: see HANDOFF).

What the alignments show (small-against-all)
- H1_scaffold_44, 1.96 Mb: 0-0.55 Mb aligns to itself only; 0.55-1.95 Mb aligns loosely
  (0.3-0.45 of the bases in matches, mapq 0) to H1_scaffold_36, H2_scaffold_40 and itself:
  the array. The last 262 kb aligns densely (0.83, mapq 60) to H2_scaffold_40 0-262 kb
  and has no comparable match in H1_scaffold_36.
- The five real duplicates align over their whole length to both copies at 0.91-1.00.
- The H2-against-H1 alignment misses the 262 kb match: H2_scaffold_40's first 262 kb is
  already taken by its long chain to H1_scaffold_36, and one target per stretch is reported.

Done
- `place_small`: where a small scaffold has a dense alignment (half or more of the bases
  in matches) over at least 50 kb and 10% of its length but under 90%, that part is its
  non-repetitive sequence and the chromosome, the copy and the duplicate call rest on it
  alone; another copy "has the same sequence" only if it aligns there at 0.7 of the best
  density. The evidence text gives the size of the non-repetitive part.
- Unique alignments of small scaffolds to the other haplotype are added to the block
  index used for order and orientation when the main alignment does not cover them.
- Owl result: H1_scaffold_44 is no longer a false duplicate. It is joined after
  H1_scaffold_36, reversed, parallel to the start of H2_scaffold_40 (in the H1 copy: only
  H2_scaffold_40 has that sequence). 5 false duplicates remain. New: H1_scaffold_31
  (4.98 Mb, 972 kb non-repetitive) attached as Unloc to H1_scaffold_33.

## 2026-11-14: H1_scaffold_44 into the H2 chain; H1_scaffold_31 to HAP2; v1.0

Curator
- H1_scaffold_44 belongs in the H2 superscaffold, next to H2_scaffold_40. This breaks the
  Hi-C phasing rule, but hifiasm's Hi-C phasing is not perfect, and it is more
  parsimonious to partition homologous sequence between the haplotypes than to have a
  segmental duplication in one against a complete deletion in the other.
- H1_scaffold_31 (hard even for the curator) seems to belong in HAP2, maybe as an unloc,
  judging by repeat patterns and the balance of haplotype content.
- Release a stable v1.0 for the team to pre-curate with; then Red mask, telo_scan.py,
  minimap2 against FASTGA, and a comparison with hapscotch (Chenxi Zhou). The curator
  trusts the work of Richard Durbin's group and anyone who has worked with him.

Done
- `place_small._extends`: a piece whose non-repetitive (dense) alignment lies at one of its
  own ends and at an end of the one copy that has that sequence, with the strand that lets
  it run on past that end, continues that copy: it is joined there (overlapping ends
  meeting), whatever its prefix (not under `--phasing trio`). H1_scaffold_44: its last
  262 kb overlaps the first 262 kb of H2_scaffold_40; laid out as H1_scaffold_44,
  H2_scaffold_40 (the curator's layout reversed as a whole). Pieces whose match is inside
  the copy still fill the hole in the other copy (H1_scaffold_168).
- A small scaffold whose two homologs both carry its own prefix goes to the shorter copy
  (balance). H1_scaffold_31 now attaches as Unloc to H1_scaffold_40, tagged HAP2.
- v1.0.0: `--mode auto` (default; allhaps when two prefixes each hold a fifth of the
  sequence), `--hap-regex`, `autocurate.review.tsv` (the items to look at, by importance),
  `autocurate.run.txt` (version, command, options), clear errors, warnings for missing
  tracks; README rewritten as a user guide; CHANGELOG.md.
- Whole test suite: 57 tests pass (single mode included, first full run in two days).

## 2026-11-14: tracks (telomere strands, repeats); minimap2 against FASTGA; hapscotch

Installed (conda): FASTGA 1.5.20260729 into `pretext-autocurate` (native arm64); Red
2015.05.22 and PretextGraph 0.1.0 into a new env `pretext-autocurate-x86` (bioconda has
them for Intel Macs only; they run under Rosetta).

Tracks (`tracks.py`)
- Bedgraph / BED files or map tracks become per-pixel arrays: `telomere_5p`,
  `telomere_3p` (and their sum as `telomere`), `repeat` (fraction of each pixel).
- telo_scan.py labels strands empirically: 5p = the strand enriched at scaffold starts,
  3p = at scaffold ends. Rule at a telomere site: 3p just left of the gap and 5p just right
  (both facing it) = two chromosome ends joined, cut; only the long side facing it = cut;
  only the end piece facing it = invert the end piece. Tested on the owl's H1_scaffold_2
  site with artificial peaks (`tests/test_tracks.py`).
- Repeat track: alignment breaks inside repeats (80% of the pixels around) are ignored;
  a small scaffold 90% repeat is not moved.
- CLI: `--telomere-motif` (runs telo_scan.py from `--telo-scan PATH`), `--telomere-5p`,
  `--telomere-3p`, `--repeats BED`, `--red`. The tool prints the PretextGraph command for
  adding computed tracks to a copy of the map.
- Not run here: telo_scan.py itself. Executing the downloaded script was refused by this
  session's permission settings; the curator can run it (or allow it). Its arguments and
  output names were checked against its source.

minimap2 against FASTGA on the owl (H2 against H1, 8 threads, this Mac)
| | wall | peak memory | records | aligned | in records >= 10 kb |
|---|---|---|---|---|---|
| minimap2 asm5 | 2.1 min | 8.1 GB | 25,772 | 1,260 Mb | 1,203 Mb |
| FastGA, PAF | 6.6 min | 4.8 GB | 400,618 | 16,364 Mb | 15,493 Mb |
| FastGA .1aln + ALNchain | 5.7 min + 2 s | 4.8 GB | 32,172 | 1,696 Mb | 1,617 Mb |
- FastGA reports every local alignment (repeats many times over) and no mapping quality
  (always 255). Chained (as hapscotch does for its plots) it is comparable in size.
- Calls on the owl with each alignment (map + PAF, nothing else): cuts, tip inversions,
  kept sites and the main pairs identical. With FastGA: 3 more small pairs and an extra
  unequal pair, 32 false duplicates instead of 2, 11 joins instead of 2. All from repeats:
  with no mapq the tool's uniqueness filter (mapq >= 30) passes everything.
- Troglobisium (2.2 Gb per haplotype): FastGA aborted after 4.3 min when the disk filled.
  Its k-mer index (FastK tables) for one haplotype alone was 25 GB; both plus temporary
  files need over 50 GB of free disk, which this Mac does not have. minimap2 needed 94 min
  and 10 GB of memory but little disk. On a cluster with scratch space FastGA is likely the
  faster of the two for large genomes; not measured here.
- To use FastGA, the tool needs its own uniqueness filter (e.g. best chain per stretch of
  the query, or drop stretches aligned to more than one place), or hapscotch's approach.

hapscotch: read and compared in `docs/RELATED.md` (all-against-all FastGA, synteny groups,
MILP phasing against inter-haplotype Hi-C, a shared Hi-C backbone per group, hapcount,
palindrome correction).

Curator's decision (2026-11-14): FastGA will not work on the curator's MacBook; stay with
minimap2 for the haplotype alignments.

## 2026-11-14: Red on the owl; v1.0.1, license, repository moved

- Released: v1.0.0 (main, tag), then v1.0.1 = the same code under the PolyForm
  Noncommercial License 1.0.0 with a CNAG notice (commercial use needs a license from
  CNAG). The repository moved to https://github.com/cnag-aat/pretext-autocurate (private).
- Red on the owl, each haplotype separately: 23.5 min, 8.3 GB peak, 706 Mb of 2.7 Gb
  called repeat (26%). `run_red` now deletes its per-haplotype FASTA copies afterwards.
- First run with the track changed four placements: H1_scaffold_44, H1_scaffold_31,
  H1_scaffold_310 and H2_scaffold_61 fell back to candidates because they are 90% repeat.
  Wrong by the curator's rule (repeat content must not decide; judge the non-repetitive
  part). The repeat guard now applies only when Hi-C alone places a scaffold. With that,
  the owl output with and without Red is identical: the track adds no decision on the owl,
  and is mainly for display (PretextGraph `repeat_density`) and the Hi-C-only cases.
- `CLAUDE.md`: always run commands under bash, not zsh (curator, 2026-11-14).
- `--add-tracks` (curator: Red may be redundant for decisions but the tracks are worth
  seeing; on a high-MQ map a repeat track tells repeat from contaminant; of BUSCO only the
  genes missing in one haplotype and the duplicates are worth showing, paired so the
  curator can tag a haplotig): a copy of each map with `repeat_density`, `5p_telomere`,
  `3p_telomere`, `busco_dup_pairs` (label per pair of scaffolds: d1, d2 ...) and
  `busco_missing_hap`, via PretextGraph. Owl: d1 H1_scaffold_14 / H1_scaffold_192 (9
  genes), d2 H2_scaffold_3 / H2_scaffold_43 (5), d3 H1_scaffold_12 / H1_scaffold_209 (4),
  d4 within H2_scaffold_1 (2); 188 gene copies missing from the other haplotype. The copies
  keep the header hash, so the same savestate loads on both. 60 tests pass.

## 2026-11-14: telo_scan.py on the owl (strand tracks)

- Permission: `.claude/settings.local.json` allows running
  `/Users/talioto/repositories/scripts/telo_scan.py` (curator's clone of cnag-aat/scripts).
- Toy check: a chromosome with telomeres at both ends gives 5p at its start and 3p at its
  end; two chromosomes joined end to end give both at the junction. At map resolution
  both fall in the gap pixel, so the strand check now includes the pixel across the gap.
- Owl run (`out_bird_telo/`, `--telomere-motif TTAGGG`, Red, BUSCO, `--add-tracks`):
  telo_scan.py plus the whole run took 2 min 54 s.
  - H1_scaffold_1 134.6 Mb: cut, telomeres from both sides face the gap (confirms).
  - H1_scaffold_2 166.3 Mb: invert_tip, the end piece's telomere faces the gap (confirms).
  - H1_scaffold_5 5.3 Mb: keep -> invert_tip. The original curation (savestate_6) inverts
    that 5.3 Mb piece: now matched.
  - H2_scaffold_2 5.1 Mb: keep -> invert_tip. The curator said to keep it (only the H1 tip
    of scaffold_2 needed inverting). Its telomere faces the gap as H1_scaffold_2's does,
    and the alignment says both are colinear as assembled. To ask the curator.
  - HAP1 choice changed for a few pairs (telomeric ends now from telo_scan).
- Troglobisium with telo_scan.py, motif TTAGG (`out_tro_telo/`, map + PAF, 2 min 8 s):
  the arachnid motif is clearly there (forward strand skews to scaffold ends 33,536 to
  92; 73 5p and 101 3p windows at score 100 or more), unlike the map's own telomere track.
  Four flagged sites become cuts from the strands: H1_scaffold_1 12.4 Mb, H1_scaffold_10
  37.8 Mb, H2_scaffold_4 45.4 Mb, H2_scaffold_6 37.75 Mb, all places the curator broke;
  H2_scaffold_17 0.47 Mb becomes a tip inversion (not a curator site). Against the final
  curation: 21 of 29 new junctions (was 22 of 32), 27 reproduced (unchanged).
- Curator (2026-11-14) on H2_scaffold_2 5.1 Mb: invert it. "If a curator wants to change
  it, then they have all the evidence explained before a decision is made." The tool
  already does this with telo_scan tracks; without them it keeps the site and flags it.
- telo_scan.py output in PretextGraph (tested on the owl MQ10 map): the names
  `5p_telomere` / `3p_telomere` are recognised as the telomere type (extension type 5),
  the same as `telomere`. The signal is sparse: 0.16% of pixels reach 10% of the maximum,
  so peaks stand out without post-processing. PretextGraph itself keeps 185 of the 399
  pixels that have any 3p window; dropping windows below telo_scan's own threshold (100)
  leaves 126 and removes weak interstitial runs. Optional, not needed.
- telo_scan.py bundled (curator's own script, cnag-aat/scripts commit 666994c) as
  `src/pretext_autocurate/vendor/telo_scan.py`, unchanged apart from a provenance header;
  `--telomere-motif` uses it by default (`--telo-scan` for another copy). Test: the bundled
  copy on a toy assembly gives 5p at a start, 3p at an end, both at an end-to-end join.

## 2026-11-14: v1.1.0 released; code commented; suggestions recorded

- v1.1.0 (evidence tracks: telo_scan strands, repeats, BUSCO pairs, `--add-tracks`): all
  61 tests pass; main fast-forwarded and tagged `v1.1.0` at 60e7665, pushed to
  cnag-aat/pretext-autocurate.
- Commenting pass over every module in `src/pretext_autocurate/` (not `vendor/`):
  docstrings for the functions that had none, inline comments at each step of the long
  functions (why a threshold or guard exists, units, which decision row it produces),
  stale docstrings corrected (allhaps.py still said joins were not done). Code unchanged:
  checked by comparing each module's syntax tree, docstrings removed, with the v1.1.0
  version.
- `docs/HANDOFF.md` brought to v1.1.0 (tracks.py and vendor/ in the module table,
  `--phasing` state, the curator's list with Red, telo_scan and hapscotch done) and a new
  section "Suggestions for further development" (18 items: contig-led reference copy,
  hemizygous chains from coverage, a shared Hi-C backbone, alignment-checked merges, a
  scoring script over all fixtures, fixture-free tests for CI, splitting allhaps.py ...).
- Correction to the `--add-tracks` entry above: the track copies keep the compressed
  header, but PretextView's header hash also covers the file name, so a savestate for
  `X.pretext` does not load on `X.tracks.pretext`. That is why `run` writes a separate
  savestate for each copy (the code always did; the note was wrong).

## 2026-11-15: contigs outrank scaffold joins (reference copy per pair)

Curator's rule (2026-11-14): when the two copies of a chromosome disagree, the copy that is
contiguous sequence there is the reference, not the copy in fewer pieces.

- `_contig_conflicts(raw, idx, pieces, layout)`: each contig (gap-free stretch, from the
  map's `gap` track) of one copy is followed through its alignments (stretches of 50 kb or
  more) into the other copy's layout. A contig that runs across two layout pieces must meet
  them in layout order and orientation; otherwise the aligned bp count as contradicted.
  Rearrangements inside one layout piece are not counted.
- Phase 2 (`curate_allhaps`): both directions are measured for each pair. The reference is
  the copy the other's contigs contradict less, when that is at least `contig_min`
  (`Level`: 1 / 0.5 / 0.2 Mb) and twice the reverse; when both exceed `contig_min` the
  fewer-pieces copy stays and the pair is flagged `arrangement_differs`. Every row of the
  pair says which contigs led.
- Troglobisium (map + PAF + TTAGG): contradictions in 10 pairs, up to 12.6 Mb (one 37.8 Mb
  contig of H2_scaffold_3 runs across H1_scaffold_68 and H1_scaffold_5 the other way round
  from the assembly). In every case the contradicting contigs belong to the copy in fewer
  pieces, which was already the reference: decisions and savestate unchanged, 21 of 29 tool
  junctions in the final curation as before. The rows now carry the contig evidence
  (rearrangements of H1_scaffold_5, H2_scaffold_47, H1_scaffold_10, H1_scaffold_63).
- The swap path was exercised by forcing the measurement the other way on H1_scaffold_10 /
  H2_scaffold_27: H1_scaffold_10 became the reference and the row explains why. No fixture
  genome needs it yet; it is for YaHS joins in the copy with fewer pieces that a contig of
  the other copy contradicts.
- Test `test_contig_conflicts`: H2_scaffold_3 contradicts H1_scaffold_5 + H1_scaffold_68 as
  assembled by over 5 Mb, the arrangement the tool applies by under 0.5 Mb; the
  H1_scaffold_5 rearrange row names H2_scaffold_3 as leading.

## 2026-11-15: sex chromosomes from coverage and alignment (`sexchrom.py`)

Curator (2026-11-15): hemizygous chromosomes from coverage, using the PAF when there is
one; PARs have normal coverage and behave like recombining autosomes; X/Y or Z/W needs
knowledge of the species; if two sex chromosomes are in one haplotype, align them to each
other to find PARs; tracks for PARs and sex chromosome candidates. Also: systems such as
X1X2Y, X0 (and the ZW equivalents); homomorphic sex chromosomes, or a sample of the
homogametic sex, look like autosomes.

Measured first (coverage track divided by its median over the large scaffolds; alignment
H2 against H1, mapq >= 30)
- Troglobisium: nearly every autosome has a 1-4 Mb stretch at about 0.4-0.5 coverage, in
  both copies at homologous places and aligned to each other (divergent heterozygous
  stretches). Half coverage alone is therefore not hemizygous. The sex chromosomes are
  what is left: H1_scaffold_1 58-119.5 Mb and H2_scaffold_46 (half coverage, not aligned,
  61 and 23 Mb), and H1_scaffold_1 25-58 Mb with H2_scaffold_4 0-33 Mb (half coverage in
  both, aligned to each other). The curator's final state has exactly this: X =
  H1_scaffold_1 25-119.5 Mb + H2_scaffold_46, Y = H2_scaffold_4 0-33 Mb. The
  normal-coverage ends H1_scaffold_1 0-25 Mb and H2_scaffold_4 33-58 Mb align to each other
  like PARs would, but have no contact with the sex chromosome (0.000 against 0.1-0.7 within
  one chromosome): an autosome joined to it in both haplotypes, as the curator said.
- Owl (ZW female): Z = H1_scaffold_3 6-96 Mb and W = H2_scaffold_8 0-28.4 Mb at 0.51, not
  aligned; their normal-coverage ends H1_scaffold_3 0-6.05 Mb and H2_scaffold_8 28.5-34.5 Mb
  align to each other and have chromosome-level contact with their own sex-specific part
  (0.69, 0.23): the PAR. The two PARs share MQ0 reads, so their MQ0 contact equals that
  within each.

Done
- `sexchrom.analyse`: per scaffold, smoothed coverage (0.25 Mb) and aligned fraction, runs
  of half-unaligned / half-aligned / other (runs under 0.5 Mb absorbed), half-coverage
  stretches joined across normal blips up to 1.5 Mb. Sex-linked: a stretch with at least
  `core_min` hemizygous bp (8 / 5 / 4 Mb by level, doubled without alignments), median
  coverage 0.38-0.68, and at least 30% of its scaffold (or 3x core_min); or a half-coverage
  stretch aligned mostly into one (the other sex chromosome's homologous part). Next to each
  sex-linked stretch: `joined` (contact under 0.08 / 0.05 / 0.03), `par` (aligned to the
  normal-coverage end of another sex-linked scaffold, or without alignment MQ0 contact with
  such an end of 0.5 or more), `fused` (aligned to an ordinary chromosome: neo-sex
  chromosome or misjoin), else `par_candidate`.
- Guards found on the way: at a 3 Mb core the owl's low-coverage chromosome ends (0.33-0.38,
  subtelomeric) and a 4 Mb heterozygous insertion in H2_scaffold_5 were called, and whole
  autosomes then looked like PARs. With the coverage band, the share of the scaffold and the
  4 Mb floor, levels 0, 0.5 and 1 find the same stretches on both genomes, with or without
  alignments.
- `homology.align_within`: minimap2 -X among a few scaffolds; `run` uses it with `--fasta`
  when one haplotype holds two sex-linked scaffolds (Troglobisium: H2_scaffold_4 and
  H2_scaffold_46, 1 min; they share only repeats, as a Y and an X piece should).
- `curate_allhaps`: rows `sex_linked_region`, `par`, `par_candidate`, `sex_boundary`
  (uncut joins and fusions), `sex_summary` (also when nothing is found, saying that a
  homogametic sample or homomorphic sex chromosomes look like autosomes); cut rows at a
  coverage step say "a sex chromosome joined to another chromosome". Hemizygous scaffolds
  in no pair join the chain of the sex chromosome copy they contact (`sex_chromosome_join`;
  contact 0.405 against 0.002 for H2_scaffold_46). A pair whose copies are both mostly
  sex-linked is a sex pair even when the medians alone would not say so.
- `--sex-system XY|ZW|X0|Z0|X1X2Y|XY1Y2|Z1Z2W|ZW1W2` (`sexchrom.sex_labels`): XY/ZW the
  larger of the homologous pair X (Z); X0/Z0 every copy X (Z); three-chromosome systems: the
  copy aligned to both others is the single one, the others numbered by size. Without it no
  X/Y/Z/W tags.
- Tracks `sex_chromosome` (2 hemizygous, 1 rest of a stretch) and `par` (2 PAR, 1
  candidate), written to `OUT_DIR/tracks/` and added to the map copies with `--add-tracks`.

Results
- Troglobisium (map + PAF + TTAGG, `--sex-system XY`): X tagged on H1_scaffold_1
  25-119.5 Mb and H2_scaffold_46, Y on H2_scaffold_4 0-33 Mb, both within 3% of the
  curator's tagged bp; the cuts at 24.96 and 33.29 Mb carry the coverage step. Junctions: 21
  of 30 in the final curation (21 of 29 before): the new one puts H2_scaffold_46 after the
  X's end, while the curator split H2_scaffold_46 at its gap (8.59 Mb) and put the first
  part into H1_scaffold_1's gap at 109.78 Mb.
- Owl (`--sex-system ZW`): both PARs, Z and W tagged (the W chain includes H2_scaffold_34,
  which end contact already put in front of the W); layout unchanged.
- 14 tests in `tests/test_sexchrom.py`.

Next
- Interleaving hemizygous pieces: split a chained piece at its gaps and place the parts by
  end contact (the Troglobisium X above). Done, next entry.
- Homomorphic sex chromosomes are not found from coverage; that would need another signal
  (k-mers from reads of both sexes, as in SCINKD), which is not among this tool's inputs.

## 2026-11-15: extra pieces split at their kept gaps (interleaving)

Curator's X in Troglobisium (savestate_final): H1_scaffold_1 25.0-58.5 Mb reversed,
58.5-109.78, H2_scaffold_46 0-8.65, H1_scaffold_1 109.78-119.55, H2_scaffold_46 8.65-23.07.
Both gaps carry a `keep` site (sides belong together, contact across 0.51 and 0.48 of
expectation). Phase 1 split a copy's own units at kept sites, but put extra pieces into the
search whole: here both H1_scaffold_1 58.4-119.5 Mb (an extra of the copy) and
H2_scaffold_46 (chained by contact), so neither gap could take anything.

Change: extra pieces (`extra_of`, unlocs above `chain_px`) are split at their kept sites
too; one that does not fit under the 9-piece limit split goes in whole, as before.

- Troglobisium: phase 1 now searches five pieces and picks the curator's order
  exactly (end contact 0.98/1.23/1.27/1.28, sum 4.75, against 2.66 as assembled; before:
  three pieces, 2.09). The H1_scaffold_1 boundary is at 109.84 Mb, one pixel from the gap.
  `tools/compare_layouts.py` against savestate_final: 24 of 32 tool junctions in the
  curation (21 of 30 before), 30 of the curator's 53 reproduced (27 before); all three X
  junctions reproduced, and the wrong one (X end to H2_scaffold_46) gone. Other decisions
  unchanged.
- Owl: savestate, AGP and decisions identical.
- Test: the X's layout order in `TestTroglobisium.test_curation`.

## 2026-11-16: shared Hi-C backbone measured (HANDOFF suggestion 3): no gain

Question: does pooling the Hi-C of both copies of a chromosome through the alignment
(hapscotch / YaHS "backbone") tell the gaps a curator cut from those kept better than one
copy? Measured before building anything: `tools/backbone_gaps.py`. For every gap inside a
scaffold, contact across it as observed / expected (windows of up to 60 px on each side,
2 px skipped); in the own copy, in the homologous windows of the other copy (when both map
colinearly into one scaffold), the weaker of the two, pooled (sums of observed and
expected), and pooled with the two cross blocks. Curator's call from the final savestate:
kept if the two sides are still neighbours the same way round. AUC = how well the measure
ranks cut gaps below kept ones.

| Map | gaps (homolog mapped) | own | other | min | pooled | pooled + cross |
|---|---|---|---|---|---|---|
| Troglobisium MQ0 (savestate_final) | 154 (137 kept, 17 cut) | 0.962 | 0.969 | 0.966 | 0.969 | 0.967 |
| owl MQ0 (curator_1_on_autocurate) | 15 (11 kept, 4 cut) | 0.955 | 0.955 | 0.955 | 0.955 | 0.886 |
| owl MQ10, same gaps | 15 | 0.795 | 0.727 | 0.773 | 0.773 | 0.773 |

- The two copies give nearly the same value at each gap: correlation 0.97 (Troglobisium
  MQ0), 0.99 (owl MQ0), 0.91 (owl MQ10); median difference 0.03, 0.01, 0.1. At this depth
  counting noise is not what makes a gap hard, and at MQ0 the reads of identical sequence
  are shared between the copies anyway, so the other copy is not independent evidence.
- Best single threshold on Troglobisium: 5 errors of 154 with the own copy, 5 pooled.
- Where the other copy is also a gap (scaffold join in both, the case the backbone was
  meant for): own already 0.995 (Troglobisium, 57 gaps) and 1.0 (owl, 7).
- The hard gaps are those with no mapped homolog: AUC 0.69 (Troglobisium, 23 gaps, 9 cut)
  and 0.60 (owl, 11 gaps, 6 cut). The backbone cannot help there: by definition there is
  nothing to pool.
- Limits: the map holds 8-bit display values, not read counts; owl numbers are small; the
  curator's call is read from adjacency in the final layout.

Decision: not built. Effort is better spent on gaps with no mapped homolog (divergent or
hemizygous sequence, unpaired scaffolds), where single-copy contact is weakest.

## 2026-11-16: merges checked by alignment (HANDOFF suggestion 4): not built

Question: with `--paf`, should a merge of two homolog pairs (P = a/b, Q = c/d; a joined to
c in one haplotype, b to d in the other) require the merged copies to be colinear across
the join? Measured on every merge `merge_pairs` makes on Troglobisium at level 1 and with
relaxed thresholds (link_min 0.2, specific 1.5), against savestate_final: where the last
3 Mb at each joined end align in the homolog, and whether any sequence of one pair aligns
into the other.

| Merge | Curator | Joined ends, aligned position from the homolog's end |
|---|---|---|
| H2_scaffold_8 + H2_scaffold_59 (H1_scaffold_2 halves) | one chromosome | 1.5-4.2 Mb inside |
| H1_scaffold_10 halves (H2_scaffold_27 + H2_scaffold_62) | one chromosome | 0.1 and 3.3 Mb |
| H2_scaffold_18 + H2_scaffold_66 | one chromosome | up to 1.2 Mb |
| H2_scaffold_9 + H2_scaffold_57 | one chromosome | up to 3.9 Mb; one end unaligned |
| H1_scaffold_24 + H1_scaffold_48 | undone | 0.0-0.4 Mb (end to end) |
| H1_scaffold_16 + H1_scaffold_54 | undone | 0.0-0.1 Mb (end to end) |

(H2_scaffold_60 + H2_scaffold_64 is one chromosome in the curation with mixed prefixes;
its ends are 1.0-2.5 Mb inside.) Cross alignment (a to d, b to c) is at most 0.5 Mb in
every case: no sequence spans a merge join, so the alignment cannot confirm one.

- The false merges align perfectly end to end, as two real chromosome ends assembled to
  the telomere in both haplotypes do; in the real ones the haplotypes broke at different
  places. A colinearity rule would keep the false merges and reject the real.
- Owl: no merge at any of these thresholds.
- Decision: not built. Possible guard, six cases only: joined ends that coincide in both
  haplotypes (within ~0.5 Mb) are likely chromosome ends; would need checking on more
  merges before use.

## 2026-11-16: hifiasm graphs and YaHS AGPs as evidence (exploratory)

Data (curator): `seed/fixtures/qqTroRaco_hifiasm_gfas_yahs_agps/`: hifiasm noseq GFAs
(hap1/hap2 p_ctg, p_ctg, p_utg, r_utg) and the YaHS AGPs of each haplotype. Script:
`tools/graph_junctions.py` (13 s).

How contigs reach the graph without sequence
- YaHS used all 546 hap1 and 440 hap2 contigs unbroken; the AGPs (scaffold_N =
  H1_/H2_scaffold_N; all 765 lengths agree with the map) put every contig in its scaffold,
  so every scaffold gap is a junction between two contig ends.
- The GFAs' `A` lines list the reads in each contig and each unitig; the outermost reads of a
  contig give the unitig, and its orientation, where the contig starts and ends.
- A junction (contig A's exit, contig B's entry, as laid out) is checked for: a direct link
  in hifiasm's haplotype contig graph (rare: 4 of 197), both ends in one connected
  component of the unitig graph (components are chromosome-sized, 400-5,000 unitigs), and
  the shortest oriented path between them (bp of unitigs walked, up to 2 Mb).
- p_utg and r_utg give the same results.

Junctions YaHS made, by the curator's verdict (savestate_final)
| | n | same component | oriented path <= 2 Mb | <= 200 kb | A's end a dead end |
|---|---|---|---|---|---|
| kept by the curator | 197 | 164 | 107 | 81 | 49 |
| broken by the curator | 24 | 17 | 1 | 0 | 15 |
- A path is near-certain support for keeping a gap (107 to 1). The exception, H2_scaffold_6
  at 46.67 Mb, is a haplotype switch: the curator moved h2tg000106l to H1_scaffold_78's
  contig, which has a path of the same length to it (a bubble).
- No path says little: 90 kept gaps have none, mostly at dead ends of the graph.
- Same component says little: 17 of the 24 broken gaps are within one component.
- The tool (out of `out_tro_split`): of the YaHS gaps it broke, 5 the curator kept, 1 of them
  with a path (H1_scaffold_16 0.45 Mb, an end piece the tool inverted); 15 the curator
  broke, none with a path. A path as a veto on cuts would fix that one and cost nothing.

Junctions the curator made that YaHS did not: 102, 23 with a path (H1_scaffold_10's
interleaving, H1_scaffold_72 into H1_scaffold_2, the small scaffolds between H2_scaffold_59
and H2_scaffold_77, H2_scaffold_70 + H2_scaffold_67 with a direct contig link, ...).

The graph as a proposer: from each contig end, the nearest same-haplotype contig end
reachable in the graph
- Contigs of 1 Mb or more at both ends: 46 proposals; 27 are YaHS junctions the curator
  kept, 7 junctions only the curator made, 2 YaHS junctions the curator broke, 10 neither
  (several within one scaffold: the graph would orient a contig differently; not checked).
- With small contigs: 614 of 790 proposals are in neither layout, nearly all contigs under
  0.3 Mb (below 5 px) that the curator left as standalone scaffolds.
- The tool's own joins: a path for 6 of its 18 right ones and 8 of its 29 wrong ones (some
  of the wrong ones put a small scaffold between two pieces where the graph allows it and
  the curator did not).

Conclusion
- Worth adding as supporting evidence for keeping gaps (and against cutting them), in the
  same role as telomeres: inputs `--gfa` (noseq is enough) and `--yahs-agp`.
- As a source of joins: review hints only. Sequence GFAs are not needed for any of this;
  they would only matter to align unitigs that hifiasm left out of both haplotypes.

## 2026-11-16: AGP as PretextView exports it; pretext-to-asm

The curator's post-processing is pretext-to-asm (sanger-tol/agp-tpf-utils): input FASTA +
PretextView's AGP export -> curated assembly files. It is wrapped, not reimplemented
(`asm.py`; `pretext-autocurate asm`, `run --pretext-to-asm`): its naming rules are
maintained there and it streams the sequence through a .fai index.

What pretext-to-asm needs, from its source and PretextView's GenerateAGP (commit in
docs/pretextview_commit.txt)
- `# HiC MAP RESOLUTION: <bp> bp/texel` in the header (it snaps breaks to the assembly's
  gaps within that tolerance; without it, it exits). Our AGP lacked it.
- Tags in columns 10+, every unknown tag read as a haplotype name. Our AGP carried
  `vertPaint` (PretextView leaves its paint-display tags out), which would have made a
  "vertPaints" assembly.
- `layout.write_agp` now mirrors GenerateAGP: pixel-derived coordinates, painted
  scaffolds first (Scaffold_<paint id>), "-" only inside painted scaffolds of several
  pieces, 100 bp U gaps.
- Scaffolds under one pixel are in no PretextView export, and pretext-to-asm groups such
  sequences by name prefix (h1, h2: 91 and 43 scaffolds, 4.4 Mb, on Troglobisium). They are
  now appended to our AGP tagged HAP1/HAP2 by prefix.

Results (pretext-to-asm 1.4.1, installed in the pretext-autocurate env with Python 3.12)
- Troglobisium, the tool's result (out_tro_split): hap1 2.284 Gb (68 autosomes, X 117.6 Mb,
  Y 33.4 Mb, 362 unplaced), hap2 2.076 Gb (68 autosomes, 240 unplaced), 5 false
  duplicates 2.8 Mb; with the 44,200 bp of gaps, exactly the input's 4.363 Gb; 97.0% in
  chromosomes; 41 s.
- Owl, the tool's result: hap1 with Z 96.0 Mb and W 38.4 Mb, 40 autosomes per haplotype,
  93.6% in chromosomes.
- Troglobisium, the curator's savestate_final: pretext-to-asm stops ("Expecting 2 abutting
  sub fragments but got 1") at a 1-pixel unpainted piece inside a contig (H2_scaffold_6
  25.83-25.90 Mb). A PretextView export would hold the same piece, so this is not ours to
  fix; to confirm with a real export.

## 2026-11-16: hifiasm contig graph in run and asm (`graph.py`)

Curator's design: sequence the graph shows adjacent needs no gap, and a small piece that is
the only way across a gap fills it. hifiasm's GFAs here have S, L and A lines only (no W, P
or J); links between contigs carry the read overlap (7-21 kb, median 16 kb), so closing a
gap means dropping the overlap from one side, not just the N's.

- Inputs: `--gfa` (hifiasm `*.hapN.p_ctg.gfa`, with or without sequence) and `--yahs-agp`
  (each haplotype; scaffold_N -> HN_scaffold_N). The AGPs say where each contig sits; the
  assembly FASTA (`--fasta`, or `asm`'s) gives the bases.
- At a YaHS gap: `link` if the graph links the two contig ends; `bridge` if the only link out
  of the left contig and the only link into the right one go to the same contig and that
  contig is a scaffold of its own. Overlaps are checked as an exact suffix-prefix match (the
  graph's length first, then within +-2 kb).
- `run`: such a gap is never cut, tip-inverted or used as a split point; rows `graph_link`
  and `graph_bridge` (the latter in the review list).
- `asm`: bridges go into their gaps in the PretextView-style AGP before pretext-to-asm (the
  host line is split at the gap; the bridge gets the host's tags); after it, every gap in
  its output AGPs whose two sides are graph-linked contig ends is closed (the right-hand
  component loses the overlap), the edited records are rebuilt from the input assembly and
  the rest copied; `<name>.graph_edits.tsv` lists every edit.

Troglobisium (noseq GFAs; tool's result at level 0.5)
- 4 links (H1_scaffold_18, _38, _44, H2_scaffold_23; overlaps 13.3-17.8 kb) and 2 bridges
  (h1tg274, H1_scaffold_229, 71 kb, into H1_scaffold_54; h1tg492, H1_scaffold_373, 29 kb and
  under one pixel, into H1_scaffold_3). Every overlap is exact at the graph's length.
- None of the tool's cuts fell on these gaps, so the decisions are otherwise unchanged.
- `asm`: 2 bridges inserted, 9 gaps closed: the 4 links, both sides of each bridge, and
  h2tg8-h2tg78, a join the tool made (and the curator too) that the graph confirms.
- Check (scratchpad `checkfa.py`, independent of the code): untouched records byte-identical
  to pretext-to-asm's without the graph; each of the 7 edited chromosomes equal to the old
  record with its gaps and overlaps removed, bridges equal to the input FASTA.
- Tests: `tests/test_graph.py` (synthetic, both orientations; Troglobisium junctions).

## 2026-11-16: hifiasm unitig graph: keep evidence and gap filling (`--utg-gfa`, `--fill-gaps`)

Data: `qqTroRaco.hfsm.asm.bp.p_utg.gfa` / `r_utg.gfa` with sequence (5.6 GB each; curator).
(The first files copied, `hap1.gfa` / `hap2.gfa`, were the haplotype contig graphs with
sequence, identical to `hapN.p_ctg.gfa`; they add nothing to the FASTA. The graph code reads
them all the same.)

Contig ends on the unitig graph: the outermost reads of each haplotype contig (p_ctg `A`
lines) and the unitig of every read (p_utg `A` lines). Unitig lengths in the sequence file
agree with the noseq file (417 of 417 checked).

Across the 221 YaHS gaps (p_utg; scratchpad `utgfill.py`)
- 106 have a unitig path between the two contig ends within 1 Mb: 77 a single path, 29
  with alternatives (bubbles, up to 20 paths). The 4 contig-graph links and 2 bridges are
  among them; 100 are new.
- In 105 of 106 the last 1 kb of the left contig and the first 1 kb of the right one each
  occur once in the shortest path's sequence: 97 have sequence between them (209 bp to
  391 kb), 8 overlap (the 4 links, the 2 bridges' sides as one, and 2 more: H1_scaffold_43
  8.95 Mb, 131 bp; H1_scaffold_50 22.68 Mb, 8.8 kb).
- Where the fill comes from: in 93 of 97 the reads of the inner unitigs are in other
  contigs, mostly the other haplotype's (e.g. H1_scaffold_1's gaps: reads in hap2 contigs);
  in 4 mostly in no contig. hifiasm's HG tags are all 'a' in this run, so no phase labels.

Decisions (`run --utg-gfa`, `_unitig_evidence`)
- Weaker than a contig link (107 of the curator's 197 kept YaHS gaps have a path, 1 of the 24
  broken: H2_scaffold_6 46.67 Mb, a haplotype switch where the alignment breaks). So a path
  only stops the gap from being a split point; a tip inversion or a contact-only cut there
  becomes a keep; a cut backed by the alignment or a telomere stands (it does at
  H2_scaffold_6 46.67). The homologous kept site in the other haplotype is held too.
- Troglobisium: H1_scaffold_16's 0.5 Mb end piece is no longer inverted (the curator kept
  it), nor H2_scaffold_17 rearranged to match; five kept sites leave the review list
  (backed by a path). But H2_scaffold_17's own site is an invert_tip on telomere strand +
  alignment, which now applies (the curator kept it as assembled): the count of YaHS joins
  the tool breaks and the curator kept stays 5. Junction score unchanged (24 of 32, 30 of
  the curator's 53).

Assembly (`asm --utg-gfa --fill-gaps`, tool's result)
- 2 bridges, 9 gaps closed by contig links, 2 more closed by an overlap along a unitig path,
  68 filled (3.2 Mb, `qqTroRaco.gapfills.fa`), 41 kept with the reason (bubbles, or no
  single path in the final layout's junctions). 4 min, 1.1 GB RAM.
- Check (scratchpad `checkfa2.py`): untouched records byte-identical; each of the 30
  edited chromosomes equal to pretext-to-asm's own records plus the fills, minus the closed
  gaps and overlaps.
- Fills are an option, not the default: most bring sequence hifiasm put in the other
  haplotype, which is right for a homozygous stretch and brings the other allele for a
  heterozygous one.

## 2026-11-16: gap fills without duplicating the other haplotype (read depth)

Curator: do not copy sequence that is in the other haplotype; if the only path is collapsed
homozygous sequence, duplicating it is right; a shorter alternative path or a direct join
should be looked for first; coverage should tell the cases apart.

- hifiasm's unitig read depth (rd:i:) is cleanly bimodal on Troglobisium (unitigs >= 50 kb,
  by length): haploid peak ~26 (2.6 Gb), diploid peak ~56 (0.5 Gb, collapsed homozygous).
  `UnitigGraph.haploid` takes the lower peak (25).
- Every path across each YaHS gap (up to 20, within 150 kb of the shortest; scratchpad
  `utgdepth.py`): 7 gaps have a direct join; of the rest, 21 have a collapsed path, 53 only
  haploid paths whose reads are in the other haplotype's contigs, 20 haploid paths in this
  haplotype's contigs or in none, 5 mixed.
- The map's Hi-C coverage track is no help here: about 0.9 at the other haplotype's copy of
  collapsed fills, where 2 would be expected for the assembly's only copy.
- Ground truth by alignment: the 72 candidate fills of >= 1 kb as the minimap2 reference,
  the whole assembly as query (asm5, identity >= 0.95; scratchpad `fillcopies.py`). 12 of
  the 15 collapsed fills are only in the other haplotype; the 3 others are also in their own
  haplotype (H1_scaffold_9 24.03 Mb and H2_scaffold_2 21.67 Mb, the same 124 kb; H1_scaffold_46
  17.73 Mb), all next to a circular contig (hifiasm 'c': a tandem repeat unit).

Rule (`graph.choose_fill`), per gap, over all paths:
1. the contigs overlap or abut along a path: close;
2. next to a circular contig: not filled (number of repeat units unknown);
3. a path whose inner unitigs have reads in this haplotype's other contigs (the two flanking
   contigs excluded, their reads reach into the next unitigs): this haplotype already has
   the sequence; the contig is named for placement, nothing copied;
4. the shortest path at collapsed depth (1.6-2.6x haploid) with at most 10% of its reads in
   this haplotype's contigs: fill;
5. otherwise not filled (haploid: one haplotype's allele; repeat; ends not found).

Result (asm on the tool's Troglobisium result): 2 bridges, 11 gaps closed, 15 filled
(1.4 Mb; was 68 and 3.2 Mb), 94 kept and listed, 18 of them with a contig of the same
haplotype named for placement. Each of the 13 fills of >= 1 kb matches an aligned candidate
present only in the other haplotype (own haplotype <= 4%, short repeats). `checkfa2.py`:
edited chromosomes equal to pretext-to-asm's records plus fills, minus gaps and overlaps;
the rest byte-identical.

## 2026-11-16: contigs of the gap's own haplotype placed along a unitig path

The 18 contigs named "place that contig here": a contig is put into the gap (like a
bridge, before pretext-to-asm; its two new junctions are then closed or filled like any
other) when it lies wholly on the path (its first and last unitig both on it, in order:
`graph._along`), is a scaffold of its own, is named for this gap only, and the gap is not a
repeat tangle (fewer than 20 paths). Junctions between two pieces of the layout (joins the
tool or the curator made) are included (`graph.place_between`, piece ends snapped to contig
ends within 1.5 px).

- Most named contigs extend beyond the path (e.g. h2tg416, 42 kb: its first unitig is on
  the path into the gap, its last is not): placing them whole would add sequence that is
  not in the gap. Left as suggestions.
- Placed on Troglobisium (tool's result): H2_scaffold_139 (152 kb) into H2_scaffold_40's
  gap at 18.65 Mb; H2_scaffold_276 (38 kb, under one pixel) between h2tg95 and h2tg129 (a
  join the tool made in H2_scaffold_6), both junctions then closed by overlaps checked in the
  sequence (15.7 and 13.3 kb). The curator had left both unplaced (2-3 px).
- A first version also put H1_scaffold_171 into H1_scaffold_25 at 12.65 Mb: a repeat tangle
  (20 paths, Hi-C coverage ~4); its junctions then made no sense. Hence the tangle guard.
- Overlaps longer than the 1 kb anchor (the right contig starts before the left one ends
  along the path) are now measured and checked in the sequence: 4 more gaps closed.
- Totals: 2 bridges, 2 contigs placed, 16 gaps closed, 15 filled (1.4 Mb), 91 kept and
  listed; `checkfa2.py` all OK.

## 2026-11-16: v1.2.0

Released: sex chromosomes, hifiasm's graphs (`--gfa`, `--yahs-agp`, `--utg-gfa`,
`--fill-gaps`), PretextView-exact AGP, `asm` / `run --pretext-to-asm`. 80 tests. What is left
is in `docs/HANDOFF.md` (suggestions 5-18, the open points of 9, known gaps).

## 2026-11-16: one scoring script over all genomes (baseline, v1.2.0)

`tools/score_all.py` runs the tool on every case of `tools/score_cases.json` (a sweep at
levels 0, 0.5, 1 for two haplotypes; reuses existing outputs in `score_out/` unless
`--rerun`) and scores each savestate against each reference curation with
`pretext_autocurate.score.compare`: junctions as in `tools/compare_layouts.py` (pieces of
>= 15 px next to each other inside a chromosome, Unloc / Haplotig / FalseDuplicate aside),
breaks (where neighbouring pixels of an input scaffold stop being neighbours; matched within
3 px), the share of reference chromosomes with >= 90% in one tool chromosome, and tag
agreement by bp. 9 min for the three genomes with cached alignments.

| genome | level | reference | junctions supported / reproduced | breaks P / R | chromosomes >=90% in one | HAP1 agree | sex tags agree | Unloc agree | edits tool / ref | review |
|---|---|---|---|---|---|---|---|---|---|---|
| xgTylRafi | - | curator | 25/29 / 41/68 | 94% / 47% | 89% of 18 | - | - | 83% | 241 / 357 | 69 |
| bAegFun | 0.00 | curator | 2/4 / 2/7 | 86% / 67% | 98% of 82 | 53% | 100% | 0% | 74 / 264 | 74 |
| bAegFun | 0.50 | curator | 3/5 / 2/7 | 75% / 67% | 98% of 82 | 53% | 100% | 11% | 76 / 264 | 74 |
| bAegFun | 1.00 | curator | 5/8 / 2/7 | 67% / 67% | 98% of 82 | 52% | 100% | 11% | 79 / 264 | 72 |
| bAegFun | 0.00 | corrections | 0/4 / 0/0 | 100% / 58% | 100% of 74 | 91% | - | - | 74 / 80 | 74 |
| bAegFun | 0.50 | corrections | 1/5 / 0/0 | 88% / 58% | 100% of 74 | 91% | - | - | 76 / 80 | 74 |
| bAegFun | 1.00 | corrections | 2/8 / 0/0 | 78% / 58% | 100% of 74 | 89% | - | - | 79 / 80 | 72 |
| qqTroRaco | 0.00 | curator | 19/25 / 25/53 | 82% / 47% | 88% of 128 | 78% | 41% | 0% | 178 / 335 | 36 |
| qqTroRaco | 0.50 | curator | 24/32 / 30/53 | 85% / 61% | 89% of 128 | 78% | 100% | 0% | 191 / 335 | 38 |
| qqTroRaco | 1.00 | curator | 29/42 / 34/53 | 81% / 68% | 91% of 128 | 79% | 100% | 0% | 218 / 335 | 46 |

References: xgTylRafi savestate_4; owl savestate_6 (curator, from scratch) and
savestate_curator_1_on_autocurate (the curator's corrections of an earlier tool output,
biased towards it); Troglobisium savestate_final.

- The trade-off the level controls shows on Troglobisium: from level 0 to 1, junction
  recall 25 -> 34 of 53 and break recall 47% -> 68%, junction precision 76% -> 69%.
- The owl needs few new junctions (7 in the curation); breaks and grouping are high.
- Unloc tags of the tool and the curators hardly coincide on the two-haplotype genomes
  (0-11%), with similar totals. See HANDOFF 12 (a)-(c).
- HAP1 against the owl's savestate_6 (53%) is expected, not an error: that curation kept
  each haplotype as assembled, while the tool tags the better copy of each chromosome HAP1
  (the curator's rule, not yet used in their own curations; curator, 2026-11-16). Against
  their corrections of the tool, which follow the rule: 91%.

## 2026-11-16: Unloc disagreement (HANDOFF 12a)

Curator: the curations are human work done by zooming and scrolling; they miss things and
predate some rules. The tool should follow the stated rules consistently.

What the two sides tag (scratchpad `unlocs.py`, level 0.5)
- Troglobisium: the curator has 47 Unlocs, all small (0.07-2 Mb, 1-30 px); the tool left 36
  of them unpainted, never weighed: without the small-scaffold alignment (`--no-small-align`)
  only Hi-C can place them, and Hi-C alone was only tried from 20 px (1.3 Mb). The tool's 4
  Unlocs are large (2-10 Mb); the curator joined all of them into chromosomes.
- Owl: of the curator's 6 Unlocs (savestate_6), 3 the tool tags FalseDuplicate (100% of each
  aligns to sequence its haplotype already has; BUSCO genes doubled), 1 agrees; the tool's
  one large Unloc, H1_scaffold_31 (5 Mb), the curator left unplaced.

Small pieces by Hi-C (scratchpad `smallhic.py`, `smallhi.py`)
- MQ0 contact ratio does not separate the curator's small Unlocs from the small scaffolds
  they left unplaced (Troglobisium, < 20 px: medians 2.3 and 2.2; coverage 1.7 and 2.2): at
  1-3 px and ratio >= 2.5, 164 would be placed, 7 in the curator's chromosome.
- Curator's rule (2026-11-16): a small piece is an Unloc when its high-MQ signal is unique
  to one homolog pair; unique at MQ0 as well is stronger evidence. On the owl the curator's
  small Unlocs are not unique at MQ10 (ratio 1.0-5; they rest on sequence), and many unique
  ones were left unplaced: the curation predates the rule. Extreme ratios (1e5) come with
  a handful of reads: uniqueness needs a floor on the signal.

Changes (`allhaps.place_small`, `curate_allhaps` phase 1)
- Pieces under `hic_alone_px` without sequence evidence: Unloc when the high-MQ contact of the
  best homolog pair is `hi_unique` (10 / 5 / 3 by level) times the next pair and at least
  0.5 mean contact (MQ0 agreeing if it is unique too); without a high-MQ map, MQ0 uniqueness
  (10 / 6 / 4). The repeat and coverage guards apply; if the rule does not reach an Unloc it
  says nothing (no candidate rows; other routes, e.g. alignment into a copy, stay open:
  without that, H2_scaffold_78 lost its alignment-based join).
- Attached pieces large enough for the arrangement search, when the copy itself is not
  rearranged: placed where their ends touch the copy best (copy's own pieces as assembled),
  if the best new end contact is >= `attach_link` (0.5 / 0.3 / 0.15), the place beats any
  other by `attach_margin` (0.3 / 0.15 / 0.05), and the orientation stays as assembled when
  turning it gains less than the usual orientation margin. Troglobisium: H1_scaffold_67
  placed after H1_scaffold_57 as the curator did (1.26 against 0.47 before it; 1.26 / 1.22
  reversed / forward: as assembled); H2_scaffold_65 stays an Unloc (0.79 after against 0.77
  before H2_scaffold_66: too close). Owl: H1_scaffold_31 placed after H1_scaffold_40 (0.44,
  nothing comparable elsewhere; the curator left it unplaced: worth a look).

Scores (`tools/score_all.py`), against v1.2.0
| genome | level | junctions supported | reproduced | Unlocs (tool) |
|---|---|---|---|---|
| qqTroRaco | 0 / 0.5 / 1 | 20/26, 25/33, 30/43 (was 19/25, 24/32, 29/42) | 25, 30, 34 of 53 (unchanged) | 10 at 0.5 (was 4) |
| bAegFun | 0 / 0.5 / 1 | unchanged | unchanged | 4 at 0.5 (was 3) |

New small Unlocs: Troglobisium 7 (2-15 px, MQ0 unique 6-16x; no high-MQ map; 17 in a
map-only run, where fewer pieces are paired), owl 2 (1-2 px, MQ10 unique 7-8x). Tests
updated: pieces the curator left unplaced are counted, not failed; Unloc/FalseDuplicate
pieces may sit between the two copies of a pair. The curator had left them unplaced: possible misses to check, or the rule
needs tuning. Unloc agreement with the curations stays low (0% / 11%) because those
curations do not follow the rule.

## 2026-11-16: tests without fixtures, GitHub Actions (HANDOFF 16)

- `pretext_autocurate/synth.py`: `write_pretext` writes a real .pretext (deflated header,
  BC4 textures with every mip level, `psgh` tracks; the BC4 encoder is the inverse of
  `decode_bc4`, exact on flat blocks, within 2 on smooth ones); `toy_map` builds the map of
  a toy genome (scaffolds made of chromosome pieces, gaps between them, contact decaying as
  1/(1+d) within a copy, the two copies of a chromosome sharing 35% of that along the
  homology line, as at MQ0). A 512-pixel map is written and read in 0.06 s.
- `tests/test_synthetic.py` (8 tests, 0.5 s): format round trip; on a toy diploid with A and
  B misjoined in H1, allhaps cuts the misjoin ("different homologs"), finds the three pairs,
  writes a savestate that is a permutation with every piece painted and HAP-tagged and each
  chromosome's copies side by side; the command line writes a savestate bound to the map,
  the decisions, review, run info and a PretextView-style AGP, also with `--sweep`; single
  mode cuts the misjoin; `score.compare` runs on the result.
- Fixture-dependent tests in `test_phase1.py` now skip without fixtures; its savestate quirk
  tests use the packaged blank template and run everywhere.
- Without fixtures: 92 tests, 25 run (1.2 s), 67 skip. `.github/workflows/tests.yml` runs them
  on Python 3.10 and 3.12 after `pip install .` (checks the package data and the console
  script too).
- Found on the way: `Signal.decay` assumed textures of at least `max_dist` (200) pixels and
  crashed on smaller ones; real maps have 1024. Two empty-mean warnings on tiny maps
  silenced (`allhaps._mean0`).

## 2026-11-17: joins only the graph sees (HANDOFF 9a)

Joins between the outer ends of two YaHS scaffolds (ends no YaHS join uses) that hifiasm's
graphs support, against the curator's and the tool's layouts (scratchpad `graphjoins.py`;
Troglobisium, 765 scaffolds, 1530 outer ends).

- Contig-graph links between two outer ends: 54; the curator made 1, the tool 0. Nearly all
  join scaffolds under 0.5 Mb to each other, and many are ambiguous (two links at one end:
  a bubble).
- Mutually nearest outer ends along the unitig graph (within 1 Mb): 144; the curator made 2,
  the tool 4; again mostly small pieces with each other.
- Those that extend a chromosome (one side >= 1 Mb) and are unique (a contig link alone at
  both ends, or no other end within 3x the distance): 18. Hi-C end contact (`_adjacency`
  on MQ0) against the best other large scaffold end for the small piece: the graph's end
  is the best for 3 (e.g. H1_scaffold_57 + H1_scaffold_156: 1.89 against 1.39), 6 pieces are
  under 2 px (no Hi-C), the rest touch another end more (small pieces touch many ends
  alike). The curator made 1 (H2_scaffold_71 + H2_scaffold_55, which the tool makes too).
  H2_scaffold_91 (0.66 Mb) is proposed at the ends of both H2_scaffold_27 and H2_scaffold_62:
  the two halves the curator joined into one chromosome, with small pieces between them.

Done: review rows `graph_join_candidate` (`graph.end_joins`, `graph.join_rows`; run with
`--gfa` / `--yahs-agp`, unitig paths with `--utg-gfa`), one per join the level's layout does
not make, with the graph evidence, the Hi-C contact against the best other end, and where the
smaller scaffold went. Troglobisium at level 0.5: 15 rows (the tool makes 3 of the 18). Not
acted on: the curator's verdict on some of them should decide whether graph-unique + Hi-C
best, or graph alone for pieces under 2 px, is enough to join. Scores unchanged.

## 2026-11-17/18: icLarUrsu against the curated primary; loose pieces inside a copy

When only the curated FASTA is kept (no curator savestate), score by alignment:
`tools/compare_assembly.py` projects every map pixel onto the curated assembly and reports
purity, junctions (right; right across sequence the haplotype lacks; unloc placed; wrong
chromosome; wrong order), curated chromosomes split, missed joins, a dot plot, and D-GENIES
input (`--dgenies`: a PAF lifted onto each layout plus query/target .idx files; a local
D-GENIES is in conda env `dgenies`). `tools/align_chunked.py` makes the alignment:
whole-chromosome `minimap2 -x asm5` on this weevil ran to 73 GB on the 24 GB MacBook (swap on
a nearly full disk: the earlier kernel panic). With 10 Mb query chunks, `-K 20M --cap-kalloc 1g
-f 0.002` and a 12 GB watchdog: diploid vs curated 25 min, H2 vs H1 12.5 min, 5.3 GB peak.
Chunk boundaries split alignments, and a jump can show at a boundary (H2_scaffold_8 at
30.00 Mb; a whole-chromosome alignment of that pair puts it at 29.65 Mb).

icLarUrsu (weevil, 310 YaHS scaffolds; curated primary 1.74 Gb): as assembled purity 97.0%,
32 curated chromosome x haplotype split, 67 missed joins. Tool at levels 0 / 0.5 / 1: purity
98.6 / 98.5 / 98.6%, split 27 / 27 / 23, missed joins 57 / 57 / 50. Most joins first called
wrong at every level were not: right order across sequence one haplotype lacks (SUPER_4), a
curated unplaced bit at a scaffold end (H1_scaffold_9 + 111 + 72 on SUPER_5), unlocs placed
next to their chromosome (SUPER_11, SUPER_22). The scorer now reports these separately.

One real error: H2_scaffold_53 and H2_scaffold_65 fill a hole inside H2_scaffold_8 (SUPER_8
~30-34 Mb), but the loose-scaffold chain adds only at ends, so both went to its start, 30 Mb
off. Now a loose piece whose stretch lies inside one uncut piece is left out with an
`insert_candidate` row. `allhaps/insert.py` cuts at a gap where the alignment jumps over the
hole so the piece can go in, but only when the gap is within 250 kb of the flanks and the copy
has little sequence of its own between them. It fires on no fixture yet: bAegFun
H1_scaffold_5 (gap 0.85 Mb outside overlapping flanks) and H1_scaffold_6 (1.76 Mb of its own
sequence for a 0.88 Mb hole) were wrong to cut, and icLarUrsu H2_scaffold_8's gap is 1 Mb from
the jump. Wrong-order joins, levels 0 / 0.25 / 0.5 / 0.75 / 1: 5 / 4 / 6 / 9 / 9 (were
7 / 7 / 8 / 10 / 10). H2_scaffold_53's row said "not moved" then "joined by end contact"; it is
placed by sequence alignment (repeat content does not count then) and now says so.

Still wrong at level 0.5: H2_scaffold_99 at the wrong end of H2_scaffold_1 (44 Mb),
H2_scaffold_103 (curated unplaced) on SUPER_7, H2_scaffold_27 + 32 (overlapping duplicate),
two on SUPER_18 around H2_scaffold_92, H1_scaffold_12 + 34 (SUPER_13). SUPER_5's first 34 Mb
disagree between both haplotypes and the curation (H1_scaffold_9's tip inverted;
H2_scaffold_4 reversed with a nested inversion); the cut at H1_scaffold_9 21.93 Mb follows
H2_scaffold_4's own misassembly. The graph's "0 gaps between linked contig ends" is real: no
contig link spans a YaHS gap in this assembly.

A decision map of the whole pass (rules, overrides, every route to each edit, holes) was
published as an artifact. Among the holes: a mirror cut survives when the hifiasm graph turns
its source site into a keep (sites.py:289-293 runs before evidence.py).

## 2026-11-18: adjacency solver, first version, next to the rules pass

`solver.py`, `run --mode allhaps --engine adjacency` (curator's request, after the decision
map: one decision, which contig ends sit together, instead of ~75 rules and 13 placement
routes). Nodes: contigs between gaps and sharp in-contig breaks (`misjoin.detect_misjoins` at
the level's `sharp_max`, as the rules pass). Contigs under two 8-px bins are merged into their
neighbour (they cannot be judged by contact, and alone were squeezed out by a join between
their neighbours). Edge score = yahs prior (2.5 / 1.5 / 0.75 at levels 0 / 0.5 / 1) + 1.5 x
log2(end contact / 0.25) + 1.5 x log2(contact / next best end, clipped -2..1) + 2 x homology
(+1 next to each other on the other haplotype, -1 far apart or the wrong way) + 3 graph - 1
per telomeric end; a hifiasm link at a gap is fixed. Joins on Hi-C alone must be specific
(2 / 1.5 / 1.25x the next end) and not touch a contig at >= 1.5x coverage. Greedy path cover,
paths paired by aligned bp, small leftovers to a chromosome by contact (2x the next) or by
alignment. Runs in 35 s on icLarUrsu (rules pass 2.5 min). Tuning on icLarUrsu found: at MQ0
each end's strongest contact is the same end of its homolog (1.2-1.85 on the toy), so edges are
within one haplotype only; misjoins sit at background (0.06-0.08), right joins near 1.

`tools/compare_assembly.py`: junctions with an Unloc-tagged piece are skipped (an unloc goes at
the chromosome's end by convention), and joins involving sequence the curation left unplaced
(not SUPER_*) are a note, not "wrong chromosome". `tools/score_all.py --extra` passes more run
arguments (`--extra "--engine adjacency" --out score_adj`).

| genome, level 0.5 | engine | junctions P / R | breaks P / R | chromosomes >=90% in one |
|---|---|---|---|---|
| qqTroRaco vs curator | rules | 30/43 / 34/53 | 86% / 68% | 96% of 128 |
| | adjacency | 15/20 / 24/53 | 78% / 39% | 88% |
| bAegFun vs curator | rules | 3/6 / 2/7 | 75% / 67% | 98% of 82 |
| | adjacency | 0/4 / 1/7 | 100% / 11% | 95% |
| icLarUrsu vs curated FASTA | rules | right 16, wrong order 6, wrong chr 0 | purity 98.6% | split 29 |
| | adjacency | right 10, wrong order 1, wrong chr 0 | purity 98.4% | split 33 |

So far the solver is about as precise and finds much less. What it lacks that the rules pass
uses: the high-MQ map (the bird's breaks: only 26 of its YaHS joins are at gaps), the
small-against-all alignment and BUSCO for small scaffolds, alignment breaks at gaps (the
other haplotype's own scaffolding errors make homology near-circular: it backs H1_scaffold_9
4.39 Mb, a curated misjoin, because H2_scaffold_4 has one too), pairing without alignment,
false duplicates, sex chromosomes, HAP choice by quality.

Second round (same day). The high-MQ map was measured first and left out: at the bird's YaHS
gaps MQ10 contact across the junction is 0.00-0.44 at the five curator misjoins and 0.00 at
eight joins the curator kept (reads near gaps rarely map uniquely), so it does not separate
them; the rules pass's bird cuts rest on telomere strands, alignment breaks and different
homologs, not contact. Added instead:
- telomere strand: an end whose telomere faces outwards (`telomere_5p` at a head, `_3p` at a
  tail, the pixel across the gap included) is a chromosome end, -2.5 per end (two ends facing
  the gap: -5); without strand tracks the old unstranded -1 per end.
- contig-aware homology: support or contradiction counts fully only when it rests on the other
  haplotype's contigs (both ends land on one contig; or the contig an end lands on goes on),
  and 0.25 when it rests on that haplotype's YaHS joins, which can share the misjoin.
- ends project through the largest alignment block near the end, carried on to the end (the
  outermost was often a subtelomeric repeat matching another chromosome's end: bAegFun
  H1_scaffold_1's chromosome B start landed on the end of A's homolog).

| level 0.5 | rules | adjacency v1 | adjacency now |
|---|---|---|---|
| bAegFun breaks P / R | 75% / 67% | 100% / 11% | 100% / 67% |
| bAegFun junctions | 3/6 / 2/7 | 0/4 / 1/7 | 2/5 / 2/7 |
| qqTroRaco breaks P / R | 86% / 68% | 78% / 39% | 85% / 47% |
| qqTroRaco junctions | 30/43 / 34/53 | 15/20 / 24/53 | 15/19 / 24/53 |
| icLarUrsu right / wrong order / missed | 16 / 6 / 58 | 10 / 1 / 57 | 11 / 4 / 53 |

Third round (same day): exact search and a no-YaHS test.
- exact search: after the greedy path cover, each path is re-arranged by a windowed Held-Karp
  (8 contigs, links to the neighbours included, 3 passes) on contact + specificity + telomere +
  graph, without the YaHS prior; a re-arrangement is kept only if it gains `rearrange_gain`
  (2 / 1 / 0.5 at levels 0 / 0.5 / 1). Only the collapsed-repeat veto applies there (the
  specificity veto blocked true turns). Offline it found the curator's order for qqTroRaco
  H1_scaffold_26 (0 1 2r) and the turn in H2_scaffold_1. `solve(trace=dict)` returns the
  nodes, edges and paths for such checks.
- `--no-yahs`: no prior, every contig its own node (no merging of small ones), YaHS edges face
  the same veto and bar as any other.

| level 0.5 | adjacency + exact search | same, `--no-yahs` |
|---|---|---|
| qqTroRaco junctions P / R | 28/34 / 28/53 | 37/43 / 22/53 |
| qqTroRaco breaks P / R | 66% / 63% | 28% / 74% |
| qqTroRaco chromosomes >=90% in one | 88% | 82% |
| bAegFun breaks P / R | 67% / 67% | 62% / 89% |
| bAegFun chromosomes >=90% in one | 95% | 91% |
| icLarUrsu split / wrong order / missed | 33 / 4 / 53 | 43 / 4 / 52 |

Why YaHS joins were cut without YaHS: icLarUrsu (197 joins) 107 kept, 56 specificity veto,
29 contig under two 8-px bins, 0 contact below adjacent, 5 other; qqTroRaco (186) 55 kept,
90 contig too small, 27 specificity, 2 contact below adjacent, 12 other. The contact almost
never contradicts YaHS; YaHS supplies what the binned map cannot resolve. Next: score small
contigs at full map resolution and revisit the specificity bar for neighbours, then rerun.

Fourth round (same day): small contigs at full resolution (`--fine-px N`, off by default).
Contigs under two 8-px bins are no longer merged into a neighbour (down to N px) and their end
contact is read from the full-resolution map: each small node's end rows (up to half the node,
gap pixels left out) against every end's window (up to 64 px), observed / expected from the
full-resolution decay. One strip per small node; the solver takes ~75-85 s instead of ~35 s.

Without YaHS, the YaHS joins kept rise from 107 to 123 of 197 (icLarUrsu) and from 55 to 120
of 186 (qqTroRaco); the specificity veto is now the main reason a YaHS join is lost (63 / 36).

| level 0.5 | with YaHS | + fine 2 px | without YaHS | + fine 2 px |
|---|---|---|---|---|
| qqTroRaco junctions P / R | 28/34 / 28/53 | 35/43 / 25/53 | 37/43 / 22/53 | 30/37 / 22/53 |
| qqTroRaco breaks P / R | 66% / 63% | 48% / 66% | 28% / 74% | 41% / 74% |
| qqTroRaco chromosomes >=90% in one | 88% | 88% | 82% | 81% |
| bAegFun junctions | 2/5 / 2/7 | 1/7 / 1/7 | - | 1/4 / 1/7 |
| bAegFun breaks P / R | 67% / 67% | 67% / 67% | 62% / 89% | 62% / 89% |
| icLarUrsu right / wrong order / split | 11 / 4 / 33 | 16 / 4 / 34 | 13 / 4 / 43 | 17 / 5 / 44 |

A floor of 4 or 8 px instead of 2 changes little (qqTroRaco with YaHS: breaks 52% / 53%).
With YaHS, the 15 new qqTroRaco breaks the curator did not make all sit next to pieces of
2-36 px (H1_scaffold_43 4.39-4.93 and 8.65-8.79 Mb, H2_scaffold_40 18.71-19.84 Mb,
H2_scaffold_55, _61): small pieces taken out of place, which the curator left; not yet checked
in PretextView whether they are misplaced.

Fifth round (same day): mutual-best rule. A Hi-C-alone join that misses the specificity margin
is still allowed when its two ends are each other's best same-haplotype end and their contact
ratio is at least `mutual_min` (1.0 / 0.85 / 0.7 at levels 0 / 0.5 / 1; ~1 = adjacent
sequence). Measured first on qqTroRaco without YaHS, on the YaHS joins the specificity veto
blocked: the 8 mutual-best ones at contact >= 0.97 were all kept by the curator; the 2 the
curator broke had contact 0.70 and 0.77.

| level 0.5, with mutual-best | with YaHS | + fine 2 px | without YaHS | + fine 2 px |
|---|---|---|---|---|
| qqTroRaco junctions P / R | 27/35 / 28/53 | 34/45 / 25/53 | 40/48 / 24/53 | 31/41 / 23/53 |
| qqTroRaco breaks P / R | 66% / 63% | 48% / 66% | 30% / 74% | 44% / 76% |
| qqTroRaco chromosomes >=90% in one | 88% | 88% | 84% | 82% |
| bAegFun junctions | 3/23 / 3/7 | 2/27 / 2/7 | 3/22 / 3/7 | 2/27 / 2/7 |
| bAegFun breaks P / R | 75% / 67% | 67% / 67% | 62% / 89% | 57% / 89% |
| icLarUrsu right / wrong order / split | 12 / 4 / 33 | 17 / 4 / 34 | 14 / 4 / 43 | 18 / 5 / 44 |

On the bird the rule makes ~18 more joins, nearly all wrong. The ones seen so far join
chromosome ends (H1_scaffold_1 0.25 Mb to its own 264.4 Mb end; two scaffold starts):
chromosome ends cluster at MQ0, so two ends can be each other's best. Next: keep the rule away
from scaffold ends / telomeric ends.

`tools/mark_disagreements.py MAP REF OUT_DIR TOOL...`: a copy of the map (same file name, so
every savestate still opens) with one `vs_ref_<label>` track per tool savestate marking extra
breaks (4), wrong joins (3), missed breaks (2), missed joins (1) as `score.compare` counts them,
the savestates with that track shown, numbered waypoints when a savestate had none, and a TSV.
For the curator's review: work/review/qqTroRaco and work/review/bAegFun (rules, adjacency,
mutual-best at level 0.5).

Sixth round (same day): fixes from the curator's review of the bird (work/review/bAegFun).
- H1_scaffold_8 / H2_scaffold_7 unpaired: mutual-best joined the starts of H2_scaffold_6 and
  H2_scaffold_7 (contact 0.93, 1.06x the next end; the curator: subtelomeric repeat contact), and
  that path paired with H1_scaffold_7. The starts carry a telomere in the map's own track (192,
  128) that both telo_scan strand tracks miss (0), and are 100% Red repeat for 12+ px.
  `tracks.load` now keeps the map's track as `telomere_map`, and `telomere_ends` uses it where
  the strand tracks see nothing. Mutual-best does not apply at a telomeric end, and no join is
  made on Hi-C alone at an end whose outer 8 px are mostly repeat (`repeat_ends`).
- painting only for putative chromosomes (curator): paired paths, and unpaired ones only at
  0.6-1.5x coverage and under half repeat (`unpaired_cov`). On the bird this unpaints ~50
  repeat-rich pieces of 0.8-3.7 Mb at ~0.35x coverage. Unloc tags only on painted pieces.
- pairs in parallel, never an X (curator): the second copy is turned by alignment strand
  (`_projections` returns `sense`, signed bp per node pair); the old test compared node
  positions and needed two nodes per path, so one-contig pairs stayed crossed. Crossed homolog
  scaffold pairs (>= 5 Mb aligned): bAegFun 10 of 29 -> 0, qqTroRaco 18 of 73 -> 4 (rules 1 / 4).

| level 0.5 | mutual-best | + review fixes | without YaHS, + fixes |
|---|---|---|---|
| bAegFun junctions P / R | 3/23 / 3/7 | 1/4 / 1/7 | 1/2 / 1/7 |
| bAegFun breaks P / R | 75% / 67% | 67% / 67% | 62% / 89% |
| qqTroRaco junctions P / R | 27/35 / 28/53 | 27/35 / 27/53 | 40/48 / 23/53 |
| qqTroRaco breaks P / R | 66% / 63% | 66% / 63% | 29% / 74% |
| icLarUrsu right / wrong order / split | 12 / 4 / 33 | 11 / 4 / 33 | 13 / 3 / 43 |

With the X fixed the solver now also cuts the H2_scaffold_2 tip at its 5.10 Mb gap and turns it
to match its homolog, as the curator asked. Not done: turning a contig to match its homolog in
general (review 5-6, 10-11), an all-vs-all PAF for mis-phased pieces (12-14: H1_scaffold_40
belongs in HAP2), keeping cuts out of repetitive shrapnel. New review copies:
work/review/bAegFun_v2, work/review/qqTroRaco_v2.

Correction (same day): the map's own telomere track (tidk) is not a fallback for telo_scan.
Every bAegFun scaffold start where it peaks and telo_scan does not (51 of 51) has no
TTAGGG/CCCTAA in its first 2 kb (most start with a ~190 bp satellite array); the scaffold before
it in the map ends in ~330 motifs per 2 kb. The peak is that telomere, in the pixel the two
scaffolds share. `telomere_map` and the fallback are removed again; the H2_scaffold_6 / _7 join
stays out on the repeat-end guard alone (bAegFun level 0.5: junctions 1/5 / 1/7, breaks 67% /
67%, no pair crossed). Listings: work/review/bAegFun_v2/telomere_*.txt. The small-vs-all PAF
already shows the mis-phasing the curator found (H1_scaffold_40 aligns 3.0 Mb to H1_scaffold_33,
same haplotype); the solver reads only the H2-vs-H1 PAF so far.

All-vs-all alignment of the bird (curator's choice over H2-vs-H1 + small-vs-all): `tools/align_chunked.py
--no-chunk --secondary --mm-args "-D -N 5 -p 0.5 -K 20M --cap-kalloc 1g"`, the diploid FASTA
against itself (minimap2 -D skips a sequence's own diagonal only for a query of the target's
name and length, hence whole sequences). 19 min, minimap2 peak 7.9 GB (watchdog, top: 12.0),
65,084 alignments: work/bAegFun/all_vs_all.asm5.paf. Across haplotypes 2.40 Gb primary; same
haplotype, other scaffold 110 Mb primary + 42 Mb secondary; same scaffold off the diagonal
218 + 111 Mb. Largest same-haplotype scaffold pairs: H1_scaffold_33 / _40 5.15 Mb and
H1_scaffold_36 / _44 2.59 Mb, the two mis-phasings the curator found by eye; then
H1_scaffold_31 / _33 1.76 Mb, H1_scaffold_33 / _62 1.35 Mb, and below 1 Mb mostly shrapnel
against shrapnel. Not used by the solver yet.

All-vs-all in the solver (same day). `--paf` may now be an all-vs-all PAF: `solve` keeps one
direction of each alignment (`_one_direction`, only when the file has both) and sorts blocks by
haplotype into across (projections, pairing, as before), within (other scaffolds of the same
haplotype) and self. An H2-vs-H1 PAF gives the same result as before.
- `rephase`: a scaffold moves to the other haplotype when its best homolog is a larger scaffold of
  its own haplotype (>= 300 kb and 30% of it, 2x the second, more than it aligns across), at
  0.6-1.5x coverage, and < 30% of the partner region it copies has an other-haplotype homolog (it is
  the missing copy; with one there it is a third copy and stays). On the bird: H1_scaffold_40
  only (4.69 Mb on H1_scaffold_33, 1% of that region with an H2 homolog), the curator's call.
  H1_scaffold_44 (on H1_scaffold_36, whose region H2_scaffold_40 covers 100%) is not moved.
  Row `rephase`; in the review list.
- no new join between two scaffolds of one haplotype that align >= 300 kb ("two copies of one
  region in one haplotype").
- repeats: without a Red track, `copies_track` (>= 3 alignments over a pixel) stands in for it.
- no cuts inside shrapnel: a YaHS join between two contigs that are mostly repeat is kept as
  assembled unless an outward telomere (telo_scan) marks a chromosome end there (bAegFun
  H2_scaffold_30 2.84 Mb, H2_scaffold_31 2.55 Mb are curator cuts of that kind).

bAegFun level 0.5, all-vs-all vs H2-vs-H1: the same scores (junctions 1/5 / 1/7, breaks 67% / 67%,
no pair crossed) and the same disagreements with the curator; pairs 39 vs 40. After the move,
H1_scaffold_33 has two H2 homologs of equal weight (H1_scaffold_40 1.74 Mb, H2_scaffold_37 1.61 Mb
at MAPQ >= 30; 9.5 Mb aligned within itself), so the 2x pairing rule leaves all three unpaired;
H1_scaffold_40 and H2_scaffold_37 are likely two H2 pieces of one chromosome. Next: pair one path
with several, and join pieces of one haplotype that tile one homolog. Review copies:
work/review/bAegFun_v3.

All-vs-all on the other genomes (2026-11-18/19). Whole-sequence all-vs-all (`--no-chunk`, `-D`)
passed the 16 GB cap on qqTroRaco and icLarUrsu within 6 min. With 10 Mb chunks and
`--drop-self` (each chunk's alignment to its own place dropped after the fact): qqTroRaco 31 min,
9.3 GB; icLarUrsu 58 min, 11 GB. But a chunk's best hit is then itself, so the other haplotype's
copy is a secondary, and minimap2 gives secondaries MAPQ 0: the first qqTroRaco run had no pair at
all (HAP1 agreement 0%). align_chunked.py now works MAPQ out again without the self hits, from the
chain scores as minimap2 does (`--requalify` for an existing file). After it: qqTroRaco 3.9 Gb
across the haplotypes at MAPQ >= 30 (both directions; the H2-vs-H1 file 2.1 Gb in one);
icLarUrsu only 1.95 Gb (both directions) against 1.6 Gb in its H2-vs-H1 file: on this repeat-rich
genome `-p 0.5` against the chunk's own score drops many of the other haplotype's copies.

Level 0.5, adjacency engine:
- qqTroRaco all-vs-all: junctions 27/34 / 28/53, breaks 66% / 63%, HAP1 agree 78%, review 12 (H2-vs-H1,
  earlier code: 27/35 / 27/53, 66% / 63%, 76%, 11). Rephased 2 (H2_scaffold_72, H2_scaffold_172).
- icLarUrsu: H2-vs-H1 with the current code as before (right 11, wrong order 4, split 33, missed
  joins 53); all-vs-all fewer joins (13 vs 28): right 3, wrong order 0, split 34, missed joins 64.
  The weaker cross-haplotype signal costs homology-backed joins. Next: across-haplotype blocks from
  the H2-vs-H1 file with the all-vs-all for within-haplotype and self blocks, or realign this one
  whole-sequence on a larger machine.
- bAegFun all-vs-all, current code: the same junctions and breaks (1/5 / 1/7, 67% / 67%); HAP1
  agreement 30% -> 40% (curator), 66% -> 76% (corrections), because unplaced shrapnel is now tagged
  by its haplotype (below).

Haplotype tags on shrapnel (solver). The solver tagged only painted paths; pretext-to-asm reads an
untagged H1_scaffold_N as haplotype "H1", which is not "HAP1" (it compares names ignoring case only:
h1 vs hap1), so unplaced pieces would form separate H1 / H2 assemblies, and a re-phased scaffold
would go back by its name. Now every contig without a pair's tag gets its haplotype's HAP tag after
re-phasing, as the rules engine does.

Non-target sequence: `--contaminants FILE` and `--tiara` (contam.py). A list (name in the first
column; comma, space or tab; optional class column) or Tiara's calls on the scaffolds up to
`--tiara-max-mb` (10). They are tagged Contaminant and nothing else (no HAP tag, no paint; sub-pixel
ones in the tool's AGP too), never joined, cut, paired, placed as unlocs (solver `exclude`);
`OUT.contaminants.tsv` keeps the class, and after pretext-to-asm the organelles are moved out of its
contaminants FASTA into `*.mito.fa` / `*.plastid.fa` (also `asm --contaminants`). PretextView has no
organelle tag, and pretext-to-asm reads an unknown tag as a haplotype name.
Tiara (bioconda; on Apple silicon an Intel env, pytorch 1.7, libcxx 14) alone calls nuclear
sequence non-eukaryotic, so a call needs other evidence (`contam.confirm`):
- qqTroRaco: 29 calls. Two real foreign genomes: GC 0.22, 2.6 Mb, single copy (H1_scaffold_76,
  H2_scaffold_92, H1_scaffold_91, H1_scaffold_101, H2_scaffold_152), and GC 0.35 (genome 0.315),
  2.2 Mb, aligned only to each other (H1_scaffold_85 ~ H2_scaffold_96 ...). Their pieces contact
  each other at 105-255 per bin pair, any chromosome at 2-66. The curator had four of the first
  as Unlocs of H1_scaffold_45 / H2_scaffold_42 (their strongest chromosome contact, ~1/4 of the
  contact within the group): to review. False calls align to chromosomes (H2_scaffold_100 413 kb
  on H1_scaffold_3) and contact them.
- icLarUrsu: no bacteria; 22 organelle calls: 11 "organelle, unknown" at GC 0.27 and 4-10x (a
  high-copy AT-rich repeat), 7 "mitochondrion" at 0.5-3x (NUMT-like), 4 "plastid"; the curated
  icLarUrsu1.2 keeps all but the 4 smallest as unplaced scaffolds. bAegFun: one 55 kb
  mitochondrion (H2_scaffold_97).
Rules: prokaryote / plastid only if its second-best contact with another call (bins next to its
own left out: one neighbour in map order shares contact) is >= 2x its best chromosome contact; a
mitochondrion only at >= 5x coverage; an organelle of unknown kind never; none that aligns to a
chromosome along >= 20% of its length. Result: qqTroRaco 15 tagged (both genomes), 12 flagged
(`flag_contaminant`; 9 align to chromosomes); icLarUrsu 0 tagged, 22 flagged. The single-linkage
Hi-C rule alone passed the icLarUrsu repeat: pieces of one collapsed repeat contact each other too.
