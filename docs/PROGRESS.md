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
