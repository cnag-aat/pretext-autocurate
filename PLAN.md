# pretext-autocurate: project description and plan

Written for a Claude Code session to execute. Read `README.md`, `docs/FORMATS.md` and
`docs/CAE71_scaffold1_case.md` first; they hold what is already known and verified.

## 0. Curator's priorities (2026-11-12; these override the ordering in sections 6 and 7)

Automate the easiest work first, in this order. Language: **Python first**; port to Rust only
once the methods settle.

1. **Haplotypic duplications (haplotigs).** Signs: off-diagonal lines of contact parallel to
   the diagonal, lower long-read coverage, and position at gaps.
2. **Placement.** Put sequences into chromosomes: localise them where order and orientation
   are clear, otherwise attach them as `Unloc`.
3. **Obvious cuts and joins.** Cut at gap positions by preference (pretext-to-asm snaps to the
   exact gap downstream). A cut inside a contig needs a clear, reported reason. A/B
   compartment and other chromatin checkerboard patterns are normal and must not trigger
   cuts. Compare MQ0 with the high-MQ map to recognise repetitive, poorly mappable regions.
4. **Contaminants.** True lack of contact with everything else plus different coverage; also
   the contact pattern of a circular bacterial genome.

5. **Chromosome orientation.** Orient every chromosome short arm (p) to long arm (q) where
   the evidence allows: inter-chromosomal centromere-centromere and telomere-telomere contact
   (Rabl-like clustering) to locate centromeres, plus the telomere track(s). This must degrade
   gracefully: holocentric genomes (e.g. butterflies) have no centromere signal and some
   genomes (e.g. flies) have no conserved telomeric repeat. When neither signal is clear,
   leave the orientation alone and say so in the report.

**Use genome biology as a prior.** The tool should use what is known about genome structure
to steer calls toward a biologically sensible result (expected karyotype shape, telomeres at
chromosome ends, one centromere per monocentric chromosome, sex-chromosome behaviour), and
report when the map disagrees with that expectation instead of forcing it.

**Sex chromosomes.**
- Pseudoautosomal regions (PARs) are shared by X and Y (or Z and W): they look like a
  homologous pair, have autosome-like coverage, and contact both sex chromosomes. Never call
  a PAR copy a haplotig or false duplication, and never cut a sex chromosome at the PAR
  boundary because coverage or homology changes there.
- Hemizygous sex-chromosome sequence has about half the autosomal coverage in the
  heterogametic sex. That is not a haplotig signal and not a contaminant signal.
- Pseudohaploid (single) assembly of a heterogametic individual: keep the full X or Z, plus
  only the non-recombining part of Y or W.
- Diploid (allhaps) assembly: put the whole X and Y (or Z and W) into HAP1. HAP2 holds the
  alternate autosomes. (The bAegFun fixture follows this: Z and W are both tagged HAP1.)

**MQ0 maps, repeats and telomeres: curator's notes (2026-11-13).**
- Telomeric and subtelomeric contacts between different chromosomes are common, especially
  in MQ0 maps. Strong end-to-end contact at a telomeric end is not evidence for a join.
- A high-MQ map is not always available (it can mean re-downloading reads). When there is
  one, comparing it with MQ0 separates unique-mapping (real) contacts from multi-mapping.
  Without it, coverage and the Hi-C pattern itself (saturation) can flag repetitive regions.
- MQ0 is preferred for curation because it distinguishes different types of repeat arrays,
  and joins can sometimes be made from matching patterns at junctions. (Not used by the
  tool yet: a possible later method for joins at repeat-rich ends.)

**allhaps first step: curator's spec (2026-11-13).**
- First, put the haplotypic sequences of each chromosome next to each other, oriented the
  same way (homolog line parallel to the diagonal, not an X pattern).
- Tag sequences prefixed `H1_` / `H2_` with HAP1 / HAP2. Most assemblies have no parental
  data, so the "best" haplotype of a pair can be HAP1 (fewer gaps, more complete, telomeres).
- Identify potential sex chromosomes if possible; without homology information this is
  hard and not required.
- One pass is preferred. A second pass (write each haplotype, remap, make haploid maps,
  run single mode on each) is acceptable as an option for complicated cases.

**Phasing and haplotype assignment: curator's spec (2026-11-13).**
- New option `--phasing none|hic|trio` (no hyphen in `hic`):
  - `trio`: the assembly is fully phased. No reassignment between haplotypes at all.
  - `hic`: phased per chromosome only (Hi-C phasing does not work across chromosomes).
    Never mix `H1_` and `H2_` prefixes within one chromosomal superscaffold. Whole
    superscaffolds may be swapped between HAP1 and HAP2, so the better copy of each pair
    can be HAP1. This is the setting for the bAegFun bird.
  - `none`: unphased. Mix and match freely to build the best HAP1, as long as it stays
    reasonably balanced with HAP2.
- Exception for `hic` (curator, 2026-11-13): a small piece (an unloc) may go to the other
  prefix's copy of its chromosome when the evidence points there. Not under `trio`.
- In every mode, avoid introducing haplotype switches inside a superscaffold.
- With low heterozygosity the assembler sometimes collapses a chromosome (one copy where
  two are expected). Do not force a pair; note it in the report.

**Misjoins and orientation in allhaps: curator's note (2026-11-13).**
- Use telomere signal as evidence: a telomere inside a scaffold suggests two chromosomes
  joined end to end; telomeres at ends help orient and confirm chromosome ends.
- Use the other haplotype's contiguous sequence as evidence (as hifiasm's dual scaffolding
  does): where one haplotype is contiguous across a point, a break or different arrangement
  in the other at the homologous point is suspect, and vice versa.

**Bird microchromosomes ("dot chromosomes"): curator's notes (2026-11-13).**
- They need more matching up between haplotypes than the tool does so far. They are hard
  to work with: gene rich, but also repetitive, with non-B-DNA structures, and they
  cluster in the nucleus (high contact between different microchromosomes).
- They are conserved and well known, so expect a set of small paired chromosomes.
- The curator's `savestate_curator_1_on_autocurate` did not finish them; the original
  `bAegFun*.savestate_6` is the more complete reference for microchromosome pairs.
- scaffold_2 (both haplotypes) is NOT a misjoin at the telomere-bearing gap: keep joined.
  The H1 tip needs inverting. Swapping the tips between haplotypes looks better in the map
  but would mix haplotypes within a superscaffold, which `--phasing hic` forbids: do not.

**Sequence as input: curator's proposal (2026-11-13).**
- The assembly FASTA (with `H1_`/`H2_` prefixes) can be given as an argument. The tool can
  then run minimap2 internally (H1 against H2, or all against all) and read homology
  directly, instead of inferring it from MQ0 Hi-C signal.
- If a BUSCO lineage is given, its genes can be aligned to the assembly with miniprot.
  Uses: homology, false duplications, and choosing the most complete / least duplicated
  haplotype as HAP1.
- Design: `--fasta` (runs minimap2), `--paf` (use an existing alignment), `--busco-lineage`
  or a miniprot GFF/PAF. Hi-C homology stays as the fallback and as a cross-check. Tools
  are external: detect them, say clearly when missing, never download silently.

**Telomere detection: curator's note for a future feature (2026-11-13).**
- `telo_scan.py` (https://github.com/cnag-aat/scripts/blob/main/telo_scan.py) is the lab's
  production telomere detector for a given telomeric motif. It gives fewer false positives
  than the current map track and reports the forward and reverse-complement motif
  separately (so the strand tells which way a chromosome end points).
- Future PretextMap files will carry its output as extension tracks. It could also be run
  from this tool on `--fasta` with an option (e.g. `--telomere-motif`). Not implemented.

**BUSCO: curator's notes (2026-11-13).**
- BUSCO evidence must be entirely optional; nothing may depend on it (the search is slow).
- If BUSCO was already run, feed its `full_table.tsv` (`--busco-table`). It must have been
  run on the same assembly as the map (the diploid one for allhaps). In a diploid assembly
  most genes are reported "Duplicated"; count the copies per haplotype from the rows.
- Useful in single mode too (duplicated genes as haplotig / false-duplication evidence);
  not wired into single mode yet.
- BUSCO and OrthoDB are not perfect: some species have real duplications and really
  missing genes. Use BUSCO as a guide or tie-breaker only. No hard rule may rest on a
  single duplicated or missing gene.
- Per-haplotype tables are fine: `--busco-table H1=... --busco-table H2=...` adds the
  prefix. The Complete/Duplicated label is ignored; copies are counted from the rows.

**Inversions: curator's notes (2026-11-13).**
- If a putative inversion scores about the same inverted or not, do not break and re-join:
  it is probably a heterozygous inversion. Flag it.
- Inversions can be nested; an iterative approach (fix the outer one, re-test the inside)
  may be needed.

**Haplotigs: curator's notes (2026-11-12, to come back to after placement).**
- Cutting haplotigs off the ends of contigs is a valid approach. Pixel-level cuts are a
  little imprecise; that is acceptable for now. Later, a separate program can help
  pretext-to-asm make precise edits from long-read alignments (BAM).
- allhaps with more than two haplotypes present in a "diploid" assembly (as in the CAE71
  protist): the job is to separate sequence into HAP1, HAP2, and either `Haplotig` or
  `FalseDuplicate`. This separation matters.
- allhaps with exactly two haplotypes: for each chromosome, place the two haplotypes side
  by side in PretextView, and label the better-assembled one HAP1. "Better" means fewer
  gaps, more complete sequence representation, and presence of telomeres. (Sex chromosomes
  are the exception: both go to HAP1, see above.)

Inversions, SV flagging and haplotype-aware ordering come after these.

Ground truth for these tasks: `seed/fixtures/xgTylRafi*` (132 haplotigs, 59 unlocs, 72 cuts,
1 contaminant; single haplotype, MQ0 only) and `seed/fixtures/bAegFun*` (well-behaved allhaps
bird, MQ0 and MQ10). Both maps are 65536 px, so a full-resolution matrix is 4.3 GB: work
from mip levels and decode full-resolution tiles on demand.

## 1. Goal

A command-line tool that reads PretextMap `.pretext` files and does most of a first curation
pass automatically, handing the curator:

- a PretextView **savestate for each input map** (e.g. MQ0 and MQ10/MQ30 of the same assembly),
- an **AGP** of the proposed assembly and a TSV of every decision,
- a **report** (HTML with figures, plus JSON) explaining each edit and each thing it flagged
  but did not touch.

The curator opens the savestate, reviews waypoints and the report, and fixes what is left.
YaHS-style scaffolding is not the aim; the aim is the manual work that comes after it.

### Modes

- **allhaps** (primary): two or more haplotypes scaffolded separately and concatenated into one
  map. Needs the MQ0 map (homology shows up as off-diagonal signal) and benefits from a
  high-MQ map (haplotype-specific signal).
- **single**: one haplotype or a collapsed assembly.

### What it should do

| Task | allhaps | single |
|---|---|---|
| Group scaffolds into homologous sets (chromosomes) | yes | n/a |
| Order and orient pieces within a set, haplotypes side by side | yes | order/orient only |
| Detect misjoins (at gaps first, then anywhere) and cut | yes | yes |
| Detect inversions and fix them | heterozygous-looking and assembly errors | homozygous / assembly errors |
| Flag structural variants between haplotypes without "fixing" them | yes | n/a |
| Place unlocalised pieces (unlocs) | yes | yes |
| Find and tag haplotigs / false duplications | yes | yes (needs coverage) |
| Tag likely contaminants | yes | yes |
| Paint scaffolds, tag HAP1/HAP2/..., drop waypoints at things to review | yes | yes |

## 2. What already exists (session 1)

- A Python reader for `.pretext` (header, BC4 textures, graph tracks), a savestate
  reader/writer that round-trips the curator's file byte for byte, and a pixel-level simulator
  of PretextView's move/invert/break edits: `seed/scripts/pretext.py`.
- The savestate header hash, verified against the curator's file: `seed/scripts/header_hash.py`.
- Exploratory analyses that worked on one hard region (see the case document):
  - `homology.py`: forward/reverse diagonal tracing in off-diagonal blocks of the MQ0 map.
  - `junction2.py`, `junctionR.py`, `junctionL.py`, `distcurve.py`: "arm" contact tests. Sum
    contacts over all assembled copies of each arm, compare by distance band against cis
    references and background.
  - `scan.py`, `scan2.py`: genome-wide scan for where a loose end makes contact.
  - `write_state.py`: turning decisions into edits, paint, tags and waypoints.
- Not yet done: anything automatic. Every piece, breakpoint and arm in those scripts was
  chosen by hand.

## 3. Constraints to design around

- Map values are 8-bit, BC4-compressed and normalised by PretextMap. Use ratios and ranks
  within a map, never absolute thresholds. Every threshold must be derived from the map's own
  distance-decay and background, or be a tunable parameter with a data-derived default.
- Resolution is one pixel (genome length / 32768 at the usual setting). Scaffolds shorter than
  a pixel are invisible; pieces of a few pixels carry little signal. Report these as
  "too small to assess" instead of guessing.
- At MQ0, reads in identical regions are split at random among copies, so evidence must be
  summed across copies ("collapsed") before testing a join.
- No sequence is available by default. Optional inputs should improve calls when given:
  FASTA or a per-scaffold table (GC, length), telomere/gap/coverage/repeat tracks (already in
  the map as graphs), a taxonomy or blob table for contaminants, PAF alignments between
  haplotypes.
- The savestate is tied to the map's file name and header (see FORMATS.md). The edit list is
  shared between maps of the same assembly; only the hash differs. The tool must check that all
  input maps have identical contig layout and refuse otherwise.
- A curator must be able to trust the output. Prefer flagging to editing when evidence is
  marginal; never hide an edit.

## 4. Architecture

Recommendation: **Rust** for the tool, keeping the Python as a reference implementation and
test oracle. Decoding 528 textures and working on a ~1 GB matrix per map is where speed matters,
and a single static binary is easy to hand to other curators. If Rust slows early progress,
phases 1-3 can be done in Python first and ported once the algorithms settle; decide after
phase 1 and record the decision.

    crates/
      pretext-io      .pretext reader (parallel texture decode), graph tracks, savestate
                      read/write, header hash, AGP/TSV writers
      mapstate        pixel permutation model: break / move / invert, contig list,
                      conversion of a target layout into a minimal PretextView edit list
      signal          distance-decay P(s), background, per-pixel bias, collapsed-copy sums,
                      cis reference curves
      detect          homology, grouping, misjoin, inversion, unloc, haplotig, contaminant
      plan            turn detections into an ordered, non-conflicting edit plan with confidence
      report          HTML + JSON + figures
      cli             the `pretext-autocurate` binary
    python/           reference implementation (from seed/), used in tests
    tests/            fixtures, simulated scrambles, regression cases

Memory: a full u8 matrix is n^2 bytes (1.07 GB at 32768). Either hold one map at a time, or
work from mip levels and decode full-resolution tiles on demand. Start simple (full matrix),
measure, then optimise.

## 5. Command line

    pretext-autocurate run  --map mq0.pretext [--map mq30.pretext ...]
                            --mode allhaps|single
                            [--hap-regex '^(H\d+)_']        # how to read the haplotype from a name
                            [--start-state in.savestate]     # continue from a curator's edits
                            [--out-dir DIR] [--prefix NAME]
                            [--apply high|medium|none]       # which confidence tier becomes edits
                            [--tasks misjoin,inversion,group,unloc,haplotig,contam]
                            [--fasta F | --table T.tsv] [--paf P] [--blob B.tsv]
                            [--threads N] [--config params.toml]
    pretext-autocurate info     MAP            # header, contigs, tracks, hash
    pretext-autocurate dump     STATE --map M  # human-readable savestate, resulting layout
    pretext-autocurate convert  --map M (--agp A | --state S) --to (state|agp|tsv)
    pretext-autocurate simulate --map M --seed N --ops K   # scramble a curated map for testing

Outputs per run: `<prefix>.<mapname>.savestate` for each map, `<prefix>.agp`,
`<prefix>.decisions.tsv`, `<prefix>.report.html`, `<prefix>.report.json`, `<prefix>.params.toml`
(the parameters actually used, for reproducibility).

Tunable parameters (all in the TOML, the common ones also as flags): diagonal window and minimum
length for homology traces; score thresholds expressed as multiples of background or fractions
of the cis reference; minimum arm length for junction tests; distance bands; confidence cut-offs
for edit vs flag; minimum piece size to act on; coverage ratio bounds for haplotig calls;
whether to trust gap positions as preferred cut sites; maximum number of waypoints (255 limit).

## 6. Methods to implement (and how each is tested)

1. **Signal model.** Per map: distance-decay curve from long scaffolds, background from
   distant inter-scaffold pairs, per-pixel coverage bias. Everything downstream is expressed
   relative to these.
2. **Homology map (allhaps).** Generalise `homology.py`: for all scaffold pairs, find
   off-diagonal forward and reverse traces in MQ0, ideally using MQ0 minus scaled high-MQ signal
   to suppress true cis contacts. Use coarse mip levels to find candidate pairs, full resolution
   to refine. Output a set of colinear blocks (a PAF-like table). If a real PAF is supplied,
   use it instead and keep Hi-C homology as a cross-check.
3. **Grouping.** Build a graph of scaffolds linked by homology blocks and by cis contact;
   connected components are chromosome groups. Expect and handle groups where pieces
   share only part of their length (as in the CAE71 case).
4. **Junction test.** The core primitive from session 1. Given a candidate breakpoint and the
   homologous positions on other copies, define arms, sum over copies, and compare each arm pair
   with the cis reference at matching distance and with background. Returns a likelihood ratio
   for "adjacent" vs "not adjacent" and an estimated separation.
5. **Misjoin detection.** Run the junction test at every gap, then at every point where a
   homology block ends or changes partner, then at local minima of a sliding insulation score.
   Cut when adjacency is rejected with high confidence; flag otherwise.
6. **Inversions.** Candidates come from reverse-oriented homology blocks (allhaps) and from
   the anti-diagonal pattern in cis (single). Decide with junction tests on both breakpoints in
   both orientations. If both orientations are supported, call it a heterozygous SV and flag,
   do not edit.
7. **Ordering and orientation.** Within a group, lay out one haplotype path at a time using
   homology to a reference path plus cis contacts; place the other haplotypes alongside in the
   same orientation. Unresolvable multi-way junctions are reported as such.
8. **Unlocs.** Generalise `scan.py`: for each small or unplaced piece, rank contact with every
   window in the genome, corrected for coverage and repeats; place when one location clearly
   wins and its orientation can be told, tag `Unloc` and attach to the scaffold end when only
   the chromosome is clear.
9. **Haplotigs / false duplications.** Combine: homology to a longer scaffold over most of the
   piece's length, coverage near half (or 1/n) of the local norm, strong MQ0 but weak high-MQ
   signal, and position (ends of scaffolds, small pieces). Tag `Haplotig` or `FalseDuplicate`.
10. **Contaminants.** Pieces with inter-scaffold contact near zero to everything, outlier
    coverage, no telomere/homology partner, and, when supplied, outlier GC or a blob/taxonomy
    hit. Tag `Contaminant`; never delete.
11. **Edit planning.** Convert the target layout into PretextView edits by diffing pixel
    permutations (the simulator is the oracle), minimising the number of edits, appended after
    any `--start-state` edits. Add paint, tags and waypoints. Verify by replaying the written
    savestate through the simulator and comparing with the target layout.
12. **Report.** For each decision: what, where (scaffold and bp coordinates), evidence numbers,
    confidence, a thumbnail before/after. A summary table, a list of things flagged for review,
    and a list of things it could not assess.

## 7. Phases

Each phase ends with tests passing and a short entry in `docs/PROGRESS.md`.

- **Phase 0, verify the foundation.**
  Ask the curator to load `seed/fixtures/*claude_scaffold1` in PretextView for MQ0 and MQ30 and
  confirm it loads, the layout matches the case document, and paint, tags and waypoints show.
  Fix the writer if not. This is the only step that cannot be automated; do it first.
- **Phase 1, I/O and state.** Package the reader, savestate reader/writer, hash, simulator;
  AGP and TSV out; `info`, `dump`, `convert`. Tests: byte-identical round trip of the curator's
  savestate; hash known-answer tests; replay of both fixtures gives the expected layout; the
  nEdits-multiple-of-8 quirk; identical-layout check across maps.
- **Phase 2, test harness.** `simulate`: take a well-curated map, apply random breaks,
  inversions and moves through the simulator, and score how well later phases recover the
  original (edit distance on the permutation, per-task precision and recall). This is the main
  way to tune parameters without hand labels.
- **Phase 3, single-haplotype tasks.** Signal model, junction test, misjoins, inversions,
  unlocs, edit planning, report. Target: the butterfly genome, original and simulated scrambles.
- **Phase 4, allhaps.** Homology map, grouping, haplotype-aware ordering, SV flagging,
  haplotigs. Targets: CAE71 scaffold 1 (must reproduce the case document from the unedited
  map), then the whole CAE71 map compared with the curator's 97 edits, then CAE72, MED, PAR.
- **Phase 5, haplotigs in single mode, contaminants, optional inputs** (FASTA/table, PAF, blob).
- **Phase 6, polish.** Speed and memory, parameter defaults, documentation, a run on a genome
  the tool has never seen with the curator scoring the output.

## 8. Test data

- **CAE71 hap-merged** (hard case): MQ0 and MQ30 maps, curator's savestate, session-1 result.
  Related maps in the same folder: per-haplotype maps, `nhap_4/6/8` runs, CAE72, MED, PAR.
- **A well-behaved butterfly genome** (easy case): to be supplied by the curator, ideally the
  pre-curation map, the curated AGP or savestate, and MQ0 plus high-MQ maps for a dual-haplotype
  assembly. With the curated result available it doubles as ground truth for phases 3 and 4.
- **Simulated scrambles** of any curated map (phase 2).

## 9. Questions for the curator (ask at the start of the next session)

1. Did the two fixture savestates load correctly in PretextView?
2. Which butterfly genome, and where are its maps and curated result?
3. Rust from the start, or Python first and port later?
4. Naming convention for haplotypes in scaffold names (is `H1_`, `H2_` always used)?
5. Which tag conventions should the output follow (the default tag list in the savestate, or
   the lab's own)?
6. Are FASTA files and haplotype-to-haplotype alignments usually at hand when curating, or
   should the tool assume map-only input?
7. How conservative should the default be: apply only high-confidence edits, or apply
   medium-confidence ones too and rely on the report?

## 10. Risks

- PretextView's savestate format is undocumented and has changed between versions; pin the
  commit, keep AGP output as the fallback hand-over.
- Thresholds that work on one genome may not transfer. The simulated-scramble harness and
  map-relative statistics are the mitigation.
- Non-clonal samples and real structural variation can look like misassembly. The tool must be
  able to say "both arrangements are supported" and leave the map alone.
- The session-1 findings on CAE71 are a single unreviewed analysis; if the curator disagrees
  with any of them, update the case document before using it as a regression target.
