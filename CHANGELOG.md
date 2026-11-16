# Changelog

## Unreleased

- Contigs outrank scaffold joins: with sequence alignments, the copy of a chromosome that
  leads the other's order and orientation is the one whose arrangement the other copy's
  contigs contradict less (before: always the copy in fewer pieces). When the contigs of
  each copy contradict the other, the pair is flagged `arrangement_differs`. The rows give
  the contigs and the contradicted Mb. No change on the test genomes, where the contigs
  agree with the old choice.
- Sex chromosomes from coverage and alignment (`sexchrom.py`). Hemizygous stretches (half
  coverage, no alignment to the other haplotype) mark the sex-specific part of a sex
  chromosome; a half-coverage stretch aligned into one is the other sex chromosome's
  homologous part. Divergent heterozygous stretches of autosomes, also at half coverage but
  aligned to each other, are no longer mistaken for them. Normal-coverage sequence next to
  a sex-linked stretch is called a PAR when it has chromosome-level contact with it and
  aligns to (or, without alignments, shares MQ0 reads with) the other sex chromosome's
  end; with no contact it is another chromosome joined there, and the cut row says so.
  Two sex-linked scaffolds in one haplotype are aligned to each other with `--fasta`
  (`homology.align_within`). Hemizygous scaffolds in no pair join the chain of the sex
  chromosome they contact. New rows `sex_linked_region`, `par`, `par_candidate`,
  `sex_boundary`, `sex_chromosome_join`; new tracks `sex_chromosome` and `par`
  (`OUT_DIR/tracks/`, added to the map with `--add-tracks`); `--sex-system` (XY, ZW, X0,
  Z0, X1X2Y, XY1Y2, Z1Z2W, ZW1W2) tags the copies (without it they are not labelled: that
  needs species knowledge). A `sex_summary` row reports what was found, and when nothing
  was, that a homogametic sample or homomorphic sex chromosomes look like autosomes.
- Extra pieces of a chromosome copy (scaffolds aligned or chained into it) are split at
  their kept gaps before order and orientation are searched, so a piece can go into
  another's gap. Troglobisium: the X is laid out as the curator did (H2_scaffold_46 split
  at 8.6 Mb, its first part in H1_scaffold_1's gap at 109.8 Mb); 24 of 32 tool junctions
  in the curation (21 of 30 before). No change on the owl.
- Every module commented inline; `docs/HANDOFF.md` lists suggestions for further
  development.

## 1.1.0 (2026-11-14)

Evidence tracks: telomeres by strand, repeats, BUSCO pairs.

- Telomeres by strand from telo_scan.py (CNAG), now bundled in `pretext_autocurate/vendor/`:
  `--telomere-motif MOTIF` with `--fasta` runs it (a few minutes for a 4 Gb assembly);
  `--telomere-5p` / `--telomere-3p` take existing telo_scan bedgraphs; maps that already
  carry its `5p_telomere` / `3p_telomere` tracks are read directly. At a gap: telomeres
  facing it from both sides, or from the long side, mean two chromosome ends joined (cut);
  from the end piece only, an inverted end piece (inverted, with the evidence in the
  decisions table). The check includes the gap's own pixel, since telomeres are shorter
  than a bin.
- Repeats: `--repeats BED`, `--red` (runs Red on each haplotype of `--fasta`), or a
  `repeat_density` track in the map. Alignment breaks inside repeats are ignored; a
  scaffold placed by Hi-C alone is not moved when it is 90% repeat. With sequence evidence
  the decision rests on unique alignments and repeat content does not count.
- `--add-tracks`: a copy of each map with the evidence tracks for PretextView
  (`5p_telomere`, `3p_telomere`, `repeat_density`, `busco_dup_pairs`, `busco_missing_hap`)
  and a savestate for each copy (needs PretextGraph). On a high-MQ map a repeat track tells
  repeat from contaminant. BUSCO duplicate pairs share one colour per pair of scaffolds, so
  both places of a possible haplotig light up together; genes missing from one haplotype
  are shown where the other has them.
- Results on the test genomes: on the boreal owl the strands confirm the scaffold_1 cut and
  the H1_scaffold_2 tip, and invert the end pieces of H1_scaffold_5 (as in the original
  curation) and H2_scaffold_2; on Troglobisium (TTAGG) four flagged sites become cuts, all
  at places the curator broke.
- Tools and notes: `tools/bench_aligners.py` (minimap2 against FastGA; minimap2 stays the
  aligner), `docs/RELATED.md` (hapscotch), `CLAUDE.md` (run commands under bash).

## 1.0.1 (2026-11-14)

- License: PolyForm Noncommercial License 1.0.0 (`LICENSE.md`). Free for noncommercial use,
  including by educational institutions and public research organizations; commercial use
  requires a license from CNAG (`NOTICE`).
- Repository moved to https://github.com/cnag-aat/pretext-autocurate.
- No change to the code or its results.

## 1.0.0 (2026-11-14)

First release for testing by curators ("pre-curation").

Two-haplotype mode (`allhaps`)
- Homolog pairing from sequence alignment (minimap2, `--fasta` / `--paf`) or MQ0 Hi-C;
  pairs laid out side by side and parallel, better copy tagged HAP1 (`--phasing`).
- Cuts at gaps with reduced contact, at alignment breaks and at interstitial telomeres,
  carried to the matching gap of the other haplotype; sharp breaks away from gaps.
  Telomere signal is supporting evidence only.
- Order and orientation of pieces within each chromosome copy, along the other haplotype
  when aligned, by end contact otherwise; merging of pairs that are parts of one
  chromosome; scaffolds joined where they align.
- Small scaffolds: unlocs where sequence and Hi-C agree, joins when a piece continues a
  copy, `FalseDuplicate` for second copies, decisions on non-repetitive sequence only.
- Sex chromosome candidates from half coverage.
- BUSCO tables as an optional tie-breaker and for false-duplication candidates.
- `--level` (conservative to aggressive) and `--sweep` (one savestate per level).

Both modes
- `--mode auto` (default) from the scaffold names; `--hap-regex`.
- New outputs: `autocurate.review.tsv` (what to look at) and `autocurate.run.txt`
  (version, command and options).
- Clear errors for missing maps, maps of different assemblies and missing prefixes.

Known limitations
- Hemizygous chromosomes built from several scaffolds are not chained.
- Two haplotypes of one chromosome inside one scaffold are not separated.
- No repeat track yet (Red / PretextGraph planned); repeats are inferred from alignments.
- The AGP has not been checked against PretextView's own export.
