# Changelog

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
