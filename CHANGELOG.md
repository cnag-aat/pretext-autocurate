# Changelog

## Unreleased

- Telomeres by strand: `--telomere-motif` runs telo_scan.py (CNAG) on `--fasta`
  (`--telo-scan PATH`), or give its bedgraphs with `--telomere-5p` / `--telomere-3p`, or
  carry them in the map as `5p_telomere` / `3p_telomere` tracks. At a gap: telomeres facing
  the gap from both sides, or from the long side, mean two chromosome ends (cut); from the
  end piece only, an inverted end piece.
- Repeats: `--repeats BED`, `--red` (runs Red on each haplotype of `--fasta`), or a
  `repeat_density` track in the map. Alignment breaks inside repeats are ignored; a small
  scaffold that is 90% repeat is not moved. The tool prints the PretextGraph command to add
  computed tracks to a copy of the map.
- `--add-tracks`: writes a copy of each map with the evidence tracks for PretextView
  (repeat_density, 5p/3p_telomere, busco_dup_pairs, busco_missing_hap) and a savestate for
  each copy. Useful on the high-MQ map too, where a repeat track tells repeat from
  contaminant. BUSCO duplicate pairs share a colour per pair of scaffolds (both places of a
  possible haplotig); genes missing from one haplotype are shown where the other has them.
  Needs PretextGraph.
- The repeat track blocks only placements that rest on Hi-C alone; with sequence evidence
  the decision rests on unique alignments and repeat content does not count.
- `tools/bench_aligners.py`: minimap2 against FASTGA for the haplotype alignment.
- `docs/RELATED.md`: hapscotch (Chenxi Zhou) read and compared.

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
