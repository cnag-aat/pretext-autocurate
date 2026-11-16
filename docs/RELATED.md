# Related work: hapscotch (Chenxi Zhou) and what to take from it

Read 2026-11-14 from https://github.com/c-zhou/hapscotch (version 0.2; Zenodo
10.5281/zenodo.23573025). Chenxi Zhou also wrote YaHS and, with Gene Myers and Richard
Durbin, FastGA. The curator regards the work of Durbin's group as a reference to learn from.

## What hapscotch does

A pipeline that phases and scaffolds the haplotypes of a polyploid or multi-haplotype
assembly, starting from **contigs**:

1. `seqtools idx`: index of all contigs from all haplotypes together.
2. **All-against-all self-alignment** with FastGA (`FastGA -1:genome.1aln genome.fa`;
   the .1aln format is preferred over PAF).
3. Optional Hi-C alignment (bwa-mem2 or minibwa), converted to a binary contact file.
4. `hapcure`: contig error correction from Hi-C, including palindromes (with the
   self-alignment).
5. `hapscotch`: bins contigs into haplotype groups from synteny in the self-alignment
   ("scaffold groups constructed from synteny"), estimates the haplotype count (`hapcount`),
   builds longer sequences from overlaps, and with Hi-C phases the haplotypes by
   minimising Hi-C contact between them: a mixed-integer program (HiGHS) solved on a
   rolling window of about 300 sequences, with Leiden clustering of Hi-C linkage groups.
6. YaHS scaffolds one **backbone** per group: the Hi-C contacts of all haplotypes of the
   group are moved onto a single backbone coordinate system and scaffolded together, then
   split back into haplotypes (`seqtools hap -p N`).

Outputs: per-haplotype AGP and FASTA, haplotype counts, group assignments, alignment and
Hi-C plots.

## Comparison with pretext-autocurate

| | hapscotch | pretext-autocurate |
|---|---|---|
| Starts from | contigs (before scaffolding) | scaffolds and a PretextMap map (after YaHS), for curation |
| Homology | FastGA all-against-all | minimap2 H2 against H1, plus small scaffolds against all |
| Haplotype count | estimated (polyploids) | two prefixes assumed; more is open (CAE71) |
| Phasing | MILP minimising inter-haplotype Hi-C | keeps hifiasm's phasing; whole swaps, small pieces move |
| Scaffolding | YaHS on a shared backbone per group | order and orientation along the other haplotype, end contact as a check |
| Contig errors | hapcure (Hi-C, palindromes) | cuts at gaps; inside contigs only at sharp breaks |
| Output | new assembly | savestate for PretextView, decisions, review list |

They are complementary: hapscotch builds the haplotype scaffolds, pretext-autocurate
pre-curates whatever scaffolds a curator is given.

## Ideas to adopt

1. **All-against-all alignment with FastGA** in place of two minimap2 runs. One run gives
   the H1-H2 homology, same-haplotype homologs and duplicates together, and hapscotch
   shows how to use it. See the benchmark in PROGRESS.md.
2. **A shared backbone per chromosome**: the tool already orders one copy along the other
   (`_project`). Pooling the Hi-C contacts of both copies on one coordinate system would
   give twice the signal for joins and cuts, which matters for sparse or MQ0-only maps.
   Measured 2026-11-16 on curated gaps: no gain, the copies give nearly the same contact
   (PROGRESS.md).
3. **Phasing as an optimisation** (minimise inter-haplotype Hi-C within a group): a
   principled form of `--phasing none`, and a way to decide hard cases such as
   H1_scaffold_31 or H1_scaffold_44 on the owl.
4. **Haplotype count** (`hapcount`): for assemblies with more than two copies (CAE71).
5. **Palindrome correction** (`hapcure -a`): a contig error type the tool does not look for.
