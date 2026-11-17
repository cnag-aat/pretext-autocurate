"""allhaps mode: an assembly with both haplotypes in one map (scaffolds named H1_ / H2_).
Homologous scaffolds of the two haplotypes are paired and laid side by side in the same
orientation, tagged HAP1 / HAP2, with cuts, joins and small scaffolds decided on the way.

Homology comes from sequence alignments between the haplotypes when there are any
(minimap2 PAF, see homology.py), else from the MQ0 map: reads from sequence shared by two
haplotypes are split between the copies, so a scaffold shows a line of contact with its
homolog, along the diagonal of their block if both are assembled in the same direction and
along the anti-diagonal if not. Unrelated scaffolds show only background.

Pipeline (`curate_allhaps`):
  pair_homologs      homolog pairs (alignment blocks, or the MQ0 homology line)
  telomere_sites     break sites: gaps with reduced contact, alignment breaks, sharp breaks,
                     telomere strands; each decided cut / invert_tip / keep, and a cut in
                     one haplotype carried to the matching gap of the other
  same_prefix_pairs  homologs that ended up in the same haplotype file
  place_small        small scaffolds: unloc, join when it extends a copy, FalseDuplicate when
                     it is a second copy; judged on non-repetitive aligned sequence
  merge_pairs        pairs that are pieces of one chromosome (contact at both haplotypes' ends)
  arrange_pieces     order and orientation of a copy's pieces by end-to-end contact; then,
                     with alignments, the other copy follows the reference (`_project`)
  layout             parallel pairs, HAP tags, paint, savestate edits
Every act-or-flag threshold goes through `Level` (0 conservative .. 1 aggressive).
Open items are listed in docs/HANDOFF.md.

Two copies of one region in a scaffold (`selfdup`): review rows; a haplotypic copy (by read
depth) with a gap between the copies is cut out.

Modules: common, pairing, sites, evidence, selfdup, small, arrange, driver, state.
A diagnostic that monkeypatches a step (e.g. `place_small`) must patch it where it is used (`allhaps.driver`), not here.
"""
from .common import Level, Pair, THRESHOLDS, Thresholds, _mean0, _signal, _level, _union, write_rows  # noqa: F401
from .pairing import _quality, pair_homologs, same_prefix_pairs, merge_pairs  # noqa: F401
from .sites import telomere_sites, alignment_breaks  # noqa: F401
from .evidence import _hold_mate, _unitig_evidence, _graph_evidence  # noqa: F401
from .small import _small_alignments, place_small  # noqa: F401
from .arrange import _adjacency, _homolog_dir, _block_index, _aligned_dir, _aligned_bp, _project, _contig_conflicts, arrange_pieces  # noqa: F401
from .driver import curate_allhaps  # noqa: F401
