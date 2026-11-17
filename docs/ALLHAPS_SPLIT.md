# Plan: split `allhaps.py` into steps that can be tested on their own (HANDOFF 15)

Prepared 2026-11-17. Goal: the same results, byte for byte, from code in small modules,
each step a function with its own tests on synthetic maps.

## Status (2026-11-17, end of session)

- Step 1 done: `tools/golden.py` written; `record` ran (47 files from 5 runs, 8.5 min, in
  `golden_out/golden.json`, git-ignored). The determinism check (`check` on unchanged code)
  was started and stopped before it finished: run it first next time; it must say
  "identical" before any move.
- Step 2 prepared, not run: `tools/split_allhaps_moves.py` generates the package from
  `allhaps.py` (dry run prints the plan; `--write` writes it and deletes allhaps.py). No
  import cycles. Before `--write`, fix one thing in it: the moved functions import siblings
  inside their bodies (`from .misjoin`, `from .tracks`, `from . import sexchrom`,
  `from .busco` at allhaps.py lines 352, 469, 1559, 2154, 2169); in the package these must
  become `from ..` (rewrite `from .` to `from ..` inside each moved segment). Then
  `golden.py check`, the full test suite, and delete the script.
- `allhaps.py` itself is unchanged.

## Where it stands

`src/pretext_autocurate/allhaps.py`: 2,457 lines.

| Part | Lines | What |
|---|---|---|
| helpers | ~40 | `Level`, `Pair`, `_mean0`, `_signal`, `_level`, `_quality`, `write_rows` |
| `pair_homologs`, `same_prefix_pairs`, `merge_pairs` | ~315 | homolog pairs; pairs that are parts of one chromosome |
| `telomere_sites`, `alignment_breaks` | ~345 | cut / invert_tip / keep sites at gaps |
| `_graph_evidence`, `_unitig_evidence`, `_hold_mate` | ~100 | hifiasm graph evidence at sites |
| `_small_alignments`, `place_small` | ~325 | small scaffolds: unloc, candidate, false duplicate |
| `_adjacency`, `_homolog_dir`, `_aligned_*`, `_project`, `_contig_conflicts`, `arrange_pieces` | ~320 | order and orientation inside a chromosome copy |
| `curate_allhaps` | **907** | the driver: eight phases in one function |

The eight phases of `curate_allhaps` (lines as of today):

| # | Phase | Lines | Notes |
|---|---|---|---|
| 0 | setup: sex chromosomes, first pairing, sites, graph evidence, cut map, pairing | 35 | |
| 1 | small scaffolds (`place_small`), state for arranging | 53 | defines most shared state |
| 2 | merge pairs that are parts of one chromosome | 37 | |
| 3 | hemizygous pieces chained to a sex chromosome | 43 | |
| 4 | phase 1: each copy arranged by end contact; large attached pieces placed | 125 | |
| 5 | reference copy per pair; loose scaffolds into chains | 96 | |
| 6 | phase 2: the other copy follows the reference (alignment) | 187 | |
| 7 | layout order, decision rows, paint, tags, savestate, sex summary | 331 | largest; two jobs |

The phases share about 30 variables through the function's locals (the cut map,
`pairs`, `info`, `hap`, `attach`, `mate_of`, `extra_of`, `host_of`, `members`, `arranged`,
`joined`, `arr_rows`, `sol`, `soft`, `tips`, `small_calls`, `sx`, `sex_forced`, ... and the
thresholds). A dependency scan also shows single-letter loop variables (`x`, `v`, `c`, `copy`,
`s`, `why`) read in a later phase than the one that set them: either deliberate re-use or
latent bugs; the split forces each one to be made explicit, and the byte-identical check
below shows which.

There are 30 level-dependent thresholds (`lv.at(...)`), 16 in `curate_allhaps`, the rest
passed into the helpers; they are the real meaning of `--level` and are scattered.

## Target

A package with the same public names, so callers and tests keep working:

    src/pretext_autocurate/allhaps/
        __init__.py      re-exports: curate_allhaps, write_rows, Level, Pair, pair_homologs,
                         telomere_sites, place_small, merge_pairs, arrange_pieces, and the
                         private helpers other modules use (_signal, _adjacency, _block_index,
                         _aligned_dir, _contig_conflicts)
        common.py        Level, Pair, Thresholds, _signal, _level, _mean0, _union, write_rows
        pairing.py       _quality, pair_homologs, same_prefix_pairs, merge_pairs
        sites.py         telomere_sites, alignment_breaks
        evidence.py      _graph_evidence, _unitig_evidence, _hold_mate
        small.py         _small_alignments, place_small
        arrange.py       _adjacency, _homolog_dir, _block_index, _aligned_*, _project,
                         _contig_conflicts, arrange_pieces
        state.py         Run: the shared state, one dataclass with every field documented
        steps.py         one function per phase: step_sites(run), step_small(run),
                         step_merge(run), step_sex_chain(run), step_arrange(run),
                         step_reference(run), step_follow(run), each adding its rows to run
        layout.py        step_layout(run) and step_rows(run): phase 7 split in two
        driver.py        curate_allhaps: builds Run, calls the steps in order, returns
                         (rows, savestate, info) as today

- `Thresholds(level)`: one dataclass holding all 30 level-dependent values with a comment
  each, built from `Level`. `--level` then has one place to read and document, and a
  sweep prints its table.
- `Run`: fields grouped by the step that fills them; steps read only what earlier steps
  filled (checked by a test that runs each step on a fresh Run).

## Safety net: byte-identical output

Before the first move, record golden outputs, and compare after every move:

- the three fixture genomes (`tools/score_cases.json`) at levels 0, 0.5, 1:
  `autocurate.decisions.tsv` and every savestate, compared byte for byte;
- the synthetic maps (`tests/test_synthetic.py`) likewise, so the check also runs on GitHub.

A small script (`tools/golden.py record|check`) does it; every step of the work below must
end with `check` passing. Any difference is either a latent bug found (decide and record it)
or a mistake in the move.

## Order of work

1. `tools/golden.py`; record.
2. Pure moves: helpers into `common`, `pairing`, `sites`, `evidence`, `small`, `arrange`;
   `allhaps/__init__.py` re-exports. No edits inside functions. Check.
3. `Thresholds`: replace every `lv.at(...)` by a named field. Check.
4. `Run` and the steps, one phase at a time, from the top (setup, small, merge, ...). Each
   loop variable that crossed phases is made explicit. Check after each phase.
5. Split phase 7 into layout and rows. Check.
6. Tests per step on synthetic toy genomes (`synth.toy_map`), no fixtures:
   - sites: a misjoin at a gap is cut; a gap between two pieces of one chromosome is kept;
   - pairing: copies of the same chromosome pair; a chromosome in two scaffolds in both
     haplotypes is merged;
   - small: a small piece with signal unique to one pair becomes an Unloc, one with signal
     everywhere does not;
   - arrange: a piece inverted in the assembly is turned back by end contact; a large
     attached piece goes where its ends touch best;
   - follow: with synthetic alignment blocks, the other copy takes the reference's order;
   - layout: the savestate is a permutation; each pair's copies sit side by side; every
     piece is painted and tagged.
7. Docs: HANDOFF module table, README development notes.

## Costs and risks

- About a day of work, most of it steps 4 and 6; the golden check makes it safe to stop
  between any two steps.
- Scratch scripts that monkeypatch `allhaps.place_small` or `allhaps.arrange_pieces` (used
  for diagnostics in this project) must patch the step module instead; documented in
  `allhaps/__init__.py`.
- No behaviour change is intended; if the golden check exposes a latent bug, it is fixed in
  a separate commit after the split, with its own test.
