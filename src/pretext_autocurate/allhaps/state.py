"""`Run`: the state the steps of `curate_allhaps` share (allhaps.driver, `STEPS`).

Each step reads the attributes it needs and sets the ones later steps read. In the order the
steps set them:

  parameters   pt (the map; replaced by the cut map in _setup), raw (the map as given), th
               (Thresholds at the run's level), hi, blocks, small_blocks, busco_hits,
               phasing, hap_regex, sex_system, same_hap_blocks, graph_junctions,
               unitig_paths, do_cuts, do_place, do_arrange, level
  _setup       sx (sex chromosomes, sexchrom.analyse), sites, cuts, flips, graph_rows,
               pairs, info
  _small       small_calls, attach (unlocs per copy), hap, npx_, start_, sg_arr, soft, tips,
               chain_px, gain_min, attach_link, attach_margin, orient_gain_, mate_of,
               extra_of, arranged, joined, arr_rows, sol, block_idx
  _merge       host_of (member unit -> host unit), members
  _sex_chain   sexfrac, sex_forced, sex_rows, span
  _arrange     (changes arranged, joined, arr_rows in place)
  _reference   pair_of, along, loose, false_dup, conflict_min, contig_min, contig_ratio,
               g230, join_rows
  _follow      arr_rows
  _rows        rows (decision rows so far), genes, dup_of (BUSCO), npx, start, mb, bpp
  _layout      returns (rows, savestate, info)

`curate_allhaps(..., stop_after=step)` returns the Run after that step (tests/test_steps.py).
"""


class Run:
    """Attribute container for one allhaps run; see the module docstring for its fields."""

    def __init__(self, **kw):
        self.__dict__.update(kw)

    def __repr__(self):
        return f"Run({', '.join(sorted(self.__dict__))})"
