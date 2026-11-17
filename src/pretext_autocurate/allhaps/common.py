"""Shared pieces of the allhaps pass: `Level` (the conservative-to-aggressive setting), `Pair`,
the per-map `Signal` cache and small helpers."""
from dataclasses import dataclass
from ..signal import Signal


@dataclass
class Level:
    """One setting for how readily the tool acts: 0 = conservative (act only on the clearest
    evidence, flag the rest), 0.5 = default, 1 = aggressive (act on weaker evidence too).
    Every threshold that decides between acting and flagging is read through `at`."""
    a: float = 0.5

    def at(self, conservative, default, aggressive):
        """The threshold at this level: piecewise linear through the three values given for
        levels 0, 0.5 and 1 (the level is clamped to 0..1)."""
        a = min(max(self.a, 0.0), 1.0)
        return conservative + (default - conservative) * a / 0.5 if a <= 0.5 else \
            default + (aggressive - default) * (a - 0.5) / 0.5


# Every level-dependent threshold of the allhaps pass: name -> (value at level 0, 0.5, 1,
# meaning). The level moves them all together (`--level`); this table is where they live.
THRESHOLDS = {
    "pair_margin": (1.5, 1.3, 1.15, "how far the best homolog must beat the next (pair_homologs)"),
    "gap_max": (0.45, 0.6, 0.75, "contact across a gap below this (of expectation) makes it a site"),
    "sharp_max": (0.05, 0.1, 0.2, "contact across a point with no gap below this: a sharp break"),
    "tip_max": (0.35, 0.5, 0.65, "end piece contact ratio deciding a tip inversion"),
    "seq_min": (0.6, 0.5, 0.4, "share of a small scaffold aligned to a chromosome to place it"),
    "hic_agree": (2.0, 1.5, 1.2, "MQ0 contact ratio best / next chromosome confirming sequence"),
    "hic_alone": (3.5, 2.5, 2.0, "MQ0 contact ratio to place a small scaffold on Hi-C alone"),
    "hi_unique": (10.0, 5.0, 3.0, "high-MQ contact ratio best / next pair: a small piece's signal is unique"),
    "mq0_unique": (10.0, 6.0, 4.0, "the same on the MQ0 map (without a high-MQ map)"),
    "chain_px": (60, 40, 20, "smallest attached piece (px) put into a copy's arrangement search"),
    "gain_min": (0.3, 0.15, 0.05, "least gain in end contact to rearrange a copy"),
    "attach_link": (0.5, 0.3, 0.15, "least end contact for an attached piece to be placed in its copy"),
    "attach_margin": (0.3, 0.15, 0.05, "how much better that place must be than any other"),
    "orient_gain": (0.2, 0.12, 0.05, "least gain for a piece to be turned round"),
    "merge_link_min": (0.7, 0.5, 0.3, "end contact for two pairs to be parts of one chromosome"),
    "merge_specific": (4.0, 3.0, 2.0, "that contact over the pairs' third-best partner"),
    "merge_link_weak": (0.25, 0.15, 0.1, "weaker end contact accepted with chromosome-level contact (merge_block_min)"),
    "merge_block_min": (8.0, 5.0, 3.0, "whole-block contact of two copies over their trans background, both haplotypes"),
    "sex_chain_min": (0.08, 0.05, 0.03, "chromosome-level contact chaining a hemizygous piece"),
    "sex_chain_ratio": (4.0, 3.0, 2.0, "that contact over the next chromosome copy"),
    "conflict_min": (0.3, 0.5, 0.8, "end-contact gain that keeps an arrangement against the alignment"),
    "contig_min": (1e6, 5e5, 2e5, "bp of contradicting contigs to make the other copy the reference"),
    "min_loose": (8, 4, 2, "smallest loose scaffold (px) put into a chain where it aligns"),
}


class Thresholds:
    """The thresholds of `THRESHOLDS` at one level, as attributes (`Thresholds(0.5).gap_max`)."""

    def __init__(self, level=0.5):
        lv = Level(level)
        self.level = level
        for name, (c, d, a, _) in THRESHOLDS.items():
            setattr(self, name, lv.at(c, d, a))

    def table(self):
        """[(name, value, meaning)] in table order, for documentation and run reports."""
        return [(n, getattr(self, n), v[3]) for n, v in THRESHOLDS.items()]


@dataclass
class Pair:
    """Two homologous units (scaffolds or cut pieces) of different haplotypes, or a further
    piece of one homolog (kind 'extra'). `merge_pairs` later sets kind 'merged' on a pair
    that became part of another pair's chromosome."""
    a: int                # scaffold index, placed first (HAP1)
    b: int                # its homolog (HAP2)
    score: float          # mean of the homology line (map units at the working level)
    second: float         # best line of either scaffold with any other scaffold
    same_strand: bool
    kind: str             # pair | unequal (possible sex chromosomes) | extra
    note: str


def _mean0(v):
    """Mean of v, 0 when empty (small maps: no pair of pixels at some distance)."""
    return float(v.mean()) if len(v) else 0.0


def _signal(m):
    """One Signal per map, kept on the map object: its decay curve, insulation profiles,
    binned matrix and decoded tiles are reused between steps and between sweep levels."""
    m = getattr(m, "pt", m)
    if "_signal" not in m.__dict__:
        m.__dict__["_signal"] = Signal(m)
    return m.__dict__["_signal"]


def _level(pt, max_side=8192):
    """The finest mip level whose matrix side is at most `max_side` bins (pairing works on a
    whole-map matrix, so 65,536 px is too large)."""
    lv = 0
    while (pt.npix >> lv) > max_side and lv < pt.n_mips - 1:
        lv += 1
    return lv


def _union(iv):
    """Total length covered by half-open intervals (a, b), overlaps counted once (same units
    as the intervals: bp or pixels)."""
    tot, end = 0, -1
    for a, b in sorted(iv):
        if b > end:
            tot += b - max(a, end)
            end = b
    return tot


def write_rows(rows, path):
    """Write the decision rows of `curate_allhaps` to `path` as a TSV with a header line."""
    cols = ["action", "scaffold", "length_mb", "partner", "strand", "tag", "score", "next_best", "evidence"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) for c in cols) + "\n")
