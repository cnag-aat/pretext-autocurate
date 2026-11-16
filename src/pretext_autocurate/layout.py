"""Turn a map state plus savestate annotations into scaffolds; write AGP and TSV."""
from typing import NamedTuple

from .mapstate import MapState


class Piece(NamedTuple):
    """One contig of the edited map, placed in an output scaffold (a row of the TSV)."""
    index: int          # contig index along the map (what savestate paint/tags refer to)
    out_scaffold: int   # 1-based output scaffold number
    painted_id: int     # PretextView scaffold id, 0 if unpainted
    name: str           # original scaffold name
    start_bp: int       # 1-based, inclusive, in the original scaffold
    end_bp: int         # inclusive
    strand: str
    start_px: int       # position in the edited map
    length_px: int
    tags: tuple         # PretextView tag names set on this contig


def build_layout(pt, state=None, ms=None):
    """Replay a savestate on a map and return the list of pieces in map order.

    Painted contigs that are adjacent and share a scaffold id form one output scaffold;
    every unpainted contig is its own. bp coordinates are pixel-derived (rounded to the
    pixel grid, clipped to the scaffold length), so cut points are approximate: downstream
    tools are expected to snap them to real gaps.
    """
    if ms is None:
        ms = MapState(pt)
        if state is not None:
            ms.replay(state.edits)
    paint = dict(state.scaffs) if state else {}
    flags = dict(state.meta_flags) if state else {}
    bpp = pt.bp_per_pix
    lens = pt.scaffold_len_bp()
    pieces, out_id, prev_paint = [], 0, None
    for i, c in enumerate(ms.contigs):
        pid = paint.get(i, 0)
        if not pid or pid != prev_paint:
            out_id += 1
        prev_paint = pid or None        # None, so two unpainted neighbours never merge
        last_px = pt.scaffold_npix[c.scaffold] - 1
        ln = int(lens[c.scaffold])
        # pixel p covers bp [p*bpp, (p+1)*bpp); a piece reaching the scaffold's last pixel
        # ends at its full length (pixel rounding would otherwise drop or overshoot the tail)
        s = min(int(round(c.lo * bpp)), ln - 1) + 1
        e = ln if c.hi >= last_px else max(min(int(round((c.hi + 1) * bpp)), ln), s)
        tags = tuple(state.tag_names(flags.get(i, 0))) if state else ()
        pieces.append(Piece(i, out_id, pid, pt.names[c.scaffold], s, e,
                            "-" if c.inverted else "+", c.start, c.length, tags))
    return pieces


def invisible_scaffolds(pt):
    """Scaffolds shorter than a pixel: absent from the map, so absent from any layout."""
    lens = pt.scaffold_len_bp()
    return [(pt.names[i], int(lens[i])) for i in range(pt.n_scaffolds) if pt.scaffold_npix[i] == 0]


# tags PretextView uses to paint the map display; its AGP export leaves them out
PAINT_TAGS = {"vertPaint", "horzPaint", "crossPaint"}


def unplaced_hap_tags(pt, state, hap_regex=r"^H(\d+)_"):
    """Scaffolds too short for a pixel (absent from the map, so from PretextView's AGP),
    each with the haplotype tag of its name prefix (H1_ -> HAP1, as the tool tags unplaced
    scaffolds), when the savestate has such tags: [(name, length bp, [tag])]. Without them
    pretext-to-asm sorts these sequences by name prefix into extra assemblies (h1, h2)."""
    import re
    have = set(state.tags.values())
    out = []
    for name, ln in invisible_scaffolds(pt):
        m = re.match(hap_regex, name)
        tag = f"HAP{m.group(1)}" if m else None
        out.append((name, ln, [tag] if tag in have else []))
    return out


def write_agp(pt, state, path, ms=None, unplaced=()):
    """The AGP PretextView would export after loading `state` on map `pt` (GenerateAGP in
    PretextView.cpp at the commit in docs/pretextview_commit.txt, default options: no
    reverse-complemented singletons, painted scaffolds first). Downstream tools, notably
    pretext-to-asm (sanger-tol/agp-tpf-utils), are built on that export, so it is mirrored
    exactly:
      - header with "# HiC MAP RESOLUTION: <bp> bp/texel" (pretext-to-asm needs it to snap
        breaks to the assembly's gaps);
      - bp coordinates from pixels: start = int(lo / N * G) + 1, end = int((lo + n) / N * G)
        (not clipped to the scaffold: pretext-to-asm snaps them);
      - painted contigs first, in map order; contigs painted alike and next to each other
        (unpainted ones skipped) form Scaffold_<paint id>, joined by 100 bp U gaps; then
        each unpainted contig as Scaffold_<last painted id + 1, + 2, ...>;
      - "-" only for contigs in a painted scaffold of more than one piece;
      - "Painted" and the tag names (bit order) as columns 10+, without the paint tags.
    `unplaced`: [(name, length, tags)] appended after everything else, one scaffold each
    (the one departure from PretextView: sequences too short for a pixel, see
    `unplaced_hap_tags`)."""
    from . import __version__
    if ms is None:
        ms = MapState(pt).replay(state.edits)
    paint = dict(state.scaffs)
    flags = dict(state.meta_flags)
    N, G = float(pt.npix), float(pt.total_len)
    cs = ms.contigs
    n = len(cs)
    sid = [paint.get(i, 0) for i in range(n)]

    def strand(i):
        c = cs[i]
        lone = (not sid[i]
                or (0 < i < n - 1 and sid[i] != sid[i + 1] and sid[i] != sid[i - 1])
                or (i == 0 and n > 1 and sid[i] != sid[i + 1])
                or (i == n - 1 and n > 1 and sid[i] != sid[i - 1]))
        return "-" if c.inverted and not lone else "+"

    def coords(c):
        s = int(c.lo / N * G) + 1
        e = int((c.lo + c.length) / N * G)
        return s, e

    def tags(i):
        return [t for t in state.tag_names(flags.get(i, 0)) if t not in PAINT_TAGS]

    out = ["##agp-version\t2.1\n",
           f"# DESCRIPTION: Generated by pretext-autocurate {__version__}\n",
           f"# HiC MAP RESOLUTION: {G / N:f} bp/texel\n"]
    cur, pos, part = 0, 1, 0
    for i, c in enumerate(cs):                      # painted scaffolds
        if not sid[i]:
            continue
        if sid[i] != cur:
            cur, pos, part = sid[i], 1, 0
        if part:
            part += 1
            out.append(f"Scaffold_{cur}\t{pos}\t{pos + 99}\t{part}\tU\t100\tscaffold\tyes\tproximity_ligation\n")
            pos += 100
        s, e = coords(c)
        part += 1
        out.append("\t".join(map(str, [f"Scaffold_{cur}", pos, pos + e - s, part, "W", pt.names[c.scaffold], s, e,
                                       strand(i), "Painted"] + tags(i))) + "\n")
        pos += e - s + 1
    for i, c in enumerate(cs):                      # unpainted contigs, one scaffold each
        if sid[i]:
            continue
        cur += 1
        s, e = coords(c)
        out.append("\t".join(map(str, [f"Scaffold_{cur}", 1, e - s + 1, 1, "W", pt.names[c.scaffold], s, e,
                                       strand(i)] + tags(i))) + "\n")
    for name, ln, tg in unplaced:                   # not in the map at all (see docstring)
        cur += 1
        out.append("\t".join(map(str, [f"Scaffold_{cur}", 1, ln, 1, "W", name, 1, ln, "+"] + list(tg))) + "\n")
    with open(path, "w") as f:
        f.writelines(out)


def write_tsv(pieces, path):
    """One row per Piece; tags comma-separated in the last column."""
    cols = ["contig_index", "out_scaffold", "painted_id", "scaffold", "start_bp", "end_bp",
            "strand", "map_start_px", "length_px", "tags"]
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for p in pieces:
            f.write("\t".join(map(str, p[:-1])) + "\t" + ",".join(p.tags) + "\n")
