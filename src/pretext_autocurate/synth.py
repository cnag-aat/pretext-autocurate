"""Small synthetic .pretext maps, for tests that run without the (private) fixtures.

`write_pretext` writes a real .pretext file (header, BC4 textures with every mip level,
graph tracks; docs/FORMATS.md) from a full-resolution contact matrix. `toy_map` builds one
from a toy genome: chromosomes of each haplotype, assembled into scaffolds that may join
pieces of different chromosomes (misjoins) or leave a chromosome in several scaffolds.
Contact decays with distance within a chromosome copy; at MQ0 the two copies of a chromosome
also share contact along the homology line (reads of identical sequence map to either copy).
`toy_blocks` gives the alignment blocks between the haplotypes of the same toy genome.
"""
import struct

import numpy as np

from . import homology
from .pretextio import deflate


def _pixel_owner(fracs, npix):
    """Pixel -> scaffold index, as PretextView (and `Pretext._init_pixels`) compute it."""
    ends = np.floor(np.float64(npix) * np.cumsum(np.asarray(fracs, dtype=np.float32), dtype=np.float32)
                    .astype(np.float64))
    ends = np.minimum(np.maximum.accumulate(ends.astype(np.int64)), npix)
    starts = np.concatenate(([0], ends[:-1]))
    own = np.repeat(np.arange(len(fracs)), ends - starts)
    if len(own) < npix:
        own = np.concatenate((own, np.full(npix - len(own), len(fracs) - 1)))
    return own


def encode_bc4(img):
    """BC4 (RGTC1) blocks for a res x res u8 image (res a multiple of 4): bytes, block rows
    first, the inverse of `pretextio.decode_bc4` (endpoints max / min, nearest palette entry)."""
    res = img.shape[0]
    n4 = res // 4
    blk = img.reshape(n4, 4, n4, 4).transpose(0, 2, 1, 3).reshape(-1, 16).astype(np.float32)
    a0, a1 = blk.max(1), blk.min(1)
    pal = np.zeros((len(blk), 8), dtype=np.float32)
    pal[:, 0], pal[:, 1] = a0, a1
    for k in range(2, 8):
        pal[:, k] = ((8 - k) * a0 + (k - 1) * a1) / 7.0
    eq = a0 == a1                       # a single value: decoder's other mode, index 0 is a0
    codes = np.abs(blk[:, :, None] - pal[:, None, :]).argmin(2)
    codes[eq] = 0
    out = np.zeros((len(blk), 8), dtype=np.uint8)
    out[:, 0], out[:, 1] = a0.astype(np.uint8), a1.astype(np.uint8)
    bits = np.zeros(len(blk), dtype=np.uint64)
    for i in range(16):
        bits |= codes[:, i].astype(np.uint64) << np.uint64(3 * i)
    for i in range(6):
        out[:, 2 + i] = ((bits >> np.uint64(8 * i)) & np.uint64(255)).astype(np.uint8)
    return out.tobytes()


def write_pretext(path, names, lengths_bp, matrix, tex_res_log=8, n_mips=4, graphs=None):
    """Write `matrix` (npix x npix u8, symmetric) as a .pretext file for scaffolds `names` of
    `lengths_bp`. npix must be a power of two times 2**tex_res_log. `graphs`: {name: int32
    array of npix} written as extension tracks."""
    npix = matrix.shape[0]
    n_tex = npix >> tex_res_log
    n_tex_log = int(np.log2(n_tex))
    assert (1 << n_tex_log) == n_tex and (n_tex << tex_res_log) == npix
    total = int(sum(lengths_bp))
    hdr = struct.pack("<QI", total, len(names))
    for n, L in zip(names, lengths_bp):
        hdr += struct.pack("<f", L / total) + n.encode()[:63].ljust(64, b"\0")
    hdr += bytes([tex_res_log, n_tex_log, n_mips])
    comp = deflate(hdr)
    out = [b"pstm", struct.pack("<II", len(comp), len(hdr)), comp]
    res = 1 << tex_res_log
    m = matrix.astype(np.float32)
    for x in range(n_tex):
        for y in range(x, n_tex):
            tile = m[x * res:(x + 1) * res, y * res:(y + 1) * res]
            data = b""
            for lv in range(n_mips):
                r = res >> lv
                t = tile.reshape(r, 1 << lv, r, 1 << lv).mean((1, 3)) if lv else tile
                data += encode_bc4(np.clip(np.rint(t), 0, 255).astype(np.uint8))
            c = deflate(data)
            out += [struct.pack("<I", len(c)), c]
    for g, v in (graphs or {}).items():
        d = deflate(g.encode()[:63].ljust(64, b"\0") + np.asarray(v, dtype="<i4").tobytes())
        out += [b"psgh", struct.pack("<I", len(d)), d]
    with open(path, "wb") as f:
        f.write(b"".join(out))


def toy_map(path, scaffolds, chrom_px, npix=512, bp_per_px=10000, homolog=0.35, peak=220.0,
            background=1.0, seed=0, coverage=None):
    """A synthetic two-haplotype (or single) map.

    `chrom_px`: {chromosome: length in px}. `scaffolds`: [(name, [(chromosome, hap, start px,
    end px, reversed), ...])], pieces separated by a 1-px gap (in the `gap` track). Contact
    within a copy decays as peak / (1 + d); between the two haplotypes of a chromosome,
    `homolog` times that at the homologous positions (MQ0). `coverage`: {(chromosome, hap):
    factor} for the `coverage` track (1 = normal). Returns {pixel: (chromosome, hap, pos)}."""
    rng = np.random.default_rng(seed)
    names, lengths, truth = [], [], []
    for name, pieces in scaffolds:
        cells = []
        for k, (ch, hap, a, b, rev) in enumerate(pieces):
            if k:
                cells.append(None)                 # a gap pixel between pieces
            pos = list(range(a, b))
            cells += [(ch, hap, p) for p in (pos[::-1] if rev else pos)]
        names.append(name)
        lengths.append(len(cells) * bp_per_px)
        truth.append(cells)
    flat = [c for t in truth for c in t]
    own = _pixel_owner(np.array(lengths) / sum(lengths), npix)
    # pixel -> its cell (scaffold pixels map onto the scaffold's cells by relative position)
    px_cell, gap = [], np.zeros(npix, dtype=np.int32)
    first = np.concatenate(([True], own[1:] != own[:-1]))
    start = np.maximum.accumulate(np.where(first, np.arange(npix), 0))
    count = np.bincount(own, minlength=len(names))
    for p in range(npix):
        s = own[p]
        cells = truth[s]
        c = cells[min(len(cells) - 1, int((p - start[p]) * len(cells) / max(count[s], 1)))]
        px_cell.append(c)
        if c is None:
            gap[p] = 1
    ch = np.array([c[0] if c else "" for c in px_cell])
    hp = np.array([c[1] if c else -1 for c in px_cell])
    ps = np.array([c[2] if c else -10 ** 6 for c in px_cell], dtype=np.float64)
    d = np.abs(ps[:, None] - ps[None, :])
    same_ch = (ch[:, None] == ch[None, :]) & (ch[:, None] != "")
    same_hp = hp[:, None] == hp[None, :]
    m = np.where(same_ch & same_hp, peak / (1.0 + d), 0.0)
    m += np.where(same_ch & ~same_hp, homolog * peak / (1.0 + d), 0.0)
    m += background * rng.random((npix, npix))
    m = np.triu(m) + np.triu(m, 1).T                 # symmetric
    m[gap > 0, :] *= 0.2                             # gaps carry few reads
    m[:, gap > 0] *= 0.2
    cov = np.full(npix, 100, dtype=np.int32)
    for (c_, h_), f in (coverage or {}).items():
        cov[(ch == c_) & (hp == h_)] = int(100 * f)
    write_pretext(path, names, lengths, np.clip(m, 0, 255), graphs={"gap": gap, "coverage": cov})
    return {p: px_cell[p] for p in range(npix)}


def toy_blocks(scaffolds, bp_per_px=10000, query_hap=2):
    """Blocks (homology.Block) of every `query_hap` scaffold piece against the pieces of the
    same chromosome in the other haplotypes, in each scaffold's own bp (a gap pixel between
    pieces counts, as in `toy_map`): what minimap2 would give for the toy genome."""
    pieces = []                       # (scaffold, chromosome, hap, a, b, rev, first cell)
    for name, ps in scaffolds:
        k = 0
        for i, (ch, hap, a, b, rev) in enumerate(ps):
            if i:
                k += 1                # the gap pixel
            pieces.append((name, ch, hap, a, b, rev, k))
            k += b - a
    cell = lambda p, pos: p[6] + ((p[4] - 1 - pos) if p[5] else (pos - p[3]))
    out = []
    for q in pieces:
        if q[2] != query_hap:
            continue
        for t in pieces:
            if t[2] == query_hap or t[1] != q[1]:
                continue
            a, b = max(q[3], t[3]), min(q[4], t[4])
            if b <= a:
                continue
            qc = sorted((cell(q, a), cell(q, b - 1)))
            tc = sorted((cell(t, a), cell(t, b - 1)))
            n = (b - a) * bp_per_px
            out.append(homology.Block(q[0], qc[0] * bp_per_px, (qc[1] + 1) * bp_per_px, t[0], tc[0] * bp_per_px,
                             (tc[1] + 1) * bp_per_px, q[5] == t[5], n, 60))
    return out
