"""Reader for PretextMap .pretext files (see docs/FORMATS.md)."""
import os
import struct
import zlib

import numpy as np

_M64 = (1 << 64) - 1
_FH_M = 0x880355F21E6D1965
_HASH_SEED = 0xBAFD06832DE619C2


def inflate(buf, expected=None):
    out = zlib.decompress(buf, -15)
    if expected is not None and len(out) != expected:
        raise ValueError(f"inflate: got {len(out)} bytes, expected {expected}")
    return out


def deflate(buf):
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return c.compress(buf) + c.flush()


def _mix(h):
    h ^= h >> 23
    h = (h * 0x2127599BF4325C37) & _M64
    h ^= h >> 47
    return h


def fasthash64(buf, seed):
    n = len(buf)
    h = seed ^ ((n * _FH_M) & _M64)
    for i in range(n // 8):
        (v,) = struct.unpack_from("<Q", buf, 8 * i)
        h ^= _mix(v)
        h = (h * _FH_M) & _M64
    tail = buf[8 * (n // 8):]
    if tail:
        h ^= _mix(int.from_bytes(tail, "little"))
        h = (h * _FH_M) & _M64
    return _mix(h)


def header_hash(path):
    """8-byte hash PretextView stores in a savestate to tie it to one .pretext file.

    Depends on the file name (first 64 bytes) and the compressed header.
    """
    name = os.path.basename(path).encode()[:64].ljust(64, b"\0")
    with open(path, "rb") as f:
        if f.read(4) != b"pstm":
            raise ValueError(f"{path}: not a .pretext file")
        nc, _ = struct.unpack("<II", f.read(8))
        comp = f.read(nc)
    return struct.pack("<Q", fasthash64(comp, fasthash64(name, _HASH_SEED)))


class Pretext:
    """Header, texture atlas and graph tracks of a .pretext file.

    Scaffolds of the input assembly are called "scaffolds" here; PretextView calls the
    pieces it shows "contigs", and so does `MapState`.
    """

    def __init__(self, path, load_graphs=True):
        self.path = path
        self.f = f = open(path, "rb")
        if f.read(4) != b"pstm":
            raise ValueError(f"{path}: not a .pretext file")
        nc, nh = struct.unpack("<II", f.read(8))
        self.header_comp = f.read(nc)
        hdr = inflate(self.header_comp, nh)
        self.total_len, self.n_scaffolds = struct.unpack("<QI", hdr[:12])
        off = 12
        self.names, fracs = [], []
        for _ in range(self.n_scaffolds):
            fracs.append(struct.unpack_from("<f", hdr, off)[0])
            self.names.append(hdr[off + 4:off + 68].split(b"\0")[0].decode())
            off += 68
        self.fracs = np.array(fracs, dtype=np.float32)
        self.tex_res_log, self.n_tex_log, self.n_mips = hdr[off], hdr[off + 1], hdr[off + 2]
        self.tex_res = 1 << self.tex_res_log
        self.n_tex = 1 << self.n_tex_log
        self.npix = self.tex_res * self.n_tex
        self.bp_per_pix = self.total_len / self.npix
        self.bytes_per_tex = sum((self.tex_res >> l) ** 2 for l in range(self.n_mips)) >> 1
        self.atlas = []
        for _ in range((self.n_tex + 1) * self.n_tex // 2):
            (n,) = struct.unpack("<I", f.read(4))
            self.atlas.append((f.tell(), n))
            f.seek(n, 1)
        self.graphs = {}
        while load_graphs:
            if f.read(4) != b"psgh":
                break
            (n,) = struct.unpack("<I", f.read(4))
            d = inflate(f.read(n), 64 + 4 * self.npix)
            self.graphs[d[:64].split(b"\0")[0].decode()] = np.frombuffer(d[64:], dtype="<i4").copy()
        self._init_pixels()

    def close(self):
        if not self.f.closed:
            self.f.close()

    def __del__(self):
        try:
            self.close()
        except Exception:
            pass

    @property
    def hash(self):
        return header_hash(self.path)

    def scaffold_len_bp(self):
        """Approximate scaffold lengths (the header stores them as f32 fractions)."""
        return np.rint(self.fracs.astype(np.float64) * self.total_len).astype(np.int64)

    def _init_pixels(self):
        """Pixel -> (scaffold, coordinate), replicating PretextView's f32 accumulation."""
        n = self.npix
        ends = np.floor(np.float64(n) * np.cumsum(self.fracs, dtype=np.float32).astype(np.float64))
        ends = np.minimum(np.maximum.accumulate(ends.astype(np.int64)), n)
        starts = np.concatenate(([0], ends[:-1]))
        self.orig_scaffold = np.repeat(np.arange(self.n_scaffolds), ends - starts)
        if len(self.orig_scaffold) < n:      # rounding leftovers go to the last scaffold
            pad = n - len(self.orig_scaffold)
            self.orig_scaffold = np.concatenate(
                (self.orig_scaffold, np.full(pad, self.n_scaffolds - 1)))
        o = self.orig_scaffold
        first = np.concatenate(([True], o[1:] != o[:-1]))
        run_start = np.maximum.accumulate(np.where(first, np.arange(n), 0))
        self.rel_coord = np.arange(n) - run_start
        self.scaffold_start = np.full(self.n_scaffolds, -1, dtype=np.int64)
        self.scaffold_start[o[first]] = np.nonzero(first)[0]
        self.scaffold_npix = np.bincount(o, minlength=self.n_scaffolds)

    # ---- contact matrix ----
    def tex_index(self, x, y):
        """Linear index of texture (x, y), x <= y (row-major upper triangle)."""
        return x * self.n_tex - (x * (x - 1)) // 2 + (y - x)

    def read_tex(self, x, y, level=0):
        """Texture (x, y) at a mip level, as rows x*res.., columns y*res.. (x <= y)."""
        base, nb = self.atlas[self.tex_index(x, y)]
        self.f.seek(base)
        raw = inflate(self.f.read(nb), self.bytes_per_tex)
        off = sum(((self.tex_res >> l) ** 2) >> 1 for l in range(level))
        res = self.tex_res >> level
        return decode_bc4(raw[off:off + ((res * res) >> 1)], res)

    def matrix(self, level):
        """Whole symmetric map at a mip level (side = npix >> level), u8."""
        if not 0 <= level < self.n_mips:
            raise ValueError(f"level must be 0..{self.n_mips - 1}")
        res = self.tex_res >> level
        side = res * self.n_tex
        if side * side > 2 ** 31:
            raise MemoryError(f"level {level} needs {side * side / 2**30:.1f} GB; use a higher level or read_tex")
        m = np.zeros((side, side), dtype=np.uint8)
        for x in range(self.n_tex):
            for y in range(x, self.n_tex):
                t = self.read_tex(x, y, level)
                m[x * res:(x + 1) * res, y * res:(y + 1) * res] = t
                if x != y:
                    m[y * res:(y + 1) * res, x * res:(x + 1) * res] = t.T
        return m


class SegmentedMap:
    """View of a map in which scaffolds are split at cut points (pixels are not moved).

    Behaves like `Pretext` for the detectors: each visible piece is a "scaffold" with its
    own name (`parent:k` for pieces of a cut scaffold), start and length. `parent_of`
    translates back to original scaffold names and coordinates.
    """

    def __init__(self, pt, cut_after):
        self.pt = pt
        n = pt.npix
        first = np.concatenate(([True], pt.orig_scaffold[1:] != pt.orig_scaffold[:-1]))
        cut_after = sorted({int(q) for q in cut_after if 0 <= q < n - 1})
        first[np.array(cut_after, dtype=int) + 1] = True
        starts = np.nonzero(first)[0]
        self.scaffold_start = starts
        self.scaffold_npix = np.diff(np.append(starts, n))
        self.n_scaffolds = len(starts)
        self.orig_scaffold = np.repeat(np.arange(self.n_scaffolds), self.scaffold_npix)
        self.rel_coord = np.arange(n) - starts[self.orig_scaffold]
        parent = pt.orig_scaffold[starts]
        count = np.bincount(parent, minlength=pt.n_scaffolds)
        self.names, self._parent, k = [], [], {}
        for i, par in enumerate(parent.tolist()):
            k[par] = k.get(par, 0) + 1
            self.names.append(pt.names[par] if count[par] == 1 else f"{pt.names[par]}:{k[par]}")
            self._parent.append((pt.names[par], int(pt.rel_coord[starts[i]])))
        self.was_cut = count[parent] > 1

    def parent_of(self, name):
        """(original scaffold name, pixel offset of this piece within it)."""
        return self._parent[self.names.index(name)]

    def scaffold_len_bp(self):
        return np.rint(self.scaffold_npix * self.pt.bp_per_pix).astype(np.int64)

    def __getattr__(self, attr):
        return getattr(self.pt, attr)


def decode_bc4(buf, res):
    """Decode BC4 (RGTC1) blocks into a res x res u8 array."""
    nb = (res // 4) ** 2
    b = np.frombuffer(buf, dtype=np.uint8).reshape(nb, 8)
    a0 = b[:, 0].astype(np.float32)
    a1 = b[:, 1].astype(np.float32)
    gt = a0 > a1
    pal = np.zeros((nb, 8), dtype=np.float32)
    pal[:, 0], pal[:, 1] = a0, a1
    for k in range(2, 8):
        pal[:, k] = np.where(gt, ((8 - k) * a0 + (k - 1) * a1) / 7.0, 0)
    for k in range(2, 6):
        pal[:, k] = np.where(gt, pal[:, k], ((6 - k) * a0 + (k - 1) * a1) / 5.0)
    pal[:, 7] = np.where(gt, pal[:, 7], 255.0)
    bits = np.zeros(nb, dtype=np.uint64)
    for i in range(6):
        bits |= b[:, 2 + i].astype(np.uint64) << np.uint64(8 * i)
    codes = np.empty((nb, 16), dtype=np.int64)
    for i in range(16):
        codes[:, i] = ((bits >> np.uint64(3 * i)) & np.uint64(7)).astype(np.int64)
    px = np.rint(np.take_along_axis(pal, codes, axis=1)).astype(np.uint8)
    n4 = res // 4
    return px.reshape(n4, n4, 4, 4).transpose(0, 2, 1, 3).reshape(res, res)


def check_same_layout(maps):
    """Raise ValueError unless all maps describe the same assembly at the same resolution."""
    ref = maps[0]
    for m in maps[1:]:
        for what, a, b in (("total length", ref.total_len, m.total_len),
                           ("pixel count", ref.npix, m.npix),
                           ("scaffold names", ref.names, m.names)):
            if a != b:
                raise ValueError(f"{m.path}: {what} differs from {ref.path}")
        if not np.array_equal(ref.orig_scaffold, m.orig_scaffold):
            raise ValueError(f"{m.path}: scaffold boundaries differ from {ref.path}")
