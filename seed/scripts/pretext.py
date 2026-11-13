"""Minimal reader for PretextMap .pretext files and PretextView savestate files.

Formats taken from the PretextView source (sanger-tol/PretextView, PretextView.cpp).
"""
import struct
import zlib
import numpy as np

BREAK_DELTA = -2147483648


def inflate(buf, expected=None):
    out = zlib.decompress(buf, -15)
    if expected is not None and len(out) != expected:
        raise ValueError(f"inflate: got {len(out)} expected {expected}")
    return out


def deflate(buf):
    c = zlib.compressobj(9, zlib.DEFLATED, -15)
    return c.compress(buf) + c.flush()


class Pretext:
    def __init__(self, path):
        self.path = path
        f = open(path, 'rb')
        self.f = f
        if f.read(4) != b'pstm':
            raise ValueError('not a pretext file')
        nc, nh = struct.unpack('<II', f.read(8))
        self.header_comp = f.read(nc)
        hdr = inflate(self.header_comp, nh)
        self.total_len, self.n_contigs = struct.unpack('<QI', hdr[:12])
        off = 12
        self.names, self.fracs = [], []
        for _ in range(self.n_contigs):
            frac, = struct.unpack('<f', hdr[off:off + 4])
            off += 4
            self.names.append(hdr[off:off + 64].split(b'\0')[0].decode())
            off += 64
            self.fracs.append(frac)
        self.tex_res_log, self.n_tex_log, self.n_mips = hdr[off], hdr[off + 1], hdr[off + 2]
        self.tex_res = 1 << self.tex_res_log
        self.n_tex = 1 << self.n_tex_log
        self.npix = self.tex_res * self.n_tex
        self.bp_per_pix = self.total_len / self.npix
        self.bytes_per_tex = sum((self.tex_res >> l) ** 2 for l in range(self.n_mips)) >> 1
        # atlas
        ntex = (self.n_tex + 1) * (self.n_tex >> 1)
        self.atlas = []
        for _ in range(ntex):
            n, = struct.unpack('<I', f.read(4))
            self.atlas.append((f.tell(), n))
            f.seek(n, 1)
        # extensions
        self.graphs = {}
        while True:
            m = f.read(4)
            if len(m) < 4:
                break
            if m == b'psgh':
                n, = struct.unpack('<I', f.read(4))
                d = inflate(f.read(n), 64 + 4 * self.npix)
                name = d[:64].split(b'\0')[0].decode()
                self.graphs[name] = np.frombuffer(d[64:], dtype='<i4').copy()
            else:
                break
        self._init_map_state()

    # ---- initial pixel -> contig mapping (mirrors PretextView LoadFile) ----
    def _init_map_state(self):
        n = self.npix
        orig = np.zeros(n, dtype=np.int64)
        rel = np.zeros(n, dtype=np.int64)
        total = np.float32(0.0)
        last = 0
        for i, fr in enumerate(self.fracs):
            total = np.float32(total + np.float32(fr))
            pixel = int(np.float64(n) * np.float64(total))
            r = 0
            while last < pixel and last < n:
                orig[last] = i
                rel[last] = r
                last += 1
                r += 1
        r_last = None
        while last < n:
            lc = self.n_contigs - 1
            orig[last] = lc
            if lc != orig[last - 1]:
                r = 0
            rel[last] = r
            r += 1
            last += 1
        self.orig_contig = orig
        self.rel_coord = rel

    def tex_index(self, x, y):
        """linear index of texture (x, y) with x <= y (row-major upper triangle)."""
        n = self.n_tex
        return x * n - (x * (x - 1)) // 2 + (y - x)

    def read_tex(self, x, y, level=0):
        base, nb = self.atlas[self.tex_index(x, y)]
        self.f.seek(base)
        raw = inflate(self.f.read(nb), self.bytes_per_tex)
        off = 0
        for l in range(level):
            off += ((self.tex_res >> l) ** 2) >> 1
        res = self.tex_res >> level
        return decode_bc4(raw[off:off + ((res * res) >> 1)], res)


def decode_bc4(buf, res):
    nb = (res // 4) ** 2
    b = np.frombuffer(buf, dtype=np.uint8).reshape(nb, 8)
    a0 = b[:, 0].astype(np.float32)
    a1 = b[:, 1].astype(np.float32)
    pal = np.zeros((nb, 8), dtype=np.float32)
    pal[:, 0] = a0
    pal[:, 1] = a1
    gt = a0 > a1
    for k in range(2, 8):
        pal[:, k] = np.where(gt, ((8 - k) * a0 + (k - 1) * a1) / 7.0, 0)
    for k in range(2, 6):
        pal[:, k] = np.where(gt, pal[:, k], ((6 - k) * a0 + (k - 1) * a1) / 5.0)
    pal[:, 6] = np.where(gt, pal[:, 6], 0.0)
    pal[:, 7] = np.where(gt, pal[:, 7], 255.0)
    bits = np.zeros(nb, dtype=np.uint64)
    for i in range(6):
        bits |= b[:, 2 + i].astype(np.uint64) << np.uint64(8 * i)
    codes = np.zeros((nb, 16), dtype=np.int64)
    for i in range(16):
        codes[:, i] = ((bits >> np.uint64(3 * i)) & np.uint64(7)).astype(np.int64)
    px = np.take_along_axis(pal, codes, axis=1)
    px = np.rint(px).astype(np.uint8)
    n4 = res // 4
    return px.reshape(n4, n4, 4, 4).transpose(0, 2, 1, 3).reshape(res, res)


# ----------------------------------------------------------------------------
# Savestate
# ----------------------------------------------------------------------------
class SaveState:
    """Parsed PretextView savestate (magic ptsx\\x02, 12-byte edits)."""

    def __init__(self, path=None):
        if path:
            self.load(path)

    def load(self, path):
        d = open(path, 'rb').read()
        assert d[:5] == b'ptsx\x02', d[:6]
        self.tail = d[5]
        off = 6
        if self.tail == 3:
            self.header_hash = d[off:off + 8]
            off += 8
        else:
            self.header_hash = None
        nc, nf = struct.unpack('<II', d[off:off + 8])
        off += 8
        raw = inflate(d[off:off + nc], nf)
        self.raw = raw
        p = 0
        self.settings = raw[p:p + 2]; p += 2
        self.colours = raw[p:p + 288]; p += 288
        self.sizes = raw[p:p + 28]; p += 28
        self.colour_map = raw[p]; p += 1
        self.gamma = raw[p:p + 12]; p += 12
        self.camera = raw[p:p + 12]; p += 12
        n_edits, = struct.unpack('<I', raw[p:p + 4]); p += 4
        self.edits = []
        flags = raw[p + 12 * n_edits: p + 12 * n_edits + ((n_edits + 7) >> 3)]
        for i in range(n_edits):
            x, y, dlt = struct.unpack('<IIi', raw[p:p + 12]); p += 12
            inv = bool(flags[(i + 1) >> 3] & (1 << ((i + 1) & 7)))
            self.edits.append((x, y, dlt, inv))
        p += (n_edits + 7) >> 3
        n_wp = raw[p]; p += 1
        self.waypoints = raw[p:p + 13 * n_wp]; p += 13 * n_wp
        self.n_wp = n_wp
        n_sc, = struct.unpack('<I', raw[p:p + 4]); p += 4
        self.scaffs = []
        for _ in range(n_sc):
            self.scaffs.append(struct.unpack('<II', raw[p:p + 8])); p += 8
        self.meta_active_tag = raw[p:p + 4]; p += 4
        self.meta_mode_data = raw[p:p + 36]; p += 36
        n_mf, = struct.unpack('<I', raw[p:p + 4]); p += 4
        self.meta_flags = []
        for _ in range(n_mf):
            self.meta_flags.append(struct.unpack('<IQ', raw[p:p + 12])); p += 12
        n_tags = raw[p]; p += 1
        self.tags = {}
        for _ in range(n_tags):
            tid = raw[p]; p += 1
            e = raw.index(b'\0', p)
            self.tags[tid] = raw[p:e].decode()
            p = e + 1
        self.ext = raw[p:]
        return self

    def dump(self, path):
        out = bytearray()
        out += self.settings + self.colours + self.sizes + bytes([self.colour_map]) + self.gamma + self.camera
        n = len(self.edits)
        out += struct.pack('<I', n)
        flags = bytearray((n + 7) >> 3)
        for i, (x, y, dlt, inv) in enumerate(self.edits):
            out += struct.pack('<IIi', x, y, dlt)
            if inv and dlt != BREAK_DELTA:
                flags[(i + 1) >> 3] |= (1 << ((i + 1) & 7))
        out += flags
        out += bytes([self.n_wp]) + self.waypoints
        out += struct.pack('<I', len(self.scaffs))
        for c, s in self.scaffs:
            out += struct.pack('<II', c, s)
        out += self.meta_active_tag + self.meta_mode_data
        out += struct.pack('<I', len(self.meta_flags))
        for c, fl in self.meta_flags:
            out += struct.pack('<IQ', c, fl)
        out += bytes([len(self.tags)])
        for tid in sorted(self.tags):
            out += bytes([tid]) + self.tags[tid].encode() + b'\0'
        out += self.ext
        comp = deflate(bytes(out))
        with open(path, 'wb') as f:
            f.write(b'ptsx\x02' + bytes([3]) + self.header_hash)
            f.write(struct.pack('<II', len(comp), len(out)))
            f.write(comp)
        return bytes(out)


class MapState:
    """Pixel-level simulation of PretextView's Map_State under edits."""

    def __init__(self, pt):
        self.n = pt.npix
        self.n_orig = pt.n_contigs
        self.pix = np.arange(self.n, dtype=np.int64)          # original pixel index at each position
        self.orig = pt.orig_contig.copy()                       # originalContigIds (may get += n_orig on breaks)
        self.rel = pt.rel_coord.copy()                          # contigRelCoords
        self.update_contigs()

    def update_contigs(self):
        o, r = self.orig, self.rel
        n = self.n
        # boundaries follow UpdateContigsFromMapState: new contig when orig id changes or coord is not +-1 continuation
        cid = np.zeros(n, dtype=np.int64)
        contigs = []  # (start_pix, length, orig_base, startCoord, inverted)
        start = 0
        inverted = r[1] < r[0]
        cp = 0
        for i in range(1, n):
            if o[i] != o[i - 1] or (inverted and r[i] != r[i - 1] - 1) or ((not inverted) and r[i] != r[i - 1] + 1):
                contigs.append((start, i - start, int(o[start] % self.n_orig), int(r[start]), bool(inverted), int(o[start])))
                cp += 1
                start = i
                if i < n - 1:
                    inverted = r[i + 1] < r[i]
            cid[i] = cp
        contigs.append((start, n - start, int(o[start] % self.n_orig), int(r[start]), bool(inverted), int(o[start])))
        self.cid = cid
        self.contigs = contigs

    def rearrange(self, a, b, delta):
        if a > b:
            a, b = b, a
        if delta == 0:
            return
        n = self.n
        for arr_name in ('pix', 'orig', 'rel'):
            arr = getattr(self, arr_name)
            seg = arr[a:b + 1].copy()
            if delta > 0:
                assert b + delta < n, 'wrap not supported'
                arr[a:a + delta] = arr[b + 1:b + 1 + delta].copy()
                arr[a + delta:b + delta + 1] = seg
            else:
                d = -delta
                assert a - d >= 0, 'wrap not supported'
                moved = arr[a - d:a].copy()
                arr[a - d:b - d + 1] = seg
                arr[b - d + 1:b + 1] = moved

    def invert(self, a, b):
        if a > b:
            a, b = b, a
        for arr_name in ('pix', 'orig', 'rel'):
            arr = getattr(self, arr_name)
            arr[a:b + 1] = arr[a:b + 1][::-1].copy()

    def break_at(self, loc):
        """BreakMap(loc, 1): cut so that loc is the last pixel of the left piece."""
        cid = self.cid[loc]
        s, ln = self.contigs[cid][0], self.contigs[cid][1]
        e = s + ln - 1
        l_len = loc - s + 1
        r_len = e - loc
        if l_len <= 1 or r_len <= 0:
            return False
        self.orig[loc + 1:e + 1] += self.n_orig
        return True

    def apply(self, edit):
        x, y, d, inv = edit
        if d == BREAK_DELTA:
            if self.break_at(x):
                self.update_contigs()
            return
        self.rearrange(x, y, d)
        if inv:
            self.invert(x + d, y + d)
        self.update_contigs()
