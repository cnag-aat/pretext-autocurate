"""PretextView savestate reader/writer (magic ptsx\\x02; see docs/FORMATS.md)."""
import struct
from dataclasses import dataclass
from typing import NamedTuple

from .pretextio import inflate, deflate, header_hash

BREAK_DELTA = -2147483648
MAX_WAYPOINTS = 255
MAX_EDITS = 32768


class Edit(NamedTuple):
    """Move pixels [x, y] by delta, then reverse them if invert; or break after pixel x."""
    x: int
    y: int
    delta: int
    invert: bool = False

    @property
    def is_break(self):
        return self.delta == BREAK_DELTA


class Waypoint(NamedTuple):
    x: float
    y: float
    z: float
    id: int


@dataclass
class SaveState:
    head: bytes = b""           # theme, colours, sizes, colour map, gamma (331 bytes)
    camera: bytes = b""         # 3 x f32
    edits: list = None
    waypoints: list = None
    scaffs: list = None         # (contig index after edits, scaffold id)
    meta_head: bytes = b""      # active tag + tag display data (40 bytes)
    meta_flags: list = None     # (contig index after edits, tag bitmask)
    tags: dict = None           # tag id -> name
    ext: bytes = b""            # per-graph display settings
    header_hash: bytes = None
    tail: int = 3

    _HEAD = 2 + 288 + 28 + 1 + 12

    @classmethod
    def load(cls, path):
        with open(path, "rb") as f:
            d = f.read()
        if d[:5] != b"ptsx\x02":
            raise ValueError(f"{path}: not a version-2 PretextView savestate")
        s = cls(tail=d[5])
        off = 6
        if s.tail == 3:
            s.header_hash = d[off:off + 8]
            off += 8
        nc, nf = struct.unpack_from("<II", d, off)
        s._parse(inflate(d[off + 8:off + 8 + nc], nf))
        return s

    def _parse(self, raw):
        p = self._HEAD
        self.head = raw[:p]
        self.camera = raw[p:p + 12]
        p += 12
        (n,) = struct.unpack_from("<I", raw, p)
        p += 4
        flag_off = p + 12 * n
        self.edits = []
        for i in range(n):
            x, y, dlt = struct.unpack_from("<IIi", raw, p + 12 * i)
            # edit i uses bit i+1; for n % 8 == 0 the last bit falls in the byte after the
            # flag array (the waypoint count), exactly as PretextView reads it
            inv = bool(raw[flag_off + ((i + 1) >> 3)] >> ((i + 1) & 7) & 1)
            self.edits.append(Edit(x, y, dlt, inv and dlt != BREAK_DELTA))
        p = flag_off + ((n + 7) >> 3)
        n_wp = raw[p]
        p += 1
        self.waypoints = [Waypoint(*struct.unpack_from("<fffB", raw, p + 13 * i)) for i in range(n_wp)]
        p += 13 * n_wp
        (n_sc,) = struct.unpack_from("<I", raw, p)
        p += 4
        self.scaffs = [struct.unpack_from("<II", raw, p + 8 * i) for i in range(n_sc)]
        p += 8 * n_sc
        self.meta_head = raw[p:p + 40]
        p += 40
        (n_mf,) = struct.unpack_from("<I", raw, p)
        p += 4
        self.meta_flags = [struct.unpack_from("<IQ", raw, p + 12 * i) for i in range(n_mf)]
        p += 12 * n_mf
        n_tags = raw[p]
        p += 1
        self.tags = {}
        for _ in range(n_tags):
            tid = raw[p]
            e = raw.index(b"\0", p + 1)
            self.tags[tid] = raw[p + 1:e].decode()
            p = e + 1
        self.ext = raw[p:]

    @classmethod
    def blank(cls):
        """Empty savestate with PretextView's default tags and display settings."""
        import importlib.resources as ir
        s = cls(tail=3)
        s._parse(ir.files(__package__).joinpath("template.bin").read_bytes())
        return s

    def tag_id(self, name):
        for tid, n in self.tags.items():
            if n == name:
                return tid
        raise KeyError(name)

    def safe_edits(self):
        """Edit list that survives the nEdits % 8 == 0 flag quirk.

        With a multiple of 8 edits, PretextView reads the last edit's invert flag from bit 0
        of the waypoint count. If that disagrees with the flag we want, add one edit so the
        count is no longer a multiple of 8, without changing the result.
        """
        edits = list(self.edits)
        n = len(edits)
        if n == 0 or n % 8 or (len(self.waypoints) & 1) == int(edits[-1].invert):
            return edits
        last = edits[-1]
        if last.invert:   # split "move and invert" into a move and an in-place inversion
            edits[-1] = Edit(last.x, last.y, last.delta, False)
            edits.append(Edit(last.x + last.delta, last.y + last.delta, 0, True))
        else:             # inverting a single pixel in place changes nothing
            edits.append(Edit(0, 0, 0, True))
        return edits

    def payload(self):
        edits = self.safe_edits()
        if len(edits) > MAX_EDITS:
            raise ValueError(f"{len(edits)} edits; PretextView allows {MAX_EDITS}")
        if len(self.waypoints) > MAX_WAYPOINTS:
            raise ValueError(f"{len(self.waypoints)} waypoints; PretextView allows {MAX_WAYPOINTS}")
        n = len(edits)
        out = bytearray(self.head + self.camera + struct.pack("<I", n))
        flags = bytearray((n + 7) >> 3)
        for i, e in enumerate(edits):
            out += struct.pack("<IIi", e.x, e.y, e.delta)
            if e.invert and not e.is_break and ((i + 1) >> 3) < len(flags):
                flags[(i + 1) >> 3] |= 1 << ((i + 1) & 7)
        out += flags + bytes([len(self.waypoints)])
        for w in self.waypoints:
            out += struct.pack("<fffB", *w)
        out += struct.pack("<I", len(self.scaffs))
        for c, sid in self.scaffs:
            out += struct.pack("<II", c, sid)
        out += self.meta_head + struct.pack("<I", len(self.meta_flags))
        for c, fl in self.meta_flags:
            out += struct.pack("<IQ", c, fl)
        out += bytes([len(self.tags)])
        for tid in sorted(self.tags):
            out += bytes([tid]) + self.tags[tid].encode() + b"\0"
        return bytes(out + self.ext)

    def save(self, path, map_path=None):
        """Write the savestate. With map_path, bind it to that .pretext file (sets the hash)."""
        if map_path is not None:
            self.header_hash = header_hash(map_path)
        if self.header_hash is None:
            raise ValueError("no header hash: pass map_path")
        raw = self.payload()
        comp = deflate(raw)
        with open(path, "wb") as f:
            f.write(b"ptsx\x02\x03" + self.header_hash + struct.pack("<II", len(comp), len(raw)) + comp)

    def tag_names(self, mask):
        return [name for tid, name in sorted(self.tags.items()) if mask >> tid & 1]
