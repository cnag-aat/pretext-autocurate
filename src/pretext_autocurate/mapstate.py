"""Pixel-level model of PretextView's map state under break / move / invert edits."""
from typing import NamedTuple

import numpy as np

from .savestate import Edit, BREAK_DELTA


class Contig(NamedTuple):
    """A piece as PretextView shows it: a run of pixels from one scaffold, in order."""
    start: int          # first pixel in the current map
    length: int         # pixels
    scaffold: int       # index of the original scaffold
    start_coord: int    # scaffold coordinate (pixels) of the first pixel
    inverted: bool
    raw_id: int         # scaffold index + k * n_scaffolds after explicit breaks

    @property
    def lo(self):
        """Lowest scaffold coordinate covered."""
        return self.start_coord - self.length + 1 if self.inverted else self.start_coord

    @property
    def hi(self):
        return self.start_coord if self.inverted else self.start_coord + self.length - 1


class MapState:
    def __init__(self, pt):
        self.n = pt.npix
        self.n_orig = pt.n_scaffolds
        self.pix = np.arange(self.n, dtype=np.int64)     # original pixel at each position
        self.orig = pt.orig_scaffold.astype(np.int64).copy()
        self.rel = pt.rel_coord.astype(np.int64).copy()
        self.update_contigs()

    def update_contigs(self):
        """Mirror of PretextView's UpdateContigsFromMapState.

        A new contig starts where the scaffold id changes or the coordinate stops running
        +1 (or -1 for an inverted contig); the direction is set by a contig's first step.
        """
        o, r, n = self.orig, self.rel, self.n
        step = np.diff(r)
        ok = (o[1:] == o[:-1]) & (np.abs(step) == 1)       # index j describes pixel j + 1
        change = np.zeros(n - 1, dtype=bool)
        change[1:] = step[1:] != step[:-1]
        starts = [0]
        for j in np.nonzero(~ok | change)[0]:
            i = int(j) + 1
            if not ok[j] or starts[-1] != i - 1:
                starts.append(i)
        starts = np.array(starts)
        lengths = np.diff(np.append(starts, n))
        contigs, inv = [], False
        for s, ln in zip(starts.tolist(), lengths.tolist()):
            if s < n - 1:
                inv = bool(r[s + 1] < r[s])
            contigs.append(Contig(s, ln, int(o[s] % self.n_orig), int(r[s]), inv, int(o[s])))
        self.contigs = contigs
        self.cid = np.repeat(np.arange(len(starts)), lengths)

    def _rearrange(self, a, b, delta):
        n = self.n
        if not (0 <= a + delta and b + delta < n):
            raise ValueError(f"move of [{a},{b}] by {delta} leaves the map")
        for arr in (self.pix, self.orig, self.rel):
            seg = arr[a:b + 1].copy()
            if delta > 0:
                arr[a:a + delta] = arr[b + 1:b + 1 + delta].copy()
            else:
                arr[b + delta + 1:b + 1] = arr[a + delta:a].copy()
            arr[a + delta:b + delta + 1] = seg

    def _invert(self, a, b):
        for arr in (self.pix, self.orig, self.rel):
            arr[a:b + 1] = arr[a:b + 1][::-1].copy()

    def _break(self, loc):
        c = self.contigs[self.cid[loc]]
        end = c.start + c.length - 1
        if loc - c.start + 1 <= 1 or end - loc <= 0:
            return False
        self.orig[loc + 1:end + 1] += self.n_orig
        return True

    def apply(self, edit, update=True):
        x, y, d, inv = edit
        if d == BREAK_DELTA:
            if self._break(x):
                self.update_contigs()
            return
        a, b = min(x, y), max(x, y)
        if d:
            self._rearrange(a, b, d)
        if inv:
            self._invert(a + d, b + d)
        if update:
            self.update_contigs()

    def replay(self, edits):
        """Apply edits in order. Breaks need the contig list, so it is refreshed lazily."""
        dirty = False
        for e in edits:
            e = Edit(*e)
            if e.is_break and dirty:
                self.update_contigs()
                dirty = False
            self.apply(e, update=False)
            dirty = dirty or not e.is_break
        if dirty:
            self.update_contigs()
        return self


def plan_edits(target, current=None):
    """Edits (moves, with inversion where needed) that turn `current` into `target`.

    Both are permutations giving the original pixel at each map position; `current`
    defaults to the unedited map. Works left to right: at the first position that differs,
    find the longest run of the target that is contiguous (forward or reversed) in the
    current map and move it into place. One edit per moved block, not guaranteed minimal.
    """
    target = np.asarray(target, dtype=np.int64)
    n = len(target)
    cur = np.arange(n, dtype=np.int64) if current is None else np.asarray(current, dtype=np.int64).copy()
    if not np.array_equal(np.sort(target), np.sort(cur)):
        raise ValueError("target is not a permutation of the current map")
    edits = []
    p = 0
    while True:
        diff = np.nonzero(cur[p:] != target[p:])[0]
        if len(diff) == 0:
            return edits
        p += int(diff[0])
        pos = np.empty(n, dtype=np.int64)
        pos[cur] = np.arange(n)
        q = pos[target[p:]]                        # where each remaining target pixel is now
        d = np.diff(q)
        fwd = int(np.argmax(d != 1)) + 1 if (d != 1).any() else len(q)
        rev = int(np.argmax(d != -1)) + 1 if (d != -1).any() else len(q)
        inv = rev > fwd
        ln = rev if inv else fwd
        a = int(q[0]) - ln + 1 if inv else int(q[0])
        b = a + ln - 1
        delta = p - a
        edits.append(Edit(a, b, delta, inv))
        seg = cur[a:b + 1].copy()
        cur[p + ln:b + 1] = cur[p:a].copy()
        cur[p:p + ln] = seg[::-1] if inv else seg
        p += ln
