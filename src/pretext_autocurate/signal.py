"""Map-relative signal: distance decay, row strips at full resolution, coverage norm."""
import functools

import numpy as np


class Signal:
    def __init__(self, pt, max_dist=200, big_px=400, cache=80):
        self.pt = pt
        self.max_dist = max_dist
        self._tex = functools.lru_cache(maxsize=cache)(pt.read_tex)
        self.big = pt.scaffold_npix[pt.orig_scaffold] > big_px      # pixels in large scaffolds
        if not self.big.any():
            self.big = np.ones(pt.npix, dtype=bool)
        self._decay = None
        self._long = {}
        self._ins = {}
        self._binned = {}

    @property
    def decay(self):
        """Mean map value by pixel distance, within scaffolds (index = distance)."""
        if self._decay is None:
            pt, R, D = self.pt, self.pt.tex_res, self.max_dist
            tot, cnt = np.zeros(D), np.zeros(D)
            for t in range(0, pt.n_tex, max(1, pt.n_tex // 11)):
                m = self._tex(t, t).astype(np.float32)
                o = pt.orig_scaffold[t * R:(t + 1) * R]
                for d in range(D):
                    v = np.diagonal(m, d)[o[:R - d] == o[d:]]
                    tot[d] += v.sum()
                    cnt[d] += len(v)
            self._decay = tot / np.maximum(cnt, 1)
        return self._decay

    def adjacency_level(self, length):
        """Expected mean of the best off-diagonal line for a piece of this length that is
        merely adjacent to the partner locus (not a copy of it), relative to distance 0."""
        g = self.decay
        idx = np.minimum(2 * np.arange(length) + 1, len(g) - 1)
        return float(g[idx].mean() / g[0])

    def strip(self, a, b):
        """Rows a..b (inclusive) of the full-resolution map against every column."""
        pt, R = self.pt, self.pt.tex_res
        out = np.zeros((b - a + 1, pt.npix), dtype=np.uint8)
        for tx in range(a // R, b // R + 1):
            r0, r1 = max(a, tx * R), min(b, tx * R + R - 1)
            for ty in range(pt.n_tex):
                if tx <= ty:
                    t = self._tex(tx, ty)[r0 - tx * R:r1 - tx * R + 1, :]
                else:
                    t = self._tex(ty, tx)[:, r0 - tx * R:r1 - tx * R + 1].T
                out[r0 - a:r1 - a + 1, ty * R:(ty + 1) * R] = t
        return out

    def block(self, r0, r1, c0, c1):
        """Full-resolution rows r0..r1, columns c0..c1 (inclusive), clipped to the map."""
        pt, R = self.pt, self.pt.tex_res
        r0, c0 = max(r0, 0), max(c0, 0)
        r1, c1 = min(r1, pt.npix - 1), min(c1, pt.npix - 1)
        out = np.zeros((r1 - r0 + 1, c1 - c0 + 1), dtype=np.uint8)
        for tx in range(r0 // R, r1 // R + 1):
            ra, rb = max(r0, tx * R), min(r1, tx * R + R - 1)
            for ty in range(c0 // R, c1 // R + 1):
                ca, cb = max(c0, ty * R), min(c1, ty * R + R - 1)
                if tx <= ty:
                    t = self._tex(tx, ty)[ra - tx * R:rb - tx * R + 1, ca - ty * R:cb - ty * R + 1]
                else:
                    t = self._tex(ty, tx)[ca - ty * R:cb - ty * R + 1, ra - tx * R:rb - tx * R + 1].T
                out[ra - r0:rb - r0 + 1, ca - c0:cb - c0 + 1] = t
        return out

    def insulation(self, w, skip=0):
        """Contact across every boundary q|q+1, relative to the decay expectation.

        Sums the square of rows q-skip-w+1..q-skip against columns q+1+skip..q+skip+w and
        divides by what an intact scaffold gives at those distances. NaN where the square
        would leave the scaffold that holds q. Index q describes the boundary after pixel q.
        """
        if (w, skip) in self._ins:
            return self._ins[(w, skip)]
        pt, n, R = self.pt, self.pt.npix, self.pt.tex_res
        g = self.decay_long(2 * (w + skip) + 2)
        i = np.arange(w)
        expect = g[(i[:, None] + i[None, :] + 1 + 2 * skip)].sum()
        out = np.full(n, np.nan)
        pad = w + skip
        for s in range(0, n, R):
            lo, hi = max(0, s - pad), min(n - 1, s + R - 1 + pad)
            m = self.block(lo, hi, lo, hi).astype(np.float64)
            ii = np.zeros((m.shape[0] + 1, m.shape[1] + 1))
            ii[1:, 1:] = m.cumsum(0).cumsum(1)
            q = np.arange(max(s, pad - 1), min(s + R, n - pad))
            ra, rb = q - skip - w + 1 - lo, q - skip + 1 - lo        # row range [ra, rb)
            ca, cb = q + 1 + skip - lo, q + 1 + skip + w - lo
            out[q] = (ii[rb, cb] - ii[ra, cb] - ii[rb, ca] + ii[ra, ca]) / expect
        sc = pt.orig_scaffold
        q = np.arange(n)
        ok = (q - pad + 1 >= 0) & (q + pad < n)
        ok[ok] &= (sc[q[ok] - pad + 1] == sc[q[ok]]) & (sc[q[ok] + pad] == sc[q[ok]])
        out[~ok] = np.nan
        self._ins[(w, skip)] = out
        return out

    def decay_long(self, max_dist):
        """Decay curve out to max_dist (computed separately from the default 200 px)."""
        if max_dist <= self.max_dist:
            return self.decay
        key = int(max_dist)
        if key not in self._long:
            pt, R = self.pt, self.pt.tex_res
            tot, cnt = np.zeros(key), np.zeros(key)
            for t in range(0, pt.n_tex, max(1, pt.n_tex // 11)):
                m = self._tex(t, t).astype(np.float32)
                o = pt.orig_scaffold[t * R:(t + 1) * R]
                for d in range(min(key, R)):
                    v = np.diagonal(m, d)[o[:R - d] == o[d:]]
                    tot[d] += v.sum()
                    cnt[d] += len(v)
            self._long[key] = tot / np.maximum(cnt, 1)
        return self._long[key]

    def binned(self, f):
        """Whole map averaged from full resolution into f x f pixel bins (float32).

        The file's own mip levels are re-normalised and saturate, so their values cannot
        be compared with a decay curve; this average can. Needed for sparse maps (high MQ),
        where full-resolution windows of tens of pixels hold too little signal.
        """
        if f in self._binned:
            return self._binned[f]
        pt, R = self.pt, self.pt.tex_res
        if R % f:
            raise ValueError("bin size must divide the texture size")
        n, k = pt.npix // f, R // f
        out = np.zeros((n, n), dtype=np.float32)
        for x in range(pt.n_tex):
            for y in range(x, pt.n_tex):
                t = pt.read_tex(x, y).astype(np.float32).reshape(k, f, k, f).mean(axis=(1, 3))
                out[x * k:(x + 1) * k, y * k:(y + 1) * k] = t
                if x != y:
                    out[y * k:(y + 1) * k, x * k:(x + 1) * k] = t.T
        self._binned[f] = out
        return out

    def insulation_binned(self, f, w, skip=0):
        """Like `insulation`, on the map binned by f. Returns (values per bin, decay per bin).
        Index b describes the boundary after bin b; NaN where the window leaves the scaffold."""
        M = self.binned(f).astype(np.float64)
        n = M.shape[0]
        sc = self.pt.orig_scaffold[np.minimum(np.arange(n) * f + f // 2, self.pt.npix - 1)]
        D = 2 * (w + skip) + 2
        g = np.zeros(D)
        for d in range(D):
            v = np.diagonal(M, d)
            same = sc[:n - d] == sc[d:]
            g[d] = v[same].mean() if same.any() else 0.0
        ii = np.zeros((n + 1, n + 1))
        ii[1:, 1:] = M.cumsum(0).cumsum(1)
        i = np.arange(w)
        expect = g[i[:, None] + i[None, :] + 1 + 2 * skip].sum()
        out = np.full(n, np.nan)
        q = np.arange(w + skip - 1, n - w - skip)
        ra, rb, ca, cb = q - skip - w + 1, q - skip + 1, q + 1 + skip, q + 1 + skip + w
        out[q] = (ii[rb, cb] - ii[ra, cb] - ii[rb, ca] + ii[ra, ca]) / max(expect, 1e-9)
        ok = np.zeros(n, dtype=bool)
        ok[q] = (sc[q - w - skip + 1] == sc[q]) & (sc[q + w + skip] == sc[q])
        out[~ok] = np.nan
        return out, g

    def coverage_ratio(self):
        """Coverage track divided by its median over large scaffolds; None if no track."""
        cov = self.pt.graphs.get("coverage")
        if cov is None:
            return None
        norm = np.median(cov[self.big])
        return cov / norm if norm > 0 else None
