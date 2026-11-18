#!/usr/bin/env python3
"""Mark where tool savestates and a reference curation disagree, for viewing in PretextView.

    python3 tools/mark_disagreements.py MAP REF_SAVESTATE OUT_DIR TOOL_SAVESTATE [TOOL_SAVESTATE ...]
        [--label NAME ...] [--width 2]

Writes to OUT_DIR:
  - MAP copied under its own file name (the savestate hash covers the name and header, not the
    tracks, so every savestate of MAP still opens on the copy), plus one graph track per tool
    savestate, `vs_ref_<label>`, marking in the pixels of the input assembly (tracks move with
    the edits):
        4  extra break    the tool cut or moved here, the reference did not
        3  wrong join     the tool put these two pieces together, the reference does not
        2  missed break   the reference cut here, the tool did not
        1  missed join    the reference put these together, the tool did not
    Definitions and tolerances are those of `score.compare` (what score_all.py counts).
  - each tool savestate as <name>_<label>, its track shown in magenta; waypoints numbered from 1 at the
    disagreements in layout order if it had none (PretextView holds 255);
  - the reference savestate, with the first track shown;
  - <label>.disagreements.tsv: waypoint, kind, scaffold and Mb on both sides, layout pixel.
"""
import argparse
import os
import shutil
import struct
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "src"))
from pretext_autocurate import Pretext, SaveState                      # noqa: E402
from pretext_autocurate.pretextio import deflate                        # noqa: E402
from pretext_autocurate.savestate import MAX_WAYPOINTS, Waypoint        # noqa: E402
from pretext_autocurate.score import Layout, _breaks, _junctions, _match, _supported   # noqa: E402

KINDS = {4: "extra break", 3: "wrong join", 2: "missed break", 1: "missed join"}
# per-graph display settings (32 bytes): f32 scale, f32 base, f32 line width, u8 shown, 3 bytes
# as PretextView writes them, f32 r, g, b, a
MAGENTA = struct.pack("<fff4B4f", 0.05, 32.0, 2.0, 1, 0, 1, 0, 1.0, 0.0, 1.0, 1.0)


def n_extensions(path):
    pt = Pretext(path, load_graphs=False)
    f, n = pt.f, 0
    while f.read(4) == b"psgh":
        (c,) = struct.unpack("<I", f.read(4))
        f.seek(c, 1)
        n += 1
    pt.close()
    return n


def disagreements(pt, T, R, min_px=15, near_px=80, tol_px=3):
    """[(kind, pixel, partner pixel or None)] in input-assembly pixels."""
    sc = pt.orig_scaffold
    out = []
    bt, br = _breaks(pt, T), _breaks(pt, R)
    for a, b, kind in ((bt, br, 4), (br, bt, 2)):
        for p in a:
            if not _match(np.array([p]), b, sc, tol_px):
                out.append((kind, int(p), None))
    jt, jr = _junctions(T, min_px), _junctions(R, min_px)
    out += [(3, int(a), int(b)) for a, b in jt if not _supported((a, b), R, near_px)]
    out += [(1, int(a), int(b)) for a, b in jr if not _supported((a, b), T, near_px)]
    return out


def where(pt, p):
    s = int(pt.orig_scaffold[p])
    return pt.names[s], (p - int(pt.scaffold_start[s])) * pt.bp_per_pix / 1e6


def with_track(st, n_graphs, k, show):
    """Copy of st whose display settings cover n_graphs + k tracks; new track i shown if show[i]."""
    ext = bytearray(st.ext)
    old = [bytes(ext[i:i + 32]) for i in range(0, len(ext), 32)]
    head, tail = old[:n_graphs], old[n_graphs:]          # PretextView writes one more than it has
    hidden = MAGENTA[:12] + b"\0" + MAGENTA[13:]
    st.ext = b"".join(head + [MAGENTA if s else hidden for s in show] + tail)
    return st


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("map")
    ap.add_argument("ref")
    ap.add_argument("out_dir")
    ap.add_argument("tool", nargs="+")
    ap.add_argument("--label", nargs="*", help="one name per tool savestate (default: 1, 2, ...)")
    ap.add_argument("--width", type=int, default=2, help="pixels marked either side [2]")
    a = ap.parse_args()
    labels = a.label or [str(i + 1) for i in range(len(a.tool))]
    if len(labels) != len(a.tool):
        sys.exit("one --label per tool savestate")
    os.makedirs(a.out_dir, exist_ok=True)
    name = os.path.basename(a.map)
    out_map = os.path.join(a.out_dir, name)
    if os.path.abspath(out_map) == os.path.abspath(a.map):
        sys.exit("OUT_DIR must not be the map's own directory")
    pt = Pretext(a.map)
    n_graphs = n_extensions(a.map)
    ref = SaveState.load(a.ref)
    R = Layout(pt, ref)
    tracks = []
    for path, lab in zip(a.tool, labels):
        st = SaveState.load(path)
        T = Layout(pt, st)
        dis = disagreements(pt, T, R)
        v = np.zeros(pt.npix, dtype="<i4")
        for kind, p, q in dis:
            for x in (p, q) if q is not None else (p, p + 1):
                lo, hi = max(0, x - a.width), min(pt.npix, x + a.width + 1)
                v[lo:hi] = np.maximum(v[lo:hi], kind)
        tracks.append((f"vs_ref_{lab}", v))
        # waypoints and the table, in tool layout order
        dis.sort(key=lambda d: int(T.pos[d[1]]))
        add_wp = not st.waypoints
        with open(os.path.join(a.out_dir, f"{lab}.disagreements.tsv"), "w") as f:
            f.write("waypoint\tkind\tscaffold\tMb\tpartner\tpartner_Mb\tlayout_px\n")
            for k, (kind, p, q) in enumerate(dis):
                wp = k + 1 if add_wp and k < MAX_WAYPOINTS else ""
                if wp:
                    x = int(T.pos[p]) / pt.npix - 0.5
                    st.waypoints.append(Waypoint(x, x, 30.0, wp))
                s, mb = where(pt, p)
                s2, mb2 = where(pt, q) if q is not None else ("", None)
                f.write(f"{wp}\t{KINDS[kind]}\t{s}\t{mb:.2f}\t{s2}\t{'' if mb2 is None else f'{mb2:.2f}'}\t{int(T.pos[p])}\n")
        n = {KINDS[k]: sum(d[0] == k for d in dis) for k in KINDS}
        print(f"{lab}: " + ", ".join(f"{c} {kname}" for kname, c in n.items()))
        show = [t == len(tracks) - 1 for t in range(len(a.tool))]
        st = with_track(st, n_graphs, len(a.tool), show)
        st.save(os.path.join(a.out_dir, f"{os.path.basename(path)}_{lab}"), map_path=a.map)
    ref = with_track(ref, n_graphs, len(a.tool), [t == 0 for t in range(len(a.tool))])
    ref.save(os.path.join(a.out_dir, os.path.basename(a.ref)), map_path=a.map)
    pt.close()
    shutil.copyfile(a.map, out_map)
    with open(out_map, "ab") as f:
        for tname, v in tracks:
            body = deflate(tname.encode()[:63].ljust(64, b"\0") + v.tobytes())
            f.write(b"psgh" + struct.pack("<I", len(body)) + body)
    print(f"map with {len(tracks)} more track(s): {out_map}")


if __name__ == "__main__":
    main()
