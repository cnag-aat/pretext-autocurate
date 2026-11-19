"""Breakpoints placed to the base, and the assembly edits that use them.

A savestate holds cuts at pixel resolution (tens of kb). When the sequence says exactly where a
cut belongs (the end of a haplotypic copy from an alignment; later inversion and move breakpoints
from the homolog, the edge of a repeat array, a telomere array), `run` writes it to
OUT.breaks.tsv, one row per break:

  kind        cut: one input scaffold is cut, the left piece ending at `bp` (1-based, included);
              switch: a gapless join: the output runs along `scaffold` to `bp` (its last base,
              1-based) and goes on in `partner` after `partner_bp`, the base identical to it
              (the two are copies there), with no gap between them
  scaffold    input scaffold
  px_bp       where the pixel cut puts it in PretextView's AGP (cut: the left line's end;
              switch: the scaffold's end coordinate that faces the join)
  bp, partner, partner_bp, method, evidence

`asm` uses the file (its own --breaks, or the run's): cuts move the matching line ends of the
PretextView AGP before pretext-to-asm (`refine_cuts`); switches remove the gap between the two
pieces in pretext-to-asm's output and rebuild the sequence (`apply_switches`). A row applies only
where the savestate still has that cut or join: a cut the curator moved or removed in PretextView
is left at its pixel position, and nothing else changes.
"""
import os

from .graph import AgpRow, read_agp, rebuild, write_agp_rows

COLS = ["kind", "scaffold", "px_bp", "bp", "partner", "partner_bp", "method", "evidence"]


def write(rows, path):
    with open(path, "w") as f:
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r.get(c, "")) for c in COLS) + "\n")


def read(path):
    out = []
    with open(path) as f:
        cols = f.readline().rstrip("\n").split("\t")
        for line in f:
            c = dict(zip(cols, line.rstrip("\n").split("\t")))
            for k in ("px_bp", "bp", "partner_bp"):
                c[k] = int(c[k] or 0)
            out.append(c)
    return out


def refine_cuts(agp_in, agp_out, rows, tol_bp):
    """Move cut line ends in a PretextView-style AGP from the pixel position to the exact base.
    Each cut row needs a line of `scaffold` ending within `tol_bp` of px_bp and another starting
    right after it. Returns [(row, done, note)]."""
    head, objs = read_agp(agp_in)
    lines = [(o, i) for o, rows_ in objs.items() for i, r in enumerate(rows_) if r.kind == "W"]
    report = []
    for c in (r for r in rows if r["kind"] == "cut"):
        left = [(o, i) for o, i in lines if objs[o][i].name == c["scaffold"] and abs(objs[o][i].end - c["px_bp"]) <= tol_bp]
        right = [(o, i) for o, i in lines if objs[o][i].name == c["scaffold"] and left
                 and objs[o][i].beg == objs[left[0][0]][left[0][1]].end + 1]
        if len(left) != 1 or len(right) != 1:
            report.append((c, False, "no cut there in the savestate (moved or removed): pixel position kept"))
            continue
        (lo, li), (ro, ri) = left[0], right[0]
        a, b = objs[lo][li], objs[ro][ri]
        if not (a.beg < c["bp"] < b.end):
            report.append((c, False, "exact position outside the two pieces: pixel position kept"))
            continue
        objs[lo][li] = a._replace(end=c["bp"])
        objs[ro][ri] = b._replace(beg=c["bp"] + 1)
        report.append((c, True, f"moved {c['bp'] - a.end:+,} bp"))
    write_agp_rows(agp_out, head, objs)
    return report


def _facing(r, right):
    """Coordinate of row r at its right (or left) side in the object."""
    return (r.end if r.orient == "+" else r.beg) if right else (r.beg if r.orient == "+" else r.end)


def _set_facing(r, right, v):
    if (r.orient == "+") == right:
        return r._replace(end=v)
    return r._replace(beg=v)


def _new_haplotig(files, parts, out_dir):
    """Add an object H_<n+1> made of `parts` to the haplotigs AGP/FASTA among `files` (created
    next to the primary files if pretext-to-asm wrote none). Returns its name."""
    import re
    n = max([int(m.group(1)) for fi in files for o in fi[3] for m in [re.fullmatch(r"H_(\d+)", o)] if m] or [0])
    name = f"H_{n + 1}"
    hap = next((fi for fi in files if "haplotig" in os.path.basename(fi[0])), None)
    if hap is None:
        prim = next((fi for fi in files if "primary" in os.path.basename(fi[0])), files[0])
        agp = prim[0].replace("primary", "haplotigs") if "primary" in prim[0] else prim[0][:-4] + ".haplotigs.agp"
        fa = agp[:-4] + ".fa"
        open(fa, "w").close()
        hap = [agp, fa, ["##agp-version\t2.1\n"], {}, {}]
        files.append(hap)
    rows = [p._replace(obj=name) for p in parts if p.beg <= p.end]
    hap[3][name] = rows
    hap[4][name] = rows
    return name


def apply_switches(out_dir, fa_in, rows, report_rows, tol_bp=2_000_000):
    """Gapless joins in pretext-to-asm's output (every NAME*.agp with its FASTA in `out_dir`).

    For a switch row: an object with `scaffold`'s line and `partner`'s line on either side of
    one gap, `scaffold` facing the gap with its end (within `tol_bp` of px_bp) and `partner`
    with its cut. `scaffold` is trimmed to `bp`, `partner` starts after `partner_bp` (it gains
    the part of its copy beyond the switch), and the gap goes. Then the object holding the cut
    copy of `partner` (one line, the Haplotig) ends at partner_bp and is followed by the rest of
    `scaffold` past bp, so no sequence is lost: the haplotig carries the other copy of each side.
    `fa_in` is the input assembly (graph.Faidx). Appends to `report_rows`; returns the number
    of joins made."""
    files = []
    for f in sorted(os.listdir(out_dir)):
        if not f.endswith(".agp") or ".pretext." in f:
            continue
        stem = os.path.join(out_dir, f[:-4])
        fa = next((stem + e for e in (".fa", ".fa.gz") if os.path.exists(stem + e)), None)
        if fa:
            head, objs = read_agp(os.path.join(out_dir, f))
            files.append([os.path.join(out_dir, f), fa, head, objs, {}])
    made = []
    for sw in (r for r in rows if r["kind"] == "switch"):
        K, T = sw["scaffold"], sw["partner"]
        Lk = fa_in.length(K)
        k_tail = sw["px_bp"] >= Lk // 2            # the kept end: its tail (or its head)
        gap_mode = K == T                           # a gap inside one scaffold: px_bp is its left side
        done = None
        for fi in files:
            for obj, rs in fi[3].items():
                for i in range(1, len(rs) - 1):
                    r1, g, r2 = rs[i - 1], rs[i], rs[i + 1]
                    if r1.kind != "W" or g.kind == "W" or r2.kind != "W" or {r1.name, r2.name} != {K, T}:
                        continue
                    if gap_mode:                    # the kept row ends exactly at the gap
                        if _facing(r1, True) == sw["px_bp"]:
                            kr, kright, tr, tright = r1, True, r2, False
                        elif _facing(r2, False) == sw["px_bp"]:
                            kr, kright, tr, tright = r2, False, r1, True
                        else:
                            continue
                    else:
                        kr, kright = (r1, True) if r1.name == K else (r2, False)
                        tr, tright = (r2, False) if r1.name == K else (r1, True)
                    kf = _facing(kr, kright)
                    if not gap_mode and abs(kf - (Lk if k_tail else 1)) > tol_bp:
                        continue
                    # the kept row ends at the switch; the partner row starts one base past its
                    # match, going away from the gap (rightwards in the object, or leftwards)
                    away = 1 if (tr.orient == "+") != tright else -1
                    t_in = sw["partner_bp"] + away
                    kr2, tr2 = _set_facing(kr, kright, sw["bp"]), _set_facing(tr, tright, t_in)
                    if kr2.beg > kr2.end or tr2.beg > tr2.end:
                        continue
                    old_t = _facing(tr, tright)
                    rs[i - 1:i + 2] = [kr2, tr2] if r1.name == K else [tr2, kr2]
                    fi[4][obj] = rs
                    done = (obj, kf, old_t, t_in)
                    break
                if done:
                    break
            if done:
                break
        if not done:
            report_rows.append(dict(file="", object="", left=K, right=T, graph_overlap="", overlap="",
                                    action="switch_skipped: the two pieces are not across one gap in the output"))
            continue
        obj, kf, old_t, t_in = done
        made.append(sw)
        # the haplotig: the cut copy of T, now ending at the switch, then K past the switch
        tail = None
        if gap_mode:
            # a new haplotig: T from the gap to its switch base, then K past its switch base up to
            # the gap (the other path through the bubble)
            tpart = AgpRow("", "W", T, min(old_t, sw["partner_bp"]), max(old_t, sw["partner_bp"]),
                           "+" if sw["partner_bp"] >= old_t else "-", ())
            kpart = AgpRow("", "W", K, sw["bp"] + 1, kf, "+", ()) if kf > sw["bp"] else \
                AgpRow("", "W", K, kf, sw["bp"] - 1, "-", ())
            tail = _new_haplotig(files, [tpart, kpart], out_dir)
        for fi in files if not gap_mode else []:
            for hobj, rs in fi[3].items():
                if len(rs) == 1 and rs[0].kind == "W" and rs[0].name == T and \
                        (rs[0].beg <= sw["partner_bp"] <= rs[0].end or abs(rs[0].end - old_t) <= 1 or abs(rs[0].beg - old_t) <= 1):
                    r = rs[0]
                    t_head = r.beg == 1                      # the copy is T's start
                    r = r._replace(beg=1, end=sw["partner_bp"], orient="+") if t_head else \
                        r._replace(beg=sw["partner_bp"], end=fa_in.length(T), orient="-")
                    kpart = AgpRow(hobj, "W", K, sw["bp"] + 1, Lk, "+", ()) if k_tail else \
                        AgpRow(hobj, "W", K, 1, sw["bp"] - 1, "-", ())
                    fi[3][hobj] = [r, kpart] if kpart.beg <= kpart.end else [r]
                    fi[4][hobj] = fi[3][hobj]
                    tail = hobj
                    break
            if tail:
                break
        report_rows.append(dict(file="", object=obj, left=f"{K}:{sw['bp']}", right=f"{T}:{t_in}", graph_overlap="",
                                overlap=f"{abs(kf - sw['bp'])}+{abs(old_t - t_in)}",
                                action="switch: gap removed, joined where the two copies are identical; " +
                                       (f"haplotig {tail} holds the other copy of each side" if tail else
                                        f"no haplotig object found for the copy of {T}: {abs(kf - sw['bp']):,} bp of {K} left out")))
    for path, fa, head, objs, changed in files:
        if changed:
            rebuild(path, fa, fa_in, head, objs, changed)
    return len(made)
