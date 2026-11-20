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


def _delta(old, new):
    """Bases a row gained and lost when its range went from `old` to `new` (AgpRows of one
    component): ([(lo, hi)] gained, [(lo, hi)] lost), 1-based inclusive."""
    gain, loss = [], []
    if new.beg < old.beg:
        gain.append((new.beg, old.beg - 1))
    if new.beg > old.beg:
        loss.append((old.beg, new.beg - 1))
    if new.end > old.end:
        gain.append((old.end + 1, new.end))
    if new.end < old.end:
        loss.append((new.end + 1, old.end))
    return gain, loss


def apply_switches(out_dir, fa_in, rows, report_rows, tol_bp=2_000_000):
    """Gapless joins in pretext-to-asm's output (every NAME*.agp with its FASTA in `out_dir`).

    For a switch row: an object with a line of `scaffold` and a line of `partner` on either side
    of one gap, `scaffold`'s facing the gap at px_bp (the nearest such line, within `tol_bp`).
    That line now ends at `bp`, `partner`'s starts one base past `partner_bp` (its base identical
    to bp's), and the gap goes. Each line may gain bases (they were in a haplotig, cut off as a
    copy) or lose them; nothing is lost: bases a line gains leave the one-line haplotig object
    that held them, and bases a line loses go to the haplotig (that same object, else a new one),
    so the haplotig holds the other path through the bubble. The cases: two scaffold ends that
    overlap (one line gains, the other loses), a gap with a copy on both sides (both lose), a
    piece filling a gap (one junction loses on both sides; at the other the piece may gain back
    part of a copy cut off it). `fa_in` is the input assembly (graph.Faidx). Appends to
    `report_rows`; returns the number of joins made."""
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
        best = None                                  # (distance, file, object, index, k on the left)
        for fi in files:
            for obj, rs in fi[3].items():
                for i in range(1, len(rs) - 1):
                    r1, g, r2 = rs[i - 1], rs[i], rs[i + 1]
                    if r1.kind != "W" or g.kind == "W" or r2.kind != "W" or {r1.name, r2.name} != {K, T}:
                        continue
                    for k_left in (True, False):
                        kr = r1 if k_left else r2
                        if kr.name != K:
                            continue
                        d = abs(_facing(kr, k_left) - sw["px_bp"])
                        if d <= tol_bp and (best is None or d < best[0]):
                            best = (d, fi, obj, i, k_left)
        if best is None:
            report_rows.append(dict(file="", object="", left=K, right=T, graph_overlap="", overlap="",
                                    action="switch_skipped: the two pieces are not across one gap in the output"))
            continue
        _, fi, obj, i, k_left = best
        rs = fi[3][obj]
        r1, r2 = rs[i - 1], rs[i + 1]
        kr, tr, kright, tright = (r1, r2, True, False) if k_left else (r2, r1, False, True)
        # the kept line ends at the switch; the partner's starts one base past its match, going
        # away from the gap (rightwards in the object, or leftwards)
        away = 1 if (tr.orient == "+") != tright else -1
        t_in = sw["partner_bp"] + away
        kr2, tr2 = _set_facing(kr, kright, sw["bp"]), _set_facing(tr, tright, t_in)
        if kr2.beg > kr2.end or tr2.beg > tr2.end:
            report_rows.append(dict(file="", object=obj, left=K, right=T, graph_overlap="", overlap="",
                                    action="switch_skipped: the switch lies outside the lines it joins"))
            continue
        rs[i - 1:i + 2] = [kr2, tr2] if k_left else [tr2, kr2]
        fi[4][obj] = rs
        made.append(sw)
        gained, lost = [], []
        for old, new in ((tr, tr2), (kr, kr2)):      # the partner's first: the bubble's other path
            g_, l_ = _delta(old, new)
            gained += [(old.name, a, b) for a, b in g_]
            lost += [(old.name, a, b) for a, b in l_]
        # a line that gained bases took them from the one-line haplotig that held them
        hobj = None
        for name, a, b in gained:
            for fj in files:
                for o, hrs in fj[3].items():
                    if o == obj or len(hrs) != 1 or hrs[0].kind != "W" or hrs[0].name != name:
                        continue
                    h = hrs[0]
                    if not (h.beg <= a and b <= h.end):
                        continue
                    if a == h.beg:
                        h = h._replace(beg=b + 1)
                    elif b == h.end:
                        h = h._replace(end=a - 1)
                    else:
                        continue                     # the gained stretch is not at either end
                    fj[3][o] = [h] if h.beg <= h.end else []
                    fj[4][o] = fj[3][o]
                    hobj = (fj, o)
        parts = [AgpRow("", "W", n, a, b, "+", ()) for n, a, b in lost if a <= b]
        if hobj is not None:
            fj, o = hobj
            fj[3][o] = fj[3][o] + [p._replace(obj=o) for p in parts]
            fj[4][o] = fj[3][o]
            tail = o
        elif parts:
            tail = _new_haplotig(files, parts, out_dir)
        else:
            tail = None
        moved_bp = sum(b - a + 1 for _, a, b in lost)
        report_rows.append(dict(file="", object=obj, left=f"{K}:{sw['bp']}", right=f"{T}:{t_in}", graph_overlap="",
                                overlap=f"{moved_bp}",
                                action="switch: gap removed, joined where the two copies are identical" +
                                       (f"; haplotig {tail} holds the other path ({moved_bp:,} bp moved there)"
                                        if tail else "")))
    for path, fa, head, objs, changed in files:
        if changed:
            rebuild(path, fa, fa_in, head, objs, changed)
    return len(made)
