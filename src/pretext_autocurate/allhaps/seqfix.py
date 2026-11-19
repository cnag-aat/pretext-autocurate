"""Edits the sequence places to the base, in allhaps (as in the single rules pass, see overlap.py
and breaks.py): contig-end overlaps within one haplotype, and gaps with a copy on both sides.

An overlap here is two ends of one haplotype's scaffolds that are copies of each other: one of
them is a second copy of sequence the haplotype already has. The copy on the longer scaffold is
cut off at the base the alignment gives and kept out of the chain as a FalseDuplicate at the end
of the chromosome copy that holds the other end (the curator's rule for second copies,
2026-11-14). The curators re-phased such copies when the other haplotype lacked the region
(bAegFun H1_scaffold_44 and 40), so every applied overlap is also a review row.

Gaps whose two sides are copies change nothing on the map; asm closes them (breaks.tsv)."""
import re

from ..overlap import gap_overlaps, overlap_evidence, resolve_overlaps
from ..pretextio import SegmentedMap


def _row(action, scaffold, mb, partner, ev, score=""):
    return dict(action=action, scaffold=scaffold, length_mb=mb, partner=partner, strand="", tag="", score=score,
                next_best="", evidence=ev)


def find(raw, cut_qs, fasta, hap_regex, threads=4, minimap2=None):
    """Overlaps among the pieces after `cut_qs` (ends of one haplotype only) and gaps with a copy
    on both sides. Returns (new cuts [(q, s, reason)], rows, breaks rows, applied overlap records)."""
    pre = lambda n: (re.match(hap_regex, n).group(1) if re.match(hap_regex, n) else "")
    vt = SegmentedMap(raw, sorted(cut_qs))
    recs = resolve_overlaps(raw, vt, fasta, minimap2, threads, group=pre)
    cuts, rows, exact, applied = [], [], [], []
    for rec in recs:
        ts, ks, al, o = rec["tside"], rec["kside"], rec["al"], rec["ov"]
        ev = overlap_evidence(rec, raw.bp_per_pix)
        mb = round((al.copy_bp if al is not None else o.overlap_px * raw.bp_per_pix) / 1e6, 2)
        if rec["apply"]:
            q = rec["q"]
            cuts.append((q, int(raw.orig_scaffold[q]),
                         f"end of a copy of {ks.name} {'start' if ks.head else 'end'} in the same haplotype, placed by "
                         f"the sequence at {rec['cut_bp']:,} bp (see the overlap row)"))
            exact += rec["exact"]
            # where the kept end is on the raw map (its piece end pixel)
            u, e = rec["keep"]
            rec["keep_px"] = int(vt.scaffold_start[u]) if e == "h" else int(vt.scaffold_start[u] + vt.scaffold_npix[u] - 1)
            applied.append(rec)
            rows.append(_row("overlap", ts.name, mb, ks.name, ev + f"; cut at {rec['cut_bp']:,} bp: the copy is a "
                             "second copy of this haplotype's sequence (FalseDuplicate); check whether the other "
                             "haplotype lacks it (re-phase it then)", round(o.excess, 2)))
        else:
            rows.append(_row("flag_overlap", ts.name, mb, ks.name, ev + "; " + rec["why"], round(o.excess, 2)))
    if raw.graphs.get("gap") is not None:
        for g in gap_overlaps(raw, fasta, minimap2, threads):
            al = g["al"]
            if al is None or al.copy_bp < 10_000:
                continue
            ev = (f"gap at {g['gap'][0] + 1:,} bp: the {al.copy_bp / 1e3:.0f} kb after it aligns in order with the "
                  f"sequence before it ({al.aligned_bp / 1e3:.0f} kb aligned, {al.divergence:.2%} divergence, reaching "
                  f"{al.keep_reach_bp / 1e3:.0f} kb from the gap)")
            if g["ok"]:
                exact.append(g["switch"])
                rows.append(_row("close_gap", g["scaffold"], round(al.copy_bp / 1e6, 2), "",
                                 ev + f"; asm closes it where {al.match_bp:,} bases are identical; the other copy of "
                                      "each side becomes a haplotig"))
            elif al.aligned_bp >= 0.5 * al.copy_bp and al.divergence <= 0.05:
                rows.append(_row("gap_copy", g["scaffold"], round(al.copy_bp / 1e6, 2), "",
                                 ev + "; not closed: the copy stops short of the gap or has no identical stretch to "
                                      "switch at"))
    return cuts, rows, [{k: v for k, v in r.items() if k != "src"} for r in exact], applied


def tag_copies(run_):
    """In _follow: each applied overlap's copy piece leaves its chain as a FalseDuplicate of the
    chromosome copy holding the kept end. Returns rows for copies that could not be tagged."""
    pt, joined, attach = run_.pt, run_.joined, run_.attach
    host_of, mate_of = getattr(run_, "host_of", {}), run_.mate_of
    out = []
    for rec in run_.seq_overlaps:
        ts, q = rec["tside"], rec["q"]
        v = int(pt.orig_scaffold[q if ts.head else q + 1])           # the copy piece
        K = int(pt.orig_scaffold[rec["keep_px"]])
        c = joined.get(K, K)
        c = host_of.get(c, c)
        if c not in mate_of:
            out.append(_row("flag_overlap", ts.name, "", rec["kside"].name,
                            f"{pt.names[v]} was cut off as a copy, but {pt.names[K]} is in no homolog pair: left as it is"))
            continue
        if v in joined:
            del joined[v]
            for lay in run_.arranged.values():
                lay[:] = [x for x in lay if not (int(pt.scaffold_start[v]) <= x[0] <= int(pt.scaffold_start[v] + pt.scaffold_npix[v] - 1))]
        cs = [x for x in run_.small_calls if x["u"] == v]
        why = (f"the end of {ts.name} that is another copy of {rec['kside'].name} {'start' if rec['kside'].head else 'end'} "
               f"(see the overlap row): a second copy of this haplotype's sequence, kept out of the chain")
        if cs:
            cs[0]["third_copy"], cs[0]["evidence"] = True, cs[0]["evidence"] + "; " + why
            cs[0]["copy"] = c
        else:
            x = dict(u=v, copy=c, action="unloc", third_copy=True, evidence=why, score="", ratio="")
            run_.small_calls.append(x)
            attach.setdefault(c, []).append(x)
    return out
