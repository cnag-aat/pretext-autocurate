"""Joins only hifiasm's graph sees, between scaffold ends (`graph.end_joins`), made in the
allhaps pass (curator, Troglobisium review 2026-11-17: all 15 such joins should be made).

  links       the copy ends the graph joins, directly or through one small scaffold (both its
              ends joined): evidence for `merge_pairs` that two pairs are one chromosome;
  place       each remaining join puts the smaller scaffold at the end of the chain of the
              chromosome copy the larger one ends, oriented as the graph says, taking it from
              wherever the tool had put it (unplaced, unloc, or a piece of another chain), unless
              it is a chromosome copy of its own.

Joins are given in scaffold terms (`graph.scaffold_end_joins`): dict(kind, sa, ea, sb, eb,
dist, overlap, size_a, size_b), ea / eb 'h' or 't'."""


def unit_end(pt, raw, name, end):
    """(unit of the cut map, its end 'h' / 't') at end `end` of input scaffold `name`, or None."""
    idx = {n: i for i, n in enumerate(raw.names)}
    s = idx.get(name)
    if s is None or raw.scaffold_npix[s] == 0:
        return None
    p = int(raw.scaffold_start[s]) + (0 if end == "h" else int(raw.scaffold_npix[s]) - 1)
    u = int(pt.orig_scaffold[p])
    a, b = int(pt.scaffold_start[u]), int(pt.scaffold_start[u] + pt.scaffold_npix[u] - 1)
    return (u, "h" if p == a else "t") if p in (a, b) else None


def _chain(run_, c):
    """Layout of copy c as [(first px, last px, reversed)]: its arrangement, else the unit as
    assembled with any inverted end piece turned (as `_layout`'s `pixels`)."""
    if c in run_.arranged:
        return list(run_.arranged[c])
    a, b = int(run_.start_[c]), int(run_.start_[c] + run_.npx_[c] - 1)
    out, lo = [], a
    for x, y in sorted(run_.flips):
        if a <= x and y <= b:
            if x > lo:
                out.append((lo, x - 1, False))
            out.append((x, y, True))
            lo = y + 1
    if lo <= b:
        out.append((lo, b, False))
    return out


def place(run_):
    """Each join (larger scaffold's end, smaller scaffold) puts the smaller one at that end of
    the chain of the copy the larger one ends. Returns decision rows (`graph_join`, or
    `graph_join_candidate` with why it was not made)."""
    pt, raw = run_.pt, run_.raw
    first = lambda ch: ch[0][1] if ch[0][2] else ch[0][0]
    last = lambda ch: ch[-1][0] if ch[-1][2] else ch[-1][1]
    beg = lambda x: x[1] if x[2] else x[0]            # first / last pixel of a chain piece
    fin = lambda x: x[0] if x[2] else x[1]
    idx = {n: i for i, n in enumerate(raw.names)}
    end_px = lambda name, e: int(raw.scaffold_start[idx[name]]) + (0 if e == "h" else int(raw.scaffold_npix[idx[name]]) - 1)
    # every pair of scaffold-end pixels the graph joins (for a small scaffold between two pieces)
    joinpx = {frozenset({end_px(j["sa"], j["ea"]), end_px(j["sb"], j["eb"])}) for j in run_.end_joins
              if j["sa"] in idx and j["sb"] in idx}
    rows = []
    # where each scaffold's graph joins would put it: V unit -> roots of the copies it joins
    target = {}
    for j in run_.end_joins:
        big, small = ("a", "b") if j["size_a"] >= j["size_b"] else ("b", "a")
        ue, ve = unit_end(pt, raw, j["s" + big], j["e" + big]), unit_end(pt, raw, j["s" + small], j["e" + small])
        if ue is not None and ve is not None:
            target.setdefault(ve[0], set()).add(run_.joined.get(ue[0], ue[0]))
    for j in run_.end_joins:
        big, small = ("a", "b") if j["size_a"] >= j["size_b"] else ("b", "a")
        S, eS, V, eV = j["s" + big], j["e" + big], j["s" + small], j["e" + small]
        how = (f"hifiasm's contig graph links the {'tail' if eS == 't' else 'head'} of {S} to the "
               f"{'tail' if eV == 't' else 'head'} of {V} (the only link at both ends)" if j["kind"] == "link" else
               f"a unitig path of {j['dist'] / 1e3:.0f} kb joins the {'tail' if eS == 't' else 'head'} of {S} to the "
               f"{'tail' if eV == 't' else 'head'} of {V}, each the other's nearest scaffold end")
        row = lambda action, text: rows.append(dict(
            action=action, scaffold=V, length_mb=round(j["size_" + small] / 1e6, 2), partner=S, strand="", tag="",
            score="", next_best="", evidence=f"{how}: {text}"))
        if S not in idx or V not in idx or 0 in (int(raw.scaffold_npix[idx[S]]), int(raw.scaffold_npix[idx[V]])):
            row("graph_join_candidate", "under one pixel at this map's resolution: a savestate cannot place it "
                                        "(PretextView cannot either); not made")
            continue
        ue, ve = unit_end(pt, raw, S, eS), unit_end(pt, raw, V, eV)
        if ue is None or ve is None:
            row("graph_join_candidate", "an end lies inside a cut piece; not made")
            continue
        u0, v = ue[0], ve[0]
        if int(pt.scaffold_npix[v]) != int(raw.scaffold_npix[idx[V]]):
            row("graph_join_candidate", f"{V} was cut; not made")
            continue
        root = run_.joined.get(u0, u0)
        if root not in run_.mate_of or root in run_.host_of:
            what = ("an unloc" if any(c["u"] == u0 for lst in run_.attach.values() for c in lst) else
                    "not in a chromosome")
            row("graph_join_candidate", f"{S} is {what}, so there is no chain to extend (PretextView cannot join "
                                        "unlocs either); not made")
            continue
        if run_.joined.get(v) == root:
            if v in getattr(run_, "graph_placed", set()):
                row("graph_join", "already placed in this chromosome by its other graph join")
            continue                                  # the layout already has it there
        mate_note, pair_ = "", None
        if v in run_.host_of:
            row("graph_join_candidate", f"{V} is part of a merged chromosome; not made")
            continue
        if v in run_.mate_of:
            # a small homolog pair the graph joins into the two copies of one chromosome: both go
            w = run_.mate_of[v]
            if run_.mate_of.get(root) not in target.get(w, set()):
                row("graph_join_candidate", f"{V} is a chromosome copy of its own (paired with {pt.names[w]}, which "
                                            "the graph does not join into the homologous copy); not made")
                continue
            pair_ = next(q for q in run_.pairs if q.kind != "extra" and v in (q.a, q.b))
            mate_note = (f"; its homolog {pt.names[w]} goes into the other copy by its own graph join, so the small "
                         "pair is no longer laid out on its own")
        ch = _chain(run_, root)
        p = end_px(S, eS)
        va, vb = int(pt.scaffold_start[v]), int(pt.scaffold_start[v] + pt.scaffold_npix[v] - 1)
        other = vb if eV == "h" else va               # V's other end
        at = next((i for i, x in enumerate(ch) if p in (beg(x), fin(x))), None)
        between = ""
        if p == last(ch):
            ch = ch + [(va, vb, eV == "t")]           # V's joined end comes first
        elif p == first(ch):
            ch = [(va, vb, eV == "h")] + ch           # V's joined end comes last
        elif at is not None and p in (beg(ch[at]), fin(ch[at])):
            # S's end is inside the chain: V goes right next to it and what the tool had put
            # there follows V (Troglobisium: H2_scaffold_23, H2_scaffold_126, H2_scaffold_63)
            after = p == fin(ch[at])
            nb = ch[at + 1] if after else ch[at - 1]
            joined_too = frozenset({other, beg(nb) if after else fin(nb)}) in joinpx
            ch = (ch[:at + 1] + [(va, vb, eV == "t")] + ch[at + 1:] if after else
                  ch[:at] + [(va, vb, eV == "h")] + ch[at:])
            between = (" between the two pieces whose ends the graph joins to its two ends" if joined_too else
                       f", between {S} and {pt.names[int(pt.orig_scaffold[nb[0]])]}, which the tool had joined there")
        else:
            row("graph_join_candidate", f"that end of {S} is not an end of a piece of its chromosome in this layout; "
                                        "not made")
            continue
        was = ""
        host = run_.joined.get(v)
        if host is not None:                          # taken out of another chain
            run_.arranged[host] = [x for x in run_.arranged.get(host, []) if not (va <= x[0] and x[1] <= vb)]
            was = f"; taken from the chromosome of {pt.names[host]}, where the tool had put it"
        for k, lst in run_.attach.items():
            if any(c["u"] == v for c in lst):
                run_.attach[k] = [c for c in lst if c["u"] != v]
                was = "; it was an unloc"
        for c in run_.small_calls:
            if c["u"] == v:
                c["action"] = "join"
                c["join_note"] = f"placed by hifiasm's graph at an end of {S} (see the graph_join row)"
        run_.arranged[root] = ch
        run_.joined[v] = root
        run_.graph_placed = getattr(run_, "graph_placed", set()) | {v}
        if pair_ is not None:
            pair_.kind = "merged"
        row("graph_join", f"placed there{between}, oriented as the graph says{was}{mate_note}")
    return rows


def links(pt, raw, joins):
    """{frozenset({(unit, end), (unit, end)})}: unit ends the graph joins, directly, or through
    one scaffold whose two ends it joins to two different units."""
    out, via = set(), {}
    for j in joins:
        x, y = unit_end(pt, raw, j["sa"], j["ea"]), unit_end(pt, raw, j["sb"], j["eb"])
        if x is None or y is None or x[0] == y[0]:
            continue
        out.add(frozenset({x, y}))
        via.setdefault(x[0], []).append((x[1], y))
        via.setdefault(y[0], []).append((y[1], x))
    for v, lst in via.items():
        ends = {e for e, _ in lst}
        if len(lst) == 2 and ends == {"h", "t"} and lst[0][1][0] != lst[1][1][0]:
            out.add(frozenset({lst[0][1], lst[1][1]}))
    return out
