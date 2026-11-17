"""Step 2 of docs/ALLHAPS_SPLIT.md: move allhaps.py's top-level definitions, unchanged, into
an allhaps/ package; each module imports what its functions use. Prints the plan and checks
for import cycles; with --write, writes the package and removes allhaps.py."""
import ast
import collections
import os
import sys

SRC = "src/pretext_autocurate/allhaps.py"
PKG = "src/pretext_autocurate/allhaps"
MODULES = {
    "common": ["Level", "Pair", "_mean0", "_signal", "_level", "_union", "write_rows"],
    "pairing": ["_quality", "pair_homologs", "same_prefix_pairs", "merge_pairs"],
    "sites": ["telomere_sites", "alignment_breaks"],
    "evidence": ["_hold_mate", "_unitig_evidence", "_graph_evidence"],
    "small": ["_small_alignments", "place_small"],
    "arrange": ["_adjacency", "_homolog_dir", "_block_index", "_aligned_dir", "_aligned_bp", "_project",
                "_contig_conflicts", "arrange_pieces"],
    "driver": ["curate_allhaps"],
}
DOCS = {
    "common": "Shared pieces of the allhaps pass: `Level` (the conservative-to-aggressive setting), `Pair`,\nthe per-map `Signal` cache and small helpers.",
    "pairing": "Homolog pairs: `pair_homologs` (alignment blocks or the MQ0 homology line),\n`same_prefix_pairs` (homologs in one haplotype file), `merge_pairs` (pairs that are\npieces of one chromosome).",
    "sites": "Break sites at gaps and sharp breaks (`telomere_sites`), each decided cut / invert_tip /\nkeep, and breaks of colinearity with the other haplotype (`alignment_breaks`).",
    "evidence": "hifiasm graph evidence at break sites: contig-graph links and bridges\n(`_graph_evidence`) and unitig paths (`_unitig_evidence`).",
    "small": "Small scaffolds: which chromosome and copy (`place_small`): unloc, candidate, join when\nit extends a copy, FalseDuplicate when it is a second copy.",
    "arrange": "Order and orientation of the pieces of a chromosome copy: end contact\n(`arrange_pieces`, `_adjacency`) and the other haplotype's alignment (`_project`,\n`_aligned_dir`, `_contig_conflicts`).",
    "driver": "The allhaps pass (`curate_allhaps`): sites and cuts, pairing, small scaffolds, merges,\narrangement, layout, rows and savestate.",
}

src = open(SRC).read()
lines = src.splitlines(keepends=True)
tree = ast.parse(src)
top = {n.name: n for n in tree.body if isinstance(n, (ast.FunctionDef, ast.ClassDef))}
imports = [n for n in tree.body if isinstance(n, (ast.Import, ast.ImportFrom))]
assigned = {name for m, names in MODULES.items() for name in names}
assert assigned == set(top), (set(top) - assigned, assigned - set(top))
home = {name: m for m, names in MODULES.items() for name in names}
imported_names = {}                              # name -> import statement source
for n in imports:
    for a in n.names:
        imported_names[a.asname or a.name.split(".")[0]] = (n, a)


def used(node):
    return {x.id for x in ast.walk(node) if isinstance(x, ast.Name)} | \
           {x.value.id for x in ast.walk(node) if isinstance(x, ast.Attribute) and isinstance(x.value, ast.Name)}


def segment(name):
    """Source of a definition with the comment lines right above it."""
    n = top[name]
    start = (n.decorator_list[0].lineno if getattr(n, "decorator_list", None) else n.lineno) - 1
    while start > 0 and lines[start - 1].lstrip().startswith("#"):
        start -= 1
    return "".join(lines[start:n.end_lineno])


deps = collections.defaultdict(set)
plan = {}
for m, names in MODULES.items():
    u = set().union(*(used(top[n]) for n in names))
    local = {x for x in u if x in top and home[x] != m}
    ext = sorted(x for x in u if x in imported_names)
    for x in local:
        deps[m].add(home[x])
    plan[m] = (local, ext)
# cycle check
def cyc(m, stack):
    for d in deps[m]:
        if d in stack:
            return stack + [d]
        r = cyc(d, stack + [d])
        if r:
            return r
    return None
for m in MODULES:
    c = cyc(m, [m])
    if c:
        print("CYCLE:", " -> ".join(c))
for m, (local, ext) in plan.items():
    print(f"{m:9} {sum(top[n].end_lineno - top[n].lineno + 1 for n in MODULES[m]):5} lines; "
          f"from siblings: {sorted(local)}; external: {ext}")

if "--write" in sys.argv:
    os.makedirs(PKG, exist_ok=True)
    for m, names in MODULES.items():
        local, ext = plan[m]
        out = [f'"""{DOCS[m]}"""\n']
        std = sorted({imported_names[x][0].lineno for x in ext})
        for ln in std:
            out.append(lines[ln - 1].replace("from .", "from .."))
        by_mod = collections.defaultdict(list)
        for x in sorted(local):
            by_mod[home[x]].append(x)
        for mm in sorted(by_mod):
            out.append(f"from .{mm} import {', '.join(by_mod[mm])}\n")
        out.append("\n\n")
        out.append("\n\n".join(segment(n).rstrip("\n") + "\n" for n in names))
        with open(os.path.join(PKG, m + ".py"), "w") as f:
            f.write("".join(out))
    doc = ast.get_docstring(tree)
    names_all = [n for m in MODULES for n in MODULES[m]]
    init = [f'"""{doc}\n\nModules: ' + ", ".join(MODULES) + ".\nA diagnostic that monkeypatches a step "
            "(e.g. `place_small`) must patch it where it is used (`allhaps.driver`), not here.\n\"\"\"\n"]
    for m in MODULES:
        init.append(f"from .{m} import {', '.join(MODULES[m])}  # noqa: F401\n")
    with open(os.path.join(PKG, "__init__.py"), "w") as f:
        f.write("".join(init))
    os.remove(SRC)
    print("written", PKG)
