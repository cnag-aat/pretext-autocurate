"""BUSCO genes as evidence: align a lineage's proteins to the assembly with miniprot and
count, per haplotype and per scaffold, which genes are present once, more than once, or
missing. Uses: choosing the more complete / less duplicated copy as HAP1, spotting false
duplications, and homology between scaffolds that share genes (also within one haplotype).

BUSCO sets are a guide, not truth: species have real duplications and real losses. Nothing
here acts on a single gene; calls need several genes agreeing, and are flagged, not applied.

A lineage is given as a directory (an unpacked BUSCO lineage, e.g. .../aves_odb12) or as a
name. A name is looked up in the cache directory and downloaded from the BUSCO server only
when the caller allows it; the file name and size are printed before the download starts.
"""
import os
import re
import subprocess
import sys
import tarfile
import urllib.request
from collections import defaultdict

from .homology import find_tool, _open

BUSCO_URL = "https://busco-data.ezlab.org/v6/data/lineages/"
CACHE = os.path.join(os.path.expanduser("~"), ".cache", "pretext-autocurate", "busco")


def list_remote(name):
    """[(file name, size text)] of the server's archives for a lineage name, newest last."""
    html = urllib.request.urlopen(BUSCO_URL, timeout=60).read().decode()
    pat = re.compile(r'href="(%s\.\d{4}-\d{2}-\d{2}\.tar\.gz)".*?\s(\d+[KMG]?)\s*$' % re.escape(name), re.M)
    return sorted(pat.findall(html))


def resolve_lineage(lineage, cache=CACHE, download=False, log=sys.stderr):
    """Directory of an unpacked lineage. `lineage` is a directory or a name like aves_odb12."""
    if os.path.isdir(lineage):
        return lineage
    local = os.path.join(cache, lineage)
    if os.path.isdir(local):
        return local
    if not download:
        raise RuntimeError(f"BUSCO lineage '{lineage}' is not a directory and is not in {cache}. "
                           "Give the path to an unpacked lineage, or allow the download (--busco-download).")
    found = list_remote(lineage)
    if not found:
        raise RuntimeError(f"no archive named {lineage}.<date>.tar.gz at {BUSCO_URL}")
    fname, size = found[-1]
    os.makedirs(cache, exist_ok=True)
    archive = os.path.join(cache, fname)
    print(f"downloading {BUSCO_URL}{fname} ({size}) to {archive}", file=log)
    urllib.request.urlretrieve(BUSCO_URL + fname, archive)
    with tarfile.open(archive) as t:
        t.extractall(cache, filter="data")
    os.remove(archive)
    if not os.path.isdir(local):
        raise RuntimeError(f"{fname} did not unpack to {local}")
    return local


def lineage_proteins(lineage_dir, variants=False):
    """Protein FASTA to align. Default is one ancestral sequence per gene: for aves_odb12
    that is about 8 thousand proteins against 50 thousand variants, and miniprot on a 2.7 Gb
    assembly takes roughly half an hour instead of several hours. `variants=True` uses all."""
    for name in (("ancestral_variants", "ancestral") if variants else ("ancestral", "ancestral_variants")):
        p = os.path.join(lineage_dir, name)
        if os.path.exists(p):
            return p
    raise RuntimeError(f"no ancestral or ancestral_variants file in {lineage_dir}")


def run_miniprot(genome, proteins, out_paf, threads=8, miniprot=None):
    """Protein-to-genome alignment (PAF). Cached: an existing non-empty output is reused."""
    if os.path.exists(out_paf) and os.path.getsize(out_paf) > 0:
        return out_paf
    exe = find_tool("miniprot", miniprot)
    if exe is None:
        raise RuntimeError("miniprot not found (PATH, --miniprot, or conda env 'pretext-autocurate').")
    os.makedirs(os.path.dirname(out_paf) or ".", exist_ok=True)
    tmp = out_paf + ".tmp"
    with open(tmp, "w") as out, open(out_paf + ".log", "w") as log:
        subprocess.run([exe, "-t", str(threads), "-I", "--outn", "4", genome, proteins], stdout=out, stderr=log, check=True)
    os.replace(tmp, out_paf)
    return out_paf


def read_hits(paf, min_cov=0.8):
    """Complete gene copies: {busco id: [(scaffold, start, end, strand, score)]}, one per locus.

    A hit is complete when it covers at least `min_cov` of the protein. Variants of one gene
    hitting the same place are merged (the best score is kept)."""
    loci = defaultdict(list)
    with _open(paf) as f:
        for line in f:
            if line.startswith("#"):
                continue
            c = line.split("\t")
            if len(c) < 12:
                continue
            qlen, qs, qe = int(c[1]), int(c[2]), int(c[3])
            if (qe - qs) / max(qlen, 1) < min_cov:
                continue
            gene = re.sub(r"_\d+$", "", c[0])
            score = next((int(x[5:]) for x in c[12:] if x.startswith("AS:i:")), int(c[9]))
            loci[gene].append((c[5], int(c[7]), int(c[8]), c[4], score))
    out = {}
    for gene, hits in loci.items():
        kept = []
        for h in sorted(hits, key=lambda h: -h[4]):
            if not any(k[0] == h[0] and h[1] < k[2] and k[1] < h[2] for k in kept):
                kept.append(h)
        # a second locus counts only if it scores close to the best one (paralogs score lower)
        out[gene] = [k for k in kept if k[4] >= 0.9 * kept[0][4]]
    return out


def read_full_tables(specs):
    """Several tables merged. Each spec is a path, or `PREFIX=path` for a table made on one
    haplotype alone: the prefix (e.g. H1) is put in front of every sequence name, with an
    underscore, so `scaffold_14` becomes `H1_scaffold_14` as in the diploid map."""
    out = defaultdict(list)
    for spec in specs:
        prefix, path = spec.split("=", 1) if "=" in spec and not os.path.exists(spec) else ("", spec)
        if prefix and not prefix.endswith("_"):
            prefix += "_"
        for gene, loci in read_full_table(path).items():
            out[gene] += [(prefix + sc, a, b, strand, score) for sc, a, b, strand, score in loci]
    return dict(out)


def false_duplicates(hits, hap_of, min_genes=3, min_frac=0.8):
    """Scaffolds whose genes nearly all have a second copy on one other scaffold of the same
    haplotype: candidates for a false duplication (or haplotig) of that region.

    Returns [(scaffold, other scaffold, start, end, shared genes, genes on scaffold)], the
    smaller-gene-count scaffold first; start/end delimit the copies on the other scaffold."""
    on = defaultdict(dict)
    for gene, loci in hits.items():
        for sc, a, b, _, _ in loci:
            on[sc].setdefault(gene, []).append((a, b))
    out = []
    for sc, genes in on.items():
        h = hap_of(sc)
        if not h or len(genes) < min_genes:
            continue
        other = defaultdict(list)
        for gene in genes:
            for sc2, a, b, _, _ in hits[gene]:
                if sc2 != sc and hap_of(sc2) == h:
                    other[sc2].append((gene, a, b))
        for sc2, lst in other.items():
            shared = {g for g, _, _ in lst}
            if len(shared) >= min_genes and len(shared) >= min_frac * len(genes) and len(on[sc2]) > len(genes):
                out.append((sc, sc2, min(a for _, a, _ in lst), max(b for _, _, b in lst), len(shared), len(genes)))
    return sorted(out, key=lambda r: -r[4])


def read_full_table(path):
    """Gene copies from a BUSCO `full_table.tsv` (BUSCO already run on this assembly).

    Same structure as `read_hits`. Rows with status Complete or Duplicated are kept, one
    entry per row, so a gene present once in each haplotype of a diploid assembly (which
    BUSCO reports as Duplicated) simply has two entries; who holds which copy is worked out
    from the scaffold names. Fragmented and Missing rows are ignored."""
    out = defaultdict(list)
    with _open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.rstrip("\n").split("\t")
            if len(c) < 7 or c[1] not in ("Complete", "Duplicated"):
                continue
            a, b = int(c[3]), int(c[4])
            out[c[0]].append((c[2], min(a, b), max(a, b), c[5], float(c[6])))
    return dict(out)


def summarise(hits, hap_of, n_genes=None):
    """Per haplotype: complete, single, duplicated, missing. Per scaffold: gene counts.

    `hap_of(scaffold name)` gives the haplotype label ('' to ignore)."""
    haps = defaultdict(lambda: dict(single=0, duplicated=0))
    per_scaffold = defaultdict(lambda: dict(genes=0, duplicated=0))
    genes_on = defaultdict(set)
    for gene, loci in hits.items():
        by = defaultdict(list)
        for sc, s, e, strand, score in loci:
            by[hap_of(sc)].append(sc)
            genes_on[sc].add(gene)
        for h, scs in by.items():
            if not h:
                continue
            haps[h]["duplicated" if len(scs) > 1 else "single"] += 1
            for sc in scs:
                per_scaffold[sc]["genes"] += 1
                per_scaffold[sc]["duplicated"] += len(scs) > 1
    total = n_genes or len(hits)
    for h in haps.values():
        h["complete"] = h["single"] + h["duplicated"]
        h["missing"] = total - h["complete"]
    return dict(haps), dict(per_scaffold), dict(genes_on)
