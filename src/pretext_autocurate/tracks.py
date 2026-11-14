"""Extra evidence tracks as per-pixel arrays: telomeres by strand (telo_scan.py) and repeats (Red).

Each track can come from the map itself (added with PretextGraph), from a file given on the
command line (bedgraph or BED), or be computed from `--fasta` by running the external tool.
Nothing is downloaded; tools are looked for on PATH and in the conda environments
`pretext-autocurate` and `pretext-autocurate-x86` (Red and PretextGraph are built for Intel
Macs only and run under Rosetta).

Strand convention (telo_scan.py, PretextGraph names "5p_telomere" / "3p_telomere"):
telo_scan labels empirically, per motif, the strand that is enriched at scaffold starts as
"5p" and the one enriched at scaffold ends as "3p". So a 5p peak marks the left end of a
chromosome as assembled and a 3p peak its right end, whatever the motif.
"""
import gzip
import os
import re
import shutil
import subprocess
import sys

import numpy as np

from .homology import find_tool, split_by_haplotype

X86_ENV = "pretext-autocurate-x86"


def find_x86_tool(name, explicit=None):
    exe = find_tool(name, explicit)
    if exe:
        return exe
    for r in (os.environ.get("CONDA_PREFIX", ""), os.path.expanduser("~/miniconda3"), "/opt/miniconda3",
              os.path.expanduser("~/miniforge3"), os.path.expanduser("~/mambaforge")):
        cand = os.path.join(r, "envs", X86_ENV, "bin", name)
        if r and os.path.exists(cand):
            return cand
    return None


def _run(cmd, **kw):
    """Run a command; Intel-only binaries on Apple silicon go through Rosetta."""
    if sys.platform == "darwin" and os.uname().machine == "arm64" and X86_ENV in cmd[0]:
        cmd = ["arch", "-x86_64"] + cmd
    return subprocess.run(cmd, check=True, **kw)


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def to_pixels(pt, path, value=True, how="max"):
    """A bedgraph (name, start, end, value) or BED (name, start, end) file as an array over
    the map's pixels. BED: fraction of each pixel covered. Bedgraph: the maximum (or the
    coverage-weighted mean, how="mean") of the values overlapping each pixel. Names not in
    the map are ignored (returns the array and the number of lines used)."""
    raw = getattr(pt, "pt", pt)
    ix = {n: i for i, n in enumerate(raw.names)}
    bpp = raw.bp_per_pix
    out = np.zeros(raw.npix, dtype=np.float64)
    wsum = np.zeros(raw.npix, dtype=np.float64)
    used = 0
    with _open(path) as f:
        for line in f:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            c = line.split()
            i = ix.get(c[0])
            if i is None or raw.scaffold_npix[i] == 0:
                continue
            s, e = int(c[1]), int(c[2])
            v = float(c[3]) if value and len(c) > 3 else 1.0
            a0, n = int(raw.scaffold_start[i]), int(raw.scaffold_npix[i])
            p0, p1 = min(int(s / bpp), n - 1), min(int(max(e - 1, s) / bpp), n - 1)
            used += 1
            for p in range(p0, p1 + 1):
                ov = min(e, (p + 1) * bpp) - max(s, p * bpp)
                if ov <= 0:
                    continue
                if not value or len(c) <= 3:
                    out[a0 + p] += ov / bpp
                elif how == "mean":
                    out[a0 + p] += v * ov
                    wsum[a0 + p] += ov
                else:
                    out[a0 + p] = max(out[a0 + p], v)
    if how == "mean":
        out = np.where(wsum > 0, out / np.maximum(wsum, 1e-9), 0.0)
    return np.minimum(out, 1.0) if (not value) else out, used


def from_map(pt):
    """Tracks the map already carries: {'telomere_5p', 'telomere_3p', 'repeat'} (any subset)."""
    out = {}
    for name, g in pt.graphs.items():
        low = name.lower()
        if "5p" in low and "telo" in low:
            out["telomere_5p"] = g.astype(np.float64)
        elif "3p" in low and "telo" in low:
            out["telomere_3p"] = g.astype(np.float64)
        elif low in ("repeat_density", "repeats", "repeat", "red", "repeatmask", "repeat_mask"):
            g = g.astype(np.float64)
            out["repeat"] = g / max(float(np.percentile(g, 99.5)), 1e-9)       # scale unknown: relative
    return out


def run_telo_scan(fasta, workdir, motif, script=None, threads=8):
    """telo_scan.py on the assembly. Returns {'telomere_5p': bedgraph, 'telomere_3p': bedgraph}."""
    exe = script or shutil.which("telo_scan.py")
    if not exe or not os.path.exists(exe):
        raise RuntimeError("telo_scan.py not found: give its path with --telo-scan "
                           "(https://github.com/cnag-aat/scripts/blob/main/telo_scan.py)")
    os.makedirs(workdir, exist_ok=True)
    prefix = os.path.join(workdir, "telo_scan")
    m = motif.upper()
    bg5, bg3 = f"{prefix}.{m}.5p_telomere.bg", f"{prefix}.{m}.3p_telomere.bg"
    if not (os.path.exists(bg5) and os.path.exists(bg3)):
        with open(prefix + ".log", "w") as log:
            subprocess.run([sys.executable, exe, "-i", fasta, "-m", m, "-o", prefix, "--threads", str(threads),
                            "--no-hits", "--no-karyogram"], stdout=log, stderr=log, check=True)
    return {"telomere_5p": bg5, "telomere_3p": bg3}


def run_red(fasta, workdir, hap_regex=r"^(H\d+)_", red=None):
    """Red on each haplotype separately (homologous copies would otherwise count as repeats).
    Returns one BED of repeat intervals for the whole assembly (cached)."""
    exe = find_x86_tool("Red", red)
    if exe is None:
        raise RuntimeError("Red not found (PATH, conda env 'pretext-autocurate-x86'); install it with "
                           "CONDA_SUBDIR=osx-64 conda create -n pretext-autocurate-x86 -c conda-forge -c bioconda red")
    os.makedirs(workdir, exist_ok=True)
    bed = os.path.join(workdir, "red.repeats.bed")
    if os.path.exists(bed) and os.path.getsize(bed) > 0:
        return bed
    parts = split_by_haplotype(fasta, os.path.join(workdir, "split"), hap_regex)
    tmp = bed + ".tmp"
    with open(tmp, "w") as out:
        for hap, path in sorted(parts.items()):
            gnm, rpt = os.path.join(workdir, f"red_{hap}", "gnm"), os.path.join(workdir, f"red_{hap}", "rpt")
            os.makedirs(gnm, exist_ok=True)
            os.makedirs(rpt, exist_ok=True)
            link = os.path.join(gnm, f"{hap}.fa")
            if not os.path.exists(link):
                os.symlink(os.path.abspath(path), link)
            with open(os.path.join(workdir, f"red_{hap}.log"), "w") as log:
                _run([exe, "-gnm", gnm, "-rpt", rpt, "-frm", "2"], stdout=log, stderr=log)
            for fn in sorted(os.listdir(rpt)):
                if fn.endswith(".rpt"):
                    with open(os.path.join(rpt, fn)) as f:
                        for line in f:
                            c = line.split()
                            if len(c) >= 3:
                                out.write(f"{c[0].lstrip('>')}\t{c[1]}\t{c[2]}\n")
    os.replace(tmp, bed)
    shutil.rmtree(os.path.join(workdir, "split"), ignore_errors=True)   # the per-haplotype copies are large
    return bed


def load(pt, args_tracks):
    """Merge tracks into pt.graphs for the rules: 'telomere' (both strands, if strand tracks
    exist), 'telomere_5p', 'telomere_3p', 'repeat' (fraction of each pixel). `args_tracks`
    maps a track name to a file (bedgraph for telomeres, BED for repeats). Returns the names
    of the tracks loaded and their source."""
    got = {k: (v, "map") for k, v in from_map(pt).items()}
    for name, path in (args_tracks or {}).items():
        if path:
            arr, used = to_pixels(pt, path, value=(name != "repeat"))
            if used:
                got[name] = (arr, path)
    if "telomere_5p" in got and "telomere_3p" in got:
        t5, t3 = got["telomere_5p"][0], got["telomere_3p"][0]
        pt.graphs["telomere_5p"], pt.graphs["telomere_3p"] = t5, t3
        pt.graphs["telomere"] = t5 + t3
    if "repeat" in got:
        pt.graphs["repeat"] = got["repeat"][0]
    return {k: src for k, (_, src) in got.items()}


def telomere_threshold(track):
    """Peak level for a telomere track: telo_scan scores are 0-1000 (100 = default threshold);
    a map track has an unknown scale, so it is set relative to the track itself."""
    mx = float(np.max(track)) if len(track) else 0.0
    if 0 < mx <= 1000 and float(np.percentile(track, 99)) < 100:
        return 100.0
    return max(50.0, 5 * float(np.percentile(track, 99)))


def busco_tracks(hits, hap_of, out_dir, min_pair_genes=2, max_pairs=30):
    """BUSCO evidence worth showing on the map (curator, 2026-11-14), as PretextGraph inputs:

      busco_dup_pairs   genes with two copies in one haplotype, on two different scaffolds
                        (or far apart on one). Both copies of a pair of scaffolds get the
                        same label (d1, d2, ...), so PretextView colours the two places of a
                        possible haplotig or false duplication alike. Pairs sharing fewer
                        than `min_pair_genes` genes are left out (one gene is no evidence).
      busco_missing_hap genes present in one haplotype and absent from the other: shown where
                        the present copy is (the other haplotype lacks that region, or has it
                        unassembled or in a haplotig).

    Returns {track name: bedgraph path} for the tracks that have any interval."""
    os.makedirs(out_dir, exist_ok=True)
    haps = sorted({hap_of(sc) for loci in hits.values() for sc, *_ in loci if hap_of(sc)})
    pairs = {}
    for gene, loci in hits.items():
        by = {}
        for sc, a, b, *_ in loci:
            h = hap_of(sc)
            if h:
                by.setdefault(h, []).append((sc, a, b))
        for h, lst in by.items():
            lst = sorted(lst)
            for i in range(len(lst)):
                for j in range(i + 1, len(lst)):
                    x, y = lst[i], lst[j]
                    if x[0] == y[0] and abs(x[1] - y[1]) < 1_000_000:
                        continue                      # tandem copies close together: not a pair of places
                    pairs.setdefault((x[0], y[0]), []).append((gene, x, y))
    ranked = sorted(((len({g for g, _, _ in v}), k) for k, v in pairs.items()), reverse=True)
    ranked = [(n, k) for n, k in ranked if n >= min_pair_genes][:max_pairs]
    out = {}
    if ranked:
        path = os.path.join(out_dir, "busco_dup_pairs.bedgraph")
        with open(path, "w") as f, open(os.path.join(out_dir, "busco_dup_pairs.tsv"), "w") as t:
            t.write("label\tscaffold_1\tscaffold_2\tshared_genes\n")
            for k, (n, (s1, s2)) in enumerate(ranked):
                lab = f"d{k + 1}"
                t.write(f"{lab}\t{s1}\t{s2}\t{n}\n")
                for _, x, y in pairs[(s1, s2)]:
                    for sc, a, b in (x, y):
                        f.write(f"{sc}\t{a}\t{b}\t{lab}\n")
        out["busco_dup_pairs"] = path
    if len(haps) == 2:
        path = os.path.join(out_dir, "busco_missing_hap.bedgraph")
        n = 0
        with open(path, "w") as f:
            for gene, loci in hits.items():
                present = {hap_of(sc) for sc, *_ in loci}
                if len(present & set(haps)) == 1:
                    for sc, a, b, *_ in loci:
                        if hap_of(sc):
                            f.write(f"{sc}\t{a}\t{b}\t1\n")
                            n += 1
        if n:
            out["busco_missing_hap"] = path
    return out


def add_tracks(map_path, files, out_path, pretextgraph=None):
    """Write a copy of the map with the given tracks ({graph name: bedgraph or BED file})
    using PretextGraph. BED files (three columns) get a value of 1. Returns out_path."""
    exe = find_x86_tool("PretextGraph", pretextgraph)
    if exe is None:
        raise RuntimeError("PretextGraph not found (PATH, conda env 'pretext-autocurate-x86'); install it with "
                           "CONDA_SUBDIR=osx-64 conda create -n pretext-autocurate-x86 -c conda-forge -c bioconda pretextgraph")
    args, k = [], 0
    for name, path in files.items():
        if not path or not os.path.exists(path):
            continue
        with _open(path) as f:
            first = next((l for l in f if l.strip() and not l.startswith(("#", "track"))), "")
        if len(first.split()) == 3:                 # BED: give every interval the value 1
            bg = path + ".bedgraph"
            if not os.path.exists(bg):
                with _open(path) as f, open(bg, "w") as o:
                    for line in f:
                        c = line.split()
                        if len(c) >= 3:
                            o.write(f"{c[0]}\t{c[1]}\t{c[2]}\t1\n")
            path = bg
        k += 1
        args += [f"-t{k}", path, f"-n{k}", name]
    if not args:
        raise RuntimeError("no track files to add")
    log = out_path + ".pretextgraph.log"
    with open(log, "w") as lf:
        _run([exe, "-i", map_path, "-o", out_path] + args, stdout=lf, stderr=lf)
    return out_path
