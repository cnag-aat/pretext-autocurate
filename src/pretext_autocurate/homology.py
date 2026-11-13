"""Sequence homology between haplotypes: run minimap2 on the assembly, or read a PAF.

Optional input. With it, pairing of homologous scaffolds no longer depends on the MQ0 Hi-C
signal (which is weak for small chromosomes that contact each other a lot), and blocks
come with coordinates and strand. External tools are looked for on PATH and in the conda
environment `pretext-autocurate`; nothing is downloaded or installed by this module.
"""
import gzip
import os
import re
import shutil
import subprocess
from collections import defaultdict
from dataclasses import dataclass

TOOL_ENV = "pretext-autocurate"


def find_tool(name, explicit=None):
    """Path of an external program, or None. Order: explicit path, PATH, the conda env."""
    if explicit:
        return explicit if os.path.exists(explicit) else None
    hit = shutil.which(name)
    if hit:
        return hit
    roots = [os.environ.get("CONDA_PREFIX", ""), os.path.expanduser("~/miniconda3"), "/opt/miniconda3",
             os.path.expanduser("~/mambaforge"), os.path.expanduser("~/miniforge3")]
    for r in roots:
        for cand in (os.path.join(r, "envs", TOOL_ENV, "bin", name), os.path.join(r, "..", TOOL_ENV, "bin", name)):
            if r and os.path.exists(cand):
                return os.path.abspath(cand)
    return None


def _open(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def fasta_lengths(fasta):
    """{name: length} without loading sequences."""
    out, name, n = {}, None, 0
    with _open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                if name is not None:
                    out[name] = n
                name, n = line[1:].split()[0], 0
            else:
                n += len(line) - 1
    if name is not None:
        out[name] = n
    return out


def split_by_haplotype(fasta, workdir, hap_regex=r"^(H\d+)_"):
    """Write one FASTA per haplotype prefix into workdir; returns {hap: path}. Reuses files."""
    os.makedirs(workdir, exist_ok=True)
    done = os.path.join(workdir, "split.done")
    if os.path.exists(done):
        with open(done) as f:
            return dict(line.rstrip("\n").split("\t") for line in f)
    handles, paths, cur = {}, {}, None
    with _open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                m = re.match(hap_regex, line[1:])
                hap = m.group(1) if m else "unassigned"
                if hap not in handles:
                    paths[hap] = os.path.join(workdir, f"{hap}.fa")
                    handles[hap] = open(paths[hap], "w")
                cur = handles[hap]
            cur.write(line)
    for h in handles.values():
        h.close()
    with open(done, "w") as f:
        for k, v in paths.items():
            f.write(f"{k}\t{v}\n")
    return paths


def align_haplotypes(fasta, workdir, hap_regex=r"^(H\d+)_", threads=8, minimap2=None, preset="asm5"):
    """minimap2 of the second haplotype against the first. Returns the PAF path (cached)."""
    exe = find_tool("minimap2", minimap2)
    if exe is None:
        raise RuntimeError("minimap2 not found (PATH, --minimap2, or conda env 'pretext-autocurate'). "
                           "Install it, or give an existing alignment with --paf.")
    parts = split_by_haplotype(fasta, workdir, hap_regex)
    haps = sorted(h for h in parts if h != "unassigned")
    if len(haps) < 2:
        raise RuntimeError(f"need two haplotype prefixes in {fasta}, found {haps}")
    paf = os.path.join(workdir, f"{haps[1]}_vs_{haps[0]}.{preset}.paf")
    if os.path.exists(paf) and os.path.getsize(paf) > 0:
        return paf
    tmp = paf + ".tmp"
    with open(tmp, "w") as out, open(paf + ".log", "w") as log:
        subprocess.run([exe, "-x", preset, "-t", str(threads), "--secondary=no", parts[haps[0]], parts[haps[1]]],
                       stdout=out, stderr=log, check=True)
    os.replace(tmp, paf)
    return paf


def align_small(fasta, workdir, max_bp=10_000_000, threads=8, minimap2=None, preset="asm5"):
    """minimap2 of every sequence shorter than max_bp against the whole assembly (both
    haplotypes), ignoring each sequence's alignment to itself (-D) and keeping secondary
    alignments. This shows what the H2-against-H1 alignment cannot: a small scaffold's
    copies in its own haplotype (false duplications, haplotigs, same-prefix homologs) next
    to its homolog in the other. Returns the PAF path (cached)."""
    exe = find_tool("minimap2", minimap2)
    if exe is None:
        raise RuntimeError("minimap2 not found (PATH, --minimap2, or conda env 'pretext-autocurate'). "
                           "Install it, or give an existing alignment with --small-paf.")
    os.makedirs(workdir, exist_ok=True)
    paf = os.path.join(workdir, f"small_vs_all.{preset}.paf")
    if os.path.exists(paf) and os.path.getsize(paf) > 0:
        return paf
    lens = fasta_lengths(fasta)
    small = os.path.join(workdir, "small.fa")
    with _open(fasta) as f, open(small, "w") as out:
        keep = False
        for line in f:
            if line.startswith(">"):
                keep = lens[line[1:].split()[0]] < max_bp
            if keep:
                out.write(line)
    tmp = paf + ".tmp"
    try:
        with open(tmp, "w") as out, open(paf + ".log", "w") as log:
            subprocess.run([exe, "-x", preset, "-t", str(threads), "-D", "-N", "5", "-p", "0.5", fasta, small],
                           stdout=out, stderr=log, check=True)
    finally:
        os.remove(small)
    os.replace(tmp, paf)
    return paf


@dataclass
class Block:
    q: str
    qs: int
    qe: int
    t: str
    ts: int
    te: int
    same_strand: bool
    matches: int
    mapq: int
    primary: bool = True


def read_paf(path, min_len=10000, min_mapq=30):
    """Alignment blocks. Low mapping quality is dropped: repeats (e.g. on a W chromosome)
    otherwise align to unrelated small chromosomes and compete with the true homolog."""
    out = []
    with _open(path) as f:
        for line in f:
            c = line.split("\t", 12)
            if int(c[10]) >= min_len and int(c[11]) >= min_mapq:
                out.append(Block(c[0], int(c[2]), int(c[3]), c[5], int(c[7]), int(c[8]), c[4] == "+", int(c[9]), int(c[11]),
                                 "tp:A:S" not in (c[12] if len(c) > 12 else "")))
    return out


def summarise(blocks):
    """{(q, t): dict(bp, same_bp, q_span, t_span)}: aligned bases per scaffold pair and strand."""
    acc = defaultdict(lambda: dict(bp=0, same_bp=0, qmin=10**12, qmax=0, tmin=10**12, tmax=0))
    for b in blocks:
        a = acc[(b.q, b.t)]
        n = b.qe - b.qs
        a["bp"] += n
        a["same_bp"] += n if b.same_strand else 0
        a["qmin"], a["qmax"] = min(a["qmin"], b.qs), max(a["qmax"], b.qe)
        a["tmin"], a["tmax"] = min(a["tmin"], b.ts), max(a["tmax"], b.te)
    return dict(acc)
