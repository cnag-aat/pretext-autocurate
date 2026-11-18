"""Non-target sequence (contaminants, cobionts, organelles): tagged `Contaminant`, kept apart.

Sources:
  - a list (`--contaminants FILE`): one sequence name per line in the first column, columns split
    on commas, spaces or tabs; blank lines and lines starting with '#' are skipped. A second
    column, if any, is kept as the class (e.g. "bacteria", "mitochondrion").
  - Tiara (`--tiara`, https://github.com/ibe-uw/tiara): a k-mer classifier of eukaryotic vs
    prokaryotic sequence, with organelles as a second stage. Sequences it calls bacteria,
    archaea, prokarya, mitochondrion or plastid (or organelle) are tagged; eukarya and unknown
    are not. Only scaffolds up to `--tiara-max-mb` are classified: a chromosome-sized scaffold
    is target sequence whatever its k-mers say.

Every non-target scaffold gets PretextView's `Contaminant` tag. pretext-to-asm (agp-tpf-utils)
writes those to its contaminants file. Organelles are tagged `Contaminant` too: pretext-to-asm
reads any tag it does not know as a haplotype name, so a "Mito" tag would make a second
haplotype (or a tagging error next to HAP1); organelle genomes come from their own assembly
(MitoHiFi, oatk) in pretext-to-asm's YAML. The class is in `<prefix>.contaminants.tsv` and the
decision rows.

The adjacency solver also keeps these scaffolds out of new joins, pairs and paint (`exclude`).
"""
import hashlib
import os
import re
import shutil
import subprocess
import sys

import numpy as np

from .mapstate import MapState

TIARA_ENVS = ("tiara-x86", "tiara", "pretext-autocurate-x86", "pretext-autocurate")
# Tiara classes (first stage, then second stage for organelles) that are not nuclear target
NON_TARGET = {"bacteria": "bacteria", "archaea": "archaea", "prokarya": "prokarya",
              "mitochondrion": "mitochondrion", "plastid": "plastid", "organelle": "organelle"}


def read_list(path):
    """{name: class} from a delimited list (first column the name; second, if any, the class). A
    header line with a `class` column (this tool's `.contaminants.tsv`) names the class column."""
    out, k = {}, 1
    with open(path) as f:
        for i, line in enumerate(f):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            c = [x for x in re.split(r"[,\s]+", s) if x]
            if not out and "class" in c[1:]:
                k = c.index("class")
                continue
            out[c[0]] = c[k] if len(c) > k else "listed"
    return out


def find_tiara(explicit=None):
    """Path of the tiara script, or None: explicit path, PATH, then the conda envs above (on Apple
    silicon Tiara needs an Intel env: it pins pytorch 1.7, which has no arm64 build)."""
    if explicit:
        return explicit if os.path.exists(explicit) else None
    hit = shutil.which("tiara")
    if hit:
        return hit
    for r in (os.environ.get("CONDA_PREFIX", ""), os.path.expanduser("~/miniconda3"), "/opt/miniconda3",
              os.path.expanduser("~/miniforge3"), os.path.expanduser("~/mambaforge")):
        for env in TIARA_ENVS:
            cand = os.path.join(r, "envs", env, "bin", "tiara")
            if r and os.path.exists(cand):
                return cand
    return None


def _write_subset(fasta, keep, out):
    """The sequences named in `keep` from `fasta` (plain or .gz) into `out`."""
    import gzip
    op = gzip.open if fasta.endswith(".gz") else open
    on = False
    with op(fasta, "rt") as f, open(out, "w") as w:
        for line in f:
            if line.startswith(">"):
                on = line[1:].split()[0] in keep
            if on:
                w.write(line)


def run_tiara(fasta, workdir, lengths, max_mb=10.0, threads=4, tiara=None, min_len=3000):
    """Tiara on the scaffolds of `fasta` up to `max_mb` (cached in `workdir`). Returns the path of
    Tiara's output table. `lengths`: {name: bp}."""
    exe = find_tiara(tiara)
    if exe is None:
        raise RuntimeError("tiara not found (PATH, conda envs " + ", ".join(TIARA_ENVS) + "); install it with "
                           "conda create -n tiara -c conda-forge -c bioconda tiara (on Apple silicon: "
                           "CONDA_SUBDIR=osx-64 conda create -n tiara-x86 -c pytorch -c conda-forge -c bioconda tiara)")
    keep = sorted(n for n, l in lengths.items() if min_len <= l <= max_mb * 1e6)
    key = hashlib.sha1((os.path.abspath(fasta) + f"{max_mb}" + "\n".join(keep)).encode()).hexdigest()[:10]
    os.makedirs(workdir, exist_ok=True)
    out = os.path.join(workdir, f"tiara_{key}.txt")
    if os.path.exists(out) and os.path.getsize(out) > 0:
        return out
    sub = os.path.join(workdir, f"tiara_{key}.fa")
    _write_subset(fasta, set(keep), sub)
    cmd = [exe, "-i", sub, "-o", out + ".tmp", "-t", str(threads), "-m", str(min_len)]
    if sys.platform == "darwin" and os.uname().machine == "arm64" and "x86" in exe:
        cmd = ["arch", "-x86_64"] + cmd          # an Intel env: run under Rosetta
    with open(os.path.join(workdir, f"tiara_{key}.log"), "w") as log:
        try:
            subprocess.run(cmd, check=True, stdout=log, stderr=log)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"tiara failed ({e.returncode}); see {log.name}")
    os.replace(out + ".tmp", out)
    os.remove(sub)
    return out


def read_tiara(path):
    """{name: class} for the non-target calls in a Tiara output table."""
    out = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 2 or c[0] == "sequence_id":
                continue
            first = c[1].strip().lower()
            second = c[2].strip().lower() if len(c) > 2 else ""
            cls = second if first == "organelle" and second in NON_TARGET else first
            if cls in NON_TARGET:
                out[c[0].split()[0]] = NON_TARGET[cls]
    return out


def tag(pt, st, classes):
    """Tag every contig of the savestate's final layout that comes from a listed scaffold
    `Contaminant` and nothing else: no haplotype, sex chromosome or Unloc tag, no paint (it is
    not part of either haplotype's assembly). Returns the contig indices that sit inside a painted
    scaffold that also holds target sequence (a non-target piece inside a chromosome, possible
    only with the rules engine): those keep the paint, so the chromosome is not split, and are to
    be looked at."""
    if not classes:
        return []
    bit = 1 << st.tag_id("Contaminant")
    listed = np.zeros(len(pt.names), dtype=bool)
    for i, n in enumerate(pt.names):
        listed[i] = n in classes
    ms = MapState(pt).replay(st.edits)
    sc_of = np.repeat(np.arange(len(pt.names)), pt.scaffold_npix.astype(np.int64))
    is_listed = [bool(listed[sc_of[ms.pix[c.start]]]) for c in ms.contigs]
    flags = dict(st.meta_flags)
    for i, x in enumerate(is_listed):
        if x:
            flags[i] = bit
    st.meta_flags = sorted(flags.items())
    # painted scaffolds (scaffs: contig index -> scaffold id): unpaint the non-target contigs,
    # except inside a scaffold that also holds target sequence
    sid = dict(st.scaffs)
    by = {}
    for i, k in sid.items():
        if k:
            by.setdefault(k, []).append(i)
    mixed = []
    for k, idx in by.items():
        if any(is_listed[i] for i in idx) and not all(is_listed[i] for i in idx):
            mixed += [i for i in idx if is_listed[i]]
    keep = set(mixed)
    st.scaffs = [(i, k) for i, k in st.scaffs if not is_listed[i] or i in keep]
    return sorted(mixed)


def write_table(classes, sources, lengths, path):
    """`<prefix>.contaminants.tsv`: scaffold, length (Mb), class, source."""
    with open(path, "w") as f:
        f.write("scaffold\tlength_mb\tclass\tsource\n")
        for n in sorted(classes, key=lambda n: -lengths.get(n, 0)):
            f.write(f"{n}\t{lengths.get(n, 0) / 1e6:.3f}\t{classes[n]}\t{sources[n]}\n")


ORGANELLE = {"mitochondrion": "mito", "plastid": "plastid"}


def split_organelles(asm_dir, classes):
    """After pretext-to-asm: move the organelle sequences out of its contaminants FASTA (and AGP)
    into `<stem>.mito.fa` / `<stem>.plastid.fa` (+ .agp), by the class of the scaffolds each output
    sequence is made of (its AGP). Returns {file written: number of sequences}."""
    import gzip
    out = {}
    if not os.path.isdir(asm_dir):
        return out
    for fn in sorted(os.listdir(asm_dir)):
        if "contaminant" not in fn.lower() or not fn.endswith(".agp"):
            continue
        agp = os.path.join(asm_dir, fn)
        stem = agp[:-len(".agp")]                         # <...>.fa(.gz)
        fa = next((x for x in (stem, stem + ".gz") if os.path.exists(x)), None)   # x.fa.agp: x.fa or x.fa.gz
        if fa is None:
            continue
        comps, lines = {}, []
        with open(agp) as f:
            for line in f:
                lines.append(line)
                c = line.split("\t")
                if not line.startswith("#") and len(c) > 5 and c[4] not in ("N", "U"):
                    comps.setdefault(c[0], set()).add(c[5])
        kind = {}
        for obj, cs in comps.items():
            ks = {ORGANELLE.get(classes.get(x, "")) for x in cs}
            if len(ks) == 1 and None not in ks:
                kind[obj] = ks.pop()
        if not kind:
            continue
        gz = fa.endswith(".gz")
        op = gzip.open if gz else open
        base = fa[:-3] if gz else fa                      # <...>.fa
        root = base[:-3] if base.endswith(".fa") else base
        files = {k: root + f".{k}.fa" + (".gz" if gz else "") for k in set(kind.values())}
        handles = {k: op(p, "wt") for k, p in files.items()}
        rest = op(fa + ".tmp", "wt")
        n = dict.fromkeys(files, 0)
        w = rest
        with op(fa, "rt") as f:
            for line in f:
                if line.startswith(">"):
                    k = kind.get(line[1:].split()[0])
                    w = handles[k] if k else rest
                    if k:
                        n[k] += 1
                w.write(line)
        for h in list(handles.values()) + [rest]:
            h.close()
        os.replace(fa + ".tmp", fa)
        for k, p in files.items():
            with open(root + f".{k}.fa.agp", "w") as w:
                w.writelines(l for l in lines if l.startswith("#") or kind.get(l.split("\t")[0]) == k)
            out[p] = n[k]
        with open(agp, "w") as w:
            w.writelines(l for l in lines if l.startswith("#") or l.split("\t")[0] not in kind)
    return out


def confirm(pt, classes, sources, big_bp=10e6, blocks=None, ratio=2.0, mito_cov=5.0, chrom_frac=0.2,
            max_side=32768, gap=1, window_bp=1e6):
    """Which Tiara calls to act on. Tiara alone also calls nuclear sequence non-eukaryotic: repeat-
    rich pieces prokaryotic (qqTroRaco), a high-copy AT-rich repeat "organelle" and NUMT-like pieces
    "mitochondrion" (icLarUrsu: GC 0.27 at 4-10x, and 0.5-3x; the curated assembly keeps them all).
    So each call needs other evidence; listed names (--contaminants) are taken as given.
      - organelle of unknown kind: never tagged on Tiara alone.
      - mitochondrion: organelle copy number, median coverage >= `mito_cov` x the haploid level.
      - prokaryote or plastid: Hi-C sets it apart. The pieces of one foreign genome contact each
        other, a nuclear piece its chromosome: its second-best contact with another called scaffold
        (per bin pair, bins within `gap` of its own left out; one called neighbour in map order can
        share contact, a genome's pieces contact several) is >= `ratio` x its peak contact with a
        scaffold over `big_bp` (target sequence, never classified): the best `window_bp` window along
        it, not the mean over it, so a piece at a chromosome end or inside one is caught. Contact
        among the calls only supports them: a group that contacts itself and nothing else is kept.
      - any: with `blocks` (alignments), not aligned to a scaffold over `big_bp` along >=
        `chrom_frac` of its length (then it is nuclear sequence).
    Scaffolds under one pixel cannot be judged by the map.

    Returns (kept {name: class}, unconfirmed [(name, class, partner contact, chromosome contact,
    the reason)]). qqTroRaco (65,536 px, read at 32,768): two foreign genomes (GC 0.22, 2.6 Mb,
    single copy; GC 0.35 against the genome's 0.315, 2.2 Mb, aligned only to each other), second
    partner 125-255. The GC 0.22 group peaks at 150-250 on the first 2 Mb of H1_scaffold_45 (the mean
    over the chromosome is 33-66; between chromosomes ~1): nuclear sequence at that chromosome's start
    (the curator's Unlocs), not a genome of its own; the GC 0.35 group peaks at 46-74, scattered.
    H2_scaffold_100 contacts one called neighbour (197) and no other (0)."""
    nan = float("nan")
    test = [n for n in classes if sources.get(n) == "Tiara"]
    if not test:
        return dict(classes), []
    idx = {n: i for i, n in enumerate(pt.names)}
    L = pt.scaffold_len_bp()
    on_chrom = {}
    for b in blocks or ():
        if b.q in classes and b.q != b.t and b.t in idx and L[idx[b.t]] > big_bp:
            on_chrom.setdefault(b.q, {}).setdefault(b.t, 0)
            on_chrom[b.q][b.t] += b.qe - b.qs
    from .allhaps.common import _signal
    covr = _signal(pt).coverage_ratio() if "coverage" in pt.graphs else None
    hic = [n for n in test if classes[n] in ("bacteria", "archaea", "prokarya", "plastid")]
    M = None
    bins = {}
    lvl = 0
    while (pt.npix >> lvl) > max_side:
        lvl += 1
    fct = 1 << lvl
    for i, n in enumerate(pt.names):
        a, k = int(pt.scaffold_start[i]), int(pt.scaffold_npix[i])
        if k:
            bins[n] = np.arange(a // fct, (a + k - 1) // fct + 1)
    if any(n in bins for n in hic):
        M = pt.matrix(lvl)
    big = [n for n in pt.names if L[idx[n]] > big_bp and n in bins]

    def contact(a, b):
        ba, bb = bins[a], np.setdiff1d(bins[b], bins[a])   # bins a small neighbour shares are not contact
        if not len(bb):
            return nan
        far = np.abs(ba[:, None] - bb[None, :]) > gap
        return float(M[np.ix_(ba, bb)][far].astype(np.float32).mean()) if far.any() else nan
    win = max(1, int(round(window_bp / (fct * pt.bp_per_pix))))

    def peak(a, b):
        """Best `window_bp` window of a's contact along b (bins next to a's own left out)."""
        ba, bb = bins[a], np.setdiff1d(bins[b], bins[a])
        if not len(bb):
            return nan
        sub = M[np.ix_(ba, bb)].astype(np.float32)
        sub[np.abs(ba[:, None] - bb[None, :]) <= gap] = np.nan
        if not np.isfinite(sub).any():
            return nan
        prof = np.nan_to_num(np.nanmean(sub, axis=0))
        w = min(win, len(prof))
        return float(np.convolve(prof, np.ones(w) / w, mode="valid").max())
    called = [n for n in classes if n in bins]
    kept, unconfirmed = {}, []
    for n, c in classes.items():
        if n not in test:
            kept[n] = c
            continue
        ln = L[idx[n]] if n in idx else 0
        if n in on_chrom:
            t, bp = max(on_chrom[n].items(), key=lambda kv: kv[1])
            if bp >= chrom_frac * ln:
                unconfirmed.append((n, c, nan, nan, f"aligns to {t} ({bp / 1e3:.0f} kb): nuclear sequence"))
                continue
        if c == "organelle":
            unconfirmed.append((n, c, nan, nan, "organelle of unknown kind: Tiara alone is not enough"))
            continue
        if n not in bins:
            unconfirmed.append((n, c, nan, nan, "under one pixel of the map: no coverage or Hi-C to check"))
            continue
        if c == "mitochondrion":
            i = idx[n]
            cv = (float(np.median(covr[int(pt.scaffold_start[i]):int(pt.scaffold_start[i] + pt.scaffold_npix[i])]))
                  if covr is not None else nan)
            if cv >= mito_cov:
                kept[n] = c
            else:
                unconfirmed.append((n, c, nan, nan, f"coverage {cv:.1f}x, not organelle copy number (>= {mito_cov:g}x)"))
            continue
        ps = sorted((x for x in (contact(n, m) for m in called if m != n) if x == x), reverse=True)
        partner = ps[1] if len(ps) > 1 else 0.0
        chrom, where = max(((x, b) for x, b in ((peak(n, b), b) for b in big) if x == x), default=(0.0, ""))
        if partner > 0 and partner >= ratio * chrom:
            kept[n] = c
        else:
            unconfirmed.append((n, c, partner, chrom, f"Hi-C ties it to {where or 'a chromosome'}"))
    return kept, unconfirmed
