"""Benchmark minimap2 against FASTGA for the haplotype alignment (H2 against H1).

Usage: python3 tools/bench_aligners.py ASSEMBLY.fa WORKDIR [THREADS] [EXISTING_MINIMAP2.paf]

Splits the assembly by haplotype prefix (kept in WORKDIR/split for Red and later runs),
runs each aligner under /usr/bin/time -l (wall time, peak memory), writes
WORKDIR/<aligner>.paf and prints a summary: blocks, aligned bases, bases in blocks of
10 kb or more and with mapq >= 30 (what the tool uses).
"""
import os
import re
import subprocess
import sys

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate.homology import split_by_haplotype, find_tool  # noqa: E402

fa, work = sys.argv[1], sys.argv[2]
threads = int(sys.argv[3]) if len(sys.argv) > 3 else 8
os.makedirs(work, exist_ok=True)
parts = split_by_haplotype(fa, os.path.join(work, "split"))
h1, h2 = parts["H1"], parts["H2"]


def timed(name, cmd, out):
    log = os.path.join(work, name + ".time")
    env = dict(os.environ, PATH=os.path.dirname(cmd[0]) + os.pathsep + os.environ.get("PATH", ""))
    with open(out, "w") as fo, open(log, "w") as fe:
        subprocess.run(["/usr/bin/time", "-l"] + cmd, stdout=fo, stderr=fe, check=True, env=env)


def timing(name):
    """(wall seconds, peak GB) from the /usr/bin/time -l log (macOS reports bytes)."""
    wall = rss = float("nan")
    p = os.path.join(work, name + ".time")
    if os.path.exists(p):
        for line in open(p):
            m = re.match(r"\s*([\d.]+) real", line)
            if m:
                wall = float(m.group(1))
            m = re.match(r"\s*(\d+)\s+maximum resident set size", line)
            if m:
                rss = int(m.group(1)) / 2 ** 30
    return wall, rss


def summary(paf):
    n = bp = big = good = 0
    for line in open(paf):
        c = line.split("\t")
        if len(c) < 12:
            continue
        L = int(c[3]) - int(c[2])
        n += 1
        bp += L
        if int(c[10]) >= 10000:
            big += L
            if int(c[11]) >= 30:
                good += L
    return n, bp, big, good


mm = find_tool("minimap2")
paf = os.path.join(work, "minimap2.asm5.paf")
if len(sys.argv) > 4 and os.path.exists(sys.argv[4]):       # an existing minimap2 PAF (timing not rerun)
    paf = sys.argv[4]
if not (os.path.exists(paf) and os.path.getsize(paf)):
    timed("minimap2", [mm, "-x", "asm5", "-t", str(threads), "--secondary=no", h1, h2], paf)
fg = find_tool("FastGA")
paf2 = os.path.join(work, "fastga.paf")
if not (os.path.exists(paf2) and os.path.getsize(paf2)):
    tmp = os.path.join(work, "fastga_tmp")
    os.makedirs(tmp, exist_ok=True)
    timed("fastga", [fg, f"-T{threads}", f"-P{tmp}", "-paf", h2, h1], paf2)
# chained FastGA (as hapscotch does for its plots): .1aln, ALNchain, ALNtoPAF
paf3 = os.path.join(work, "fastga.chain.paf")
if not (os.path.exists(paf3) and os.path.getsize(paf3)):
    tmp = os.path.join(work, "fastga_tmp")
    os.makedirs(tmp, exist_ok=True)
    one, ch = os.path.join(work, "fastga.1aln"), os.path.join(work, "fastga.chain.1aln")
    timed("fastga1aln", [fg, f"-T{threads}", f"-P{tmp}", f"-1:{one}", h2, h1], os.devnull)
    timed("alnchain", [find_tool("ALNchain"), f"-o{ch}", one], os.devnull)
    with open(paf3, "w") as fo:
        subprocess.run([find_tool("ALNtoPAF"), f"-T{threads}", ch], stdout=fo, check=True)
for name, p, key in (("minimap2 asm5", paf, "minimap2"), ("FastGA", paf2, "fastga"),
                     ("FastGA chained", paf3, "fastga1aln")):
    n, bp, big, good = summary(p)
    wall, gb = timing(key)
    print(f"{name}\twall {wall / 60:.1f} min\tpeak {gb:.1f} GB\tblocks {n}\taligned {bp / 1e6:.0f} Mb\t"
          f">=10 kb {big / 1e6:.0f} Mb\t>=10 kb and mapq>=30 {good / 1e6:.0f} Mb")
