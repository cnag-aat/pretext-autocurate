"""minimap2 of an assembly (query) against a reference, with the query cut into chunks and
the PAF shifted back onto whole query sequences. Like D-GENIES' split_before: whole-chromosome
queries on a repeat-rich genome can make asm5 chaining take tens of GB.

    python3 tools/align_chunked.py REF.fa(.gz) QUERY.fa OUT.paf [--chunk-mb 10] [-t 4]
                                   [--preset asm5] [--max-rss-gb 12] [--no-chunk] [--secondary]

All-vs-all (the assembly against itself, each sequence's own diagonal skipped):
    python3 tools/align_chunked.py ASM.fa ASM.fa OUT.paf --no-chunk --secondary \
        --mm-args "-D -N 5 -p 0.5 -K 20M --cap-kalloc 1g"
minimap2's -D only skips the diagonal of a query with the target's own name and length, so it
needs whole sequences (--no-chunk), and the other copies of a region are secondary alignments.
Whole chromosomes of a large repeat-rich genome can pass the memory cap (qqTroRaco, icLarUrsu
passed 16 GB); then keep the chunks and drop each chunk's own diagonal after the fact:
    python3 tools/align_chunked.py ASM.fa ASM.fa OUT.paf --secondary --drop-self \
        --mm-args "-N 5 -p 0.5 -K 20M --cap-kalloc 1g -f 0.002"
(a chunk's best hit is then itself, so a copy is kept only as a secondary at >= half its score;
minimap2 gives secondaries MAPQ 0, so MAPQ is worked out again without the self hits, from the
chain scores as minimap2 does; `--requalify PAF --chunk-mb 10` does that for an existing file).

The chunks are streamed to minimap2 on stdin (no copy on disk). A watchdog kills minimap2 if
its resident memory passes --max-rss-gb (the machine has 24 GB; swap on a full disk crashed it).
Alignments that cross a chunk boundary come out as two blocks."""
import argparse
import gzip
import os
import re
import shlex
import subprocess
import sys
import threading
import time


def chunks(path, size, keep=None):
    """(name, offset, sequence) per chunk of every sequence (whose name matches `keep`)."""
    op = gzip.open if path.endswith(".gz") else open
    name, buf, off, nbuf = None, [], 0, 0

    def flush(final):
        nonlocal buf, off, nbuf
        s = "".join(buf)
        out = []
        while len(s) >= size or (final and s):
            out.append((name, off, s[:size]))
            off += len(s[:size])
            s = s[size:]
        buf, nbuf = ([s] if s else []), len(s)
        return out

    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                if name is not None and (keep is None or keep.search(name)):
                    yield from flush(True)
                name, buf, off, nbuf = line[1:].split()[0], [], 0, 0
            elif keep is not None and not keep.search(name):
                continue
            else:
                s = line.strip()
                buf.append(s)
                nbuf += len(s)
                if nbuf >= 4 * size:
                    yield from flush(False)
        if name is not None and (keep is None or keep.search(name)):
            yield from flush(True)


def lengths(path):
    op = gzip.open if path.endswith(".gz") else open
    out, name = {}, None
    if os.path.exists(path + ".fai"):
        with open(path + ".fai") as f:
            return {c[0]: int(c[1]) for c in (l.split("\t") for l in f)}
    with op(path, "rt") as f:
        for line in f:
            if line.startswith(">"):
                name = line[1:].split()[0]
                out[name] = 0
            else:
                out[name] += len(line.strip())
    return out


def _mapq(group):
    """MAPQ again for the alignments of one query chunk once its self hit is gone: minimap2 gives
    every secondary 0, and with the chunk itself as primary the other haplotype's copy is one.
    As minimap2 does: an alignment that an overlapping one (half the shorter) outscores (s1, the
    chain score) is secondary (0); the best gets 40 (1 - s2/s1) min(1, anchors/10) ln(s1), at most 60,
    s2 the best overlapping other chain."""
    import math
    rows = []
    for c in group:
        tags = {t[:2]: t[5:] for t in c[12:]}
        rows.append((int(c[2]), int(c[3]), int(tags.get("s1", 0)), int(tags.get("cm", 0))))
    for k, (qs, qe, s1, cm) in enumerate(rows):
        s2, beaten = 0, False
        for j, (bs, be, t1, _) in enumerate(rows):
            if j == k or min(qe, be) - max(qs, bs) < 0.5 * min(qe - qs, be - bs):
                continue
            if t1 > s1 or (t1 == s1 and j < k):
                beaten = True
                break
            s2 = max(s2, t1)
        if beaten or s1 <= 0:
            q = 0
        else:
            q = int(40 * (1 - s2 / s1) * min(1.0, cm / 10) * math.log(s1))
        group[k][11] = str(max(0, min(60, q)))
    return group


def requalify(path, chunk_bp, out=None):
    """`_mapq` over a PAF written with --drop-self, grouped by query chunk (chunks start at multiples
    of chunk_bp, and an alignment starts in its own chunk). Rewrites `path` (or writes `out`)."""
    groups = {}
    with open(path) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            groups.setdefault((c[0], int(c[2]) // chunk_bp), []).append(c)
    tmp = (out or path) + ".tmp"
    n = 0
    with open(tmp, "w") as w:
        for g in groups.values():
            for c in _mapq(g):
                w.write("\t".join(c) + "\n")
                n += int(c[11]) >= 30
    os.replace(tmp, out or path)
    return n


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref", nargs="?")
    ap.add_argument("query", nargs="?")
    ap.add_argument("out", nargs="?")
    ap.add_argument("--requalify", metavar="PAF",
                    help="only redo the MAPQ of an existing --drop-self PAF (with its --chunk-mb), in place")
    ap.add_argument("--chunk-mb", type=float, default=10)
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--preset", default="asm5")
    ap.add_argument("--max-rss-gb", type=float, default=12)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--mm-args", default="-K 20M --cap-kalloc 1g -f 0.002",
                    help="extra minimap2 options (small batches, capped per-thread memory, more repeat masking)")
    ap.add_argument("--query-regex", help="align only query sequences whose name matches (e.g. ^H2_)")
    ap.add_argument("--no-chunk", action="store_true",
                    help="stream whole sequences under their own names (for -D in all-vs-all); the watchdog stays")
    ap.add_argument("--secondary", action="store_true", help="keep secondary alignments (default: --secondary=no)")
    ap.add_argument("--drop-self", action="store_true",
                    help="all-vs-all with chunks: drop each chunk's alignment to its own place (what -D does "
                         "for whole sequences)")
    a = ap.parse_args()
    if a.requalify:
        n = requalify(a.requalify, int(a.chunk_mb * 1e6))
        print(f"{n} alignments at MAPQ >= 30 in {a.requalify}", file=sys.stderr)
        return
    if not (a.ref and a.query and a.out):
        ap.error("REF QUERY OUT are needed (or --requalify PAF)")
    size = int(a.chunk_mb * 1e6) if not a.no_chunk else 1 << 40
    qlen = lengths(a.query)
    cmd = [a.minimap2, "-x", a.preset, "-t", str(a.threads), *([] if a.secondary else ["--secondary=no"]), *shlex.split(a.mm_args), a.ref, "-"]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    state = {"killed": False, "peak": 0.0}

    def feed():
        try:
            for name, off, s in chunks(a.query, size, re.compile(a.query_regex) if a.query_regex else None):
                p.stdin.write(f">{name}\n" if a.no_chunk else f">{name}::{off}\n")
                for i in range(0, len(s), 100000):
                    p.stdin.write(s[i:i + 100000] + "\n")
        except BrokenPipeError:
            pass
        finally:
            try:
                p.stdin.close()
            except BrokenPipeError:
                pass

    def watch():
        while p.poll() is None:
            # top's MEM counts compressed and swapped pages; ps RSS does not (a 73 GB job showed 1.6)
            r = subprocess.run(["top", "-l", "1", "-pid", str(p.pid), "-stats", "mem"], capture_output=True,
                               text=True).stdout.strip().splitlines()
            v = r[-1].strip().rstrip("+-") if r else ""
            if v and v[-1] in "KMGT" and v[:-1].replace(".", "").isdigit():
                gb = float(v[:-1]) * {"K": 1 / 1048576, "M": 1 / 1024, "G": 1, "T": 1024}[v[-1]]
                state["peak"] = max(state["peak"], gb)
                if gb > a.max_rss_gb:
                    state["killed"] = True
                    p.kill()
                    return
            time.sleep(5)

    threading.Thread(target=feed, daemon=True).start()
    threading.Thread(target=watch, daemon=True).start()
    t0, n = time.time(), 0
    tmp = a.out + ".tmp"
    with open(tmp, "w") as out:
        for line in p.stdout:
            c = line.rstrip("\n").split("\t")
            name, off = (c[0], 0) if a.no_chunk else c[0].rsplit("::", 1)
            off = int(off)
            c[0], c[1] = name, str(qlen[name])
            c[2], c[3] = str(int(c[2]) + off), str(int(c[3]) + off)
            if (a.drop_self and c[0] == c[5] and c[4] == "+" and abs(int(c[2]) - int(c[7])) <= 1000
                    and abs(int(c[3]) - int(c[8])) <= 1000):
                continue                          # a chunk on its own place in its own sequence
            out.write("\t".join(c) + "\n")
            n += 1
    rc = p.wait()
    print(f"{n} alignments, {time.time() - t0:.0f} s, minimap2 peak RSS {state['peak']:.1f} GB", file=sys.stderr)
    if state["killed"]:
        sys.exit(f"minimap2 killed: resident memory passed {a.max_rss_gb} GB")
    if rc != 0:
        sys.exit(f"minimap2 exited with {rc}")
    os.replace(tmp, a.out)
    if a.drop_self and a.secondary:
        n = requalify(a.out, size)
        print(f"MAPQ redone without the self hits: {n} alignments at MAPQ >= 30", file=sys.stderr)


if __name__ == "__main__":
    main()
