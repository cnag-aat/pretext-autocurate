"""minimap2 of an assembly (query) against a reference, with the query cut into chunks and
the PAF shifted back onto whole query sequences. Like D-GENIES' split_before: whole-chromosome
queries on a repeat-rich genome can make asm5 chaining take tens of GB.

    python3 tools/align_chunked.py REF.fa(.gz) QUERY.fa OUT.paf [--chunk-mb 10] [-t 4]
                                   [--preset asm5] [--max-rss-gb 12]

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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ref")
    ap.add_argument("query")
    ap.add_argument("out")
    ap.add_argument("--chunk-mb", type=float, default=10)
    ap.add_argument("-t", "--threads", type=int, default=4)
    ap.add_argument("--preset", default="asm5")
    ap.add_argument("--max-rss-gb", type=float, default=12)
    ap.add_argument("--minimap2", default="minimap2")
    ap.add_argument("--mm-args", default="-K 20M --cap-kalloc 1g -f 0.002",
                    help="extra minimap2 options (small batches, capped per-thread memory, more repeat masking)")
    ap.add_argument("--query-regex", help="align only query sequences whose name matches (e.g. ^H2_)")
    a = ap.parse_args()
    size = int(a.chunk_mb * 1e6)
    qlen = lengths(a.query)
    cmd = [a.minimap2, "-x", a.preset, "-t", str(a.threads), "--secondary=no", *shlex.split(a.mm_args), a.ref, "-"]
    p = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    state = {"killed": False, "peak": 0.0}

    def feed():
        try:
            for name, off, s in chunks(a.query, size, re.compile(a.query_regex) if a.query_regex else None):
                p.stdin.write(f">{name}::{off}\n")
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
            name, off = c[0].rsplit("::", 1)
            off = int(off)
            c[0], c[1] = name, str(qlen[name])
            c[2], c[3] = str(int(c[2]) + off), str(int(c[3]) + off)
            out.write("\t".join(c) + "\n")
            n += 1
    rc = p.wait()
    print(f"{n} alignments, {time.time() - t0:.0f} s, minimap2 peak RSS {state['peak']:.1f} GB", file=sys.stderr)
    if state["killed"]:
        sys.exit(f"minimap2 killed: resident memory passed {a.max_rss_gb} GB")
    if rc != 0:
        sys.exit(f"minimap2 exited with {rc}")
    os.replace(tmp, a.out)


if __name__ == "__main__":
    main()
