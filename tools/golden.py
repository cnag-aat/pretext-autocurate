"""Golden outputs for refactoring: record what the tool writes now, check that it still does.

Usage:
    python3 tools/golden.py record [--out golden_out]
    python3 tools/golden.py check  [--out golden_out]

Runs the tool on every case of tools/score_cases.json (levels 0, 0.5, 1 for two haplotypes)
and on the synthetic toy maps of tests/test_synthetic.py, in parallel, and hashes what each
run writes that should not change: the decisions and review tables, the AGP and every
savestate (not run.txt, which holds the date and command line). `record` stores the hashes
in OUT/golden.json; `check` reruns and lists every file that differs (exit status 1 if any).
Cases whose inputs are missing are skipped (the synthetic ones always run).
"""
import argparse
import concurrent.futures
import hashlib
import json
import os
import shutil
import subprocess
import sys

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))

LEVELS = "0,0.5,1"
KEEP = (".decisions.tsv", ".review.tsv", ".agp")


def fill(s, F):
    return s.replace("{F}", F).replace("{R}", R)


def synthetic_maps(d):
    """The toy maps of tests/test_synthetic.py: [(name, run arguments)]."""
    sys.path.insert(0, os.path.join(R, "tests"))
    from pretext_autocurate.synth import toy_map
    import test_synthetic as T
    os.makedirs(d, exist_ok=True)
    out = []
    for name, sc in (("toy_diploid", T.DIPLOID), ("toy_single", T.SINGLE)):
        p = os.path.join(d, name + ".pretext")
        toy_map(p, sc, T.CHROMS)
        out.append((name, ["--map", p], name == "toy_single"))
    return out


def run_one(name, args, single, out):
    d = os.path.join(out, "runs", name)
    shutil.rmtree(d, ignore_errors=True)
    cmd = [sys.executable, "-m", "pretext_autocurate.cli", "run", *args, "--out-dir", d]
    if not single:
        cmd += ["--sweep", LEVELS]
    env = dict(os.environ, PYTHONPATH=os.path.join(R, "src"))
    r = subprocess.run(cmd, env=env, capture_output=True, text=True)
    if r.returncode != 0:
        return name, None, (r.stderr or r.stdout)[-2000:]
    hashes = {}
    for root, _, files in os.walk(d):
        for f in sorted(files):
            if f.endswith(KEEP) or ".savestate_autocurate" in f:
                p = os.path.join(root, f)
                with open(p, "rb") as fh:
                    hashes[os.path.relpath(p, d)] = hashlib.sha256(fh.read()).hexdigest()
    return name, hashes, ""


def run_all(out, fixtures):
    F = os.path.abspath(fixtures)
    jobs = [(n, a, s) for n, a, s in synthetic_maps(os.path.join(out, "maps"))]
    for case in json.load(open(os.path.join(R, "tools", "score_cases.json")))["cases"]:
        args = [fill(x, F) for x in case["run"]]
        if any((x.startswith("/") or "=" in x) and not os.path.exists(x.split("=", 1)[-1]) for x in args):
            print(f"{case['name']}: inputs missing, skipped", file=sys.stderr)
            continue
        jobs.append((case["name"], args, case["mode"] == "single"))
    res = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as ex:
        for name, hashes, err in ex.map(lambda j: run_one(*j, out), jobs):
            if hashes is None:
                sys.exit(f"{name}: run failed:\n{err}")
            res[name] = hashes
            print(f"{name}: {len(hashes)} files", file=sys.stderr)
    return res


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("what", choices=["record", "check"])
    ap.add_argument("--out", default=os.path.join(R, "golden_out"))
    ap.add_argument("--fixtures", default=os.path.join(R, "seed", "fixtures"))
    a = ap.parse_args()
    gold = os.path.join(a.out, "golden.json")
    if a.what == "check" and not os.path.exists(gold):
        sys.exit(f"no {gold}: run record first")
    res = run_all(a.out, a.fixtures)
    if a.what == "record":
        with open(gold, "w") as f:
            json.dump(res, f, indent=1, sort_keys=True)
        print(f"recorded {sum(len(v) for v in res.values())} files in {gold}")
        return
    old = json.load(open(gold))
    diff = []
    for name, files in old.items():
        if name not in res:
            continue
        for f, h in files.items():
            if res[name].get(f) != h:
                diff.append(f"{name}/{f}: {'missing' if f not in res[name] else 'differs'}")
        diff += [f"{name}/{f}: new" for f in res[name] if f not in files]
    if diff:
        print("DIFFERENT:\n  " + "\n  ".join(diff))
        sys.exit(1)
    print(f"identical: {sum(len(v) for v in res.values())} files in {len(res)} runs")


if __name__ == "__main__":
    main()
