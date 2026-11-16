"""Run pretext-autocurate on every test genome and score it against the curators' work.

Usage:
    python3 tools/score_all.py [--cases tools/score_cases.json] [--fixtures seed/fixtures]
                               [--out score_out] [--levels 0,0.5,1] [--only NAME ...] [--rerun]

For each case in the cases file (a JSON list: name, mode, the `run` arguments, and the
reference savestates), the tool is run once with `--sweep LEVELS` (two-haplotype cases) or
once (single-haplotype cases) into OUT/<name>/, unless that output is already there (then it
is reused: --rerun to run again). Each savestate is scored against each reference with
`pretext_autocurate.score.compare`. Writes OUT/score.tsv (every metric) and prints a summary
table (also OUT/score.md), to paste into docs/PROGRESS.md.

To add a genome: put its map(s), inputs and curated savestate in the fixtures directory and
add an entry to the cases file.
"""
import argparse
import json
import math
import os
import subprocess
import sys
import time

R = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(R, "src"))
from pretext_autocurate import Pretext, SaveState            # noqa: E402
from pretext_autocurate.score import compare                # noqa: E402


def fill(s, F):
    return s.replace("{F}", F).replace("{R}", R)


def tool_states(case, out, levels):
    """[(level or '', savestate path)] written by the run in `out`."""
    base = os.path.basename(fill(case["run"][case["run"].index("--map") + 1], ""))
    if case["mode"] == "single":
        p = os.path.join(out, base + ".savestate_autocurate")
        return [("", p)] if os.path.exists(p) else []
    res = []
    for lv in levels:
        p = os.path.join(out, f"level_{lv:.2f}", base + f".savestate_autocurate_level{lv:.2f}")
        if os.path.exists(p):
            res.append((f"{lv:.2f}", p))
    return res


def review_count(path):
    try:
        with open(path) as f:
            return sum(1 for _ in f) - 1
    except OSError:
        return ""


def fmt(v, pct=False):
    if v == "" or v is None or (isinstance(v, float) and math.isnan(v)):
        return "-"
    return f"{100 * v:.0f}%" if pct else str(v)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--cases", default=os.path.join(R, "tools", "score_cases.json"))
    ap.add_argument("--fixtures", default=os.path.join(R, "seed", "fixtures"))
    ap.add_argument("--out", default=os.path.join(R, "score_out"))
    ap.add_argument("--levels", default="0,0.5,1")
    ap.add_argument("--only", nargs="*", help="case names to run (default: all)")
    ap.add_argument("--rerun", action="store_true", help="run the tool again even if its output is there")
    a = ap.parse_args()
    F = os.path.abspath(a.fixtures)
    levels = [float(x) for x in a.levels.split(",")]
    cases = json.load(open(a.cases))["cases"]
    rows = []
    for case in cases:
        if a.only and case["name"] not in a.only:
            continue
        args = [fill(x, F) for x in case["run"]]
        missing = [x for x in args if (x.startswith("/") or "=" in x) and not os.path.exists(x.split("=", 1)[-1])]
        if missing:
            print(f"{case['name']}: skipped, missing {missing[0]}", file=sys.stderr)
            continue
        out = os.path.join(a.out, case["name"])
        states = tool_states(case, out, levels)
        took = ""
        if a.rerun or not states:
            cmd = [sys.executable, "-m", "pretext_autocurate.cli", "run", *args, "--out-dir", out]
            if case["mode"] != "single":
                cmd += ["--sweep", ",".join(f"{x:g}" for x in levels)]
            print(f"{case['name']}: running ...", file=sys.stderr)
            t = time.time()
            env = dict(os.environ, PYTHONPATH=os.path.join(R, "src"))
            r = subprocess.run(cmd, env=env, capture_output=True, text=True)
            took = round((time.time() - t) / 60, 1)
            with open(os.path.join(a.out, f"{case['name']}.log"), "w") as f:
                f.write(r.stdout + r.stderr)
            if r.returncode != 0:
                print(f"{case['name']}: run failed (see {a.out}/{case['name']}.log)", file=sys.stderr)
                continue
            states = tool_states(case, out, levels)
        pt = Pretext(args[args.index("--map") + 1], load_graphs=False)
        for ref_name, ref_path in case["references"].items():
            ref_path = fill(ref_path, F)
            if not os.path.exists(ref_path):
                continue
            ref = SaveState.load(ref_path)
            for lv, path in states:
                m = compare(pt, SaveState.load(path), ref)
                rev = review_count(os.path.join(os.path.dirname(path), "autocurate.review.tsv"))
                rows.append(dict(case=case["name"], level=lv, reference=ref_name, review_items=rev,
                                 run_minutes=took, **m))
                print(f"{case['name']} {lv or '-'} vs {ref_name}: junctions {fmt(m['junction_precision'], 1)} / "
                      f"{fmt(m['junction_recall'], 1)}", file=sys.stderr)
    if not rows:
        sys.exit("nothing scored")
    os.makedirs(a.out, exist_ok=True)
    cols = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(a.out, "score.tsv"), "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(round(r[c], 4)) if isinstance(r.get(c), float) else str(r.get(c, ""))
                              for c in cols) + "\n")
    # summary: the numbers that move with the rules
    head = ["genome", "level", "reference", "junctions P / R", "breaks P / R", "chromosomes >=90% in one",
            "HAP1 agree", "sex tags agree", "Unloc agree", "edits tool / ref", "review"]
    md = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]

    def agree(r, tags):
        both = sum(r.get(f"tag_{t}_both_mb", 0) for t in tags)
        ref = sum(r.get(f"tag_{t}_ref_mb", 0) for t in tags)
        return f"{100 * both / ref:.0f}%" if ref else "-"
    for r in rows:
        md.append("| " + " | ".join([
            r["case"], r["level"] or "-", r["reference"],
            f"{r['junctions_supported']}/{r['junctions_tool']} / {r['junctions_reproduced']}/{r['junctions_ref']}",
            f"{fmt(r['break_precision'], 1)} / {fmt(r['break_recall'], 1)}",
            f"{fmt(r['grouped_90'], 1)} of {r['ref_chromosomes']}",
            agree(r, ["HAP1"]), agree(r, ["X", "Y", "Z", "W"]), agree(r, ["Unloc"]),
            f"{r['edits_tool']} / {r['edits_ref']}", str(r["review_items"])]) + " |")
    text = "\n".join(md) + "\n"
    with open(os.path.join(a.out, "score.md"), "w") as f:
        f.write(text)
    print(text)
    print(f"all metrics: {os.path.join(a.out, 'score.tsv')}", file=sys.stderr)


if __name__ == "__main__":
    main()
