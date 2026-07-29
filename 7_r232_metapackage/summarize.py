#!/usr/bin/env python3
"""Summarise benchmark 7 (r232) results: runtime / RAM (from the snakemake
benchmark TSVs) and accuracy (F1 + Bray-Curtis from the OPAL reports), per tool.

singlem-regime3's wall time and peak RAM are the sum / max across its three stages
(pipe, weebill, condense), since those run as separate rules/jobs.

Usage: pixi run -e art python3 7_r232_metapackage/summarize.py   (run from repo root)
"""
import os
import glob
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
SAMPLE = "known1000"
RANK = "species"

# tool -> list of benchmark-dir stems that make up its runtime/RAM
TOOL_BENCHMARKS = {
    "singlem": ["singlem"],
    "sylph": ["sylph"],
    "singlem-regime3": [
        "singlem-regime3-pipe",
        "singlem-regime3-weebill",
        "singlem-regime3-condense",
    ],
}


def read_benchmark(stem):
    path = os.path.join(HERE, "benchmarks", stem, f"{SAMPLE}-8threads.benchmark")
    if not os.path.exists(path):
        return None
    with open(path) as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    r = rows[0]
    return {"s": float(r["s"]), "max_rss_mb": float(r["max_rss"])}


def resource_summary(tool):
    stems = [read_benchmark(s) for s in TOOL_BENCHMARKS[tool]]
    stems = [s for s in stems if s is not None]
    if not stems:
        return None, None
    wall_s = sum(s["s"] for s in stems)            # summed across stages
    peak_rss = max(s["max_rss_mb"] for s in stems)  # peak of any stage
    return wall_s, peak_rss


def opal_metrics(tool):
    path = os.path.join(HERE, f"output_{tool}", "opal", f"{SAMPLE}.opal_report")
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as f:
        for row in csv.reader(f, delimiter="\t"):
            if len(row) != 5:
                continue
            label, rank, metric, sample, value = row
            if label == "Gold standard" or rank != RANK:
                continue
            try:
                out[metric] = float(value)
            except ValueError:
                pass
    return out


def fmt_time(s):
    if s is None:
        return "n/a"
    m, sec = divmod(int(round(s)), 60)
    return f"{m}m{sec:02d}s" if m else f"{sec}s"


print(f"\nBenchmark 7 (GTDB r232 metapackage) — sample {SAMPLE}, rank: {RANK}\n")
hdr = f"{'tool':<18}{'wall':>9}{'peak RAM':>12}{'F1':>8}{'Bray-Curtis':>13}{'Purity':>9}{'Complet.':>10}"
print(hdr)
print("-" * len(hdr))
for tool in ["singlem", "singlem-regime3", "sylph"]:
    wall, rss = resource_summary(tool)
    m = opal_metrics(tool)
    f1 = m.get("F1 score")
    bc = m.get("Bray-Curtis distance")
    pur = m.get("Purity")
    comp = m.get("Completeness")
    rss_gb = f"{rss/1024:.1f} GB" if rss is not None else "n/a"
    print(f"{tool:<18}{fmt_time(wall):>9}{rss_gb:>12}"
          f"{(f'{f1:.3f}' if f1 is not None else 'n/a'):>8}"
          f"{(f'{bc:.3f}' if bc is not None else 'n/a'):>13}"
          f"{(f'{pur:.3f}' if pur is not None else 'n/a'):>9}"
          f"{(f'{comp:.3f}' if comp is not None else 'n/a'):>10}")
print("\n(F1/Purity/Completeness higher=better; Bray-Curtis lower=better. "
      "regime3 wall=sum of pipe+weebill+condense; RAM=peak stage.)\n")
