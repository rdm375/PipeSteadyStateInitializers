#!/usr/bin/env python3
"""Summarize per-run CSV output from pipe_init_benchmark."""

import argparse
import csv
from collections import defaultdict
from pathlib import Path
import statistics


METHOD_LABELS = {
    "analytical": "Analytical pressure-squared",
    "discrete_marching": "Discrete steady marching",
    "newton_from_flat": "Newton from flat initial state",
    "theta_relaxation": "Theta pseudo-time relaxation",
}


def read_runs(path):
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream)
        required = {"method", "repeat", "seconds", "steps",
                    "scaled_steady_residual"}
        if not required.issubset(reader.fieldnames or ()):
            raise ValueError(f"{path} is not a pipe_init_benchmark CSV")
        fields = set(reader.fieldnames or ())
        if "newton_iterations" in fields:
            iteration_column = "newton_iterations"
            has_inner_counts = True
        elif "iterations" in fields:
            iteration_column = "iterations"
            has_inner_counts = False
        else:
            raise ValueError(f"{path} has no Newton iteration-count column")
        groups = defaultdict(list)
        for row in reader:
            method = row["method"]
            if method not in METHOD_LABELS:
                raise ValueError(f"unknown benchmark method: {method}")
            groups[method].append({
                "seconds": float(row["seconds"]),
                "newton_iterations": (
                    int(row[iteration_column])
                    if has_inner_counts or method != "theta_relaxation"
                    else None
                ),
                "steps": int(row["steps"]),
                "residual": float(row["scaled_steady_residual"]),
            })
    if not groups:
        raise ValueError(f"no timing runs found in {path}")
    return groups


def duration(value):
    if value >= 1.0:
        return f"{value:.3f} s"
    if value >= 1.0e-3:
        return f"{value*1.0e3:.3f} ms"
    if value >= 1.0e-6:
        return f"{value*1.0e6:.3f} μs"
    return f"{value*1.0e9:.3f} ns"


def make_report(groups, csv_path):
    medians = {name: statistics.median(r["seconds"] for r in runs)
               for name, runs in groups.items()}
    fastest = min(medians.values())
    lines = [
        "# Initializer Timing Report",
        "",
        f"Source data: `{csv_path.name}`",
        "",
        "Standard case: 100 km pipe, 100 cells, 0.6 m diameter, 8 MPa inlet, "
        "100 kg/s outlet; theta relaxation uses θ = 0.5.",
        "Newton starts from uniform inlet pressure and prescribed flow. "
        "Relaxation runs until its scaled derivative tolerance is met.",
        "",
        "| Method | Runs | Median | Mean | Std. dev. | Min–max | "
        "Iterations / steps | Max scaled steady residual | vs. fastest |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for method in METHOD_LABELS:
        if method not in groups:
            continue
        runs = groups[method]
        times = [row["seconds"] for row in runs]
        work = {(row["newton_iterations"], row["steps"]) for row in runs}
        work = sorted(work, key=lambda item: (item[0] is None,
                                               item[0] or 0, item[1]))
        work_descriptions = []
        for iterations, steps in work:
            if iterations is None:
                work_descriptions.append(
                    f"{steps} pseudo-time steps; inner Newton count unavailable")
            elif steps:
                work_descriptions.append(
                    f"{iterations} Newton iterations / {steps} pseudo-time steps")
            else:
                work_descriptions.append(f"{iterations} Newton iterations")
        work_text = ", ".join(work_descriptions)
        stddev = statistics.stdev(times) if len(times) > 1 else 0.0
        lines.append(
            f"| {METHOD_LABELS[method]} | {len(runs)} | "
            f"{duration(statistics.median(times))} | "
            f"{duration(statistics.mean(times))} | {duration(stddev)} | "
            f"{duration(min(times))}–{duration(max(times))} | {work_text} | "
            f"{max(row['residual'] for row in runs):.3e} | "
            f"{medians[method]/fastest:.2f}× |"
        )
    lines.extend([
        "",
        "Timing ratios use the fastest measured median as the reference. "
        "The analytical initializer, direct nonlinear solve, and transient "
        "relaxation have different work and stopping criteria, so timing is "
        "a performance comparison rather than an accuracy ranking.",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV written by pipe_init_benchmark")
    parser.add_argument("-o", "--output", type=Path,
                        help="Markdown report path (default: <csv stem>_report.md)")
    args = parser.parse_args()
    output = args.output or args.csv.with_name(f"{args.csv.stem}_report.md")
    report = make_report(read_runs(args.csv), args.csv)
    output.write_text(report, encoding="utf-8")
    print(report)
    print(f"\nWrote {output}")


if __name__ == "__main__":
    main()
