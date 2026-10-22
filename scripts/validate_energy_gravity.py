#!/usr/bin/env python3
# Copyright (C) 2026 Richard Myers
# SPDX-License-Identifier: AGPL-3.0-only
#
# Licensed under the GNU Affero General Public License, version 3 only.
# See the LICENSE file for the complete terms.
#
# A separate commercial license for proprietary use may be available by
# written agreement with the relevant rights holders. This notice does not
# grant that license, guarantee its availability, or set fees, royalties,
# or revenue shares.

"""Check steady energy balances and gravity contribution from profile CSVs."""

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt


ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
METHODS = ("analytical", "marching", "newton", "relaxation")
COLORS = {"analytical": "#8a8a8a", "marching": "#1769aa",
          "newton": "#dd8b00", "relaxation": "#7b61a8"}
R = 518.3
GAMMA = 1.30
CP = GAMMA * R / (GAMMA - 1.0)
TIN = 288.15
MOUT = 100.0
GRAVITY = 9.80665
DIAMETER = 0.6
AREA = math.pi * DIAMETER**2 / 4.0
INNER_CONDUCTANCE = 10.0 * math.pi * DIAMETER


def read_rows(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"No profile data in {path}")
    return rows


def populated(rows, column):
    return [row for row in rows if row.get(column, "").strip()]


def stats(values):
    return max(map(abs, values)), math.sqrt(sum(v * v for v in values) / len(values))


def adiabatic(rows, method):
    nodes = populated(rows, "x_node_m")
    faces = populated(rows, "x_face_m")
    if len(faces) != len(nodes) - 1:
        raise ValueError("Adiabatic profile has inconsistent node and face counts")
    x = [float(r["x_node_m"]) for r in nodes]
    z = [float(r["elevation_m"]) for r in nodes]
    p = [float(r[f"p_{method}_Pa"]) for r in nodes]
    t = [float(r[f"T_{method}_K"]) for r in nodes]
    flow = [float(r[f"m_{method}_kg_s"]) for r in faces]
    recovered = [z[0]]
    residual_on, residual_off = [None], [None]
    scale = CP * TIN
    for i in range(len(flow)):
        rho_left = p[i] / (R * t[i])
        rho_right = p[i + 1] / (R * t[i + 1])
        u_left = flow[i] / (AREA * rho_left)
        u_right = flow[i] / (AREA * rho_right)
        delta_h = CP * (t[i + 1] - t[i]) + 0.5 * (u_right**2 - u_left**2)
        dz = -delta_h / GRAVITY
        recovered.append(recovered[-1] + dz)
        residual_off.append(delta_h / scale)
        residual_on.append((delta_h + GRAVITY * (z[i + 1] - z[i])) / scale)
    misfit = [a - b for a, b in zip(recovered, z)]
    return x, z, recovered, misfit, residual_on, residual_off


def thermal(rows, method):
    nodes = populated(rows, "x_node_m")
    faces = populated(rows, "x_face_m")
    if len(faces) != len(nodes) - 1:
        raise ValueError("Thermal profile has inconsistent node and face counts")
    x = [float(r["x_node_m"]) for r in nodes]
    z = [float(r["elevation_m"]) for r in nodes]
    tg = [float(r[f"Tg_{method}_K"]) for r in nodes]
    tw = [float(r[f"Tw_{method}_K"]) for r in nodes]
    flow = [float(r[f"m_{method}_kg_s"]) for r in faces]
    cells = len(nodes) - 1
    dx = (x[-1] - x[0]) / cells
    recovered = [z[0]]
    residual_on, residual_off = [None], [None]
    scale = MOUT * CP * TIN
    for i in range(1, len(nodes)):
        ml = flow[i - 1]
        mr = flow[i] if i < len(faces) else MOUT
        ell = dx * (0.5 if i == len(faces) else 1.0)
        heat = INNER_CONDUCTANCE * (tw[i] - tg[i])
        sensible_flux = CP * (mr * tg[i] - ml * tg[i - 1])
        grav_flux_known = GRAVITY * (mr * z[i] - ml * z[i - 1])
        # Invert the same discrete energy residual, anchoring z at the inlet.
        zi = (heat * ell - sensible_flux + GRAVITY * ml * recovered[-1]) / (GRAVITY * mr)
        recovered.append(zi)
        residual_on.append((sensible_flux + grav_flux_known - heat * ell) / scale)
        residual_off.append((sensible_flux - heat * ell) / scale)
    misfit = [a - b for a, b in zip(recovered, z)]
    return x, z, recovered, misfit, residual_on, residual_off


def write_and_plot(title, csv_name, plot_name, source, evaluator, methods=METHODS):
    rows = read_rows(ROOT / source)
    results = {m: evaluator(rows, m) for m in methods}
    output = RESULTS / csv_name
    plot = RESULTS / plot_name
    output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["x_node_m", "elevation_known_m"]
    for method in methods:
        columns.extend((f"elevation_recovered_{method}_m",
                        f"elevation_misfit_{method}_m",
                        f"energy_residual_with_gravity_{method}_scaled",
                        f"energy_residual_without_gravity_{method}_scaled"))
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        n = len(results["newton"][0])
        for i in range(n):
            row = [results["newton"][0][i], results["newton"][1][i]]
            for method in methods:
                _, _, recovered, misfit, on, off = results[method]
                row.extend((recovered[i], misfit[i],
                            "" if on[i] is None else on[i],
                            "" if off[i] is None else off[i]))
            writer.writerow(row)

    print(f"\n{title} — {source}")
    print("initializer     elevation max/RMS (m)   energy residual max/RMS with g   without g")
    for method, data in results.items():
        _, _, _, misfit, on, off = data
        on_vals = [v for v in on if v is not None]
        off_vals = [v for v in off if v is not None]
        em, er = stats(misfit)
        om, ore = stats(on_vals)
        fm, fre = stats(off_vals)
        print(f"{method:12s} {em:.3e}/{er:.3e}       {om:.3e}/{ore:.3e}       {fm:.3e}/{fre:.3e}")
    print(f"CSV:  {output.relative_to(ROOT)}")

    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for method, data in results.items():
        x, _, _, misfit, on, off = data
        axes[0].plot([v / 1000.0 for v in x], misfit,
                     color=COLORS[method], label=method)
        midpoint = [(x[i] + x[i + 1]) / 2000.0 for i in range(len(x) - 1)]
        axes[1].plot(midpoint, on[1:], color=COLORS[method], label=f"{method}, gravity included")
        if method == "newton":
            axes[1].plot(midpoint, off[1:], color="#c23b22", linestyle="--",
                         label="Newton, gravity omitted")
    axes[0].set_title(f"{title}: elevation recovered from the energy balance")
    axes[0].set_ylabel("Recovered − known (m)")
    axes[1].set_title("Scaled cell energy residual, with and without gravity")
    axes[1].set_ylabel("Residual / (ṁ cp Tin)")
    axes[1].set_xlabel("Distance along pipe (km)")
    for axis in axes:
        axis.axhline(0.0, color="black", linewidth=0.9)
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=2, loc="best")
    fig.tight_layout()
    fig.savefig(plot, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Plot: {plot.relative_to(ROOT)}")


def main():
    write_and_plot("Adiabatic ideal gas",
                   "adiabatic_energy_gravity_validation.csv",
                   "adiabatic_energy_gravity_validation.png",
                   "adiabatic_profiles.csv", adiabatic,
                   methods=("marching", "newton", "relaxation"))
    write_and_plot("Thermal gas and lumped wall",
                   "thermal_wall_energy_gravity_validation.csv",
                   "thermal_wall_energy_gravity_validation.png",
                   "thermal_wall_profiles.csv", thermal)
    write_and_plot("Fully thermal gas, wall, and soil-shell",
                   "soil_shell_energy_gravity_validation.csv",
                   "soil_shell_energy_gravity_validation.png",
                   "results/soil_shell_profiles.csv", thermal)


if __name__ == "__main__":
    main()
