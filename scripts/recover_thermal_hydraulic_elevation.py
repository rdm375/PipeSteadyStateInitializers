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

"""Reconstruct pipe elevation from thermal-model hydraulic state profiles.

The thermal-wall and soil-shell profile exporters provide nodal pressure and
gas temperature, plus staggered cell mass flow. This diagnostic inverts the
same steady momentum residual used by both models to recover each cell's
elevation change. It anchors the inlet elevation at the exported value (zero
for the README hill) and compares the cumulative reconstruction and local
elevation increments against the known profile.

Run at the repository root after exporting both profiles:

    python3 scripts/recover_thermal_hydraulic_elevation.py

It writes a recovery CSV and signed-misfit plot for each model.
"""

import argparse
import csv
import math
from pathlib import Path


METHODS = ("analytical", "marching", "newton", "relaxation")
MODELS = {
    "thermal_wall": {
        "label": "Thermal gas and lumped wall",
        "input": Path("thermal_wall_profiles.csv"),
        "output": Path("results/thermal_wall_elevation_recovery.csv"),
        "plot": Path("results/thermal_wall_elevation_recovery.png"),
    },
    "soil_shell": {
        "label": "Fully thermal gas, wall, and soil-shell",
        "input": Path("results/soil_shell_profiles.csv"),
        "output": Path("results/soil_shell_elevation_recovery.csv"),
        "plot": Path("results/soil_shell_elevation_recovery.png"),
    },
}


def read_profile(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no profile rows in {path}")
    nodes = [row for row in rows if row.get("x_node_m", "").strip()]
    faces = [row for row in rows if row.get("x_face_m", "").strip()]
    if len(faces) != len(nodes) - 1:
        raise ValueError(
            f"{path}: expected {len(nodes) - 1} populated face rows; found {len(faces)}"
        )
    return nodes, faces


def friction_factor(mass_flow, diameter, area, roughness, viscosity):
    reynolds = abs(mass_flow) * diameter / (area * viscosity)
    if reynolds <= 0.0:
        return 0.0
    if reynolds < 2300.0:
        return 64.0 / reynolds
    term = roughness / (3.7 * diameter) + 5.74 / reynolds**0.9
    return 0.25 / math.log10(term) ** 2


def metrics(values):
    if not values:
        return 0.0, 0.0
    return max(abs(value) for value in values), math.sqrt(
        sum(value * value for value in values) / len(values)
    )


def reconstruct(nodes, faces, method, parameters):
    diameter, roughness, viscosity, gas_constant, gravity = parameters
    area = math.pi * diameter**2 / 4.0
    x = [float(row["x_node_m"]) for row in nodes]
    known = [float(row["elevation_m"]) for row in nodes]
    pressure = [float(row[f"p_{method}_Pa"]) for row in nodes]
    temperature = [float(row[f"Tg_{method}_K"]) for row in nodes]
    flow = [float(row[f"m_{method}_kg_s"]) for row in faces]

    if any(right <= left for left, right in zip(x, x[1:])):
        raise ValueError("node positions must increase strictly")
    if any(p <= 0.0 for p in pressure) or any(t <= 0.0 for t in temperature):
        raise ValueError(f"{method}: pressure and temperature must be positive")

    density = [p / (gas_constant * t) for p, t in zip(pressure, temperature)]
    dz = []
    for i, mass_flow in enumerate(flow):
        dx = x[i + 1] - x[i]
        rho_bar = 0.5 * (density[i] + density[i + 1])
        velocity = mass_flow / (area * rho_bar)
        f = friction_factor(mass_flow, diameter, area, roughness, viscosity)
        drag = f / (2.0 * diameter) * rho_bar * velocity * abs(velocity)

        # Invert the source model's cell residual:
        # dp/dx + rho_bar*g*dz/dx + f/(2D)*rho_bar*u*|u| = 0.
        dz.append(-((pressure[i + 1] - pressure[i]) + dx * drag)
                  / (rho_bar * gravity))

    recovered = [known[0]]
    for increment in dz:
        recovered.append(recovered[-1] + increment)
    increment_misfit = [
        dz[i] - (known[i + 1] - known[i]) for i in range(len(dz))
    ]
    cumulative_misfit = [recovered[i] - known[i] for i in range(len(x))]
    return x, known, recovered, increment_misfit, cumulative_misfit


def process_model(name, settings, args, parameters):
    source = getattr(args, f"{name}_csv") or settings["input"]
    output = getattr(args, f"{name}_output") or settings["output"]
    plot = getattr(args, f"{name}_plot") or settings["plot"]
    nodes, faces = read_profile(source)
    results = {}
    for method in METHODS:
        results[method] = reconstruct(nodes, faces, method, parameters)

    output.parent.mkdir(parents=True, exist_ok=True)
    headers = ["x_node_m", "elevation_known_m"]
    for method in METHODS:
        headers.extend((f"elevation_recovered_{method}_m",
                        f"cumulative_misfit_{method}_m"))
    with output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(headers)
        for i in range(len(nodes)):
            row = [results["newton"][0][i], results["newton"][1][i]]
            for method in METHODS:
                x, known, recovered, _, misfit = results[method]
                row.extend((recovered[i], misfit[i]))
            writer.writerow(row)

    print(f"\n{settings['label']} ({source})")
    print("initializer   cumulative max/RMS (m)   cell-increment max/RMS (m)")
    for method in METHODS:
        _, _, _, increments, cumulative = results[method]
        cmax, crms = metrics(cumulative)
        imax, irms = metrics(increments)
        print(f"{method:11s}  {cmax: .6e}/{crms: .6e}     {imax: .6e}/{irms: .6e}")
    print(f"CSV:  {output}")

    import matplotlib.pyplot as plt

    colors = {"analytical": "#8a8a8a", "marching": "#1769aa",
              "newton": "#dd8b00", "relaxation": "#7b61a8"}
    fig, axes = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
    for method in METHODS:
        x, _, _, increments, cumulative = results[method]
        distance_km = [value / 1000.0 for value in x]
        cell_midpoints_km = [0.5 * (x[i] + x[i + 1]) / 1000.0
                             for i in range(len(x) - 1)]
        axes[0].plot(distance_km, cumulative, color=colors[method],
                     linewidth=1.6, label=method)
        axes[1].plot(cell_midpoints_km, increments, color=colors[method],
                     linewidth=1.4, marker=".", markersize=4, label=method)
    axes[0].set_title(f"{settings['label']}: reconstructed elevation misfit")
    axes[1].set_title("Cellwise momentum-balance elevation-increment misfit")
    axes[0].set_ylabel("Recovered − known (m)")
    axes[1].set_ylabel("Δz recovered − known (m)")
    axes[1].set_xlabel("Distance along pipe (km)")
    for axis in axes:
        axis.axhline(0.0, color="black", linewidth=1.0)
        axis.grid(True, alpha=0.3)
        axis.legend(ncol=4, loc="best")
    fig.tight_layout()
    plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(plot, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Plot: {plot}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name, settings in MODELS.items():
        parser.add_argument(f"--{name.replace('_', '-')}-csv", type=Path,
                            help=f"{name} exporter CSV (default: {settings['input']})")
        parser.add_argument(f"--{name.replace('_', '-')}-output", type=Path,
                            help=f"{name} recovery CSV")
        parser.add_argument(f"--{name.replace('_', '-')}-plot", type=Path,
                            help=f"{name} recovery plot")
    parser.add_argument("--diameter", type=float, default=0.6, help="pipe diameter (m)")
    parser.add_argument("--roughness", type=float, default=4.5e-5,
                        help="absolute roughness (m)")
    parser.add_argument("--viscosity", type=float, default=1.1e-5,
                        help="dynamic viscosity (Pa s)")
    parser.add_argument("--gas-constant", type=float, default=518.3,
                        help="specific gas constant (J/(kg K))")
    parser.add_argument("--gravity", type=float, default=9.80665,
                        help="gravitational acceleration (m/s^2)")
    args = parser.parse_args()
    parameters = (args.diameter, args.roughness, args.viscosity,
                  args.gas_constant, args.gravity)
    for name, settings in MODELS.items():
        process_model(name, settings, args, parameters)


if __name__ == "__main__":
    main()
