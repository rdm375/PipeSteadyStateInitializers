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

"""Recover elevation from the isothermal model's steady p and mass-flow data.

This inverts the same low-Mach discrete momentum residual used by the
isothermal initializer. The exported pressure nodes and staggered flows are
the only state variables used; the exported elevation is used only as a
reference after reconstruction. The inlet elevation anchors the datum.
"""

import argparse
import csv
import math
from pathlib import Path

METHODS = ("analytical", "marching", "newton", "relaxation")
COLORS = {"analytical": "#8a8a8a", "marching": "#1769aa",
          "newton": "#dd8b00", "relaxation": "#7b61a8"}


def records(path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def nonempty_rows(rows, column):
    return [row for row in rows if row.get(column, "").strip()]


def friction_factor(outlet_flow, diameter, area, roughness, viscosity):
    if outlet_flow == 0.0:
        return 0.0
    reynolds = abs(outlet_flow) * diameter / (area * viscosity)
    argument = roughness / (3.7 * diameter) + 5.74 / reynolds**0.9
    if argument <= 0.0:
        raise ValueError("invalid Swamee-Jain logarithm argument")
    return 0.25 / math.log10(argument) ** 2


def stats(values):
    return max(map(abs, values)), math.sqrt(sum(v * v for v in values) / len(values))


def reconstruct(nodes, faces, method, parameters):
    (diameter, roughness, viscosity, compressibility, gas_constant,
     temperature, gravity, outlet_flow, inlet_pressure) = parameters
    area = math.pi * diameter**2 / 4.0
    gas_product = compressibility * gas_constant * temperature
    friction = friction_factor(outlet_flow, diameter, area, roughness, viscosity)
    x = [float(row["x_pressure_m"]) for row in nodes]
    truth = [float(row["elevation_m"]) for row in nodes]
    pressure = [float(row[f"p_{method}_Pa"]) for row in nodes]
    flow = [float(row[f"m_{method}_kg_s"]) for row in faces]
    if len(flow) != len(nodes) - 1:
        raise ValueError(f"{method}: need one internal face flow per pressure cell")
    if any(p <= 0.0 for p in pressure):
        raise ValueError(f"{method}: pressure must be positive")

    recovered = [truth[0]]
    residual_with, residual_without, dz_errors = [], [], []
    for i, mass_flow in enumerate(flow):
        dx = x[i + 1] - x[i]
        p_bar = 0.5 * (pressure[i] + pressure[i + 1])
        density = p_bar / gas_product
        velocity = mass_flow / (area * density)
        drag = (friction / (2.0 * diameter) * density
                * velocity * abs(velocity))
        delta_p = pressure[i + 1] - pressure[i]
        delta_z = -(delta_p + dx * drag) / (density * gravity)
        recovered.append(recovered[-1] + delta_z)
        dz_errors.append(delta_z - (truth[i + 1] - truth[i]))

        # Cell-form residual, normalized by inlet pressure, matching the
        # source residual's (Pa/m)/(p_in/dx) scaling.
        residual_without.append((delta_p + dx * drag) / inlet_pressure)
        residual_with.append((delta_p + dx * drag
                              + density * gravity * (truth[i + 1] - truth[i]))
                             / inlet_pressure)
    misfit = [recovered[i] - truth[i] for i in range(len(truth))]
    return x, truth, recovered, misfit, dz_errors, residual_with, residual_without


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", nargs="?", type=Path, default=Path("initializer_profiles.csv"))
    parser.add_argument("--output", type=Path,
                        default=Path("results/isothermal_elevation_recovery.csv"))
    parser.add_argument("--plot", type=Path,
                        default=Path("results/isothermal_elevation_recovery.png"))
    parser.add_argument("--diameter", type=float, default=0.6)
    parser.add_argument("--roughness", type=float, default=4.5e-5)
    parser.add_argument("--viscosity", type=float, default=1.1e-5)
    parser.add_argument("--compressibility", type=float, default=0.9)
    parser.add_argument("--gas-constant", type=float, default=518.3)
    parser.add_argument("--temperature", type=float, default=288.15)
    parser.add_argument("--gravity", type=float, default=9.80665)
    parser.add_argument("--inlet-pressure", type=float, default=8.0e6)
    parser.add_argument("--outlet-flow", type=float, default=100.0)
    args = parser.parse_args()

    all_rows = records(args.csv)
    if not all_rows:
        raise ValueError(f"no profile data in {args.csv}")
    nodes = nonempty_rows(all_rows, "x_pressure_m")
    elevation_rows = nonempty_rows(all_rows, "x_elevation_m")
    if len(nodes) != len(elevation_rows):
        raise ValueError("pressure and elevation node arrays differ in length")
    for node, elev in zip(nodes, elevation_rows):
        if float(node["x_pressure_m"]) != float(elev["x_elevation_m"]):
            raise ValueError("pressure and elevation nodes are not aligned")
        node["elevation_m"] = elev["elevation_m"]

    flow_rows = nonempty_rows(all_rows, "x_flow_m")
    length = float(nodes[-1]["x_pressure_m"])
    faces = [row for row in flow_rows
             if 0.0 < float(row["x_flow_m"]) < length]
    faces.sort(key=lambda row: float(row["x_flow_m"]))
    if len(faces) != len(nodes) - 1:
        raise ValueError(f"expected {len(nodes)-1} interior flows; found {len(faces)}")
    # Confirm each flow lies at the midpoint of its pressure interval.
    for i, row in enumerate(faces):
        midpoint = 0.5 * (float(nodes[i]["x_pressure_m"])
                          + float(nodes[i + 1]["x_pressure_m"]))
        if not math.isclose(float(row["x_flow_m"]), midpoint, abs_tol=1.0e-8):
            raise ValueError("flow faces are not aligned with pressure-cell midpoints")

    params = (args.diameter, args.roughness, args.viscosity,
              args.compressibility, args.gas_constant, args.temperature,
              args.gravity, args.outlet_flow, args.inlet_pressure)
    results = {method: reconstruct(nodes, faces, method, params) for method in METHODS}

    args.output.parent.mkdir(parents=True, exist_ok=True)
    columns = ["x_node_m", "elevation_known_m"]
    for method in METHODS:
        columns.extend((f"elevation_recovered_{method}_m",
                        f"cumulative_misfit_{method}_m"))
    for method in METHODS:
        columns.extend((f"momentum_residual_with_gravity_{method}_scaled",
                        f"momentum_residual_without_gravity_{method}_scaled"))
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(columns)
        for i in range(len(nodes)):
            row = [results["newton"][0][i], results["newton"][1][i]]
            for method in METHODS:
                result = results[method]
                row.extend((result[2][i], result[3][i]))
            for method in METHODS:
                result = results[method]
                if i == 0:
                    row.extend(("", ""))
                else:
                    row.extend((result[5][i - 1], result[6][i - 1]))
            writer.writerow(row)

    print(f"Isothermal hydraulic elevation recovery — {args.csv}")
    print("initializer   cumulative max/RMS (m)   cell Δz max/RMS (m)   momentum residual max with/without g")
    for method, result in results.items():
        cmax, crms = stats(result[3])
        imax, irms = stats(result[4])
        r_on = max(map(abs, result[5]))
        r_off = max(map(abs, result[6]))
        print(f"{method:11s} {cmax:.3e}/{crms:.3e}   {imax:.3e}/{irms:.3e}"
              f"   {r_on:.3e}/{r_off:.3e}")
    print(f"CSV:  {args.output}")

    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)
    for method, result in results.items():
        distance_km = [value / 1000.0 for value in result[0]]
        midpoints_km = [(result[0][i] + result[0][i + 1]) / 2000.0
                        for i in range(len(result[0]) - 1)]
        axes[0].plot(distance_km, result[3], color=COLORS[method], label=method)
        axes[1].plot(midpoints_km, result[4], color=COLORS[method], label=method)
    newton = results["newton"]
    midpoints_km = [(newton[0][i] + newton[0][i + 1]) / 2000.0
                    for i in range(len(newton[0]) - 1)]
    axes[2].plot(midpoints_km, newton[5], color="#1769aa", label="Newton, gravity included")
    axes[2].plot(midpoints_km, newton[6], color="#c23b22", linestyle="--",
                 label="Newton, gravity omitted")
    axes[0].set_title("Elevation recovered from the isothermal hydraulic balance")
    axes[0].set_ylabel("Recovered − known (m)")
    axes[1].set_title("Cellwise elevation-increment misfit")
    axes[1].set_ylabel("Δz recovered − known (m)")
    axes[2].set_title("Scaled Newton momentum residual, with and without gravity")
    axes[2].set_ylabel("Residual / inlet pressure")
    axes[2].set_xlabel("Distance along pipe (km)")
    for axis in axes:
        axis.axhline(0.0, color="black", linewidth=0.9)
        axis.grid(True, alpha=0.3)
    axes[0].legend(ncol=4, loc="best")
    axes[1].legend(ncol=4, loc="best")
    axes[2].legend(loc="best")
    fig.tight_layout()
    args.plot.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.plot, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)
    print(f"Plot: {args.plot}")


if __name__ == "__main__":
    main()
