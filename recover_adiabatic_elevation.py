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

"""Recover the adiabatic pipe elevation from exported p, T, and flow profiles.

The exporter CSV includes the known elevation for comparison. The reconstruction
uses only pressure, temperature, and mass flow, with the inlet elevation datum
set by --datum. The energy and momentum equations provide independent estimates.
"""

import argparse
import csv
import math
from pathlib import Path

METHODS = ("marching", "newton", "relaxation")


def floats(rows, column):
    return [float(row[column]) for row in rows]


def cumulative_profile(deltas, datum):
    result = [datum]
    for delta in deltas:
        result.append(result[-1] + delta)
    return result


def swamee_jain(mass_flow, diameter, area, roughness, viscosity):
    if mass_flow == 0.0:
        return 0.0
    reynolds = abs(mass_flow) * diameter / (area * viscosity)
    argument = roughness / (3.7 * diameter) + 5.74 / reynolds**0.9
    return 0.25 / math.log10(argument)**2


def reconstruct(nodes, faces, method, config, datum):
    x = floats(nodes, "x_node_m")
    pressure = floats(nodes, f"p_{method}_Pa")
    temperature = floats(nodes, f"T_{method}_K")
    flow = floats(faces, f"m_{method}_kg_s")

    if len(x) != len(pressure) or len(x) != len(temperature):
        raise ValueError(f"{method}: inconsistent node arrays")
    if len(flow) != len(x) - 1:
        raise ValueError(
            f"{method}: expected one face flow per cell "
            f"({len(x) - 1}), found {len(flow)}"
        )
    if any(right <= left for left, right in zip(x, x[1:])):
        raise ValueError("node positions must be strictly increasing")
    if any(p <= 0.0 for p in pressure) or any(t <= 0.0 for t in temperature):
        raise ValueError("pressure and temperature must be positive")

    diameter, roughness, viscosity, gas_constant, gamma, gravity = config
    area = math.pi * diameter**2 / 4.0
    cp = gamma * gas_constant / (gamma - 1.0)
    density = [p / (gas_constant * t) for p, t in zip(pressure, temperature)]
    recovered_energy = []
    recovered_momentum = []

    for i, mass_flow in enumerate(flow):
        dx = x[i + 1] - x[i]
        rho_left = density[i]
        rho_right = density[i + 1]
        u_left = mass_flow / (area * rho_left)
        u_right = mass_flow / (area * rho_right)

        # Discrete total-enthalpy balance:
        # cp*T + u^2/2 + g*z is constant across each adiabatic cell.
        delta_z_energy = (
            cp * (temperature[i] - temperature[i + 1])
            + 0.5 * (u_left**2 - u_right**2)
        ) / gravity
        recovered_energy.append(delta_z_energy)

        # Invert the same discrete momentum residual used by the C++ model.
        rho_bar = 0.5 * (rho_left + rho_right)
        u_bar = 0.5 * (u_left + u_right)
        friction_factor = swamee_jain(
            mass_flow, diameter, area, roughness, viscosity
        )
        drag = (
            friction_factor / (2.0 * diameter)
            * rho_bar * u_bar * abs(u_bar)
        )
        delta_z_momentum = -(
            (pressure[i + 1] - pressure[i])
            + (rho_right * u_right**2 - rho_left * u_left**2)
            + dx * drag
        ) / (rho_bar * gravity)
        recovered_momentum.append(delta_z_momentum)

    return (
        x,
        cumulative_profile(recovered_energy, datum),
        cumulative_profile(recovered_momentum, datum),
    )


def error_metrics(reconstructed, truth):
    errors = [a - b for a, b in zip(reconstructed, truth)]
    return max(abs(value) for value in errors), math.sqrt(
        sum(value * value for value in errors) / len(errors)
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV from adiabatic profile exporter")
    parser.add_argument(
        "--output", type=Path, default=Path("adiabatic_elevation_recovery.csv"),
        help="reconstructed profile CSV",
    )
    parser.add_argument(
        "--plot", type=Path, default=Path("adiabatic_elevation_recovery.png"),
        help="comparison plot",
    )
    parser.add_argument("--no-plot", action="store_true", help="skip plot generation")
    parser.add_argument("--datum", type=float, default=0.0,
                        help="assumed elevation at the inlet node, in metres")
    parser.add_argument("--diameter", type=float, default=0.6, help="pipe diameter, m")
    parser.add_argument("--roughness", type=float, default=4.5e-5,
                        help="absolute roughness, m")
    parser.add_argument("--viscosity", type=float, default=1.1e-5,
                        help="dynamic viscosity, Pa s")
    parser.add_argument("--gas-constant", type=float, default=518.3,
                        help="specific gas constant, J/(kg K)")
    parser.add_argument("--gamma", type=float, default=1.30,
                        help="heat-capacity ratio")
    parser.add_argument("--gravity", type=float, default=9.80665,
                        help="gravitational acceleration, m/s^2")
    args = parser.parse_args()

    with args.csv.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no profile data found in {args.csv}")
    required = {
        "x_node_m", "elevation_m",
        *(f"p_{method}_Pa" for method in METHODS),
        *(f"T_{method}_K" for method in METHODS),
        "x_face_m",
        *(f"m_{method}_kg_s" for method in METHODS),
    }
    missing = required.difference(rows[0])
    if missing:
        raise ValueError(f"missing CSV columns: {', '.join(sorted(missing))}")

    nodes = [row for row in rows if row["x_node_m"]]
    faces = [row for row in rows if row["x_face_m"]]
    truth = floats(nodes, "elevation_m")
    config = (
        args.diameter, args.roughness, args.viscosity,
        args.gas_constant, args.gamma, args.gravity,
    )
    reconstructions = {}
    x = None
    for method in METHODS:
        x, energy_z, momentum_z = reconstruct(
            nodes, faces, method, config, args.datum
        )
        reconstructions[method] = (energy_z, momentum_z)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    headers = ["x_node_m", "elevation_true_m"]
    for method in METHODS:
        headers.extend((f"z_energy_{method}_m", f"z_momentum_{method}_m"))
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(headers)
        for i, x_i in enumerate(x):
            row = [x_i, truth[i]]
            for method in METHODS:
                energy_z, momentum_z = reconstructions[method]
                row.extend((energy_z[i], momentum_z[i]))
            writer.writerow(row)

    print(f"Recovered profile written to {args.output}")
    print("method       energy max/RMS error (m)     momentum max/RMS error (m)")
    for method in METHODS:
        energy_z, momentum_z = reconstructions[method]
        e_max, e_rms = error_metrics(energy_z, truth)
        m_max, m_rms = error_metrics(momentum_z, truth)
        print(f"{method:10s}   {e_max: .6e} / {e_rms: .6e}       "
              f"{m_max: .6e} / {m_rms: .6e}")

    if not args.no_plot:
        import matplotlib.pyplot as plt

        distance_km = [value / 1000.0 for value in x]
        fig, axes = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
        axes[0].plot(distance_km, truth, color="black", linewidth=2.5,
                     label="Known elevation")
        axes[1].plot(distance_km, truth, color="black", linewidth=2.5,
                     label="Known elevation")
        colors = {"marching": "#1769aa", "newton": "#dd8b00",
                  "relaxation": "#7b61a8"}
        for method in METHODS:
            energy_z, momentum_z = reconstructions[method]
            axes[0].plot(distance_km, energy_z, color=colors[method],
                         linestyle="--", linewidth=1.5,
                         label=f"{method}: energy")
            axes[1].plot(distance_km, momentum_z, color=colors[method],
                         linestyle="--", linewidth=1.5,
                         label=f"{method}: momentum")
        axes[0].set_title("Elevation reconstructed from the adiabatic energy balance")
        axes[1].set_title("Elevation reconstructed from the adiabatic momentum balance")
        for axis in axes:
            axis.set_ylabel("Elevation (m)")
            axis.grid(True, alpha=0.35)
            axis.legend(loc="best", ncol=2)
        axes[1].set_xlabel("Distance along pipe (km)")
        fig.tight_layout()
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.plot, dpi=180, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        print(f"Comparison plot written to {args.plot}")


if __name__ == "__main__":
    main()

