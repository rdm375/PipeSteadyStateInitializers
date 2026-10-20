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

"""Plot pressure and mass-flow profiles exported by pipe_init_profile_export."""

import argparse
import csv
from pathlib import Path
import math

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

METHODS = (
    ("Analytical pressure-squared", "analytical"),
    ("Discrete steady marching", "marching"),
    ("Newton steady solve", "newton"),
    ("Pseudo-time relaxation", "relaxation"),
)


def read_profiles(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no profile data found in {path}")
    return rows


def plot_profiles(rows, output):
    pressure_rows = [row for row in rows if row["x_pressure_m"]]
    pressure_x = [float(row["x_pressure_m"]) for row in pressure_rows]
    flow_rows = [row for row in rows if row["x_flow_m"]]
    flow_x = [float(row["x_flow_m"]) for row in flow_rows]

    fig, axes = plt.subplots(
        5, 1, figsize=(11, 15), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0, 1.0, 0.7]},
    )
    legend_handles = None
    flow_values = []
    flow_axes = []

    for axis, (title, suffix) in zip(axes[:4], METHODS):
        pressure = [float(row[f"p_{suffix}_Pa"]) / 1.0e6
                    for row in pressure_rows]
        flow = [float(row[f"m_{suffix}_kg_s"]) for row in flow_rows]
        flow_values.extend(flow)

        flow_axis = axis.twinx()
        flow_axes.append(flow_axis)
        pressure_line, = axis.plot(
            pressure_x, pressure, color="#1769aa", linewidth=2.2
        )
        flow_line, = flow_axis.plot(
            flow_x, flow, color="#c44e52", linestyle="--", linewidth=1.7,
            marker="o", markersize=2.2
        )
        if legend_handles is None:
            legend_handles = (pressure_line, flow_line)

        axis.set_title(title, loc="left", fontsize=12, fontweight="semibold")
        axis.set_ylabel("Pressure (MPa)", color="#1769aa")
        flow_axis.set_ylabel("Mass flow (kg/s)", color="#c44e52")
        axis.tick_params(axis="y", colors="#1769aa")
        flow_axis.tick_params(axis="y", colors="#c44e52")
        flow_axis.yaxis.set_major_locator(MaxNLocator(nbins=6))
        axis.grid(True, color="#d9dee5", linewidth=0.7, alpha=0.85)

    flow_min = min(flow_values)
    flow_max = max(flow_values)
    # The outlet flow is prescribed and therefore provides a stable center for
    # the zoomed flow axes across all methods.
    flow_center = float(flow_rows[-1]["m_analytical_kg_s"])
    flow_deviation = max(abs(value - flow_center) for value in flow_values)
    flow_padding = max(1.2 * flow_deviation, abs(flow_center) * 1.0e-8)
    tick_spacing = 2.0 * flow_padding / 6.0
    decimals = max(0, min(12, math.ceil(-math.log10(tick_spacing))))
    for flow_axis in flow_axes:
        flow_axis.set_ylim(flow_center - flow_padding, flow_center + flow_padding)
        flow_axis.yaxis.set_major_formatter(FormatStrFormatter(f"%.{decimals}f"))

    elevation_rows = [row for row in rows if row["x_elevation_m"]]
    elevation_x = [float(row["x_elevation_m"]) for row in elevation_rows]
    elevation = [float(row["elevation_m"]) for row in elevation_rows]
    axes[4].fill_between(elevation_x, elevation, color="#7b9e87", alpha=0.25)
    axes[4].plot(elevation_x, elevation, color="#477553", linewidth=2.0)
    axes[4].set_title("Pipe elevation", loc="left", fontsize=12,
                      fontweight="semibold")
    axes[4].set_ylabel("Elevation (m)", color="#477553")
    axes[4].tick_params(axis="y", colors="#477553")
    axes[4].set_ylim(bottom=0)
    axes[4].grid(True, color="#d9dee5", linewidth=0.7, alpha=0.85)
    axes[4].set_xlabel("Distance along pipe (m)")
    axes[4].set_xlim(pressure_x[0], pressure_x[-1])
    fig.suptitle("Isothermal pipe steady-state initializers",
                 fontsize=15, fontweight="semibold", y=0.99)
    fig.legend(legend_handles, ["Pressure", "Mass flow"],
               loc="upper center", bbox_to_anchor=(0.5, 0.965), ncol=2)
    fig.tight_layout(rect=(0.025, 0.02, 0.975, 0.93), h_pad=1.2)
    fig.savefig(output, dpi=180, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="profile CSV from the C++ exporter")
    parser.add_argument("-o", "--output", type=Path,
                        default=Path("initializer_profiles.png"),
                        help="output PNG path (default: initializer_profiles.png)")
    args = parser.parse_args()
    plot_profiles(read_profiles(args.csv), args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
