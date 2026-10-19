#!/usr/bin/env python3
"""Plot adiabatic initializer pressure, mass-flow, temperature, and elevation."""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator


METHODS = (
    ("Steady spatial marching", "marching"),
    ("Coupled steady Newton", "newton"),
    ("Implicit θ pseudo-time relaxation", "relaxation"),
)
COLORS = {"pressure": "#1769aa", "flow": "#c44e52", "temperature": "#dd8b00",
          "elevation": "#477553"}


def read_profiles(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no profile data found in {path}")
    return rows


def plot_profiles(rows, output):
    node_rows = [row for row in rows if row["x_node_m"]]
    face_rows = [row for row in rows if row["x_face_m"]]
    x_node = [float(row["x_node_m"]) / 1000.0 for row in node_rows]
    x_face = [float(row["x_face_m"]) / 1000.0 for row in face_rows]
    fig, axes = plt.subplots(
        4, 1, figsize=(12.5, 13.5), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.0, 1.0, 0.62]},
    )
    flow_values = [float(row[f"m_{key}_kg_s"])
                   for _, key in METHODS for row in face_rows]
    flow_center = sum(flow_values) / len(flow_values)
    flow_span = max(abs(value - flow_center) for value in flow_values)
    flow_pad = max(1.25 * flow_span, abs(flow_center) * 1e-7)

    for ax, (title, key) in zip(axes[:3], METHODS):
        flow_ax = ax.twinx()
        temp_ax = ax.twinx()
        temp_ax.spines["right"].set_position(("outward", 92))

        p_line, = ax.plot(
            x_node, [float(row[f"p_{key}_Pa"]) / 1e6 for row in node_rows],
            color=COLORS["pressure"], linewidth=2.0, label="Pressure",
        )
        m_line, = flow_ax.plot(
            x_face, [float(row[f"m_{key}_kg_s"]) for row in face_rows],
            color=COLORS["flow"], linestyle="--", linewidth=1.6,
            marker="o", markersize=2.4, label="Mass flow",
        )
        t_line, = temp_ax.plot(
            x_node, [float(row[f"T_{key}_K"]) for row in node_rows],
            color=COLORS["temperature"], linestyle="-.", linewidth=1.7,
            label="Temperature",
        )
        ax.set_title(title, loc="left", fontsize=11.5, fontweight="semibold")
        ax.set_ylabel("Pressure (MPa)", color=COLORS["pressure"])
        flow_ax.set_ylabel("Mass flow (kg/s)", color=COLORS["flow"])
        temp_ax.set_ylabel("Temperature (K)", color=COLORS["temperature"])
        ax.tick_params(axis="y", colors=COLORS["pressure"])
        flow_ax.tick_params(axis="y", colors=COLORS["flow"])
        temp_ax.tick_params(axis="y", colors=COLORS["temperature"])
        flow_ax.set_ylim(flow_center - flow_pad, flow_center + flow_pad)
        flow_ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        # These quantities vary only slightly about O(10^2); fixed point
        # labels avoid Matplotlib's shared 1e2/2.88e2 offset text colliding
        # across the two right-hand axes.
        flow_ax.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
        temp_ax.yaxis.set_major_locator(MaxNLocator(nbins=5))
        temp_ax.yaxis.set_major_formatter(FormatStrFormatter("%.4f"))
        ax.grid(True, color="#d9dee5", linewidth=0.7, alpha=0.85)
        if ax is axes[0]:
            handles = [p_line, m_line, t_line]

    elevation = [float(row["elevation_m"]) for row in node_rows]
    axes[3].fill_between(x_node, elevation, color=COLORS["elevation"], alpha=0.20)
    axes[3].plot(x_node, elevation, color=COLORS["elevation"], linewidth=2.0)
    axes[3].set_title("Pipe elevation", loc="left", fontsize=11.5,
                      fontweight="semibold")
    axes[3].set_ylabel("Elevation (m)", color=COLORS["elevation"])
    axes[3].tick_params(axis="y", colors=COLORS["elevation"])
    axes[3].grid(True, color="#d9dee5", linewidth=0.7, alpha=0.85)
    axes[3].set_xlabel("Distance along pipe (km)")
    axes[3].set_xlim(x_node[0], x_node[-1])
    fig.suptitle("Adiabatic pipe steady-state initializers", y=0.995,
                 fontsize=15, fontweight="semibold")
    fig.legend(handles, [handle.get_label() for handle in handles],
               loc="upper center", bbox_to_anchor=(0.5, 0.97), ncol=3,
               frameon=False)
    # Reserve a wider right margin for the two independently scaled right axes.
    fig.tight_layout(rect=(0.04, 0.02, 0.75, 0.93), h_pad=1.15)
    fig.savefig(output, dpi=190, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV from adiabatic profile exporter")
    parser.add_argument("-o", "--output", type=Path,
                        default=Path("adiabatic_initializer_profiles.png"),
                        help="output image path")
    args = parser.parse_args()
    plot_profiles(read_profiles(args.csv), args.output)
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
