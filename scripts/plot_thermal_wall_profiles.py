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

"""Plot pressure, flow, gas/wall temperature, and elevation for four methods."""
import argparse
import csv
from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter, MaxNLocator

METHODS = (
    ("Analytical approximation", "analytical"),
    ("Steady spatial marching", "marching"),
    ("Coupled steady Newton", "newton"),
    ("Implicit θ pseudo-time relaxation", "relaxation"),
)
COLORS = {"p": "#1769aa", "m": "#c44e52", "tg": "#dd8b00",
          "tw": "#7b61a8", "z": "#477553"}

def read_profiles(path):
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"no profile data found in {path}")
    return rows

def plot_profiles(rows, output):
    nodes = [r for r in rows if r["x_node_m"]]
    faces = [r for r in rows if r["x_face_m"]]
    x = [float(r["x_node_m"])/1000 for r in nodes]
    xf = [float(r["x_face_m"])/1000 for r in faces]
    fig, axes = plt.subplots(5, 1, figsize=(12.5, 16.0), sharex=True,
        gridspec_kw={"height_ratios": [1, 1, 1, 1, 0.62]})
    flows = [float(r[f"m_{key}_kg_s"]) for _, key in METHODS for r in faces]
    center = sum(flows)/len(flows)
    pad = max(1.25*max(abs(v-center) for v in flows), abs(center)*1e-7)
    handles = None
    for ax, (title, key) in zip(axes[:4], METHODS):
        maxis = ax.twinx()
        taxis = ax.twinx()
        taxis.spines["right"].set_position(("outward", 92))
        p, = ax.plot(x, [float(r[f"p_{key}_Pa"])/1e6 for r in nodes],
                     color=COLORS["p"], lw=1.9, label="Pressure")
        m, = maxis.plot(xf, [float(r[f"m_{key}_kg_s"]) for r in faces],
                        color=COLORS["m"], ls="--", marker="o", ms=2, lw=1.5,
                        label="Mass flow")
        tg, = taxis.plot(x, [float(r[f"Tg_{key}_K"]) for r in nodes],
                         color=COLORS["tg"], lw=1.8, label="Gas temperature")
        tw, = taxis.plot(x, [float(r[f"Tw_{key}_K"]) for r in nodes],
                         color=COLORS["tw"], ls=":", lw=1.8,
                         label="Wall temperature")
        if handles is None:
            handles = [p, m, tg, tw]
        ax.set_title(title, loc="left", fontsize=11, fontweight="semibold")
        ax.set_ylabel("Pressure (MPa)", color=COLORS["p"])
        maxis.set_ylabel("Mass flow (kg/s)", color=COLORS["m"])
        taxis.set_ylabel("Temperature (K)", color=COLORS["tg"])
        ax.tick_params(axis="y", colors=COLORS["p"])
        maxis.tick_params(axis="y", colors=COLORS["m"])
        taxis.tick_params(axis="y", colors=COLORS["tg"])
        maxis.set_ylim(center-pad, center+pad)
        maxis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        maxis.yaxis.set_major_formatter(FormatStrFormatter("%.6f"))
        taxis.yaxis.set_major_locator(MaxNLocator(nbins=5))
        taxis.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))
        ax.grid(True, color="#d9dee5", linewidth=0.7, alpha=0.85)
    z = [float(r["elevation_m"]) for r in nodes]
    axes[4].fill_between(x, z, color=COLORS["z"], alpha=0.2)
    axes[4].plot(x, z, color=COLORS["z"], lw=1.9)
    axes[4].set_title("Pipe elevation", loc="left", fontsize=11,
                      fontweight="semibold")
    axes[4].set_ylabel("Elevation (m)", color=COLORS["z"])
    axes[4].tick_params(axis="y", colors=COLORS["z"])
    axes[4].grid(True, color="#d9dee5", linewidth=0.7, alpha=0.85)
    axes[4].set_xlabel("Distance along pipe (km)")
    axes[4].set_xlim(x[0], x[-1])
    fig.suptitle("Thermal pipe with lumped wall: steady initializers",
                 fontsize=14, fontweight="semibold", y=0.995)
    fig.legend(handles, [h.get_label() for h in handles],
               loc="upper center", bbox_to_anchor=(0.5, 0.976), ncol=4,
               frameon=False)
    fig.tight_layout(rect=(0.04, 0.02, 0.75, 0.94), h_pad=1.0)
    fig.savefig(output, dpi=185, facecolor="white", bbox_inches="tight")
    plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="CSV from thermal profile exporter")
    parser.add_argument("-o", "--output", type=Path,
        default=Path("thermal_wall_initializer_profiles.png"))
    args = parser.parse_args()
    plot_profiles(read_profiles(args.csv), args.output)
    print(f"Wrote {args.output}")

if __name__ == "__main__":
    main()
