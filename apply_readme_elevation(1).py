#!/usr/bin/env python3
"""Set all README profile examples to the same 100 m elevation hill.

Run from anywhere, pointing to the local repository:
    python3 ~/Downloads/apply_readme_elevation.py \
        --repo ~/Sync/PipeSteadyStateInitializers
"""
import argparse
import re
from pathlib import Path


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument(
    "--repo", type=Path, default=Path.cwd(),
    help="path to the PipeSteadyStateInitializers checkout",
)
args = parser.parse_args()
repo = args.repo.expanduser().resolve()
if not (repo / "README.md").is_file() or not (repo / "src/main.cpp").is_file():
    parser.error(f"{repo} does not look like the repository root")


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    if new in text:
        return
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one matching block, found {count}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


def replace_paragraph_tail(path: Path, start: str, replacement: str) -> None:
    text = path.read_text(encoding="utf-8")
    if replacement in text:
        return
    pattern = re.compile(re.escape(start) + r"[^\n]*(?:\n(?!\n)[^\n]*)*")
    matches = list(pattern.finditer(text))
    if len(matches) != 1:
        raise RuntimeError(
            f"{path}: expected one paragraph tail beginning {start!r}, "
            f"found {len(matches)}"
        )
    match = matches[0]
    text = text[:match.start()] + replacement + text[match.end():]
    path.write_text(text, encoding="utf-8")


changes = {
    "src/main.cpp": (
        "        pipe_init::Configuration configuration;\n",
        "        pipe_init::Configuration configuration;\n"
        "        configuration.hill_height = 100.0; // 100 m cosine-shaped hill\n",
    ),
    "examples/export_initializer_profiles.cpp": (
        "        pipe_init::Configuration configuration;\n",
        "        pipe_init::Configuration configuration;\n"
        "        configuration.hill_height = 100.0; // Match the isothermal demo\n",
    ),
    "src/adiabatic_demo.cpp": (
        "        Configuration c;\n\n        auto start = Clock::now();",
        "        Configuration c;\n"
        "        c.hill_height = 100.0; // 100 m cosine-shaped hill\n\n"
        "        auto start = Clock::now();",
    ),
    "examples/export_adiabatic_profiles.cpp": (
        "        Configuration c;\n        const State marching = steady_spatial_marching(c);",
        "        Configuration c;\n"
        "        c.hill_height = 100.0; // Match the adiabatic demo\n"
        "        const State marching = steady_spatial_marching(c);",
    ),
    "src/thermal_wall_demo.cpp": (
        "        Configuration c;\n"
        "        auto t0=Clock::now(); const State analytical=analytical_approximation(c);",
        "        Configuration c;\n"
        "        c.hill_height = 100.0; // 100 m cosine-shaped hill\n"
        "        auto t0=Clock::now(); const State analytical=analytical_approximation(c);",
    ),
    "examples/export_thermal_wall_profiles.cpp": (
        "        Configuration c;\n        const State analytical=analytical_approximation(c);",
        "        Configuration c;\n"
        "        c.hill_height=100.0; // Match the thermal-wall demo\n"
        "        const State analytical=analytical_approximation(c);",
    ),
    "src/soil_shell_demo.cpp": (
        "        Configuration c;\n"
        "        auto a0=Clock::now(); const State analytical=analytical_approximation(c); auto a1=Clock::now();",
        "        Configuration c;\n"
        "        c.hill_height=100.0; // 100 m cosine-shaped hill\n"
        "        auto a0=Clock::now(); const State analytical=analytical_approximation(c); auto a1=Clock::now();",
    ),
    "examples/export_soil_shell_profiles.cpp": (
        "        Configuration c;\n"
        "        const State analytical=analytical_approximation(c),marching=steady_spatial_marching(c);",
        "        Configuration c;\n"
        "        c.hill_height=100.0; // Match the soil-shell demo\n"
        "        const State analytical=analytical_approximation(c),marching=steady_spatial_marching(c);",
    ),
}

for name, (old, new) in changes.items():
    replace_once(repo / name, old, new)

readme = repo / "README.md"
replace_once(
    readme,
    "| Fully thermal soil-shell | Gas, wall, and radial soil-shell temperatures with pressure and flow | Analytical approximation, spatial marching, coupled Newton, implicit theta relaxation |\n\n",
    "| Fully thermal soil-shell | Gas, wall, and radial soil-shell temperatures with pressure and flow | Analytical approximation, spatial marching, coupled Newton, implicit theta relaxation |\n\n"
    "All four profile examples use the same 100 m cosine-shaped hill over the 100 km pipe, with equal elevations at the two ends. Each model's `hill_height` setting controls the crest height.\n\n",
)
replace_once(
    readme,
    "![Isothermal pressure and mass-flow profiles for four initializers](initializer_profiles.png)",
    "![Isothermal pressure, mass-flow, and elevation profiles for four initializers](initializer_profiles.png)",
)
replace_once(
    readme,
    "![Thermal wall model profiles](thermal_wall_initializer_profiles.png)",
    "![Thermal wall pressure, flow, temperature, and elevation profiles](thermal_wall_initializer_profiles.png)",
)
replace_once(
    readme,
    "pseudo-time relaxation, but axial soil conduction is omitted. Default\nelevation is flat. The default nonlinear state has 133 unknowns.",
    "pseudo-time relaxation, but axial soil conduction is omitted. The README demo\nuses the shared 100 m cosine-shaped hill. The default nonlinear state has 133\nunknowns.",
)
replace_paragraph_tail(
    readme,
    "The bottom panel shows pipe",
    "The bottom panel shows the shared 100 m pipe-elevation hill.",
)

print("Set the isothermal, adiabatic, thermal-wall, and soil-shell README examples to a 100 m hill.")

