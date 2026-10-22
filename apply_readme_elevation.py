#!/usr/bin/env python3
"""Set all README profile examples to the same 100 m elevation hill.

Run from the PipeSteadyStateInitializers repository root:
    python3 apply_readme_elevation.py
"""
from pathlib import Path


def replace_once(path: Path, old: str, new: str) -> None:
    text = path.read_text(encoding="utf-8")
    count = text.count(old)
    if count != 1:
        raise RuntimeError(f"{path}: expected one matching block, found {count}")
    path.write_text(text.replace(old, new, 1), encoding="utf-8")


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
    replace_once(Path(name), old, new)

readme = Path("README.md")
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
replace_once(
    readme,
    "far-field temperature toward the outer radius. The bottom panel shows pipe\nelevation. With the default flat pipe it is a zero line; set `hill_height` in the\nconfiguration to compare elevation effects.",
    "far-field temperature toward the outer radius. The bottom panel shows the\nshared 100 m pipe-elevation hill.",
)

print("Set the isothermal, adiabatic, thermal-wall, and soil-shell README examples to a 100 m hill.")

