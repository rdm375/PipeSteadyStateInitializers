#!/usr/bin/env bash
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

set -euo pipefail

repo_root="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel --target \
  pipe_init_profile_export \
  adiabatic_pipe_init_profile_export \
  thermal_wall_profile_export \
  soil_shell_profile_export

mkdir -p results
./build/pipe_init_profile_export > initializer_profiles.csv
./build/adiabatic_pipe_init_profile_export > adiabatic_profiles.csv
./build/thermal_wall_profile_export > thermal_wall_profiles.csv
./build/soil_shell_profile_export > results/soil_shell_profiles.csv

python3 scripts/recover_isothermal_elevation.py
python3 scripts/recover_adiabatic_elevation.py adiabatic_profiles.csv \
  --output results/adiabatic_elevation_recovery.csv \
  --plot results/adiabatic_elevation_recovery.png \
  --increment-plot results/adiabatic_elevation_increment_error.png
python3 scripts/recover_thermal_hydraulic_elevation.py
python3 scripts/validate_energy_gravity.py
