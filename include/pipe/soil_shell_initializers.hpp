/*
 * Copyright (C) 2026 Richard Myers
 * SPDX-License-Identifier: AGPL-3.0-only
 *
 * Licensed under the GNU Affero General Public License, version 3 only.
 * See the LICENSE file for the complete terms.
 *
 * A separate commercial license for proprietary use may be available by
 * written agreement with the relevant rights holders. This notice does not
 * grant that license, guarantee its availability, or set fees, royalties,
 * or revenue shares.
 */

#pragma once

#include <cstddef>
#include <vector>

namespace pipe_init::soil_shell {

// Gas, locally lumped pipe wall, and concentric radial soil shells. Soil has
// radial conduction only; the far-field temperature is prescribed.
struct Configuration {
    std::size_t cells{16};
    double length{100'000.0};
    double diameter{0.6};
    double roughness{4.5e-5};
    double viscosity{1.1e-5};
    double gas_constant{518.3};
    double gamma{1.30};
    double inlet_pressure{8.0e6};
    double inlet_temperature{288.15};
    double far_field_temperature{278.15};
    double outlet_mass_flow{100.0};
    double gravity{9.80665};
    double hill_height{0.0};
    double inner_heat_transfer_coefficient{10.0};
    double wall_density{7850.0};
    double wall_heat_capacity{500.0};
    double wall_thickness{0.0127};
    double soil_conductivity{1.5};
    double soil_density{1800.0};
    double soil_heat_capacity{800.0};
    std::vector<double> soil_shell_radii{0.40, 0.58, 0.84, 1.20};
    double far_field_radius{2.0};

    [[nodiscard]] double area() const;
    [[nodiscard]] double dx() const;
    [[nodiscard]] double wall_outer_radius() const;
    [[nodiscard]] double inner_conductance_per_length() const;
    [[nodiscard]] double wall_heat_capacity_per_length() const;
    [[nodiscard]] std::vector<double> elevation_nodes() const;
    void validate() const;
};

struct State {
    std::vector<double> pressure;
    std::vector<double> gas_temperature;
    std::vector<double> wall_temperature;
    std::vector<std::vector<double>> soil_temperature; // shell-major, N+1 nodes
    std::vector<double> flow; // N staggered interior faces
};

struct SolveInfo {
    State state;
    std::size_t iterations{};
    double scaled_residual{};
};

struct RelaxationOptions {
    double tolerance{1.0e-8};
    double theta{1.0};
    double pseudo_courant{0.25};
    double max_pseudo_time{5.0e7};
    std::size_t max_steps{500};
    std::size_t consecutive_converged_steps{3};
};

struct RelaxationInfo {
    State state;
    std::size_t steps{};
    std::size_t newton_iterations{};
    double pseudo_time{};
    double dt{};
    double scaled_residual{};
};

[[nodiscard]] double swamee_jain_friction_factor(
    const Configuration& configuration, double mass_flow);
[[nodiscard]] State analytical_approximation(const Configuration& configuration);
[[nodiscard]] State steady_spatial_marching(const Configuration& configuration);
[[nodiscard]] SolveInfo coupled_steady_newton(
    const Configuration& configuration, const State& initial_guess,
    double tolerance = 1.0e-10, std::size_t max_iterations = 30);
[[nodiscard]] RelaxationInfo implicit_theta_relaxation(
    const Configuration& configuration, const RelaxationOptions& options = {});
[[nodiscard]] double scaled_steady_residual(
    const Configuration& configuration, const State& state);

} // namespace pipe_init::soil_shell
