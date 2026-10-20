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

namespace pipe_init {

struct Configuration {
    std::size_t cells{100};
    double length{100'000.0};
    double diameter{0.6};
    double roughness{4.5e-5};
    double viscosity{1.1e-5};
    double compressibility{0.9};
    double gas_constant{518.3};
    double temperature{288.15};
    double inlet_pressure{8.0e6};
    double outlet_mass_flow{100.0};
    double gravity{9.80665};
    double hill_height{0.0};

    [[nodiscard]] double area() const;
    [[nodiscard]] double dx() const;
    [[nodiscard]] double pressure_scale() const;
    [[nodiscard]] double mass_flow_scale() const;
    [[nodiscard]] double gas_constant_product() const;
    [[nodiscard]] std::vector<double> elevation_faces() const;
    void validate() const;
};

struct State {
    // Pressure is stored at N+1 nodes; flow at the N interior faces.
    // pressure[i] is at x=i*dx; flow[i] is at x=(i+1/2)*dx.
    std::vector<double> pressure;
    std::vector<double> flow;
};

struct SolveInfo {
    State state;
    std::size_t iterations{};
    double scaled_residual{};
};

struct RelaxationOptions {
    double tolerance{1.0e-9};
    double max_pseudo_time{5.0e6};
    std::size_t max_steps{2'000'000};
    std::size_t consecutive_converged_steps{8};
    double theta{0.5}; // theta in [0.5, 1]; 0.5 is trapezoidal.
    double pseudo_courant{10.0}; // dt = pseudo_courant * dx / sqrt(Z R_s T).
};

struct RelaxationInfo {
    State state;
    std::size_t steps{};
    double pseudo_time{};
    double scaled_rhs{};
    double dt{};
    std::size_t newton_iterations{};
};

[[nodiscard]] double swamee_jain_friction_factor(
    const Configuration& configuration,
    double mass_flow);

[[nodiscard]] State analytical_pressure_squared(
    const Configuration& configuration);

[[nodiscard]] State discrete_steady_marching(
    const Configuration& configuration);

[[nodiscard]] SolveInfo newton_steady_state(
    const Configuration& configuration,
    const State& initial_guess,
    double tolerance = 1.0e-10,
    std::size_t max_iterations = 30);

[[nodiscard]] RelaxationInfo relax_to_steady_state(
    const Configuration& configuration,
    const RelaxationOptions& options = {});

[[nodiscard]] double scaled_steady_residual(
    const Configuration& configuration,
    const State& state);

} // namespace pipe_init
