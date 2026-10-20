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

namespace pipe_init::thermal_wall {

// Ideal-gas gas/wall model on the simulator-aligned nodal/staggered grid.
// Soil temperature is prescribed; the wall is locally lumped through its
// thickness and has no axial-conduction term.
struct Configuration {
    std::size_t cells{20};
    double length{100'000.0};
    double diameter{0.6};
    double roughness{4.5e-5};
    double viscosity{1.1e-5};
    double gas_constant{518.3};
    double gamma{1.30};
    double inlet_pressure{8.0e6};
    double inlet_temperature{288.15};
    double soil_temperature{278.15};
    double outlet_mass_flow{100.0};
    double gravity{9.80665};
    double hill_height{0.0};
    double inner_heat_transfer_coefficient{10.0}; // W/(m^2 K)
    double soil_conductance_per_length{3.0}; // W/(m K), includes soil/insulation
    double wall_density{7850.0}; // kg/m^3
    double wall_heat_capacity{500.0}; // J/(kg K)
    double wall_thickness{0.0127};

    [[nodiscard]] double area() const;
    [[nodiscard]] double dx() const;
    [[nodiscard]] double inner_conductance_per_length() const;
    [[nodiscard]] double wall_heat_capacity_per_length() const;
    [[nodiscard]] std::vector<double> elevation_nodes() const;
    void validate() const;
};

struct State {
    std::vector<double> pressure;    // N+1 nodes
    std::vector<double> gas_temperature; // N+1 nodes
    std::vector<double> wall_temperature; // N+1 nodes
    std::vector<double> flow;        // N staggered interior faces
};

struct SolveInfo {
    State state;
    std::size_t iterations{};
    double scaled_residual{};
};

struct RelaxationOptions {
    double tolerance{1.0e-8};
    double theta{0.9};
    double pseudo_courant{20.0};
    double max_pseudo_time{2.0e6};
    std::size_t max_steps{20'000};
    std::size_t consecutive_converged_steps{4};
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

// Closed-form constant-property thermal profile plus quadrature for pressure.
// This is an approximation; its discrete residual is reported by the demo.
[[nodiscard]] State analytical_approximation(const Configuration& configuration);

// Marches the same discrete steady momentum and gas-energy balances used by
// the coupled solve. The steady wall equation is eliminated locally.
[[nodiscard]] State steady_spatial_marching(const Configuration& configuration);

// Solves continuity, momentum, gas energy, and local wall equilibrium together.
[[nodiscard]] SolveInfo coupled_steady_newton(
    const Configuration& configuration, const State& initial_guess,
    double tolerance = 1.0e-10, std::size_t max_iterations = 30);

// Advances the semidiscrete gas/wall transient with an implicit theta step
// under fixed inlet pressure/temperature and outlet mass flow.
[[nodiscard]] RelaxationInfo implicit_theta_relaxation(
    const Configuration& configuration,
    const RelaxationOptions& options = {});

[[nodiscard]] double scaled_steady_residual(
    const Configuration& configuration, const State& state);

} // namespace pipe_init::thermal_wall
