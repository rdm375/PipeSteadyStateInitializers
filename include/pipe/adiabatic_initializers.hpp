#pragma once

#include <cstddef>
#include <vector>

namespace pipe_init::adiabatic {

// Calorically perfect ideal gas in a constant-area pipe. Pressure and
// temperature are nodal (N+1); mass flow is stored at N staggered faces.
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
    double outlet_mass_flow{100.0};
    double gravity{9.80665};
    double hill_height{0.0};

    [[nodiscard]] double area() const;
    [[nodiscard]] double dx() const;
    [[nodiscard]] std::vector<double> elevation_nodes() const;
    void validate() const;
};

struct State {
    std::vector<double> pressure;
    std::vector<double> temperature;
    std::vector<double> flow;
};

struct SolveInfo {
    State state;
    std::size_t iterations{};
    double scaled_residual{};
};

struct RelaxationOptions {
    double tolerance{1.0e-8};
    double theta{0.75};
    double pseudo_courant{1.0};
    double max_pseudo_time{2.0e6};
    std::size_t max_steps{200'000};
    std::size_t consecutive_converged_steps{4};
};

struct RelaxationInfo {
    State state;
    std::size_t steps{};
    double pseudo_time{};
    double scaled_residual{};
    double dt{};
    std::size_t newton_iterations{};
};

[[nodiscard]] double swamee_jain_friction_factor(
    const Configuration& configuration,
    double mass_flow);

// Solve each cell's discrete momentum and total-enthalpy equations in order.
[[nodiscard]] State steady_spatial_marching(
    const Configuration& configuration);

// Solve the coupled mass, momentum, energy, and ideal-gas equations at once.
[[nodiscard]] SolveInfo coupled_steady_newton(
    const Configuration& configuration,
    const State& initial_guess,
    double tolerance = 1.0e-10,
    std::size_t max_iterations = 30);

// Implicit theta relaxation of the semidiscrete adiabatic pipe transient.
[[nodiscard]] RelaxationInfo implicit_theta_relaxation(
    const Configuration& configuration,
    const RelaxationOptions& options = {});

[[nodiscard]] double scaled_steady_residual(
    const Configuration& configuration,
    const State& state);

} // namespace pipe_init::adiabatic
