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

#include "pipe/steady_initializers.hpp"
#include "pipe/adiabatic_initializers.hpp"
#include "pipe/thermal_wall_initializers.hpp"
#include "pipe/soil_shell_initializers.hpp"

#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <iostream>
#include <utility>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
}

template<class Configuration, class State>
std::pair<double, double> thermal_energy_residuals(
    const Configuration& c, const State& state)
{
    const auto elevation = c.elevation_nodes();
    const double cp = c.gamma*c.gas_constant/(c.gamma - 1.0);
    const double scale = std::max(1.0, std::abs(c.outlet_mass_flow)) *
                         cp*c.inlet_temperature;
    double with_gravity = 0.0;
    double without_gravity = 0.0;
    for (std::size_t i = 1; i <= c.cells; ++i) {
        const double left_flow = state.flow[i - 1];
        const double right_flow = i < c.cells
            ? state.flow[i] : c.outlet_mass_flow;
        const double length = i == c.cells ? c.dx()/2.0 : c.dx();
        const double enthalpy_flux = cp*(
            right_flow*state.gas_temperature[i] -
            left_flow*state.gas_temperature[i - 1]);
        const double gravity_flux = c.gravity*(
            right_flow*elevation[i] - left_flow*elevation[i - 1]);
        const double heat_in = c.inner_conductance_per_length()*(
            state.wall_temperature[i] - state.gas_temperature[i])*length;
        with_gravity = std::max(with_gravity,
            std::abs(enthalpy_flux + gravity_flux - heat_in)/scale);
        without_gravity = std::max(without_gravity,
            std::abs(enthalpy_flux - heat_in)/scale);
    }
    return {with_gravity, without_gravity};
}

pipe_init::Configuration small_case()
{
    pipe_init::Configuration c;
    c.cells = 24;
    c.length = 10'000.0;
    c.inlet_pressure = 5.0e6;
    c.outlet_mass_flow = 30.0;
    c.hill_height = 100.0;
    return c;
}

void analytical_flat_pipe_matches_pressure_squared_relation()
{
    auto c = small_case();
    c.hill_height = 0.0;
    const auto state = pipe_init::analytical_pressure_squared(c);
    require(state.pressure.size() == c.cells + 1,
            "pressure state must contain N+1 simulator nodes");
    require(state.flow.size() == c.cells,
            "flow state must contain N internal staggered faces");
    require(state.pressure.front() == c.inlet_pressure,
            "inlet pressure node must equal the prescribed boundary value");
    const double friction = pipe_init::swamee_jain_friction_factor(
        c, c.outlet_mass_flow);
    const double coefficient = friction*c.gas_constant_product()
        *c.outlet_mass_flow*c.outlet_mass_flow
        /(c.diameter*c.area()*c.area());
    const double expected_outlet_squared = c.inlet_pressure*c.inlet_pressure
                                           - coefficient*c.length;
    require(expected_outlet_squared > 0.0,
            "test case must have positive outlet pressure squared");
    const double expected_outlet = std::sqrt(expected_outlet_squared);
    require(std::abs(state.pressure.back() - expected_outlet)/expected_outlet
                < 2.0e-12,
            "flat-pipe analytical pressure-squared relation mismatch");
    require(pipe_init::scaled_steady_residual(c, state) < 1.0e-3,
            "analytic profile should closely satisfy the discrete steady residual");
}

void analytical_hydrostatic_limit()
{
    auto c = small_case();
    c.outlet_mass_flow = 0.0;
    const auto state = pipe_init::analytical_pressure_squared(c);
    const auto elevation = c.elevation_faces();
    for (std::size_t i = 0; i <= c.cells; ++i) {
        const double expected = c.inlet_pressure*std::exp(
            -c.gravity*elevation[i]/c.gas_constant_product());
        require(std::abs(state.pressure[i] - expected)/expected < 2.0e-12,
                "zero-flow analytical profile should recover isothermal hydrostatics");
    }
    require(std::abs(state.pressure.back() - c.inlet_pressure)/c.inlet_pressure
                < 2.0e-12,
            "equal endpoint elevations should recover equal endpoint pressures at zero flow");
}

void newton_solves_discrete_steady_residual()
{
    const auto c = small_case();
    auto guess = pipe_init::analytical_pressure_squared(c);
    for (std::size_t i = 1; i <= c.cells; ++i) {
        const double fraction = static_cast<double>(i)/
                                static_cast<double>(c.cells);
        guess.pressure[i] *= 1.0 - 0.015*fraction;
        if (i < c.cells) {
            guess.flow[i] *= 1.0 + 0.04*fraction;
        }
    }
    const auto solved = pipe_init::newton_steady_state(c, guess);
    require(solved.state.pressure.size() == c.cells + 1,
            "Newton state must retain N+1 pressure nodes");
    require(solved.state.pressure.front() == c.inlet_pressure,
            "Newton must keep the inlet pressure fixed");
    require(solved.scaled_residual < 1.0e-10,
            "Newton steady residual did not converge tightly");
    require(pipe_init::scaled_steady_residual(c, solved.state) < 1.0e-10,
            "reported Newton state fails independent steady residual check");
    for (double flow : solved.state.flow) {
        require(std::abs(flow - c.outlet_mass_flow) < 1.0e-7,
                "steady continuity should make every segment flow equal outlet flow");
    }
}

void discrete_marching_matches_newton_and_has_small_residual()
{
    const auto c = small_case();
    const auto marched = pipe_init::discrete_steady_marching(c);
    const auto solved = pipe_init::newton_steady_state(
        c, pipe_init::analytical_pressure_squared(c), 1.0e-13);
    require(pipe_init::scaled_steady_residual(c, marched) < 2.0e-15,
            "discrete marching should satisfy each discrete steady balance");
    for (std::size_t i = 0; i <= c.cells; ++i) {
        require(std::abs(marched.pressure[i] - solved.state.pressure[i])
                    / c.inlet_pressure < 2.0e-12,
                "discrete marching and Newton pressure profiles disagree");
    }
    for (const double flow : marched.flow) {
        require(flow == c.outlet_mass_flow,
                "discrete marching should initialize constant steady flow");
    }

    const auto analytic = pipe_init::analytical_pressure_squared(c);
    require(pipe_init::scaled_steady_residual(c, analytic) >
                pipe_init::scaled_steady_residual(c, marched),
            "continuum analytical profile should retain a discretization residual");
}

void relaxation_matches_newton_steady_state()
{
    const auto c = small_case();
    const auto newton = pipe_init::newton_steady_state(
        c, pipe_init::analytical_pressure_squared(c));
    pipe_init::RelaxationOptions options;
    options.theta = 1.0;
    options.tolerance = 2.0e-8;
    options.max_pseudo_time = 200'000.0;
    options.max_steps = 500'000;
    const auto relaxed = pipe_init::relax_to_steady_state(c, options);
    require(relaxed.scaled_rhs <= options.tolerance,
            "relaxation stopped before reaching its derivative tolerance");
    require(pipe_init::scaled_steady_residual(c, relaxed.state) < 1.0e-6,
            "relaxed state fails steady residual check");
    for (std::size_t i = 0; i <= c.cells; ++i) {
        const double relative = std::abs(relaxed.state.pressure[i]
                                         - newton.state.pressure[i])
                                / c.inlet_pressure;
        require(relative < 2.0e-5,
                "relaxation and Newton pressure profiles disagree");
    }
}

void trapezoidal_relaxation_matches_newton_steady_state()
{
    const auto c = small_case();
    const auto newton = pipe_init::newton_steady_state(
        c, pipe_init::analytical_pressure_squared(c));
    pipe_init::RelaxationOptions options;
    options.theta = 0.5;
    options.tolerance = 2.0e-8;
    options.max_pseudo_time = 200'000.0;
    options.max_steps = 500'000;
    const auto relaxed = pipe_init::relax_to_steady_state(c, options);
    require(relaxed.scaled_rhs <= options.tolerance,
            "theta=0.5 relaxation stopped before reaching derivative tolerance");
    require(pipe_init::scaled_steady_residual(c, relaxed.state) < 1.0e-6,
            "theta=0.5 relaxed state fails steady residual check");
    for (std::size_t i = 0; i <= c.cells; ++i) {
        const double relative = std::abs(relaxed.state.pressure[i]
                                         - newton.state.pressure[i])
                                / c.inlet_pressure;
        require(relative < 2.0e-5,
                "theta=0.5 relaxation and Newton pressure profiles disagree");
    }
}

void adiabatic_marching_and_newton_agree()
{
    using namespace pipe_init::adiabatic;
    Configuration c;
    c.cells = 8;
    c.length = 2'000.0;
    c.inlet_pressure = 6.0e6;
    c.inlet_temperature = 300.0;
    c.outlet_mass_flow = 35.0;
    c.hill_height = 40.0;
    const State marched = steady_spatial_marching(c);
    const SolveInfo solved = coupled_steady_newton(c, marched, 1.0e-12);
    require(scaled_steady_residual(c, marched) < 2.0e-12,
            "adiabatic marching should satisfy the discrete steady residual");
    require(solved.scaled_residual < 2.0e-12,
            "adiabatic coupled Newton did not converge tightly");
    for (std::size_t i = 0; i <= c.cells; ++i) {
        require(std::abs(marched.pressure[i] - solved.state.pressure[i]) /
                    c.inlet_pressure < 2.0e-11,
                "adiabatic marching and Newton pressure profiles disagree");
        require(std::abs(marched.temperature[i] - solved.state.temperature[i]) /
                    c.inlet_temperature < 2.0e-11,
                "adiabatic marching and Newton temperature profiles disagree");
    }
}

void adiabatic_theta_relaxation_converges()
{
    using namespace pipe_init::adiabatic;
    Configuration c;
    c.cells = 6;
    c.length = 1'000.0;
    c.inlet_pressure = 6.0e6;
    c.inlet_temperature = 300.0;
    c.outlet_mass_flow = 25.0;
    c.hill_height = 10.0;
    RelaxationOptions options;
    options.theta = 1.0;
    options.pseudo_courant = 5.0;
    options.tolerance = 2.0e-7;
    options.max_pseudo_time = 1.0e6;
    options.max_steps = 200'000;
    const auto relaxed = implicit_theta_relaxation(c, options);
    require(relaxed.scaled_residual <= options.tolerance,
            "adiabatic theta relaxation did not meet its residual tolerance");
    require(scaled_steady_residual(c, relaxed.state) < 5.0e-7,
            "adiabatic theta-relaxed state fails the steady residual check");
}

void thermal_models_include_gravitational_potential_in_energy()
{
    pipe_init::thermal_wall::Configuration wall;
    wall.hill_height = 100.0;
    const auto wall_hill = pipe_init::thermal_wall::steady_spatial_marching(wall);
    const auto wall_residuals = thermal_energy_residuals(wall, wall_hill);
    require(wall_residuals.first < 2.0e-12,
            "thermal-wall marching must satisfy the gravity-inclusive energy balance");
    require(wall_residuals.second > 1.0e-5,
            "omitting gravity must produce a measurable thermal-wall energy residual");
    wall.hill_height = 0.0;
    const auto wall_flat = pipe_init::thermal_wall::steady_spatial_marching(wall);
    const std::size_t wall_crest = wall.cells/2;
    require(wall_hill.gas_temperature[wall_crest] <
                wall_flat.gas_temperature[wall_crest]-0.05,
            "thermal-wall hill should convert gravitational potential to lower gas enthalpy");

    pipe_init::soil_shell::Configuration soil;
    soil.hill_height = 100.0;
    const auto soil_hill = pipe_init::soil_shell::steady_spatial_marching(soil);
    const auto soil_residuals = thermal_energy_residuals(soil, soil_hill);
    require(soil_residuals.first < 2.0e-12,
            "soil-shell marching must satisfy the gravity-inclusive energy balance");
    require(soil_residuals.second > 1.0e-5,
            "omitting gravity must produce a measurable soil-shell energy residual");
    soil.hill_height = 0.0;
    const auto soil_flat = pipe_init::soil_shell::steady_spatial_marching(soil);
    const std::size_t soil_crest = soil.cells/2;
    require(soil_hill.gas_temperature[soil_crest] <
                soil_flat.gas_temperature[soil_crest]-0.05,
            "soil-shell hill should convert gravitational potential to lower gas enthalpy");
}

void thermal_wall_marching_newton_and_theta_agree()
{
    using namespace pipe_init::thermal_wall;
    Configuration c;
    c.cells = 5;
    c.length = 5'000.0;
    c.hill_height = 20.0;
    const State marched = steady_spatial_marching(c);
    require(scaled_steady_residual(c, marched) < 2.0e-12,
            "thermal-wall marching should satisfy discrete steady residual");
    const SolveInfo solved = coupled_steady_newton(c, analytical_approximation(c));
    require(solved.scaled_residual < 2.0e-10,
            "thermal-wall Newton did not converge");
    for (std::size_t i = 0; i <= c.cells; ++i) {
        require(std::abs(marched.pressure[i] - solved.state.pressure[i]) /
                    c.inlet_pressure < 2.0e-10,
                "thermal-wall pressure profiles disagree");
        require(std::abs(marched.gas_temperature[i] -
                         solved.state.gas_temperature[i]) / c.inlet_temperature < 2.0e-10,
                "thermal-wall gas temperatures disagree");
        require(std::abs(marched.wall_temperature[i] -
                         solved.state.wall_temperature[i]) / c.inlet_temperature < 2.0e-10,
                "thermal-wall temperatures disagree");
    }
    RelaxationOptions options;
    options.theta = 1.0;
    options.pseudo_courant = 10.0;
    options.tolerance = 2.0e-7;
    options.max_pseudo_time = 1.0e6;
    const auto relaxed = implicit_theta_relaxation(c, options);
    require(relaxed.scaled_residual <= options.tolerance,
            "thermal-wall theta relaxation missed residual tolerance");
    for (std::size_t i = 0; i <= c.cells; ++i)
        require(std::abs(relaxed.state.gas_temperature[i] -
                         marched.gas_temperature[i]) / c.inlet_temperature < 2.0e-5,
                "thermal-wall theta/marching temperatures disagree");
}

void soil_shell_marching_newton_and_theta_agree()
{
    using namespace pipe_init::soil_shell;
    Configuration c;
    const State marched=steady_spatial_marching(c);
    require(scaled_steady_residual(c,marched)<2.0e-12,
            "soil-shell marching should satisfy the discrete steady residual");
    const SolveInfo solved=coupled_steady_newton(c,analytical_approximation(c),1.0e-11,35);
    require(solved.scaled_residual<2.0e-10,"soil-shell coupled Newton did not converge");
    for(std::size_t i=0;i<=c.cells;++i) {
        require(std::abs(marched.pressure[i]-solved.state.pressure[i])/c.inlet_pressure<2.0e-9,
                "soil-shell marching/Newton pressure profiles disagree");
        require(std::abs(marched.gas_temperature[i]-solved.state.gas_temperature[i])/c.inlet_temperature<2.0e-9,
                "soil-shell marching/Newton gas profiles disagree");
        require(std::abs(marched.wall_temperature[i]-solved.state.wall_temperature[i])/c.inlet_temperature<2.0e-9,
                "soil-shell marching/Newton wall profiles disagree");
        for(std::size_t j=0;j<c.soil_shell_radii.size();++j)
            require(std::abs(marched.soil_temperature[j][i]-solved.state.soil_temperature[j][i])/c.inlet_temperature<2.0e-9,
                    "soil-shell marching/Newton soil profiles disagree");
    }
    RelaxationOptions options; options.theta=1.0; options.pseudo_courant=0.25;
    options.tolerance=3.0e-7; options.max_pseudo_time=5.0e7; options.max_steps=100;
    const RelaxationInfo relaxed=implicit_theta_relaxation(c,options);
    require(relaxed.scaled_residual<=options.tolerance,
            "soil-shell theta relaxation missed residual tolerance");
    for(std::size_t i=0;i<=c.cells;++i) {
        require(std::abs(relaxed.state.gas_temperature[i]-marched.gas_temperature[i])/c.inlet_temperature<2.0e-5,
                "soil-shell theta/marching gas temperatures disagree");
        for(std::size_t j=0;j<c.soil_shell_radii.size();++j)
            require(std::abs(relaxed.state.soil_temperature[j][i]-marched.soil_temperature[j][i])/c.inlet_temperature<2.0e-5,
                    "soil-shell theta/marching soil temperatures disagree");
    }
}

} // namespace

int main()
{
    try {
        analytical_flat_pipe_matches_pressure_squared_relation();
        analytical_hydrostatic_limit();
        newton_solves_discrete_steady_residual();
        discrete_marching_matches_newton_and_has_small_residual();
        relaxation_matches_newton_steady_state();
        trapezoidal_relaxation_matches_newton_steady_state();
        adiabatic_marching_and_newton_agree();
        adiabatic_theta_relaxation_converges();
        thermal_models_include_gravitational_potential_in_energy();
        thermal_wall_marching_newton_and_theta_agree();
        soil_shell_marching_newton_and_theta_agree();
        std::cout << "PASS: all steady-initializer checks\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
