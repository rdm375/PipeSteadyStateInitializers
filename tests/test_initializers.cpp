#include "pipe/steady_initializers.hpp"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

void require(bool condition, const std::string& message)
{
    if (!condition) {
        throw std::runtime_error(message);
    }
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
    require(std::abs(state.outlet_pressure - expected_outlet)/expected_outlet
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
    for (std::size_t i = 0; i < c.cells; ++i) {
        const double h_center = 0.5*(elevation[i] + elevation[i + 1]);
        const double expected = c.inlet_pressure*std::exp(
            -c.gravity*h_center/c.gas_constant_product());
        require(std::abs(state.pressure[i] - expected)/expected < 2.0e-12,
                "zero-flow analytical profile should recover isothermal hydrostatics");
    }
    require(std::abs(state.outlet_pressure - c.inlet_pressure)/c.inlet_pressure
                < 2.0e-12,
            "equal endpoint elevations should recover equal endpoint pressures at zero flow");
}

void newton_solves_discrete_steady_residual()
{
    const auto c = small_case();
    auto guess = pipe_init::analytical_pressure_squared(c);
    for (std::size_t i = 0; i < c.cells; ++i) {
        const double fraction = static_cast<double>(i + 1)/
                                static_cast<double>(c.cells);
        guess.pressure[i] *= 1.0 - 0.015*fraction;
        guess.flow[i] *= 1.0 + 0.04*fraction;
    }
    const auto solved = pipe_init::newton_steady_state(c, guess);
    require(solved.scaled_residual < 1.0e-10,
            "Newton steady residual did not converge tightly");
    require(pipe_init::scaled_steady_residual(c, solved.state) < 1.0e-10,
            "reported Newton state fails independent steady residual check");
    for (double flow : solved.state.flow) {
        require(std::abs(flow - c.outlet_mass_flow) < 1.0e-7,
                "steady continuity should make every segment flow equal outlet flow");
    }
}

void relaxation_matches_newton_steady_state()
{
    const auto c = small_case();
    const auto newton = pipe_init::newton_steady_state(
        c, pipe_init::analytical_pressure_squared(c));
    pipe_init::RelaxationOptions options;
    options.tolerance = 2.0e-8;
    options.max_pseudo_time = 200'000.0;
    options.max_steps = 500'000;
    const auto relaxed = pipe_init::relax_to_steady_state(c, options);
    require(relaxed.scaled_rhs <= options.tolerance,
            "relaxation stopped before reaching its derivative tolerance");
    require(pipe_init::scaled_steady_residual(c, relaxed.state) < 1.0e-6,
            "relaxed state fails steady residual check");
    for (std::size_t i = 0; i < c.cells; ++i) {
        const double relative = std::abs(relaxed.state.pressure[i]
                                         - newton.state.pressure[i])
                                / c.inlet_pressure;
        require(relative < 2.0e-5,
                "relaxation and Newton pressure profiles disagree");
    }
}

} // namespace

int main()
{
    try {
        analytical_flat_pipe_matches_pressure_squared_relation();
        analytical_hydrostatic_limit();
        newton_solves_discrete_steady_residual();
        relaxation_matches_newton_steady_state();
        std::cout << "PASS: all isothermal steady-initializer checks\n";
    } catch (const std::exception& error) {
        std::cerr << "FAIL: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
