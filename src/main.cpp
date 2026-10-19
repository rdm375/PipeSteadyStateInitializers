#include "pipe/steady_initializers.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {

void print_summary(const char* label, const pipe_init::Configuration& c,
                   const pipe_init::State& state, double residual)
{
    std::cout << std::left << std::setw(25) << label
              << " p_out=" << std::fixed << std::setprecision(3)
              << state.outlet_pressure/1.0e6 << " MPa"
              << " max|m-mdot_out|=";
    double flow_error = 0.0;
    for (double value : state.flow) {
        flow_error = std::max(flow_error,
                              std::abs(value - c.outlet_mass_flow));
    }
    std::cout << std::scientific << flow_error
              << " kg/s scaled residual=" << residual << '\n';
}

} // namespace

int main()
{
    try {
        pipe_init::Configuration configuration;
        configuration.cells = 32;
        configuration.length = 20'000.0;
        configuration.inlet_pressure = 5.0e6;
        configuration.outlet_mass_flow = 50.0;
        configuration.hill_height = 150.0;

        const auto analytical =
            pipe_init::analytical_pressure_squared(configuration);
        const auto newton = pipe_init::newton_steady_state(
            configuration, analytical);
        const auto relaxation = pipe_init::relax_to_steady_state(
            configuration);

        std::cout << "Isothermal pipeline steady initializers\n"
                  << "cells=" << configuration.cells
                  << " length=" << configuration.length << " m"
                  << " hill=" << configuration.hill_height << " m"
                  << " flow=" << configuration.outlet_mass_flow << " kg/s\n";
        print_summary("analytical pressure-squared", configuration,
                      analytical,
                      pipe_init::scaled_steady_residual(configuration, analytical));
        print_summary("Newton, d/dt = 0", configuration, newton.state,
                      newton.scaled_residual);
        print_summary("pseudo-time relaxation", configuration, relaxation.state,
                      pipe_init::scaled_steady_residual(configuration,
                                                        relaxation.state));
        std::cout << "Newton iterations=" << newton.iterations
                  << ", relaxation steps=" << relaxation.steps
                  << ", pseudo-time=" << relaxation.pseudo_time << " s\n";
    } catch (const std::exception& error) {
        std::cerr << "pipe_init_demo: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
