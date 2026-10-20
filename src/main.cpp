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
              << state.pressure.back()/1.0e6 << " MPa"
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

        const auto analytical =
            pipe_init::analytical_pressure_squared(configuration);
        const auto marching =
            pipe_init::discrete_steady_marching(configuration);
        const auto newton = pipe_init::newton_steady_state(
            configuration, marching);
        const pipe_init::RelaxationOptions relaxation_options;
        const auto relaxation = pipe_init::relax_to_steady_state(
            configuration, relaxation_options);

        std::cout << "Isothermal pipeline steady initializers\n"
                  << "cells=" << configuration.cells
                  << " length=" << configuration.length << " m"
                  << " hill=" << configuration.hill_height << " m"
                  << " flow=" << configuration.outlet_mass_flow << " kg/s\n";
        print_summary("analytical pressure-squared", configuration,
                      analytical,
                      pipe_init::scaled_steady_residual(configuration, analytical));
        print_summary("discrete steady marching", configuration, marching,
                      pipe_init::scaled_steady_residual(configuration, marching));
        print_summary("Newton, d/dt = 0", configuration, newton.state,
                      newton.scaled_residual);
        print_summary("pseudo-time relaxation", configuration, relaxation.state,
                      pipe_init::scaled_steady_residual(configuration,
                                                        relaxation.state));
        std::cout << "Newton iterations=" << newton.iterations
                  << ", relaxation theta=" << std::defaultfloat
                  << relaxation_options.theta
                  << ", relaxation steps=" << relaxation.steps
                  << ", relaxation dt=" << relaxation.dt << " s"
                  << ", pseudo-time=" << relaxation.pseudo_time << " s\n";
    } catch (const std::exception& error) {
        std::cerr << "pipe_init_demo: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
