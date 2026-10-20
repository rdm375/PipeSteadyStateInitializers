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

#include "pipe/adiabatic_initializers.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>

namespace {
using Clock = std::chrono::steady_clock;

double seconds(Clock::time_point a, Clock::time_point b)
{
    return std::chrono::duration<double>(b - a).count();
}

void report(const char* name, const pipe_init::adiabatic::Configuration& c,
            const pipe_init::adiabatic::State& s, double residual,
            double elapsed)
{
    double max_flow_error = 0.0;
    for (double m : s.flow) {
        max_flow_error = std::max(max_flow_error, std::abs(m-c.outlet_mass_flow));
    }
    std::cout << std::left << std::setw(30) << name
              << " p_out=" << std::fixed << std::setprecision(6)
              << s.pressure.back()/1.0e6 << " MPa"
              << " T_out=" << s.temperature.back() << " K"
              << " max_flow_error=" << std::scientific << max_flow_error
              << " kg/s residual=" << residual
              << " time=" << std::fixed << std::setprecision(6)
              << elapsed << " s\n";
}

void compare_to_marching(const char* name,
                         const pipe_init::adiabatic::Configuration& c,
                         const pipe_init::adiabatic::State& state,
                         const pipe_init::adiabatic::State& reference)
{
    double p_error = 0.0;
    double t_error = 0.0;
    double m_error = 0.0;
    for (std::size_t i = 0; i <= c.cells; ++i) {
        p_error = std::max(p_error, std::abs(state.pressure[i] - reference.pressure[i])
                                   / c.inlet_pressure);
        t_error = std::max(t_error, std::abs(state.temperature[i] -
                                             reference.temperature[i])
                                   / c.inlet_temperature);
    }
    for (std::size_t i = 0; i < c.cells; ++i) {
        m_error = std::max(m_error, std::abs(state.flow[i] - reference.flow[i])
                                   / c.outlet_mass_flow);
    }
    std::cout << "  " << name << " vs marching: max scaled errors p="
              << std::scientific << p_error << ", T=" << t_error
              << ", m=" << m_error << '\n';
}
} // namespace

int main()
{
    try {
        using namespace pipe_init::adiabatic;
        Configuration c;

        auto start = Clock::now();
        const State marched = steady_spatial_marching(c);
        auto finish = Clock::now();
        const double marching_time = seconds(start, finish);

        start = Clock::now();
        State flat_guess;
        flat_guess.pressure.assign(c.cells + 1, c.inlet_pressure);
        flat_guess.temperature.assign(c.cells + 1, c.inlet_temperature);
        flat_guess.flow.assign(c.cells, c.outlet_mass_flow);
        const SolveInfo newton = coupled_steady_newton(c, flat_guess);
        finish = Clock::now();
        const double newton_time = seconds(start, finish);

        RelaxationOptions options;
        options.theta = 0.9;
        options.pseudo_courant = 10.0;
        start = Clock::now();
        const RelaxationInfo relaxed = implicit_theta_relaxation(c, options);
        finish = Clock::now();
        const double relaxation_time = seconds(start, finish);

        std::cout << "Adiabatic ideal-gas pipe steady initializers\n"
                  << "cells=" << c.cells << " length=" << c.length << " m"
                  << " D=" << c.diameter << " m"
                  << " inlet=" << c.inlet_pressure/1.0e6 << " MPa, "
                  << c.inlet_temperature << " K"
                  << " outlet flow=" << c.outlet_mass_flow << " kg/s\n"
                  << "gamma=" << c.gamma << " R=" << c.gas_constant
                  << " J/(kg K); adiabatic, constant-area, forward flow\n";
        report("Steady spatial marching", c, marched,
               scaled_steady_residual(c, marched), marching_time);
        report("Coupled steady Newton", c, newton.state,
               newton.scaled_residual, newton_time);
        compare_to_marching("Newton", c, newton.state, marched);
        report("Implicit theta relaxation", c, relaxed.state,
               relaxed.scaled_residual, relaxation_time);
        compare_to_marching("Theta relaxation", c, relaxed.state, marched);
        std::cout << "Newton iterations=" << newton.iterations
                  << ", relaxation theta=" << options.theta
                  << ", relaxation steps=" << relaxed.steps
                  << ", pseudo dt=" << relaxed.dt << " s"
                  << ", pseudo-time=" << relaxed.pseudo_time << " s"
                  << ", inner Newton iterations=" << relaxed.newton_iterations
                  << '\n';
    } catch (const std::exception& error) {
        std::cerr << "adiabatic_pipe_init_demo: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
