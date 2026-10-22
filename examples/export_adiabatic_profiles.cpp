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

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main()
{
    try {
        using namespace pipe_init::adiabatic;
        Configuration c;
        c.hill_height = 100.0; // Match the adiabatic demo
        const State marching = steady_spatial_marching(c);

        State initial_guess;
        initial_guess.pressure.assign(c.cells + 1, c.inlet_pressure);
        initial_guess.temperature.assign(c.cells + 1, c.inlet_temperature);
        initial_guess.flow.assign(c.cells, c.outlet_mass_flow);
        const State newton = coupled_steady_newton(c, initial_guess).state;

        RelaxationOptions options;
        options.theta = 0.9;
        options.pseudo_courant = 10.0;
        const State relaxation = implicit_theta_relaxation(c, options).state;
        const auto elevation = c.elevation_nodes();

        std::cout << "x_node_m,elevation_m,p_marching_Pa,T_marching_K,"
                     "p_newton_Pa,T_newton_K,p_relaxation_Pa,T_relaxation_K,"
                     "x_face_m,m_marching_kg_s,m_newton_kg_s,m_relaxation_kg_s\n";
        std::cout << std::setprecision(15);
        for (std::size_t i = 0; i <= c.cells; ++i) {
            std::cout << static_cast<double>(i)*c.dx() << ',' << elevation[i]
                      << ',' << marching.pressure[i]
                      << ',' << marching.temperature[i]
                      << ',' << newton.pressure[i]
                      << ',' << newton.temperature[i]
                      << ',' << relaxation.pressure[i]
                      << ',' << relaxation.temperature[i] << ',';
            if (i < c.cells) {
                std::cout << (static_cast<double>(i) + 0.5)*c.dx()
                          << ',' << marching.flow[i]
                          << ',' << newton.flow[i]
                          << ',' << relaxation.flow[i];
            } else {
                std::cout << ",,,";
            }
            std::cout << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "adiabatic_pipe_init_profile_export: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
