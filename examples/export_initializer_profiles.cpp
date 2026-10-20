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

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main()
{
    try {
        pipe_init::Configuration configuration;

        const auto analytical =
            pipe_init::analytical_pressure_squared(configuration);
        const auto marching =
            pipe_init::discrete_steady_marching(configuration);
        const auto newton = pipe_init::newton_steady_state(
            configuration, marching).state;
        const auto relaxation =
            pipe_init::relax_to_steady_state(configuration).state;

        // The simulator grid has N+1 pressure/elevation nodes and N internal
        // flows at the staggered faces. Boundary flows are included for plots.
        const auto elevations = configuration.elevation_faces();
        std::cout << "x_pressure_m,p_analytical_Pa,p_marching_Pa,p_newton_Pa,p_relaxation_Pa,"
                     "x_flow_m,m_analytical_kg_s,m_marching_kg_s,m_newton_kg_s,m_relaxation_kg_s,"
                     "x_elevation_m,elevation_m\n";
        std::cout << std::setprecision(15);
        for (std::size_t row = 0; row <= configuration.cells + 1; ++row) {
            if (row <= configuration.cells) {
                const double x = static_cast<double>(row)*configuration.dx();
                std::cout << x << ',' << analytical.pressure[row] << ','
                          << marching.pressure[row] << ','
                          << newton.pressure[row] << ','
                          << relaxation.pressure[row] << ',';
            } else {
                std::cout << ",,,,,";
            }

            if (row == 0) {
                std::cout << 0.0 << ',' << analytical.flow.front() << ','
                          << marching.flow.front() << ','
                          << newton.flow.front() << ','
                          << relaxation.flow.front() << ',';
            } else if (row <= configuration.cells) {
                const std::size_t face = row - 1;
                const double x = (static_cast<double>(face) + 0.5)
                               * configuration.dx();
                std::cout << x << ',' << analytical.flow[face] << ','
                          << marching.flow[face] << ','
                          << newton.flow[face] << ','
                          << relaxation.flow[face] << ',';
            } else {
                std::cout << configuration.length << ','
                          << configuration.outlet_mass_flow << ','
                          << configuration.outlet_mass_flow << ','
                          << configuration.outlet_mass_flow << ','
                          << configuration.outlet_mass_flow << ',';
            }

            if (row <= configuration.cells) {
                std::cout << static_cast<double>(row)*configuration.dx()
                          << ',' << elevations[row];
            } else {
                std::cout << ",";
            }
            std::cout << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "pipe_init_profile_export: " << error.what() << '\n';
        return 1;
    }
    return 0;
}
