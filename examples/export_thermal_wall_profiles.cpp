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

#include "pipe/thermal_wall_initializers.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main()
{
    using namespace pipe_init::thermal_wall;
    try {
        Configuration c;
        c.hill_height=100.0; // Match the thermal-wall demo
        const State analytical=analytical_approximation(c);
        const State marching=steady_spatial_marching(c);
        const State newton=coupled_steady_newton(c,analytical).state;
        RelaxationOptions options; options.theta=0.9; options.pseudo_courant=20.0;
        const State relaxed=implicit_theta_relaxation(c,options).state;
        const auto z=c.elevation_nodes();
        std::cout<<"x_node_m,elevation_m,p_analytical_Pa,Tg_analytical_K,Tw_analytical_K,"
            "p_marching_Pa,Tg_marching_K,Tw_marching_K,p_newton_Pa,Tg_newton_K,Tw_newton_K,"
            "p_relaxation_Pa,Tg_relaxation_K,Tw_relaxation_K,x_face_m,"
            "m_analytical_kg_s,m_marching_kg_s,m_newton_kg_s,m_relaxation_kg_s\n";
        std::cout<<std::setprecision(15);
        for(std::size_t i=0;i<=c.cells;++i) {
            std::cout<<static_cast<double>(i)*c.dx()<<','<<z[i];
            for(const State* s:{&analytical,&marching,&newton,&relaxed})
                std::cout<<','<<s->pressure[i]<<','<<s->gas_temperature[i]
                         <<','<<s->wall_temperature[i];
            if(i<c.cells) {
                std::cout<<','<<(static_cast<double>(i)+0.5)*c.dx();
                for(const State* s:{&analytical,&marching,&newton,&relaxed})
                    std::cout<<','<<s->flow[i];
            } else std::cout<<",,,,,";
            std::cout<<'\n';
        }
    } catch(const std::exception& e) {
        std::cerr<<"thermal_wall_profile_export: "<<e.what()<<'\n';
        return 1;
    }
    return 0;
}
