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

#include "pipe/soil_shell_initializers.hpp"

#include <iomanip>
#include <iostream>
#include <stdexcept>

int main()
{
    using namespace pipe_init::soil_shell;
    try {
        Configuration c;
        c.hill_height=100.0; // Match the soil-shell demo
        const State analytical=analytical_approximation(c),marching=steady_spatial_marching(c);
        const State newton=coupled_steady_newton(c,analytical).state;
        RelaxationOptions options; const State relaxed=implicit_theta_relaxation(c,options).state;
        const State* states[]={&analytical,&marching,&newton,&relaxed};
        const char* keys[]={"analytical","marching","newton","relaxation"};
        std::cout<<"x_node_m,elevation_m";
        for(std::size_t k=0;k<4;++k) {
            std::cout<<",p_"<<keys[k]<<"_Pa,Tg_"<<keys[k]<<"_K,Tw_"<<keys[k]<<"_K";
            for(std::size_t j=0;j<c.soil_shell_radii.size();++j) std::cout<<",Ts"<<j<<'_'<<keys[k]<<"_K";
        }
        std::cout<<",x_face_m"; for(const char* key:keys) std::cout<<",m_"<<key<<"_kg_s"; std::cout<<'\n'<<std::setprecision(15);
        const auto z=c.elevation_nodes();
        for(std::size_t i=0;i<=c.cells;++i) {
            std::cout<<static_cast<double>(i)*c.dx()<<','<<z[i];
            for(const State* s:states) {
                std::cout<<','<<s->pressure[i]<<','<<s->gas_temperature[i]<<','<<s->wall_temperature[i];
                for(const auto& shell:s->soil_temperature) std::cout<<','<<shell[i];
            }
            if(i<c.cells) { std::cout<<','<<(static_cast<double>(i)+0.5)*c.dx();
                for(const State* s:states) std::cout<<','<<s->flow[i]; }
            else std::cout<<",,,,,";
            std::cout<<'\n';
        }
    } catch(const std::exception& e) { std::cerr<<"soil_shell_profile_export: "<<e.what()<<'\n'; return 1; }
    return 0;
}
