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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {
using Clock = std::chrono::steady_clock;
double elapsed(Clock::time_point a, Clock::time_point b)
{ return std::chrono::duration<double>(b-a).count(); }
void report(const char* name, const pipe_init::thermal_wall::Configuration& c,
            const pipe_init::thermal_wall::State& s, double residual, double time)
{
    double flow_error=0.0;
    for(double m:s.flow) flow_error=std::max(flow_error,std::abs(m-c.outlet_mass_flow));
    std::cout<<std::left<<std::setw(38)<<name<<std::right
        <<" p_out="<<std::fixed<<std::setprecision(6)<<s.pressure.back()/1e6<<" MPa"
        <<" Tg_out="<<s.gas_temperature.back()<<" K"
        <<" Tw_out="<<s.wall_temperature.back()<<" K"
        <<" max|m-mout|="<<std::scientific<<flow_error<<" kg/s"
        <<" residual="<<residual<<" time="<<std::fixed<<std::setprecision(6)<<time<<" s\n";
}
void compare(const char* name,const pipe_init::thermal_wall::Configuration& c,
             const pipe_init::thermal_wall::State& s,
             const pipe_init::thermal_wall::State& reference)
{
    double ep=0,eg=0,ew=0,em=0;
    for(std::size_t i=0;i<=c.cells;++i) {
        ep=std::max(ep,std::abs(s.pressure[i]-reference.pressure[i])/c.inlet_pressure);
        eg=std::max(eg,std::abs(s.gas_temperature[i]-reference.gas_temperature[i])/c.inlet_temperature);
        ew=std::max(ew,std::abs(s.wall_temperature[i]-reference.wall_temperature[i])/c.inlet_temperature);
    }
    for(std::size_t i=0;i<c.cells;++i)
        em=std::max(em,std::abs(s.flow[i]-reference.flow[i])/c.outlet_mass_flow);
    std::cout<<"  "<<name<<" vs marching: scaled max errors p="<<std::scientific<<ep
             <<", Tg="<<eg<<", Tw="<<ew<<", m="<<em<<'\n';
}
}

int main()
{
    using namespace pipe_init::thermal_wall;
    try {
        Configuration c;
        c.hill_height = 100.0; // 100 m cosine-shaped hill
        auto t0=Clock::now(); const State analytical=analytical_approximation(c);
        auto t1=Clock::now(); const double ta=elapsed(t0,t1);
        t0=Clock::now(); const State marching=steady_spatial_marching(c);
        t1=Clock::now(); const double tm=elapsed(t0,t1);
        t0=Clock::now();
        const SolveInfo newton=coupled_steady_newton(c,analytical);
        t1=Clock::now(); const double tn=elapsed(t0,t1);
        RelaxationOptions options; options.theta=0.9; options.pseudo_courant=20.0;
        t0=Clock::now(); const RelaxationInfo relaxed=implicit_theta_relaxation(c,options);
        t1=Clock::now(); const double tr=elapsed(t0,t1);

        std::cout<<"Thermal gas / lumped-wall steady initializers\n"
                 <<"cells="<<c.cells<<" length="<<c.length<<" m D="<<c.diameter<<" m"
                 <<" inlet="<<c.inlet_pressure/1e6<<" MPa, Tg_in="<<c.inlet_temperature
                 <<" K soil="<<c.soil_temperature<<" K flow="<<c.outlet_mass_flow<<" kg/s\n"
                 <<"Hi'="<<c.inner_conductance_per_length()<<" W/(m K)"
                 <<" Gsoil'="<<c.soil_conductance_per_length<<" W/(m K)"
                 <<" Cwall'="<<c.wall_heat_capacity_per_length()<<" J/(m K)\n";
        report("Analytical approximation",c,analytical,scaled_steady_residual(c,analytical),ta);
        report("Steady spatial marching",c,marching,scaled_steady_residual(c,marching),tm);
        report("Coupled steady Newton",c,newton.state,newton.scaled_residual,tn);
        compare("Newton",c,newton.state,marching);
        report("Implicit theta pseudo-time relaxation",c,relaxed.state,
               relaxed.scaled_residual,tr);
        compare("Theta relaxation",c,relaxed.state,marching);
        std::cout<<"Newton iterations="<<newton.iterations
                 <<", theta="<<options.theta<<", relaxation steps="<<relaxed.steps
                 <<", pseudo dt="<<relaxed.dt<<" s, pseudo time="<<relaxed.pseudo_time
                 <<" s, inner Newton iterations="<<relaxed.newton_iterations<<'\n';
    } catch(const std::exception& e) {
        std::cerr<<"thermal_wall_init_demo: "<<e.what()<<'\n';
        return 1;
    }
    return 0;
}
