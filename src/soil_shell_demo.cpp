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

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>

namespace {
using Clock=std::chrono::steady_clock;
double seconds(Clock::time_point a,Clock::time_point b)
{ return std::chrono::duration<double>(b-a).count(); }
void report(const char* name,const pipe_init::soil_shell::Configuration& c,
            const pipe_init::soil_shell::State& s,double residual,double runtime)
{
    double flow_error=0.0; for(double m:s.flow) flow_error=std::max(flow_error,std::abs(m-c.outlet_mass_flow));
    std::cout<<std::left<<std::setw(38)<<name<<std::right
        <<" p_out="<<std::fixed<<std::setprecision(6)<<s.pressure.back()/1e6<<" MPa"
        <<" Tg_out="<<s.gas_temperature.back()<<" K Tw_out="<<s.wall_temperature.back()<<" K"
        <<" Ts_outer_out="<<s.soil_temperature.back().back()<<" K"
        <<" max|m-mout|="<<std::scientific<<flow_error<<" kg/s residual="<<residual
        <<" time="<<std::fixed<<std::setprecision(6)<<runtime<<" s\n";
}
void compare(const char* name,const pipe_init::soil_shell::Configuration& c,
             const pipe_init::soil_shell::State& a,const pipe_init::soil_shell::State& b)
{
    double ep=0,eg=0,ew=0,em=0,es=0;
    for(std::size_t i=0;i<=c.cells;++i) {
        ep=std::max(ep,std::abs(a.pressure[i]-b.pressure[i])/c.inlet_pressure);
        eg=std::max(eg,std::abs(a.gas_temperature[i]-b.gas_temperature[i])/c.inlet_temperature);
        ew=std::max(ew,std::abs(a.wall_temperature[i]-b.wall_temperature[i])/c.inlet_temperature);
        for(std::size_t j=0;j<c.soil_shell_radii.size();++j)
            es=std::max(es,std::abs(a.soil_temperature[j][i]-b.soil_temperature[j][i])/c.inlet_temperature);
    }
    for(std::size_t i=0;i<c.cells;++i) em=std::max(em,std::abs(a.flow[i]-b.flow[i])/c.outlet_mass_flow);
    std::cout<<"  "<<name<<" vs marching: max scaled errors p="<<std::scientific<<ep
        <<", Tg="<<eg<<", Tw="<<ew<<", soil="<<es<<", m="<<em<<'\n';
}
}
int main()
{
    using namespace pipe_init::soil_shell;
    try {
        Configuration c;
        c.hill_height=100.0; // 100 m cosine-shaped hill
        auto a0=Clock::now(); const State analytical=analytical_approximation(c); auto a1=Clock::now();
        auto m0=Clock::now(); const State marching=steady_spatial_marching(c); auto m1=Clock::now();
        auto n0=Clock::now(); const SolveInfo newton=coupled_steady_newton(c,analytical); auto n1=Clock::now();
        RelaxationOptions options; options.theta=1.0; options.pseudo_courant=0.25;
        auto r0=Clock::now(); const RelaxationInfo relaxed=implicit_theta_relaxation(c,options); auto r1=Clock::now();
        std::cout<<"Fully thermal gas / pipe wall / soil-shell steady initializers\n"
            <<"cells="<<c.cells<<" length="<<c.length<<" m D="<<c.diameter<<" m inlet="<<c.inlet_pressure/1e6
            <<" MPa, Tg_in="<<c.inlet_temperature<<" K Tfar="<<c.far_field_temperature<<" K flow="<<c.outlet_mass_flow<<" kg/s\n"
            <<"radial soil shells="<<c.soil_shell_radii.size()<<", outer radius="<<c.far_field_radius
            <<" m, hill="<<c.hill_height<<" m, unknown count="<<4*c.cells+1+c.soil_shell_radii.size()*(c.cells+1)<<'\n';
        report("Analytical / semi-analytical approximation",c,analytical,scaled_steady_residual(c,analytical),seconds(a0,a1));
        report("Steady spatial marching",c,marching,scaled_steady_residual(c,marching),seconds(m0,m1));
        report("Coupled steady Newton",c,newton.state,newton.scaled_residual,seconds(n0,n1)); compare("Newton",c,newton.state,marching);
        report("Implicit theta pseudo-time relaxation",c,relaxed.state,relaxed.scaled_residual,seconds(r0,r1)); compare("Theta relaxation",c,relaxed.state,marching);
        std::cout<<"Newton iterations="<<newton.iterations<<", theta="<<options.theta
            <<", relaxation steps="<<relaxed.steps<<", pseudo dt="<<relaxed.dt<<" s, pseudo time="<<relaxed.pseudo_time
            <<" s, inner Newton iterations="<<relaxed.newton_iterations<<'\n';
    } catch(const std::exception& e) { std::cerr<<"soil_shell_init_demo: "<<e.what()<<'\n'; return 1; }
    return 0;
}
