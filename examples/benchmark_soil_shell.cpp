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

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
using Clock=std::chrono::steady_clock;
std::size_t positive_count(const char* text,const char* name)
{
    const std::string value(text);
    if(value.empty()||value.front()=='-') throw std::invalid_argument(std::string(name)+" must be a positive integer");
    std::size_t used=0;
    const auto count=std::stoull(value,&used);
    if(used!=value.size()||count==0) throw std::invalid_argument(std::string(name)+" must be a positive integer");
    return static_cast<std::size_t>(count);
}
template<class Solver,class Summarize>
void run(const char* method,std::size_t repetitions,std::size_t warmups,
         Solver solve,Summarize summarize)
{
    for(std::size_t i=0;i<warmups;++i) (void)solve();
    for(std::size_t repeat=1;repeat<=repetitions;++repeat) {
        const auto start=Clock::now();
        const auto result=solve();
        const auto stop=Clock::now();
        const double elapsed=std::chrono::duration<double>(stop-start).count();
        const auto metrics=summarize(result);
        std::cout<<method<<','<<repeat<<','<<std::setprecision(12)<<elapsed<<','
            <<metrics.steps<<','<<metrics.newton_iterations<<','<<metrics.residual<<','
            <<metrics.pressure<<','<<metrics.temperature<<'\n';
    }
}
struct TimedResult {
    std::size_t steps{};
    std::size_t newton_iterations{};
    double residual{};
    double pressure{};
    double temperature{};
};
}

int main(int argc,char** argv)
{
    try {
        const std::size_t repetitions=argc>1?positive_count(argv[1],"repetitions"):30;
        const std::size_t relaxation_repetitions=argc>2?positive_count(argv[2],"relaxation repetitions"):3;
        const std::size_t warmups=argc>3?positive_count(argv[3],"warmups"):1;
        if(argc>4) throw std::invalid_argument("usage: soil_shell_benchmark [repetitions=30] [relaxation_repetitions=3] [warmups=1]");
        using namespace pipe_init::soil_shell;
        Configuration c;
        const State analytical=analytical_approximation(c);
        RelaxationOptions options; options.theta=1.0; options.pseudo_courant=0.25;
        std::cout<<"method,repeat,seconds,steps,newton_iterations,scaled_steady_residual,p_out_Pa,Tg_out_K\n";
        run("analytical",repetitions,warmups,
            [&] { return analytical_approximation(c); },
            [&](const State& s) { return TimedResult{0,0,scaled_steady_residual(c,s),s.pressure.back(),s.gas_temperature.back()}; });
        run("marching",repetitions,warmups,
            [&] { return steady_spatial_marching(c); },
            [&](const State& s) { return TimedResult{0,0,scaled_steady_residual(c,s),s.pressure.back(),s.gas_temperature.back()}; });
        run("newton",repetitions,warmups,
            [&] { return coupled_steady_newton(c,analytical); },
            [](const SolveInfo& r) { return TimedResult{0,r.iterations,r.scaled_residual,r.state.pressure.back(),r.state.gas_temperature.back()}; });
        run("theta_relaxation",relaxation_repetitions,warmups,
            [&] { return implicit_theta_relaxation(c,options); },
            [](const RelaxationInfo& r) { return TimedResult{r.steps,r.newton_iterations,r.scaled_residual,r.state.pressure.back(),r.state.gas_temperature.back()}; });
    } catch(const std::exception& error) {
        std::cerr<<"soil_shell_benchmark: "<<error.what()<<'\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
