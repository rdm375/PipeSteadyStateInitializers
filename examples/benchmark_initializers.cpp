#include "pipe/steady_initializers.hpp"

#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {

using Clock = std::chrono::steady_clock;

double elapsed_seconds(Clock::time_point start, Clock::time_point finish)
{
    return std::chrono::duration<double>(finish - start).count();
}

pipe_init::State flat_initial_guess(const pipe_init::Configuration& c)
{
    pipe_init::State state;
    state.pressure.assign(c.cells + 1, c.inlet_pressure);
    state.flow.assign(c.cells, c.outlet_mass_flow);
    return state;
}

std::size_t parse_repetitions(const char* value, const char* label)
{
    const auto parsed = std::stoull(value);
    if (parsed == 0) {
        throw std::invalid_argument(std::string(label) + " must be positive");
    }
    return static_cast<std::size_t>(parsed);
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const std::size_t fast_repetitions = argc > 1
            ? parse_repetitions(argv[1], "fast repetition count") : 30;
        const std::size_t relaxation_repetitions = argc > 2
            ? parse_repetitions(argv[2], "relaxation repetition count") : 3;
        if (argc > 3) {
            throw std::invalid_argument(
                "usage: pipe_init_benchmark [fast_repetitions] [relaxation_repetitions]");
        }

        pipe_init::Configuration c;
        pipe_init::RelaxationOptions relaxation_options;

        // Warm each code path before collecting measurements.
        const auto analytical_warmup =
            pipe_init::analytical_pressure_squared(c);
        const auto marching_warmup = pipe_init::discrete_steady_marching(c);
        const auto newton_warmup = pipe_init::newton_steady_state(
            c, flat_initial_guess(c));
        const auto relaxation_warmup =
            pipe_init::relax_to_steady_state(c, relaxation_options);
        (void)analytical_warmup;
        (void)marching_warmup;
        (void)newton_warmup;
        (void)relaxation_warmup;

        std::cout << "method,repeat,seconds,newton_iterations,steps,scaled_steady_residual\n"
                  << std::setprecision(12);

        for (std::size_t repeat = 1; repeat <= fast_repetitions; ++repeat) {
            const auto start = Clock::now();
            const auto state = pipe_init::discrete_steady_marching(c);
            const auto finish = Clock::now();
            const double seconds = elapsed_seconds(start, finish);
            const double residual = pipe_init::scaled_steady_residual(c, state);
            std::cout << "discrete_marching," << repeat << ',' << seconds
                      << ",0,0," << residual << '\n';
        }

        for (std::size_t repeat = 1; repeat <= fast_repetitions; ++repeat) {
            auto start = Clock::now();
            const auto state = pipe_init::analytical_pressure_squared(c);
            const auto finish = Clock::now();
            const double seconds = elapsed_seconds(start, finish);
            const double residual = pipe_init::scaled_steady_residual(c, state);
            std::cout << "analytical," << repeat << ',' << seconds
                      << ",0,0," << residual << '\n';
        }

        for (std::size_t repeat = 1; repeat <= fast_repetitions; ++repeat) {
            const auto start = Clock::now();
            const auto solved = pipe_init::newton_steady_state(
                c, flat_initial_guess(c));
            const auto finish = Clock::now();
            const double seconds = elapsed_seconds(start, finish);
            const double residual =
                pipe_init::scaled_steady_residual(c, solved.state);
            std::cout << "newton_from_flat," << repeat << ',' << seconds
                      << ',' << solved.iterations << ",0," << residual << '\n';
        }

        for (std::size_t repeat = 1; repeat <= relaxation_repetitions; ++repeat) {
            const auto start = Clock::now();
            const auto solved =
                pipe_init::relax_to_steady_state(c, relaxation_options);
            const auto finish = Clock::now();
            const double seconds = elapsed_seconds(start, finish);
            const double residual =
                pipe_init::scaled_steady_residual(c, solved.state);
            std::cout << "theta_relaxation," << repeat << ',' << seconds
                      << ',' << solved.newton_iterations << ',' << solved.steps
                      << ',' << residual << '\n';
        }
    } catch (const std::exception& error) {
        std::cerr << "pipe_init_benchmark: " << error.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
