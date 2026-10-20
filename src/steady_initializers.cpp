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
#include <limits>
#include <stdexcept>
#include <utility>

namespace pipe_init {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

// Integrate the constant-Z steady pressure equation exactly over a segment
// whose elevation varies linearly between its endpoint nodes.
double pressure_advance(
    const Configuration& c,
    double friction,
    double mass_flow,
    double pressure,
    double h_left,
    double h_right,
    double distance)
{
    const double gas_product = c.gas_constant_product();
    const double alpha = 2.0*c.gravity/gas_product;
    const double slope = (h_right - h_left)/distance;
    const double exponent = alpha*slope*distance;
    const double beta = friction*gas_product*mass_flow*std::abs(mass_flow)/
                        (c.diameter*c.area()*c.area());
    const double y = pressure*pressure;

    double next_y{};
    if (std::abs(exponent) < 1.0e-10) {
        next_y = y - beta*distance;
    } else {
        const double decay = std::exp(-exponent);
        next_y = decay*y - beta*(-std::expm1(-exponent))/(alpha*slope);
    }
    if (!std::isfinite(next_y) || next_y <= 0.0) {
        throw std::runtime_error("Analytical pressure profile is inadmissible");
    }
    return std::sqrt(next_y);
}

double momentum_residual(
    const Configuration& c,
    double p_left,
    double p_right,
    double h_left,
    double h_right,
    double mass_flow,
    double friction)
{
    const double p_bar = 0.5*(p_left + p_right);
    if (!std::isfinite(p_bar) || p_bar <= 0.0) {
        throw std::runtime_error("Nonpositive pressure in momentum residual");
    }
    const double density = p_bar/c.gas_constant_product();
    const double velocity = mass_flow/(c.area()*density);
    const double dx = c.dx();
    const double pressure_slope = (p_right - p_left)/dx;
    const double elevation_slope = (h_right - h_left)/dx;
    const double friction_slope = friction/(2.0*c.diameter)*density*
                                  velocity*std::abs(velocity);
    return pressure_slope + density*c.gravity*elevation_slope + friction_slope;
}

State decode_unknowns(
    const Configuration& c,
    const std::vector<double>& scaled_unknowns)
{
    const std::size_t n = c.cells;
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();
    State state;
    state.pressure.resize(n + 1);
    state.pressure[0] = c.inlet_pressure;
    state.flow.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        state.pressure[i + 1] = scaled_unknowns[i]*p_scale;
        state.flow[i] = scaled_unknowns[n + i]*m_scale;
    }
    return state;
}

// Unknowns are p_1,...,p_N and m_{1/2},...,m_{N-1/2}. The inlet
// pressure and outlet mass flow are prescribed. The inlet half-cell balance
// determines m_in=m_{1/2} in a steady state.
std::vector<double> steady_residual_vector(
    const Configuration& c,
    const std::vector<double>& scaled_unknowns,
    double friction)
{
    const std::size_t n = c.cells;
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();
    const double momentum_scale = p_scale/c.dx();
    const auto elevations = c.elevation_faces();
    const State state = decode_unknowns(c, scaled_unknowns);
    std::vector<double> residual(2*n, 0.0);

    // Continuity at interior pressure nodes, followed by the outlet half-cell.
    for (std::size_t i = 0; i + 1 < n; ++i) {
        residual[i] = (state.flow[i + 1] - state.flow[i])/m_scale;
    }
    residual[n - 1] = (c.outlet_mass_flow - state.flow[n - 1])/m_scale;

    // Momentum equations at each interior staggered face.
    for (std::size_t i = 0; i < n; ++i) {
        const double raw = momentum_residual(
            c, state.pressure[i], state.pressure[i + 1],
            elevations[i], elevations[i + 1], state.flow[i], friction);
        residual[n + i] = raw/momentum_scale;
    }
    return residual;
}

double infinity_norm(const std::vector<double>& values)
{
    double result = 0.0;
    for (const double value : values) {
        result = std::max(result, std::abs(value));
    }
    return result;
}

std::vector<double> solve_dense(
    std::vector<std::vector<double>> matrix,
    std::vector<double> rhs)
{
    const std::size_t n = rhs.size();
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(matrix[i][k]) > std::abs(matrix[pivot][k])) {
                pivot = i;
            }
        }
        if (std::abs(matrix[pivot][k]) < 1.0e-14) {
            throw std::runtime_error("Singular Newton Jacobian");
        }
        std::swap(matrix[k], matrix[pivot]);
        std::swap(rhs[k], rhs[pivot]);

        const double diagonal = matrix[k][k];
        for (std::size_t j = k; j < n; ++j) {
            matrix[k][j] /= diagonal;
        }
        rhs[k] /= diagonal;

        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) {
                continue;
            }
            const double factor = matrix[i][k];
            if (factor == 0.0) {
                continue;
            }
            for (std::size_t j = k; j < n; ++j) {
                matrix[i][j] -= factor*matrix[k][j];
            }
            rhs[i] -= factor*rhs[k];
        }
    }
    return rhs;
}

struct Derivative {
    std::vector<double> pressure;
    std::vector<double> flow;
};

Derivative rhs(
    const Configuration& c,
    const State& state,
    double friction)
{
    const std::size_t n = c.cells;
    const double dx = c.dx();
    const double area = c.area();
    const double gas_product = c.gas_constant_product();
    const auto elevations = c.elevation_faces();

    Derivative result;
    result.pressure.assign(n + 1, 0.0); // p_0 is the prescribed inlet pressure.
    result.flow.resize(n);

    // Interior pressure control volumes.
    for (std::size_t i = 1; i < n; ++i) {
        result.pressure[i] = -gas_product/(area*dx)*
                             (state.flow[i] - state.flow[i - 1]);
    }
    // Outlet half control volume; m_out is prescribed at x=L.
    result.pressure[n] = -2.0*gas_product/(area*dx)*
                         (c.outlet_mass_flow - state.flow[n - 1]);

    // Momentum equations at the N interior staggered faces.
    for (std::size_t i = 0; i < n; ++i) {
        const double p_bar = 0.5*(state.pressure[i] + state.pressure[i + 1]);
        if (!std::isfinite(p_bar) || p_bar <= 0.0) {
            throw std::runtime_error("Relaxation produced nonpositive pressure");
        }
        const double density = p_bar/gas_product;
        const double velocity = state.flow[i]/(area*density);
        const double pressure_slope =
            (state.pressure[i + 1] - state.pressure[i])/dx;
        const double elevation_slope =
            (elevations[i + 1] - elevations[i])/dx;
        const double drag = friction/(2.0*c.diameter)*density*
                            velocity*std::abs(velocity);
        result.flow[i] = -area*(pressure_slope + density*c.gravity*
                                elevation_slope + drag);
    }
    return result;
}

std::vector<double> theta_residual_vector(
    const Configuration& c,
    const State& old_state,
    const std::vector<double>& scaled_unknowns,
    double friction,
    double theta,
    double dt)
{
    const std::size_t n = c.cells;
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();
    const double momentum_scale = p_scale/c.dx();
    const auto elevations = c.elevation_faces();
    const State state = decode_unknowns(c, scaled_unknowns);
    std::vector<double> residual(2*n, 0.0);
    const double old_weight = 1.0 - theta;
    const double gas_product = c.gas_constant_product();

    // Interior pressure control volumes, then the outlet half volume.
    for (std::size_t i = 0; i + 1 < n; ++i) {
        const double storage = c.area()*c.dx() *
            (state.pressure[i + 1] - old_state.pressure[i + 1]) /
            (gas_product*dt);
        const double right_flow = theta*state.flow[i + 1] +
                                  old_weight*old_state.flow[i + 1];
        const double left_flow = theta*state.flow[i] +
                                 old_weight*old_state.flow[i];
        residual[i] = (storage + right_flow - left_flow)/m_scale;
    }
    const double outlet_old = c.outlet_mass_flow;
    const double outlet_new = c.outlet_mass_flow;
    const double outlet_storage = 0.5*c.area()*c.dx() *
        (state.pressure[n] - old_state.pressure[n])/(gas_product*dt);
    const double outlet_flow_theta = theta*outlet_new + old_weight*outlet_old;
    const double last_flow_theta = theta*state.flow[n - 1] +
                                   old_weight*old_state.flow[n - 1];
    residual[n - 1] = (outlet_storage + outlet_flow_theta -
                       last_flow_theta)/m_scale;

    // Theta-weighted staggered momentum equations.
    for (std::size_t i = 0; i < n; ++i) {
        const double new_momentum = momentum_residual(
            c, state.pressure[i], state.pressure[i + 1],
            elevations[i], elevations[i + 1], state.flow[i], friction);
        const double old_momentum = momentum_residual(
            c, old_state.pressure[i], old_state.pressure[i + 1],
            elevations[i], elevations[i + 1], old_state.flow[i], friction);
        const double temporal = (state.flow[i] - old_state.flow[i]) /
                                (dt*c.area());
        residual[n + i] = (temporal + theta*new_momentum +
                           old_weight*old_momentum)/momentum_scale;
    }
    return residual;
}

State theta_step(
    const Configuration& c,
    const State& old_state,
    double friction,
    double theta,
    double dt,
    std::size_t& newton_iterations)
{
    const std::size_t n = c.cells;
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();
    std::vector<double> unknowns(2*n);
    for (std::size_t i = 0; i < n; ++i) {
        unknowns[i] = old_state.pressure[i + 1]/p_scale;
        unknowns[n + i] = old_state.flow[i]/m_scale;
    }

    constexpr std::size_t max_newton_iterations = 20;
    constexpr double newton_tolerance = 2.0e-11;
    for (std::size_t iteration = 0;
         iteration <= max_newton_iterations; ++iteration) {
        const auto residual = theta_residual_vector(
            c, old_state, unknowns, friction, theta, dt);
        const double norm = infinity_norm(residual);
        if (norm <= newton_tolerance) {
            return decode_unknowns(c, unknowns);
        }
        if (iteration == max_newton_iterations) {
            break;
        }

        std::vector<std::vector<double>> jacobian(
            2*n, std::vector<double>(2*n, 0.0));
        for (std::size_t column = 0; column < 2*n; ++column) {
            auto perturbed = unknowns;
            const double increment =
                1.0e-6*std::max(1.0, std::abs(unknowns[column]));
            perturbed[column] += increment;
            const auto shifted = theta_residual_vector(
                c, old_state, perturbed, friction, theta, dt);
            for (std::size_t row = 0; row < 2*n; ++row) {
                jacobian[row][column] =
                    (shifted[row] - residual[row])/increment;
            }
        }
        std::vector<double> correction(2*n);
        std::transform(residual.begin(), residual.end(), correction.begin(),
                       [](double value) { return -value; });
        correction = solve_dense(std::move(jacobian), std::move(correction));
        ++newton_iterations;

        bool accepted = false;
        double step_length = 1.0;
        for (int line_search = 0; line_search < 24; ++line_search) {
            auto candidate = unknowns;
            for (std::size_t i = 0; i < 2*n; ++i) {
                candidate[i] += step_length*correction[i];
            }
            bool positive_pressure = true;
            for (std::size_t i = 0; i < n; ++i) {
                positive_pressure = positive_pressure && candidate[i] > 0.0;
            }
            if (positive_pressure) {
                try {
                    const auto candidate_residual = theta_residual_vector(
                        c, old_state, candidate, friction, theta, dt);
                    if (infinity_norm(candidate_residual) < norm) {
                        unknowns = std::move(candidate);
                        accepted = true;
                        break;
                    }
                } catch (const std::runtime_error&) {
                    // Backtrack through inadmissible trial states.
                }
            }
            step_length *= 0.5;
        }
        if (!accepted) {
            throw std::runtime_error("Theta-step Newton line search failed");
        }
    }
    throw std::runtime_error("Theta-step Newton solve did not converge");
}

double scaled_rhs_norm(
    const Configuration& c,
    const State& state,
    double friction)
{
    const Derivative value = rhs(c, state, friction);
    double result = 0.0;
    for (std::size_t i = 1; i < value.pressure.size(); ++i) {
        result = std::max(result,
                          std::abs(value.pressure[i])/c.pressure_scale());
    }
    for (const double dm : value.flow) {
        result = std::max(result,
                          std::abs(dm)/c.mass_flow_scale());
    }
    return result;
}

} // namespace

double Configuration::area() const
{
    return kPi*diameter*diameter/4.0;
}

double Configuration::dx() const
{
    return length/static_cast<double>(cells);
}

double Configuration::pressure_scale() const
{
    return inlet_pressure;
}

double Configuration::mass_flow_scale() const
{
    return std::max(1.0, std::abs(outlet_mass_flow));
}

double Configuration::gas_constant_product() const
{
    return compressibility*gas_constant*temperature;
}

std::vector<double> Configuration::elevation_faces() const
{
    std::vector<double> result(cells + 1, 0.0);
    for (std::size_t i = 0; i <= cells; ++i) {
        const double z = length*static_cast<double>(i)/
                         static_cast<double>(cells);
        result[i] = 0.5*hill_height*(1.0 - std::cos(2.0*kPi*z/length));
    }
    return result;
}

void Configuration::validate() const
{
    const auto finite_positive = [](double value) {
        return std::isfinite(value) && value > 0.0;
    };
    if (cells < 2 || !finite_positive(length) || !finite_positive(diameter) ||
        !finite_positive(viscosity) || !finite_positive(compressibility) ||
        !finite_positive(gas_constant) || !finite_positive(temperature) ||
        !finite_positive(inlet_pressure) || !std::isfinite(outlet_mass_flow) ||
        !std::isfinite(roughness) || roughness < 0.0 ||
        !std::isfinite(gravity) || gravity < 0.0 ||
        !std::isfinite(hill_height)) {
        throw std::invalid_argument("Invalid pipe configuration");
    }
}

double swamee_jain_friction_factor(
    const Configuration& c,
    double mass_flow)
{
    c.validate();
    if (!std::isfinite(mass_flow)) {
        throw std::invalid_argument("Mass flow must be finite");
    }
    if (mass_flow == 0.0) {
        return 0.0;
    }
    const double reynolds = std::abs(mass_flow)*c.diameter/
                            (c.area()*c.viscosity);
    const double argument = c.roughness/(3.7*c.diameter) +
                            5.74/std::pow(reynolds, 0.9);
    if (!std::isfinite(argument) || argument <= 0.0 ||
        std::abs(argument - 1.0) < 1.0e-14) {
        throw std::runtime_error("Invalid Swamee-Jain logarithm argument");
    }
    const double logarithm = std::log10(argument);
    const double factor = 0.25/(logarithm*logarithm);
    if (!std::isfinite(factor) || factor <= 0.0) {
        throw std::runtime_error("Invalid Darcy friction factor");
    }
    return factor;
}

State analytical_pressure_squared(const Configuration& c)
{
    c.validate();
    const double friction =
        swamee_jain_friction_factor(c, c.outlet_mass_flow);
    const auto elevations = c.elevation_faces();
    State result;
    result.pressure.resize(c.cells + 1);
    result.pressure[0] = c.inlet_pressure;
    result.flow.assign(c.cells, c.outlet_mass_flow);

    for (std::size_t i = 0; i < c.cells; ++i) {
        result.pressure[i + 1] = pressure_advance(
            c, friction, c.outlet_mass_flow, result.pressure[i],
            elevations[i], elevations[i + 1], c.dx());
    }
    return result;
}

State discrete_steady_marching(const Configuration& c)
{
    c.validate();
    const double mass_flow = c.outlet_mass_flow;
    const double friction = swamee_jain_friction_factor(c, mass_flow);
    const auto elevations = c.elevation_faces();
    const double gas_product = c.gas_constant_product();
    const double area = c.area();
    const double dx = c.dx();

    State result;
    result.pressure.resize(c.cells + 1);
    result.pressure[0] = c.inlet_pressure;
    result.flow.assign(c.cells, mass_flow);

    // For each cell the discrete momentum residual depends on the pressure
    // sum S = p_left + p_right. Multiplying the residual by dx*S gives
    // (1+k)S^2 - 2*p_left*S + q = 0, with
    // k = g*dh/(2*ZRT) and q = dx*f*ZRT*m*|m|/(D*A^2).
    // Select the positive-pressure root continuous with p_right=p_left as
    // flow and elevation change tend to zero.
    for (std::size_t i = 0; i < c.cells; ++i) {
        const double p_left = result.pressure[i];
        const double k = c.gravity*(elevations[i + 1] - elevations[i])/
                         (2.0*gas_product);
        const double q = dx*friction*gas_product*mass_flow*std::abs(mass_flow)/
                         (c.diameter*area*area);
        const double a = 1.0 + k;
        const double discriminant =
            4.0*p_left*p_left - 4.0*a*q;
        if (!std::isfinite(discriminant) || discriminant < 0.0 ||
            std::abs(a) < 1.0e-14) {
            throw std::runtime_error(
                "Discrete steady march has no admissible pressure root");
        }
        const double root = std::sqrt(discriminant);
        const double sum_candidates[2] = {
            (2.0*p_left + root)/(2.0*a),
            (2.0*p_left - root)/(2.0*a),
        };
        double p_right = std::numeric_limits<double>::quiet_NaN();
        double smallest_change = std::numeric_limits<double>::infinity();
        for (const double pressure_sum : sum_candidates) {
            const double candidate = pressure_sum - p_left;
            if (pressure_sum > 0.0 && candidate > 0.0 &&
                std::abs(candidate - p_left) < smallest_change) {
                p_right = candidate;
                smallest_change = std::abs(candidate - p_left);
            }
        }
        if (!std::isfinite(p_right)) {
            throw std::runtime_error(
                "Discrete steady march produced nonpositive pressure");
        }
        result.pressure[i + 1] = p_right;
    }
    return result;
}

double scaled_steady_residual(
    const Configuration& c,
    const State& state)
{
    c.validate();
    if (state.pressure.size() != c.cells + 1 ||
        state.flow.size() != c.cells) {
        throw std::invalid_argument("State dimensions do not match the grid");
    }
    const double friction =
        swamee_jain_friction_factor(c, c.outlet_mass_flow);
    const auto elevations = c.elevation_faces();
    const double m_scale = c.mass_flow_scale();
    double norm = 0.0;

    for (std::size_t i = 0; i + 1 < c.cells; ++i) {
        norm = std::max(norm,
            std::abs(state.flow[i + 1] - state.flow[i])/m_scale);
    }
    norm = std::max(norm,
        std::abs(c.outlet_mass_flow - state.flow.back())/m_scale);

    for (std::size_t i = 0; i < c.cells; ++i) {
        const double momentum = momentum_residual(
            c, state.pressure[i], state.pressure[i + 1],
            elevations[i], elevations[i + 1], state.flow[i], friction);
        norm = std::max(norm,
            std::abs(momentum)*c.dx()/c.pressure_scale());
    }
    return norm;
}

SolveInfo newton_steady_state(
    const Configuration& c,
    const State& initial_guess,
    double tolerance,
    std::size_t max_iterations)
{
    c.validate();
    if (initial_guess.pressure.size() != c.cells + 1 ||
        initial_guess.flow.size() != c.cells ||
        !std::isfinite(tolerance) || tolerance <= 0.0 || max_iterations == 0) {
        throw std::invalid_argument("Invalid Newton initial guess or options");
    }
    const std::size_t n = c.cells;
    const double friction =
        swamee_jain_friction_factor(c, c.outlet_mass_flow);
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();

    std::vector<double> y(2*n);
    for (std::size_t i = 0; i < n; ++i) {
        y[i] = initial_guess.pressure[i + 1]/p_scale;
        y[n + i] = initial_guess.flow[i]/m_scale;
    }

    for (std::size_t iteration = 0; iteration <= max_iterations; ++iteration) {
        const auto residual = steady_residual_vector(c, y, friction);
        const double norm = infinity_norm(residual);
        if (norm <= tolerance) {
            return {decode_unknowns(c, y), iteration, norm};
        }
        if (iteration == max_iterations) {
            break;
        }

        std::vector<std::vector<double>> jacobian(
            2*n, std::vector<double>(2*n, 0.0));
        for (std::size_t column = 0; column < 2*n; ++column) {
            std::vector<double> perturbed = y;
            const double step = 1.0e-6*std::max(1.0, std::abs(y[column]));
            perturbed[column] += step;
            const auto shifted = steady_residual_vector(c, perturbed, friction);
            for (std::size_t row = 0; row < 2*n; ++row) {
                jacobian[row][column] = (shifted[row] - residual[row])/step;
            }
        }
        std::vector<double> rhs_value(2*n);
        std::transform(residual.begin(), residual.end(), rhs_value.begin(),
                       [](double value) { return -value; });
        const auto correction =
            solve_dense(std::move(jacobian), std::move(rhs_value));

        bool accepted = false;
        double step_length = 1.0;
        for (int line_search = 0; line_search < 24; ++line_search) {
            std::vector<double> candidate = y;
            for (std::size_t i = 0; i < 2*n; ++i) {
                candidate[i] += step_length*correction[i];
            }
            bool positive_pressure = true;
            for (std::size_t i = 0; i < n; ++i) {
                positive_pressure = positive_pressure && candidate[i] > 0.0;
            }
            if (positive_pressure) {
                try {
                    const auto candidate_residual =
                        steady_residual_vector(c, candidate, friction);
                    if (infinity_norm(candidate_residual) < norm) {
                        y = std::move(candidate);
                        accepted = true;
                        break;
                    }
                } catch (const std::runtime_error&) {
                    // Backtrack through inadmissible trial states.
                }
            }
            step_length *= 0.5;
        }
        if (!accepted) {
            throw std::runtime_error("Damped Newton line search failed");
        }
    }
    throw std::runtime_error("Newton steady solve did not converge");
}

RelaxationInfo relax_to_steady_state(
    const Configuration& c,
    const RelaxationOptions& options)
{
    c.validate();
    if (!std::isfinite(options.tolerance) || options.tolerance <= 0.0 ||
        !std::isfinite(options.max_pseudo_time) || options.max_pseudo_time <= 0.0 ||
        options.max_steps == 0 || options.consecutive_converged_steps == 0 ||
        !std::isfinite(options.theta) || options.theta < 0.5 ||
        options.theta > 1.0 || !std::isfinite(options.pseudo_courant) ||
        options.pseudo_courant <= 0.0) {
        throw std::invalid_argument("Invalid relaxation options");
    }

    const double friction =
        swamee_jain_friction_factor(c, c.outlet_mass_flow);
    State state;
    state.pressure.assign(c.cells + 1, c.inlet_pressure);
    state.flow.assign(c.cells, c.outlet_mass_flow);
    const double acoustic_speed = std::sqrt(c.gas_constant_product());
    const double dt = options.pseudo_courant*c.dx()/acoustic_speed;
    std::size_t consecutive = 0;
    double pseudo_time = 0.0;
    double scaled_rhs = scaled_rhs_norm(c, state, friction);
    std::size_t total_newton_iterations = 0;

    for (std::size_t step = 1; step <= options.max_steps; ++step) {
        if (pseudo_time + dt > options.max_pseudo_time) {
            break;
        }
        std::size_t step_newton_iterations = 0;
        state = theta_step(c, state, friction, options.theta, dt,
                           step_newton_iterations);
        total_newton_iterations += step_newton_iterations;
        pseudo_time += dt;
        scaled_rhs = scaled_rhs_norm(c, state, friction);
        if (scaled_rhs <= options.tolerance) {
            ++consecutive;
            if (consecutive >= options.consecutive_converged_steps) {
                return {std::move(state), step, pseudo_time, scaled_rhs, dt,
                        total_newton_iterations};
            }
        } else {
            consecutive = 0;
        }
    }
    throw std::runtime_error(
        "Pseudo-time relaxation did not reach steady tolerance; increase limits or inspect the case");
}

} // namespace pipe_init
