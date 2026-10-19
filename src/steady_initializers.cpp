#include "pipe/steady_initializers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>

namespace pipe_init {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

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

double face_residual(
    const Configuration& c,
    const State& state,
    std::size_t face,
    double p_left,
    double p_right,
    double h_left,
    double h_right,
    double distance,
    double friction)
{
    const double p_bar = 0.5*(p_left + p_right);
    if (!std::isfinite(p_bar) || p_bar <= 0.0) {
        throw std::runtime_error("Nonpositive face pressure in steady residual");
    }
    const double density = p_bar/c.gas_constant_product();
    const double velocity = state.flow.at(face)/(c.area()*density);
    const double elevation_slope = (h_right - h_left)/distance;
    const double pressure_slope = (p_right - p_left)/distance;
    const double friction_slope = friction/(2.0*c.diameter)*density*
                                  velocity*std::abs(velocity);
    return pressure_slope + density*c.gravity*elevation_slope + friction_slope;
}

std::vector<double> residual_vector(
    const Configuration& c,
    const std::vector<double>& scaled_unknowns,
    double friction)
{
    const std::size_t n = c.cells;
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();
    const double dx = c.dx();
    const double momentum_scale = p_scale/dx;
    const auto h_faces = c.elevation_faces();

    State state;
    state.pressure.resize(n);
    state.flow.resize(n);
    for (std::size_t i = 0; i < n; ++i) {
        state.pressure[i] = scaled_unknowns[i]*p_scale;
        state.flow[i] = scaled_unknowns[n + i]*m_scale;
    }

    std::vector<double> residual(2*n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        const double flow_right =
            (i + 1 < n) ? state.flow[i + 1] : c.outlet_mass_flow;
        residual[i] = (flow_right - state.flow[i])/m_scale;

        double p_left{};
        double p_right{};
        double h_left{};
        double h_right{};
        double distance{};
        if (i == 0) {
            p_left = c.inlet_pressure;
            p_right = state.pressure[0];
            h_left = h_faces[0];
            h_right = 0.5*(h_faces[0] + h_faces[1]);
            distance = 0.5*dx;
        } else {
            p_left = state.pressure[i - 1];
            p_right = state.pressure[i];
            h_left = 0.5*(h_faces[i - 1] + h_faces[i]);
            h_right = 0.5*(h_faces[i] + h_faces[i + 1]);
            distance = dx;
        }
        const double raw_momentum = face_residual(
            c, state, i, p_left, p_right, h_left, h_right,
            distance, friction);
        residual[n + i] = raw_momentum/momentum_scale;
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
    const auto h_faces = c.elevation_faces();

    Derivative result;
    result.pressure.resize(n);
    result.flow.resize(n);

    for (std::size_t i = 0; i < n; ++i) {
        const double flow_right =
            (i + 1 < n) ? state.flow[i + 1] : c.outlet_mass_flow;
        result.pressure[i] = -gas_product/(area*dx)*
                             (flow_right - state.flow[i]);

        double p_left{};
        double p_right{};
        double h_left{};
        double h_right{};
        double distance{};
        if (i == 0) {
            p_left = c.inlet_pressure;
            p_right = state.pressure[0];
            h_left = h_faces[0];
            h_right = 0.5*(h_faces[0] + h_faces[1]);
            distance = 0.5*dx;
        } else {
            p_left = state.pressure[i - 1];
            p_right = state.pressure[i];
            h_left = 0.5*(h_faces[i - 1] + h_faces[i]);
            h_right = 0.5*(h_faces[i] + h_faces[i + 1]);
            distance = dx;
        }
        const double p_bar = 0.5*(p_left + p_right);
        if (!std::isfinite(p_bar) || p_bar <= 0.0) {
            throw std::runtime_error("Relaxation produced nonpositive pressure");
        }
        const double density = p_bar/gas_product;
        const double velocity = state.flow[i]/(area*density);
        const double slope = (h_right - h_left)/distance;
        const double drag = friction/(2.0*c.diameter)*density*
                            velocity*std::abs(velocity);
        result.flow[i] = -area*(p_right - p_left)/distance
                         -area*density*c.gravity*slope
                         -area*drag;
    }
    return result;
}

State add_scaled(
    const State& base,
    const Derivative& k,
    double scale)
{
    State result = base;
    for (std::size_t i = 0; i < result.pressure.size(); ++i) {
        result.pressure[i] += scale*k.pressure[i];
        result.flow[i] += scale*k.flow[i];
    }
    return result;
}

State rk4_step(
    const Configuration& c,
    const State& state,
    double friction,
    double dt)
{
    const Derivative k1 = rhs(c, state, friction);
    const Derivative k2 = rhs(c, add_scaled(state, k1, 0.5*dt), friction);
    const Derivative k3 = rhs(c, add_scaled(state, k2, 0.5*dt), friction);
    const Derivative k4 = rhs(c, add_scaled(state, k3, dt), friction);

    State result = state;
    for (std::size_t i = 0; i < result.pressure.size(); ++i) {
        result.pressure[i] += dt/6.0*(k1.pressure[i] + 2.0*k2.pressure[i]
                                      + 2.0*k3.pressure[i] + k4.pressure[i]);
        result.flow[i] += dt/6.0*(k1.flow[i] + 2.0*k2.flow[i]
                                  + 2.0*k3.flow[i] + k4.flow[i]);
        if (!std::isfinite(result.pressure[i]) || result.pressure[i] <= 0.0 ||
            !std::isfinite(result.flow[i])) {
            throw std::runtime_error("Relaxation step left the admissible state set");
        }
    }
    return result;
}

double scaled_rhs_norm(
    const Configuration& c,
    const State& state,
    double friction)
{
    const Derivative value = rhs(c, state, friction);
    double result = 0.0;
    for (double dp : value.pressure) {
        result = std::max(result, std::abs(dp)/c.pressure_scale());
    }
    for (double dm : value.flow) {
        result = std::max(result, std::abs(dm)/c.mass_flow_scale());
    }
    return result;
}

double outlet_pressure_from_last_cell(
    const Configuration& c,
    double friction,
    double last_pressure)
{
    const auto h = c.elevation_faces();
    const double h_center = 0.5*(h[c.cells - 1] + h[c.cells]);
    return pressure_advance(c, friction, c.outlet_mass_flow, last_pressure,
                            h_center, h[c.cells], 0.5*c.dx());
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
    const auto h = c.elevation_faces();
    State result;
    result.pressure.resize(c.cells);
    result.flow.assign(c.cells, c.outlet_mass_flow);

    double p_left = c.inlet_pressure;
    for (std::size_t i = 0; i < c.cells; ++i) {
        const double h_mid = 0.5*(h[i] + h[i + 1]);
        const double p_mid = pressure_advance(
            c, friction, c.outlet_mass_flow, p_left,
            h[i], h_mid, 0.5*c.dx());
        result.pressure[i] = p_mid;
        p_left = pressure_advance(
            c, friction, c.outlet_mass_flow, p_mid,
            h_mid, h[i + 1], 0.5*c.dx());
    }
    result.outlet_pressure = p_left;
    return result;
}

double scaled_steady_residual(
    const Configuration& c,
    const State& state)
{
    c.validate();
    if (state.pressure.size() != c.cells || state.flow.size() != c.cells) {
        throw std::invalid_argument("State dimensions do not match the grid");
    }
    const double friction =
        swamee_jain_friction_factor(c, c.outlet_mass_flow);
    const auto h = c.elevation_faces();
    const double dx = c.dx();
    const double p_scale = c.pressure_scale();
    const double m_scale = c.mass_flow_scale();
    double norm = 0.0;

    for (std::size_t i = 0; i < c.cells; ++i) {
        const double flow_right =
            (i + 1 < c.cells) ? state.flow[i + 1] : c.outlet_mass_flow;
        norm = std::max(norm, std::abs(flow_right - state.flow[i])/m_scale);

        double p_left{};
        double p_right{};
        double h_left{};
        double h_right{};
        double distance{};
        if (i == 0) {
            p_left = c.inlet_pressure;
            p_right = state.pressure[0];
            h_left = h[0];
            h_right = 0.5*(h[0] + h[1]);
            distance = 0.5*dx;
        } else {
            p_left = state.pressure[i - 1];
            p_right = state.pressure[i];
            h_left = 0.5*(h[i - 1] + h[i]);
            h_right = 0.5*(h[i] + h[i + 1]);
            distance = dx;
        }
        const double momentum = face_residual(
            c, state, i, p_left, p_right, h_left, h_right,
            distance, friction);
        norm = std::max(norm, std::abs(momentum)*dx/p_scale);
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
    if (initial_guess.pressure.size() != c.cells ||
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
        y[i] = initial_guess.pressure[i]/p_scale;
        y[n + i] = initial_guess.flow[i]/m_scale;
    }

    for (std::size_t iteration = 0; iteration <= max_iterations; ++iteration) {
        const std::vector<double> residual = residual_vector(c, y, friction);
        const double norm = infinity_norm(residual);
        if (norm <= tolerance) {
            State state;
            state.pressure.resize(n);
            state.flow.resize(n);
            for (std::size_t i = 0; i < n; ++i) {
                state.pressure[i] = y[i]*p_scale;
                state.flow[i] = y[n + i]*m_scale;
            }
            state.outlet_pressure = outlet_pressure_from_last_cell(
                c, friction, state.pressure.back());
            return {std::move(state), iteration, norm};
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
            const auto shifted = residual_vector(c, perturbed, friction);
            for (std::size_t row = 0; row < 2*n; ++row) {
                jacobian[row][column] = (shifted[row] - residual[row])/step;
            }
        }
        std::vector<double> rhs_value(2*n);
        std::transform(residual.begin(), residual.end(), rhs_value.begin(),
                       [](double value) { return -value; });
        const std::vector<double> correction =
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
                        residual_vector(c, candidate, friction);
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
        !std::isfinite(options.cfl) || options.cfl <= 0.0 || options.cfl > 0.5) {
        throw std::invalid_argument("Invalid relaxation options");
    }

    const double friction =
        swamee_jain_friction_factor(c, c.outlet_mass_flow);
    State state;
    state.pressure.assign(c.cells, c.inlet_pressure);
    state.flow.assign(c.cells, c.outlet_mass_flow);
    const double acoustic_speed = std::sqrt(c.gas_constant_product());
    const double acoustic_dt = options.cfl*c.dx()/acoustic_speed;
    std::size_t consecutive = 0;
    double pseudo_time = 0.0;
    double scaled_rhs = scaled_rhs_norm(c, state, friction);
    double dt = acoustic_dt;

    for (std::size_t step = 1; step <= options.max_steps; ++step) {
        double maximum_speed = 0.0;
        for (std::size_t i = 0; i < c.cells; ++i) {
            const double p_left = (i == 0) ? c.inlet_pressure
                                           : state.pressure[i - 1];
            const double p_right = state.pressure[i];
            const double density = 0.5*(p_left + p_right)/
                                   c.gas_constant_product();
            maximum_speed = std::max(
                maximum_speed,
                std::abs(state.flow[i])/(c.area()*density));
        }
        const double friction_rate =
            friction*maximum_speed/c.diameter;
        const double friction_dt = friction_rate > 0.0
            ? 0.2/friction_rate
            : std::numeric_limits<double>::infinity();
        dt = std::min(acoustic_dt, friction_dt);
        if (pseudo_time + dt > options.max_pseudo_time) {
            break;
        }
        state = rk4_step(c, state, friction, dt);
        pseudo_time += dt;
        scaled_rhs = scaled_rhs_norm(c, state, friction);
        if (scaled_rhs <= options.tolerance) {
            ++consecutive;
            if (consecutive >= options.consecutive_converged_steps) {
                state.outlet_pressure = outlet_pressure_from_last_cell(
                    c, friction, state.pressure.back());
                return {std::move(state), step, pseudo_time, scaled_rhs, dt};
            }
        } else {
            consecutive = 0;
        }
    }
    throw std::runtime_error(
        "Pseudo-time relaxation did not reach steady tolerance; increase limits or inspect the case");
}

} // namespace pipe_init
