#include "pipe/adiabatic_initializers.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace pipe_init::adiabatic {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;

double density(const Configuration& c, double pressure, double temperature)
{
    return pressure/(c.gas_constant*temperature);
}

double velocity(const Configuration& c, double mass_flow,
                double pressure, double temperature)
{
    return mass_flow/(c.area()*density(c, pressure, temperature));
}

double cp(const Configuration& c)
{
    return c.gamma*c.gas_constant/(c.gamma - 1.0);
}

double total_enthalpy(const Configuration& c, double mass_flow,
                      double pressure, double temperature, double elevation)
{
    const double u = velocity(c, mass_flow, pressure, temperature);
    return cp(c)*temperature + 0.5*u*u + c.gravity*elevation;
}

double cell_momentum_residual(
    const Configuration& c, const State& state, std::size_t cell,
    double friction, const std::vector<double>& elevation)
{
    const double p_left = state.pressure[cell];
    const double p_right = state.pressure[cell + 1];
    const double t_left = state.temperature[cell];
    const double t_right = state.temperature[cell + 1];
    const double flow = state.flow[cell];
    const double rho_left = density(c, p_left, t_left);
    const double rho_right = density(c, p_right, t_right);
    const double rho_bar = 0.5*(rho_left + rho_right);
    const double u_left = flow/(c.area()*rho_left);
    const double u_right = flow/(c.area()*rho_right);
    const double pressure_gradient = (p_right - p_left)/c.dx();
    const double elevation_gradient =
        (elevation[cell + 1] - elevation[cell])/c.dx();
    const double convective_gradient =
        (rho_right*u_right*u_right - rho_left*u_left*u_left)/c.dx();
    const double drag = friction/(2.0*c.diameter)*rho_bar*
                        0.5*(u_left + u_right)*
                        std::abs(0.5*(u_left + u_right));
    return pressure_gradient + convective_gradient +
           rho_bar*c.gravity*elevation_gradient + drag;
}

std::vector<double> encode(const Configuration& c, const State& state)
{
    const std::size_t n = c.cells;
    if (state.pressure.size() != n + 1 ||
        state.temperature.size() != n + 1 || state.flow.size() != n) {
        throw std::invalid_argument("Adiabatic state dimensions do not match the grid");
    }
    const double p_scale = c.inlet_pressure;
    const double t_scale = c.inlet_temperature;
    const double m_scale = std::max(1.0, std::abs(c.outlet_mass_flow));
    std::vector<double> y(3*n);
    for (std::size_t i = 0; i < n; ++i) {
        y[i] = state.pressure[i + 1]/p_scale;
        y[n + i] = state.temperature[i + 1]/t_scale;
        y[2*n + i] = state.flow[i]/m_scale;
    }
    return y;
}

State decode(const Configuration& c, const std::vector<double>& y)
{
    const std::size_t n = c.cells;
    if (y.size() != 3*n) {
        throw std::invalid_argument("Adiabatic unknown vector has wrong size");
    }
    const double p_scale = c.inlet_pressure;
    const double t_scale = c.inlet_temperature;
    const double m_scale = std::max(1.0, std::abs(c.outlet_mass_flow));
    State state;
    state.pressure.resize(n + 1);
    state.temperature.resize(n + 1);
    state.flow.resize(n);
    state.pressure[0] = c.inlet_pressure;
    state.temperature[0] = c.inlet_temperature;
    for (std::size_t i = 0; i < n; ++i) {
        state.pressure[i + 1] = y[i]*p_scale;
        state.temperature[i + 1] = y[n + i]*t_scale;
        state.flow[i] = y[2*n + i]*m_scale;
    }
    return state;
}

std::vector<double> residual_vector(
    const Configuration& c, const std::vector<double>& y)
{
    const std::size_t n = c.cells;
    const State state = decode(c, y);
    const auto elevation = c.elevation_nodes();
    const double p_scale = c.inlet_pressure;
    const double t_scale = c.inlet_temperature;
    const double m_scale = std::max(1.0, std::abs(c.outlet_mass_flow));
    const double h_scale = cp(c)*t_scale;
    std::vector<double> r(3*n, 0.0);

    // Continuity: each internal face flow must equal the specified outlet flow.
    for (std::size_t i = 0; i < n; ++i) {
        const double right_flow = i + 1 < n
            ? state.flow[i + 1] : c.outlet_mass_flow;
        r[i] = (right_flow - state.flow[i])/m_scale;
    }
    for (std::size_t i = 0; i < n; ++i) {
        const double f = swamee_jain_friction_factor(c, state.flow[i]);
        r[n + i] = cell_momentum_residual(c, state, i, f, elevation)
                   *c.dx()/p_scale;
        const double h_left = total_enthalpy(
            c, state.flow[i], state.pressure[i], state.temperature[i], elevation[i]);
        const double h_right = total_enthalpy(
            c, state.flow[i], state.pressure[i + 1], state.temperature[i + 1],
            elevation[i + 1]);
        r[2*n + i] = (h_right - h_left)/h_scale;
    }
    return r;
}

double norm_inf(const std::vector<double>& x)
{
    double result = 0.0;
    for (double value : x) result = std::max(result, std::abs(value));
    return result;
}

std::vector<double> solve_dense(std::vector<std::vector<double>> a,
                                std::vector<double> b)
{
    const std::size_t n = b.size();
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(a[i][k]) > std::abs(a[pivot][k])) pivot = i;
        }
        if (std::abs(a[pivot][k]) < 1.0e-14) {
            throw std::runtime_error("Singular adiabatic Newton Jacobian");
        }
        std::swap(a[k], a[pivot]);
        std::swap(b[k], b[pivot]);
        const double d = a[k][k];
        for (std::size_t j = k; j < n; ++j) a[k][j] /= d;
        b[k] /= d;
        for (std::size_t i = 0; i < n; ++i) {
            if (i == k) continue;
            const double factor = a[i][k];
            for (std::size_t j = k; j < n; ++j) a[i][j] -= factor*a[k][j];
            b[i] -= factor*b[k];
        }
    }
    return b;
}

template <typename Evaluator>
std::vector<double> finite_difference_correction(
    const std::vector<double>& y, const std::vector<double>& residual,
    Evaluator&& evaluate)
{
    const std::size_t n = y.size();
    std::vector<std::vector<double>> jac(n, std::vector<double>(n, 0.0));
    for (std::size_t j = 0; j < n; ++j) {
        auto shifted_y = y;
        const double step = 1.0e-6*std::max(1.0, std::abs(y[j]));
        shifted_y[j] += step;
        const auto shifted = evaluate(shifted_y);
        for (std::size_t i = 0; i < n; ++i) {
            jac[i][j] = (shifted[i] - residual[i])/step;
        }
    }
    std::vector<double> rhs(residual.size());
    std::transform(residual.begin(), residual.end(), rhs.begin(),
                   [](double value) { return -value; });
    return solve_dense(std::move(jac), std::move(rhs));
}

std::vector<double> newton_correction(
    const Configuration& c, const std::vector<double>& y,
    const std::vector<double>& residual)
{
    return finite_difference_correction(y, residual,
        [&c](const std::vector<double>& value) {
            return residual_vector(c, value);
        });
}

bool positive_thermodynamic_state(const Configuration& c,
                                  const std::vector<double>& y)
{
    for (std::size_t i = 0; i < c.cells; ++i) {
        if (!(y[i] > 0.0) || !(y[c.cells + i] > 0.0) ||
            !std::isfinite(y[i]) || !std::isfinite(y[c.cells + i])) return false;
    }
    return true;
}

std::vector<double> physical_pseudo_rhs(const Configuration& c,
                                        const State& state)
{
    const std::size_t n = c.cells;
    const auto elevation = c.elevation_nodes();
    const double cv = c.gas_constant/(c.gamma - 1.0);
    const double p_scale = c.inlet_pressure;
    const double t_scale = c.inlet_temperature;
    const double m_scale = std::max(1.0, std::abs(c.outlet_mass_flow));

    std::vector<double> rho(n + 1);
    std::vector<double> node_flow(n + 1);
    std::vector<double> node_flow_dot(n + 1, 0.0);
    std::vector<double> u(n + 1);
    std::vector<double> mdot(n, 0.0);
    for (std::size_t i = 0; i <= n; ++i) {
        rho[i] = density(c, state.pressure[i], state.temperature[i]);
        if (i == 0) {
            node_flow[i] = state.flow[0];
        } else if (i == n) {
            node_flow[i] = 0.5*(state.flow[n - 1] + c.outlet_mass_flow);
        } else {
            node_flow[i] = 0.5*(state.flow[i - 1] + state.flow[i]);
        }
        u[i] = node_flow[i]/(c.area()*rho[i]);
    }

    for (std::size_t i = 0; i < n; ++i) {
        const double f = swamee_jain_friction_factor(c, state.flow[i]);
        mdot[i] = -c.area()*cell_momentum_residual(c, state, i, f, elevation);
    }
    for (std::size_t i = 0; i <= n; ++i) {
        if (i == 0) node_flow_dot[i] = mdot[0];
        else if (i == n) node_flow_dot[i] = 0.5*mdot[n - 1];
        else node_flow_dot[i] = 0.5*(mdot[i - 1] + mdot[i]);
    }

    std::vector<double> h_total(n + 1);
    for (std::size_t i = 0; i <= n; ++i) {
        h_total[i] = cp(c)*state.temperature[i] + 0.5*u[i]*u[i] +
                       c.gravity*elevation[i];
    }

    std::vector<double> p_dot(n + 1, 0.0);
    std::vector<double> t_dot(n + 1, 0.0);
    for (std::size_t i = 1; i <= n; ++i) {
        const double left_flow = state.flow[i - 1];
        const double right_flow = i < n ? state.flow[i] : c.outlet_mass_flow;
        const double left_energy = left_flow >= 0.0
            ? left_flow*h_total[i - 1] : left_flow*h_total[i];
        const double right_energy = i < n
            ? (right_flow >= 0.0 ? right_flow*h_total[i]
                                 : right_flow*h_total[i + 1])
            : right_flow*h_total[i];
        const double volume_length = i < n ? c.dx() : 0.5*c.dx();
        const double rho_dot = -(right_flow - left_flow)/
                               (c.area()*volume_length);
        const double energy_dot = -(right_energy - left_energy)/
                                  (c.area()*volume_length);
        const double u_dot = node_flow_dot[i]/(c.area()*rho[i]) -
                             u[i]*rho_dot/rho[i];
        const double stored_energy_slope = cv*state.temperature[i] +
            c.gravity*elevation[i] + 0.5*u[i]*u[i];
        t_dot[i] = (energy_dot - stored_energy_slope*rho_dot -
                    rho[i]*u[i]*u_dot)/(rho[i]*cv);
        p_dot[i] = c.gas_constant*(state.temperature[i]*rho_dot +
                                   rho[i]*t_dot[i]);
    }

    std::vector<double> g(3*n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        g[i] = p_dot[i + 1]/p_scale;
        g[n + i] = t_dot[i + 1]/t_scale;
        g[2*n + i] = mdot[i]/m_scale;
    }
    return g;
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

std::vector<double> Configuration::elevation_nodes() const
{
    std::vector<double> z(cells + 1);
    for (std::size_t i = 0; i <= cells; ++i) {
        const double x = length*static_cast<double>(i)/static_cast<double>(cells);
        z[i] = 0.5*hill_height*(1.0 - std::cos(2.0*kPi*x/length));
    }
    return z;
}

void Configuration::validate() const
{
    const auto positive = [](double x) { return std::isfinite(x) && x > 0.0; };
    if (cells < 2 || !positive(length) || !positive(diameter) ||
        !positive(viscosity) || !positive(gas_constant) || !positive(gamma - 1.0) ||
        !positive(inlet_pressure) || !positive(inlet_temperature) ||
        !std::isfinite(roughness) || roughness < 0.0 ||
        !positive(outlet_mass_flow) || !std::isfinite(gravity) || gravity < 0.0 ||
        !std::isfinite(hill_height)) {
        throw std::invalid_argument("Invalid adiabatic pipe configuration");
    }
}

double swamee_jain_friction_factor(const Configuration& c, double mass_flow)
{
    c.validate();
    if (!std::isfinite(mass_flow)) throw std::invalid_argument("Nonfinite mass flow");
    if (mass_flow == 0.0) return 0.0;
    const double reynolds = std::abs(mass_flow)*c.diameter/(c.area()*c.viscosity);
    const double argument = c.roughness/(3.7*c.diameter) +
                            5.74/std::pow(reynolds, 0.9);
    const double log_argument = std::log10(argument);
    const double factor = 0.25/(log_argument*log_argument);
    if (!std::isfinite(factor) || factor <= 0.0) {
        throw std::runtime_error("Invalid Swamee-Jain friction factor");
    }
    return factor;
}

double scaled_steady_residual(const Configuration& c, const State& state)
{
    c.validate();
    return norm_inf(residual_vector(c, encode(c, state)));
}

State steady_spatial_marching(const Configuration& c)
{
    c.validate();
    const double flow = c.outlet_mass_flow;
    const double friction = swamee_jain_friction_factor(c, flow);
    const auto elevation = c.elevation_nodes();
    State state;
    state.pressure.resize(c.cells + 1);
    state.temperature.resize(c.cells + 1);
    state.flow.assign(c.cells, flow);
    state.pressure[0] = c.inlet_pressure;
    state.temperature[0] = c.inlet_temperature;

    for (std::size_t cell = 0; cell < c.cells; ++cell) {
        double p = state.pressure[cell]*(1.0 - 1.0e-3);
        double t = state.temperature[cell];
        bool converged = false;
        for (std::size_t iteration = 0; iteration < 20; ++iteration) {
            State trial = state;
            trial.pressure[cell + 1] = p;
            trial.temperature[cell + 1] = t;
            const double r_m = cell_momentum_residual(c, trial, cell,
                                                       friction, elevation)*c.dx()/c.inlet_pressure;
            const double h_left = total_enthalpy(c, flow, state.pressure[cell],
                state.temperature[cell], elevation[cell]);
            const double h_right = total_enthalpy(c, flow, p, t, elevation[cell + 1]);
            const double r_e = (h_right - h_left)/(cp(c)*c.inlet_temperature);
            const double residual_norm = std::max(std::abs(r_m), std::abs(r_e));
            if (residual_norm < 2.0e-13) { converged = true; break; }

            constexpr double eps = 1.0e-6;
            State p_trial = trial;
            p_trial.pressure[cell + 1] += eps*c.inlet_pressure;
            const double rp_m = cell_momentum_residual(c, p_trial, cell,
                friction, elevation)*c.dx()/c.inlet_pressure;
            const double rp_h = (total_enthalpy(c, flow, p_trial.pressure[cell + 1], t,
                elevation[cell + 1]) - h_left)/(cp(c)*c.inlet_temperature);
            State t_trial = trial;
            t_trial.temperature[cell + 1] += eps*c.inlet_temperature;
            const double rt_m = cell_momentum_residual(c, t_trial, cell,
                friction, elevation)*c.dx()/c.inlet_pressure;
            const double rt_h = (total_enthalpy(c, flow, p, t_trial.temperature[cell + 1],
                elevation[cell + 1]) - h_left)/(cp(c)*c.inlet_temperature);
            const double j00 = (rp_m - r_m)/eps;
            const double j10 = (rp_h - r_e)/eps;
            const double j01 = (rt_m - r_m)/eps;
            const double j11 = (rt_h - r_e)/eps;
            const double det = j00*j11 - j01*j10;
            if (std::abs(det) < 1.0e-14) break;
            const double dp_scaled = (-r_m*j11 + j01*r_e)/det;
            const double dt_scaled = (-j00*r_e + j10*r_m)/det;
            double lambda = 1.0;
            bool accepted = false;
            for (int ls = 0; ls < 20; ++ls) {
                const double p_next = p + lambda*dp_scaled*c.inlet_pressure;
                const double t_next = t + lambda*dt_scaled*c.inlet_temperature;
                if (p_next > 0.0 && t_next > 0.0) {
                    State candidate = trial;
                    candidate.pressure[cell + 1] = p_next;
                    candidate.temperature[cell + 1] = t_next;
                    const double cm = cell_momentum_residual(c, candidate, cell,
                        friction, elevation)*c.dx()/c.inlet_pressure;
                    const double ce = (total_enthalpy(c, flow, p_next, t_next,
                        elevation[cell + 1]) - h_left)/(cp(c)*c.inlet_temperature);
                    if (std::max(std::abs(cm), std::abs(ce)) < residual_norm) {
                        p = p_next; t = t_next; accepted = true; break;
                    }
                }
                lambda *= 0.5;
            }
            if (!accepted) break;
        }
        if (!converged) {
            State check = state;
            check.pressure[cell + 1] = p;
            check.temperature[cell + 1] = t;
            const double rm = std::abs(cell_momentum_residual(c, check, cell,
                friction, elevation)*c.dx()/c.inlet_pressure);
            const double re = std::abs((total_enthalpy(c, flow, p, t, elevation[cell + 1])-
                total_enthalpy(c, flow, state.pressure[cell], state.temperature[cell],
                               elevation[cell]))/(cp(c)*c.inlet_temperature));
            if (std::max(rm, re) >= 2.0e-11) {
                throw std::runtime_error("Adiabatic spatial march failed in a cell");
            }
        }
        state.pressure[cell + 1] = p;
        state.temperature[cell + 1] = t;
    }
    return state;
}

SolveInfo coupled_steady_newton(const Configuration& c, const State& initial_guess,
                                double tolerance, std::size_t max_iterations)
{
    c.validate();
    if (!std::isfinite(tolerance) || tolerance <= 0.0 || max_iterations == 0) {
        throw std::invalid_argument("Invalid coupled Newton options");
    }
    auto y = encode(c, initial_guess);
    for (std::size_t k = 0; k <= max_iterations; ++k) {
        const auto r = residual_vector(c, y);
        const double norm = norm_inf(r);
        if (norm <= tolerance) return {decode(c, y), k, norm};
        if (k == max_iterations) break;
        const auto dy = newton_correction(c, y, r);
        bool accepted = false;
        double lambda = 1.0;
        for (int ls = 0; ls < 24; ++ls) {
            auto candidate = y;
            for (std::size_t j = 0; j < y.size(); ++j) candidate[j] += lambda*dy[j];
            if (positive_thermodynamic_state(c, candidate)) {
                try {
                    if (norm_inf(residual_vector(c, candidate)) < norm) {
                        y = std::move(candidate); accepted = true; break;
                    }
                } catch (const std::runtime_error&) { }
            }
            lambda *= 0.5;
        }
        if (!accepted) throw std::runtime_error("Coupled adiabatic Newton line search failed");
    }
    throw std::runtime_error("Coupled adiabatic Newton solve did not converge");
}

RelaxationInfo implicit_theta_relaxation(const Configuration& c,
                                         const RelaxationOptions& o)
{
    c.validate();
    if (!std::isfinite(o.tolerance) || o.tolerance <= 0.0 ||
        !std::isfinite(o.theta) || o.theta < 0.5 || o.theta > 1.0 ||
        !std::isfinite(o.pseudo_courant) || o.pseudo_courant <= 0.0 ||
        !std::isfinite(o.max_pseudo_time) || o.max_pseudo_time <= 0.0 ||
        o.max_steps == 0 || o.consecutive_converged_steps == 0) {
        throw std::invalid_argument("Invalid adiabatic relaxation options");
    }
    State state;
    state.pressure.assign(c.cells + 1, c.inlet_pressure);
    state.temperature.assign(c.cells + 1, c.inlet_temperature);
    state.flow.assign(c.cells, c.outlet_mass_flow);
    auto y = encode(c, state);
    const double sound_speed = std::sqrt(c.gamma*c.gas_constant*c.inlet_temperature);
    const double dt = o.pseudo_courant*c.dx()/sound_speed;
    double pseudo_time = 0.0;
    std::size_t converged_steps = 0;
    std::size_t total_newton = 0;
    double norm = norm_inf(residual_vector(c, y));

    for (std::size_t step = 1; step <= o.max_steps; ++step) {
        if (pseudo_time + dt > o.max_pseudo_time) break;
        const auto old_rhs = physical_pseudo_rhs(c, decode(c, y));
        auto candidate = y;
        bool step_converged = false;
        for (std::size_t it = 0; it < 20; ++it) {
            const auto new_rhs = physical_pseudo_rhs(c, decode(c, candidate));
            std::vector<double> step_residual(y.size());
            for (std::size_t j = 0; j < y.size(); ++j) {
                step_residual[j] = candidate[j] - y[j] - dt*
                    (o.theta*new_rhs[j] + (1.0 - o.theta)*old_rhs[j]);
            }
            const double step_norm = norm_inf(step_residual);
            if (step_norm < 2.0e-11) { step_converged = true; break; }
            const auto correction = finite_difference_correction(
                candidate, step_residual, [&](const std::vector<double>& value) {
                    const auto value_rhs = physical_pseudo_rhs(c, decode(c, value));
                    std::vector<double> result(y.size());
                    for (std::size_t j = 0; j < y.size(); ++j) {
                        result[j] = value[j] - y[j] - dt*
                            (o.theta*value_rhs[j] + (1.0-o.theta)*old_rhs[j]);
                    }
                    return result;
                });
            bool accepted = false;
            double lambda = 1.0;
            for (int ls = 0; ls < 20; ++ls) {
                auto trial = candidate;
                for (std::size_t j = 0; j < y.size(); ++j) trial[j] += lambda*correction[j];
                if (positive_thermodynamic_state(c, trial)) {
                    try {
                        const auto tr_rhs = physical_pseudo_rhs(c, decode(c, trial));
                        std::vector<double> tr_step(y.size());
                        for (std::size_t j = 0; j < y.size(); ++j) {
                            tr_step[j] = trial[j] - y[j] - dt*
                                (o.theta*tr_rhs[j] + (1.0-o.theta)*old_rhs[j]);
                        }
                        if (norm_inf(tr_step) < step_norm) {
                            candidate = std::move(trial); accepted = true; break;
                        }
                    } catch (const std::runtime_error&) { }
                }
                lambda *= 0.5;
            }
            if (!accepted) break;
            ++total_newton;
        }
        if (!step_converged) throw std::runtime_error("Adiabatic theta pseudo-step failed");
        y = std::move(candidate);
        pseudo_time += dt;
        norm = norm_inf(residual_vector(c, y));
        if (norm <= o.tolerance) {
            ++converged_steps;
            if (converged_steps >= o.consecutive_converged_steps) {
                return {decode(c, y), step, pseudo_time, norm, dt, total_newton};
            }
        } else {
            converged_steps = 0;
        }
    }
    throw std::runtime_error("Adiabatic theta relaxation did not reach steady tolerance");
}

} // namespace pipe_init::adiabatic
