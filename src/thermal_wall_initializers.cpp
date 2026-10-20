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
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pipe_init::thermal_wall {
namespace {
constexpr double pi = 3.14159265358979323846;
double cp(const Configuration& c) { return c.gamma*c.gas_constant/(c.gamma-1.0); }
double cv(const Configuration& c) { return c.gas_constant/(c.gamma-1.0); }
double rho(const Configuration& c, double p, double t) { return p/(c.gas_constant*t); }
double ueff(const Configuration& c)
{
    const double hi = c.inner_conductance_per_length();
    return hi*c.soil_conductance_per_length/(hi+c.soil_conductance_per_length);
}
double ell(const Configuration& c, std::size_t i) { return i == c.cells ? c.dx()/2.0 : c.dx(); }
double wall_rhs(const Configuration& c, double tg, double tw)
{
    return c.inner_conductance_per_length()*(tg-tw) +
           c.soil_conductance_per_length*(c.soil_temperature-tw);
}
double friction(const Configuration& c, double m)
{
    const double re = std::abs(m)*c.diameter/(c.area()*c.viscosity);
    if (re <= 0.0) return 0.0;
    if (re < 2300.0) return 64.0/re;
    const double v = c.roughness/(3.7*c.diameter)+5.74/std::pow(re,0.9);
    return 0.25/std::pow(std::log10(v),2.0);
}
double mom(const Configuration& c, const State& s, std::size_t i,
           const std::vector<double>& z)
{
    const double rl=rho(c,s.pressure[i],s.gas_temperature[i]);
    const double rr=rho(c,s.pressure[i+1],s.gas_temperature[i+1]);
    const double r=0.5*(rl+rr), v=s.flow[i]/(c.area()*r);
    return (s.pressure[i+1]-s.pressure[i])/c.dx()
        + r*c.gravity*(z[i+1]-z[i])/c.dx()
        + friction(c,s.flow[i])/(2.0*c.diameter)*r*v*std::abs(v);
}
std::vector<double> encode(const Configuration& c,const State& s)
{
    const std::size_t n=c.cells;
    if(s.pressure.size()!=n+1||s.gas_temperature.size()!=n+1||
       s.wall_temperature.size()!=n+1||s.flow.size()!=n)
        throw std::invalid_argument("thermal state dimensions do not match grid");
    std::vector<double> y(4*n+1);
    for(std::size_t i=0;i<n;++i) {
        y[i]=s.pressure[i+1]/c.inlet_pressure;
        y[n+i]=s.gas_temperature[i+1]/c.inlet_temperature;
        y[2*n+i]=s.flow[i]/c.outlet_mass_flow;
    }
    for(std::size_t i=0;i<=n;++i) y[3*n+i]=s.wall_temperature[i]/c.inlet_temperature;
    return y;
}
State decode(const Configuration& c,const std::vector<double>& y)
{
    const std::size_t n=c.cells;
    if(y.size()!=4*n+1) throw std::invalid_argument("wrong thermal unknown size");
    State s;
    s.pressure.resize(n+1); s.gas_temperature.resize(n+1);
    s.wall_temperature.resize(n+1); s.flow.resize(n);
    s.pressure[0]=c.inlet_pressure; s.gas_temperature[0]=c.inlet_temperature;
    for(std::size_t i=0;i<n;++i) {
        s.pressure[i+1]=y[i]*c.inlet_pressure;
        s.gas_temperature[i+1]=y[n+i]*c.inlet_temperature;
        s.flow[i]=y[2*n+i]*c.outlet_mass_flow;
    }
    for(std::size_t i=0;i<=n;++i) s.wall_temperature[i]=y[3*n+i]*c.inlet_temperature;
    return s;
}
double norm(const std::vector<double>& r)
{
    double v=0.0; for(double x:r) v=std::max(v,std::abs(x)); return v;
}
std::vector<double> steady_r(const Configuration& c,const std::vector<double>& y)
{
    const std::size_t n=c.cells; const State s=decode(c,y);
    const auto z=c.elevation_nodes();
    const double ms=std::max(1.0,std::abs(c.outlet_mass_flow));
    const double ps=c.inlet_pressure/c.dx(), es=ms*cp(c)*c.inlet_temperature;
    const double hs=(c.inner_conductance_per_length()+c.soil_conductance_per_length)
                    *c.inlet_temperature;
    std::vector<double> r(4*n+1);
    for(std::size_t i=1;i<n;++i) r[i-1]=(s.flow[i]-s.flow[i-1])/ms;
    r[n-1]=(c.outlet_mass_flow-s.flow[n-1])/ms;
    for(std::size_t i=0;i<n;++i) r[n+i]=mom(c,s,i,z)/ps;
    for(std::size_t i=1;i<=n;++i) {
        const double ml=s.flow[i-1], mr=(i<n?s.flow[i]:c.outlet_mass_flow);
        const double flux=cp(c)*(mr*s.gas_temperature[i]-ml*s.gas_temperature[i-1]);
        const double q=c.inner_conductance_per_length()*
                       (s.wall_temperature[i]-s.gas_temperature[i]);
        r[2*n+i-1]=(flux-q*ell(c,i))/es;
    }
    for(std::size_t i=0;i<=n;++i)
        r[3*n+i]=wall_rhs(c,s.gas_temperature[i],s.wall_temperature[i])/hs;
    return r;
}
std::vector<double> linear_solve(std::vector<std::vector<double>> a,std::vector<double> b)
{
    const std::size_t n=b.size();
    for(std::size_t k=0;k<n;++k) {
        std::size_t p=k;
        for(std::size_t i=k+1;i<n;++i) if(std::abs(a[i][k])>std::abs(a[p][k])) p=i;
        if(std::abs(a[p][k])<1e-14) throw std::runtime_error("singular thermal Jacobian");
        std::swap(a[k],a[p]); std::swap(b[k],b[p]);
        const double d=a[k][k];
        for(std::size_t j=k;j<n;++j) a[k][j]/=d;
        b[k]/=d;
        for(std::size_t i=0;i<n;++i) if(i!=k) {
            const double f=a[i][k];
            for(std::size_t j=k;j<n;++j) a[i][j]-=f*a[k][j];
            b[i]-=f*b[k];
        }
    }
    return b;
}
template<class F,class V>
std::pair<std::size_t,double> newton(std::vector<double>& y,F residual,V valid,
                                     double tol,std::size_t maxit)
{
    for(std::size_t it=0;it<=maxit;++it) {
        const auto r=residual(y); const double old=norm(r);
        if(old<=tol) return {it,old};
        if(it==maxit) break;
        const std::size_t n=y.size();
        std::vector<std::vector<double>> j(n,std::vector<double>(n));
        for(std::size_t k=0;k<n;++k) {
            auto yp=y; const double h=2e-6*std::max(1.0,std::abs(y[k]));
            yp[k]+=h; const auto rp=residual(yp);
            for(std::size_t i=0;i<n;++i) j[i][k]=(rp[i]-r[i])/h;
        }
        std::vector<double> b(r.size());
        std::transform(r.begin(),r.end(),b.begin(),[](double x){return -x;});
        const auto d=linear_solve(std::move(j),std::move(b));
        bool ok=false;
        for(double a=1.0;a>=1.0/1024.0;a*=0.5) {
            auto trial=y; for(std::size_t k=0;k<n;++k) trial[k]+=a*d[k];
            if(!valid(trial)) continue;
            try {
                if(norm(residual(trial))<old) { y=std::move(trial); ok=true; break; }
            } catch(const std::exception&) {}
        }
        if(!ok) throw std::runtime_error("thermal Newton line search failed");
    }
    throw std::runtime_error("thermal Newton did not converge");
}
bool valid(const Configuration& c,const std::vector<double>& y)
{
    const std::size_t n=c.cells;
    for(std::size_t i=0;i<n;++i) if(y[i]<=0.05||y[n+i]<=0.1) return false;
    for(std::size_t i=0;i<=n;++i) if(y[3*n+i]<=0.1) return false;
    return true;
}
double pressure_step(const Configuration& c,const State& s,std::size_t i,
                     const std::vector<double>& z)
{
    const double pl=s.pressure[i], t=0.5*(s.gas_temperature[i]+s.gas_temperature[i+1]);
    const double zs=(z[i+1]-z[i])/c.dx();
    const double a=2.0*c.gravity*zs/(c.gas_constant*t);
    const double b=friction(c,s.flow[i])*c.gas_constant*t*s.flow[i]*std::abs(s.flow[i])/
                    (c.diameter*c.area()*c.area());
    const double p2=std::abs(a)<1e-10 ? pl*pl-b*c.dx()
        : std::exp(-a*c.dx())*pl*pl-b*(-std::expm1(-a*c.dx()))/a;
    if(!(p2>0.0)) throw std::runtime_error("thermal pressure profile is inadmissible");
    return std::sqrt(p2);
}
double discrete_pressure_step(const Configuration& c, State& s, std::size_t i,
                              const std::vector<double>& z)
{
    double p = pressure_step(c, s, i, z);
    for (std::size_t it = 0; it < 20; ++it) {
        s.pressure[i+1] = p;
        const double r = mom(c, s, i, z);
        if (std::abs(r)*c.dx()/c.inlet_pressure < 1.0e-13) return p;
        const double h = std::max(1.0, 1.0e-6*p);
        s.pressure[i+1] = p + h;
        const double derivative = (mom(c, s, i, z) - r)/h;
        s.pressure[i+1] = p;
        if (std::abs(derivative) < 1.0e-20)
            throw std::runtime_error("singular local momentum solve in marching");
        const double correction = -r/derivative;
        bool accepted = false;
        for (double a = 1.0; a >= 1.0/1024.0; a *= 0.5) {
            const double trial = p + a*correction;
            if (trial <= 0.05*c.inlet_pressure) continue;
            s.pressure[i+1] = trial;
            if (std::abs(mom(c, s, i, z)) < std::abs(r)) {
                p = trial;
                accepted = true;
                break;
            }
        }
        s.pressure[i+1] = p;
        if (!accepted) throw std::runtime_error("local momentum marching failed");
    }
    throw std::runtime_error("local momentum marching did not converge");
}
std::vector<double> theta_r(const Configuration& c,const State& old,
                            const std::vector<double>& y,double dt,double theta)
{
    const std::size_t n=c.cells; const State s=decode(c,y); const auto z=c.elevation_nodes();
    const double ms=std::max(1.0,std::abs(c.outlet_mass_flow));
    const double es=ms*cp(c)*c.inlet_temperature;
    const double hs=(c.inner_conductance_per_length()+c.soil_conductance_per_length)
                    *c.inlet_temperature;
    std::vector<double> r(4*n+1);
    auto mt=[&](std::size_t i){return theta*s.flow[i]+(1.0-theta)*old.flow[i];};
    for(std::size_t i=1;i<=n;++i) {
        const double rn=rho(c,s.pressure[i],s.gas_temperature[i]);
        const double ro=rho(c,old.pressure[i],old.gas_temperature[i]);
        const double ml=mt(i-1), mr=i<n?mt(i):c.outlet_mass_flow;
        r[i-1]=c.area()*ell(c,i)*(rn-ro)/(dt*ms)+(mr-ml)/ms;
    }
    for(std::size_t i=0;i<n;++i) {
        const double mn=mom(c,s,i,z), mo=mom(c,old,i,z);
        r[n+i]=(s.flow[i]-old.flow[i])/ms+c.area()*dt*
               (theta*mn+(1.0-theta)*mo)/ms;
    }
    for(std::size_t i=1;i<=n;++i) {
        const double length=ell(c,i);
        const double rn=rho(c,s.pressure[i],s.gas_temperature[i]);
        const double ro=rho(c,old.pressure[i],old.gas_temperature[i]);
        const double storage=c.area()*length*cv(c)*(rn*s.gas_temperature[i]-
                           ro*old.gas_temperature[i])/dt;
        const std::size_t f=i-1;
        const double fl=theta*s.flow[f]*cp(c)*s.gas_temperature[f]+
                        (1.0-theta)*old.flow[f]*cp(c)*old.gas_temperature[f];
        const double fr=i<n
            ? theta*s.flow[i]*cp(c)*s.gas_temperature[i]+
              (1.0-theta)*old.flow[i]*cp(c)*old.gas_temperature[i]
            : c.outlet_mass_flow*cp(c)*(theta*s.gas_temperature[i]+
                        (1.0-theta)*old.gas_temperature[i]);
        const double qn=c.inner_conductance_per_length()*
                        (s.wall_temperature[i]-s.gas_temperature[i]);
        const double qo=c.inner_conductance_per_length()*
                        (old.wall_temperature[i]-old.gas_temperature[i]);
        r[2*n+i-1]=(storage+fr-fl-(theta*qn+(1.0-theta)*qo)*length)/es;
    }
    for(std::size_t i=0;i<=n;++i) {
        const double wn=wall_rhs(c,s.gas_temperature[i],s.wall_temperature[i]);
        const double wo=wall_rhs(c,old.gas_temperature[i],old.wall_temperature[i]);
        r[3*n+i]=(c.wall_heat_capacity_per_length()*
                   (s.wall_temperature[i]-old.wall_temperature[i])/dt-
                   theta*wn-(1.0-theta)*wo)/hs;
    }
    return r;
}
} // namespace

double Configuration::area() const { return pi*diameter*diameter/4.0; }
double Configuration::dx() const { return length/static_cast<double>(cells); }
double Configuration::inner_conductance_per_length() const
{ return inner_heat_transfer_coefficient*pi*diameter; }
double Configuration::wall_heat_capacity_per_length() const
{
    const double d=diameter+2.0*wall_thickness;
    return wall_density*wall_heat_capacity*pi*(d*d-diameter*diameter)/4.0;
}
std::vector<double> Configuration::elevation_nodes() const
{
    std::vector<double> z(cells+1);
    for(std::size_t i=0;i<=cells;++i) {
        const double x=static_cast<double>(i)*dx();
        z[i]=0.5*hill_height*(1.0-std::cos(2.0*pi*x/length));
    }
    return z;
}
void Configuration::validate() const
{
    if(cells<2||length<=0||diameter<=0||roughness<0||viscosity<=0||
       gas_constant<=0||gamma<=1||inlet_pressure<=0||inlet_temperature<=0||
       soil_temperature<=0||outlet_mass_flow<=0||
       inner_heat_transfer_coefficient<=0||soil_conductance_per_length<=0||
       wall_density<=0||wall_heat_capacity<=0||wall_thickness<=0||hill_height<0)
        throw std::invalid_argument("invalid thermal-wall configuration");
}
double swamee_jain_friction_factor(const Configuration& c,double m)
{ return friction(c,m); }

State analytical_approximation(const Configuration& c)
{
    c.validate(); const std::size_t n=c.cells; State s;
    s.pressure.resize(n+1); s.gas_temperature.resize(n+1);
    s.wall_temperature.resize(n+1); s.flow.assign(n,c.outlet_mass_flow);
    s.pressure[0]=c.inlet_pressure; s.gas_temperature[0]=c.inlet_temperature;
    const double u=ueff(c), hi=c.inner_conductance_per_length(), gs=c.soil_conductance_per_length;
    for(std::size_t i=1;i<=n;++i) {
        const double x=static_cast<double>(i)*c.dx();
        s.gas_temperature[i]=c.soil_temperature+
            (c.inlet_temperature-c.soil_temperature)*
            std::exp(-u*x/(c.outlet_mass_flow*cp(c)));
    }
    const auto z=c.elevation_nodes();
    for(std::size_t i=0;i<=n;++i)
        s.wall_temperature[i]=(hi*s.gas_temperature[i]+gs*c.soil_temperature)/(hi+gs);
    for(std::size_t i=0;i<n;++i) s.pressure[i+1]=pressure_step(c,s,i,z);
    return s;
}
State steady_spatial_marching(const Configuration& c)
{
    c.validate(); const std::size_t n=c.cells; State s;
    s.pressure.resize(n+1); s.gas_temperature.resize(n+1);
    s.wall_temperature.resize(n+1); s.flow.assign(n,c.outlet_mass_flow);
    s.pressure[0]=c.inlet_pressure; s.gas_temperature[0]=c.inlet_temperature;
    const double hi=c.inner_conductance_per_length(), gs=c.soil_conductance_per_length;
    const double u=ueff(c);
    s.wall_temperature[0]=(hi*c.inlet_temperature+gs*c.soil_temperature)/(hi+gs);
    for(std::size_t i=1;i<=n;++i) {
        const double length=ell(c,i);
        s.gas_temperature[i]=(c.outlet_mass_flow*cp(c)*s.gas_temperature[i-1]+
            u*length*c.soil_temperature)/(c.outlet_mass_flow*cp(c)+u*length);
        s.wall_temperature[i]=(hi*s.gas_temperature[i]+gs*c.soil_temperature)/(hi+gs);
    }
    const auto z=c.elevation_nodes();
    for(std::size_t i=0;i<n;++i) s.pressure[i+1]=discrete_pressure_step(c,s,i,z);
    return s;
}
SolveInfo coupled_steady_newton(const Configuration& c,const State& guess,
                                double tol,std::size_t maxit)
{
    c.validate(); auto y=encode(c,guess);
    const auto result=newton(y,[&](const auto& x){return steady_r(c,x);},
        [&](const auto& x){return valid(c,x);},tol,maxit);
    return {decode(c,y),result.first,result.second};
}
double scaled_steady_residual(const Configuration& c,const State& s)
{ return norm(steady_r(c,encode(c,s))); }
RelaxationInfo implicit_theta_relaxation(const Configuration& c,const RelaxationOptions& o)
{
    c.validate();
    if(o.theta<0.5||o.theta>1.0||o.pseudo_courant<=0||o.tolerance<=0)
        throw std::invalid_argument("invalid thermal theta options");
    const double dt=o.pseudo_courant*c.dx()/
        std::sqrt(c.gamma*c.gas_constant*c.inlet_temperature);
    const std::size_t n=c.cells; State s;
    s.pressure.assign(n+1,c.inlet_pressure);
    s.gas_temperature.assign(n+1,c.inlet_temperature);
    s.wall_temperature.assign(n+1,c.soil_temperature);
    s.flow.assign(n,c.outlet_mass_flow);
    std::size_t total=0,good=0; double time=0.0,res=scaled_steady_residual(c,s);
    for(std::size_t step=1;step<=o.max_steps;++step) {
        const State old=s; auto y=encode(c,old);
        const auto info=newton(y,[&](const auto& x){return theta_r(c,old,x,dt,o.theta);},
            [&](const auto& x){return valid(c,x);},1e-10,18);
        total+=info.first; s=decode(c,y); time+=dt; res=scaled_steady_residual(c,s);
        good=res<=o.tolerance?good+1:0;
        if(good>=o.consecutive_converged_steps)
            return {s,step,total,time,dt,res};
        if(time>=o.max_pseudo_time) break;
    }
    throw std::runtime_error("thermal theta relaxation did not reach steady state");
}
} // namespace pipe_init::thermal_wall

