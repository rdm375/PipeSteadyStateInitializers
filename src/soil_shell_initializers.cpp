#include "pipe/soil_shell_initializers.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

namespace pipe_init::soil_shell {
namespace {
constexpr double pi = 3.14159265358979323846;
double cp(const Configuration& c) { return c.gamma*c.gas_constant/(c.gamma-1.0); }
double cv(const Configuration& c) { return c.gas_constant/(c.gamma-1.0); }
double rho(const Configuration& c,double p,double t) { return p/(c.gas_constant*t); }
double ell(const Configuration& c,std::size_t i)
{ return i==c.cells ? c.dx()/2.0 : c.dx(); }
std::vector<double> radial_edges(const Configuration& c)
{
    std::vector<double> e(c.soil_shell_radii.size()+1);
    e.front()=c.wall_outer_radius(); e.back()=c.far_field_radius;
    for(std::size_t j=1;j+1<e.size();++j)
        e[j]=std::sqrt(c.soil_shell_radii[j-1]*c.soil_shell_radii[j]);
    return e;
}
double radial_g(const Configuration& c,std::size_t link)
{
    const auto& r=c.soil_shell_radii;
    if(link==0) return 2.0*pi*c.soil_conductivity/std::log(r.front()/c.wall_outer_radius());
    if(link<r.size()) return 2.0*pi*c.soil_conductivity/std::log(r[link]/r[link-1]);
    return 2.0*pi*c.soil_conductivity/std::log(c.far_field_radius/r.back());
}
double total_radial_resistance(const Configuration& c)
{
    double r=0.0;
    for(std::size_t j=0;j<=c.soil_shell_radii.size();++j) r+=1.0/radial_g(c,j);
    return r;
}
double effective_u(const Configuration& c)
{
    const double hi=c.inner_conductance_per_length();
    return 1.0/(1.0/hi+total_radial_resistance(c));
}
double friction(const Configuration& c,double m)
{
    const double re=std::abs(m)*c.diameter/(c.area()*c.viscosity);
    if(re<=0.0) return 0.0;
    if(re<2300.0) return 64.0/re;
    const double v=c.roughness/(3.7*c.diameter)+5.74/std::pow(re,0.9);
    return 0.25/std::pow(std::log10(v),2.0);
}
double mom(const Configuration& c,const State& s,std::size_t i,const std::vector<double>& z)
{
    const double rl=rho(c,s.pressure[i],s.gas_temperature[i]);
    const double rr=rho(c,s.pressure[i+1],s.gas_temperature[i+1]);
    const double r=0.5*(rl+rr),v=s.flow[i]/(c.area()*r);
    return (s.pressure[i+1]-s.pressure[i])/c.dx()+r*c.gravity*(z[i+1]-z[i])/c.dx()
        +friction(c,s.flow[i])/(2.0*c.diameter)*r*v*std::abs(v);
}
std::vector<double> encode(const Configuration& c,const State& s)
{
    const std::size_t n=c.cells,ns=c.soil_shell_radii.size();
    if(s.pressure.size()!=n+1||s.gas_temperature.size()!=n+1||s.wall_temperature.size()!=n+1||
       s.flow.size()!=n||s.soil_temperature.size()!=ns)
        throw std::invalid_argument("soil-shell state dimensions do not match grid");
    std::vector<double> y(4*n+1+ns*(n+1));
    for(std::size_t i=0;i<n;++i) {
        y[i]=s.pressure[i+1]/c.inlet_pressure;
        y[n+i]=s.gas_temperature[i+1]/c.inlet_temperature;
        y[2*n+i]=s.flow[i]/c.outlet_mass_flow;
    }
    for(std::size_t i=0;i<=n;++i) y[3*n+i]=s.wall_temperature[i]/c.inlet_temperature;
    const std::size_t offset=4*n+1;
    for(std::size_t j=0;j<ns;++j) {
        if(s.soil_temperature[j].size()!=n+1)
            throw std::invalid_argument("soil shell state must have N+1 nodes");
        for(std::size_t i=0;i<=n;++i)
            y[offset+j*(n+1)+i]=s.soil_temperature[j][i]/c.inlet_temperature;
    }
    return y;
}
State decode(const Configuration& c,const std::vector<double>& y)
{
    const std::size_t n=c.cells,ns=c.soil_shell_radii.size();
    if(y.size()!=4*n+1+ns*(n+1)) throw std::invalid_argument("wrong soil-shell unknown size");
    State s; s.pressure.resize(n+1); s.gas_temperature.resize(n+1);
    s.wall_temperature.resize(n+1); s.flow.resize(n); s.soil_temperature.assign(ns,std::vector<double>(n+1));
    s.pressure[0]=c.inlet_pressure; s.gas_temperature[0]=c.inlet_temperature;
    for(std::size_t i=0;i<n;++i) {
        s.pressure[i+1]=y[i]*c.inlet_pressure;
        s.gas_temperature[i+1]=y[n+i]*c.inlet_temperature;
        s.flow[i]=y[2*n+i]*c.outlet_mass_flow;
    }
    for(std::size_t i=0;i<=n;++i) s.wall_temperature[i]=y[3*n+i]*c.inlet_temperature;
    const std::size_t offset=4*n+1;
    for(std::size_t j=0;j<ns;++j) for(std::size_t i=0;i<=n;++i)
        s.soil_temperature[j][i]=y[offset+j*(n+1)+i]*c.inlet_temperature;
    return s;
}
double wall_rhs(const Configuration& c,const State& s,std::size_t i)
{
    return c.inner_conductance_per_length()*(s.gas_temperature[i]-s.wall_temperature[i])+
        radial_g(c,0)*(s.soil_temperature[0][i]-s.wall_temperature[i]);
}
double shell_rhs(const Configuration& c,const State& s,std::size_t j,std::size_t i)
{
    const std::size_t ns=c.soil_shell_radii.size();
    const double left=j==0 ? s.wall_temperature[i] : s.soil_temperature[j-1][i];
    const double right=j+1==ns ? c.far_field_temperature : s.soil_temperature[j+1][i];
    return radial_g(c,j)*(left-s.soil_temperature[j][i])+
        radial_g(c,j+1)*(right-s.soil_temperature[j][i]);
}
double norm(const std::vector<double>& r)
{ double v=0.0; for(double x:r) v=std::max(v,std::abs(x)); return v; }
std::vector<double> steady_r(const Configuration& c,const std::vector<double>& y)
{
    const std::size_t n=c.cells,ns=c.soil_shell_radii.size(); const State s=decode(c,y);
    const auto z=c.elevation_nodes();
    const double ms=std::max(1.0,std::abs(c.outlet_mass_flow));
    const double ps=c.inlet_pressure/c.dx(),es=ms*cp(c)*c.inlet_temperature;
    const double hs=(c.inner_conductance_per_length()+radial_g(c,0))*c.inlet_temperature;
    std::vector<double> r(y.size());
    for(std::size_t i=1;i<n;++i) r[i-1]=(s.flow[i]-s.flow[i-1])/ms;
    r[n-1]=(c.outlet_mass_flow-s.flow[n-1])/ms;
    for(std::size_t i=0;i<n;++i) r[n+i]=mom(c,s,i,z)/ps;
    for(std::size_t i=1;i<=n;++i) {
        const double ml=s.flow[i-1],mr=i<n?s.flow[i]:c.outlet_mass_flow;
        const double flux=cp(c)*(mr*s.gas_temperature[i]-ml*s.gas_temperature[i-1]);
        const double q=c.inner_conductance_per_length()*(s.wall_temperature[i]-s.gas_temperature[i]);
        r[2*n+i-1]=(flux-q*ell(c,i))/es;
    }
    for(std::size_t i=0;i<=n;++i)
        r[3*n+i]=wall_rhs(c,s,i)/hs;
    const std::size_t offset=4*n+1;
    for(std::size_t j=0;j<ns;++j) for(std::size_t i=0;i<=n;++i)
        r[offset+j*(n+1)+i]=shell_rhs(c,s,j,i)/((radial_g(c,j)+radial_g(c,j+1))*c.inlet_temperature);
    return r;
}
std::vector<double> linear_solve(std::vector<std::vector<double>> a,std::vector<double> b)
{
    const std::size_t n=b.size();
    for(std::size_t k=0;k<n;++k) {
        std::size_t p=k; for(std::size_t i=k+1;i<n;++i) if(std::abs(a[i][k])>std::abs(a[p][k])) p=i;
        if(std::abs(a[p][k])<1e-14) throw std::runtime_error("singular soil-shell Newton Jacobian");
        std::swap(a[k],a[p]); std::swap(b[k],b[p]); const double d=a[k][k];
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
template<class F,class V> std::pair<std::size_t,double> newton(std::vector<double>& y,F residual,V valid,
    double tol,std::size_t maxit)
{
    for(std::size_t it=0;it<=maxit;++it) {
        const auto r=residual(y); const double old=norm(r); if(old<=tol) return {it,old};
        if(it==maxit) break;
        const std::size_t n=y.size();
        std::vector<std::vector<double>> j(n,std::vector<double>(n));
        for(std::size_t k=0;k<n;++k) { auto yp=y; const double h=2e-6*std::max(1.0,std::abs(y[k])); yp[k]+=h;
            const auto rp=residual(yp); for(std::size_t i=0;i<n;++i) j[i][k]=(rp[i]-r[i])/h; }
        std::vector<double> b(r.size()); std::transform(r.begin(),r.end(),b.begin(),[](double x){return -x;});
        const auto d=linear_solve(std::move(j),std::move(b)); bool ok=false;
        for(double a=1.0;a>=1.0/1024.0;a*=0.5) { auto trial=y;
            for(std::size_t k=0;k<n;++k) trial[k]+=a*d[k];
            if(!valid(trial)) continue;
            try { if(norm(residual(trial))<old) { y=std::move(trial); ok=true; break; } } catch(const std::exception&) {} }
        if(!ok) throw std::runtime_error("soil-shell Newton line search failed");
    }
    throw std::runtime_error("soil-shell Newton did not converge; residual="+
        std::to_string(norm(residual(y))));
}
bool valid(const Configuration& c,const std::vector<double>& y)
{
    const std::size_t n=c.cells,ns=c.soil_shell_radii.size();
    for(std::size_t i=0;i<n;++i) if(y[i]<=0.05||y[n+i]<=0.1) return false;
    for(double v:y) if(!std::isfinite(v)||v<=0.05) return false;
    (void)ns; return true;
}
double pressure_step(const Configuration& c,const State& s,std::size_t i,const std::vector<double>& z)
{
    const double pl=s.pressure[i],t=0.5*(s.gas_temperature[i]+s.gas_temperature[i+1]);
    const double zs=(z[i+1]-z[i])/c.dx(),a=2.0*c.gravity*zs/(c.gas_constant*t);
    const double b=friction(c,s.flow[i])*c.gas_constant*t*s.flow[i]*std::abs(s.flow[i])/
        (c.diameter*c.area()*c.area());
    const double p2=std::abs(a)<1e-10?pl*pl-b*c.dx():std::exp(-a*c.dx())*pl*pl-b*(-std::expm1(-a*c.dx()))/a;
    if(!(p2>0.0)) throw std::runtime_error("soil-shell pressure profile is inadmissible");
    return std::sqrt(p2);
}
double discrete_pressure_step(const Configuration& c,State& s,std::size_t i,const std::vector<double>& z)
{
    double p=pressure_step(c,s,i,z);
    for(std::size_t it=0;it<20;++it) { s.pressure[i+1]=p; const double r=mom(c,s,i,z);
        if(std::abs(r)*c.dx()/c.inlet_pressure<1e-13) return p;
        const double h=std::max(1.0,1e-6*p); s.pressure[i+1]=p+h;
        const double derivative=(mom(c,s,i,z)-r)/h; s.pressure[i+1]=p;
        if(std::abs(derivative)<1e-20) throw std::runtime_error("singular local momentum solve");
        const double correction=-r/derivative; bool accepted=false;
        for(double a=1.0;a>=1.0/1024.0;a*=0.5) { const double trial=p+a*correction;
            if(trial<=0.05*c.inlet_pressure) continue;
            s.pressure[i+1]=trial;
            if(std::abs(mom(c,s,i,z))<std::abs(r)) { p=trial; accepted=true; break; } }
        s.pressure[i+1]=p; if(!accepted) throw std::runtime_error("local pressure marching failed"); }
    throw std::runtime_error("local pressure marching did not converge");
}
void set_radial_profile(const Configuration& c,State& s,std::size_t i)
{
    const double q=effective_u(c)*(s.gas_temperature[i]-c.far_field_temperature);
    const double hi=c.inner_conductance_per_length();
    s.wall_temperature[i]=s.gas_temperature[i]-q/hi;
    double resistance=1.0/radial_g(c,0);
    for(std::size_t j=0;j<c.soil_shell_radii.size();++j) {
        s.soil_temperature[j][i]=s.wall_temperature[i]-q*resistance;
        resistance+=1.0/radial_g(c,j+1);
    }
}
std::vector<double> theta_r(const Configuration& c,const State& old,const std::vector<double>& y,
                            double dt,double theta)
{
    const std::size_t n=c.cells,ns=c.soil_shell_radii.size(); const State s=decode(c,y);
    const auto z=c.elevation_nodes(); const double ms=std::max(1.0,std::abs(c.outlet_mass_flow));
    const double es=ms*cp(c)*c.inlet_temperature;
    const double hs=(c.inner_conductance_per_length()+radial_g(c,0))*c.inlet_temperature;
    std::vector<double> r(y.size()); auto mt=[&](std::size_t i){return theta*s.flow[i]+(1.0-theta)*old.flow[i];};
    for(std::size_t i=1;i<=n;++i) { const double rn=rho(c,s.pressure[i],s.gas_temperature[i]);
        const double ro=rho(c,old.pressure[i],old.gas_temperature[i]); const double ml=mt(i-1),mr=i<n?mt(i):c.outlet_mass_flow;
        r[i-1]=c.area()*ell(c,i)*(rn-ro)/(dt*ms)+(mr-ml)/ms; }
    for(std::size_t i=0;i<n;++i) r[n+i]=(s.flow[i]-old.flow[i])/ms+c.area()*dt*
        (theta*mom(c,s,i,z)+(1.0-theta)*mom(c,old,i,z))/ms;
    for(std::size_t i=1;i<=n;++i) { const double l=ell(c,i);
        const double rn=rho(c,s.pressure[i],s.gas_temperature[i]),ro=rho(c,old.pressure[i],old.gas_temperature[i]);
        const double storage=c.area()*l*cv(c)*(rn*s.gas_temperature[i]-ro*old.gas_temperature[i])/dt;
        const std::size_t f=i-1; const double fl=theta*s.flow[f]*cp(c)*s.gas_temperature[f]+(1.0-theta)*old.flow[f]*cp(c)*old.gas_temperature[f];
        const double fr=i<n?theta*s.flow[i]*cp(c)*s.gas_temperature[i]+(1.0-theta)*old.flow[i]*cp(c)*old.gas_temperature[i]:
            c.outlet_mass_flow*cp(c)*(theta*s.gas_temperature[i]+(1.0-theta)*old.gas_temperature[i]);
        const double qn=c.inner_conductance_per_length()*(s.wall_temperature[i]-s.gas_temperature[i]);
        const double qo=c.inner_conductance_per_length()*(old.wall_temperature[i]-old.gas_temperature[i]);
        r[2*n+i-1]=(storage+fr-fl-(theta*qn+(1.0-theta)*qo)*l)/es; }
    for(std::size_t i=0;i<=n;++i) { const double wn=wall_rhs(c,s,i),wo=wall_rhs(c,old,i);
        r[3*n+i]=(c.wall_heat_capacity_per_length()*(s.wall_temperature[i]-old.wall_temperature[i])/dt-
            theta*wn-(1.0-theta)*wo)/hs; }
    const std::size_t offset=4*n+1; const auto edges=radial_edges(c);
    for(std::size_t j=0;j<ns;++j) for(std::size_t i=0;i<=n;++i) {
        const double cap=c.soil_density*c.soil_heat_capacity*pi*(edges[j+1]*edges[j+1]-edges[j]*edges[j]);
        const double fn= shell_rhs(c,s,j,i),fo=shell_rhs(c,old,j,i);
        r[offset+j*(n+1)+i]=(cap*(s.soil_temperature[j][i]-old.soil_temperature[j][i])/dt-
            theta*fn-(1.0-theta)*fo)/((radial_g(c,j)+radial_g(c,j+1))*c.inlet_temperature);
    }
    return r;
}
} // namespace

double Configuration::area() const { return pi*diameter*diameter/4.0; }
double Configuration::dx() const { return length/static_cast<double>(cells); }
double Configuration::wall_outer_radius() const { return diameter/2.0+wall_thickness; }
double Configuration::inner_conductance_per_length() const { return inner_heat_transfer_coefficient*pi*diameter; }
double Configuration::wall_heat_capacity_per_length() const
{ const double d=diameter+2.0*wall_thickness; return wall_density*wall_heat_capacity*pi*(d*d-diameter*diameter)/4.0; }
std::vector<double> Configuration::elevation_nodes() const
{ std::vector<double> z(cells+1); for(std::size_t i=0;i<=cells;++i) z[i]=0.5*hill_height*(1.0-std::cos(2.0*pi*static_cast<double>(i)/cells)); return z; }
void Configuration::validate() const
{
    if(cells<2||length<=0||diameter<=0||roughness<0||viscosity<=0||gas_constant<=0||gamma<=1||
       inlet_pressure<=0||inlet_temperature<=0||far_field_temperature<=0||outlet_mass_flow<=0||
       gravity<0||hill_height<0||inner_heat_transfer_coefficient<=0||wall_density<=0||
       wall_heat_capacity<=0||wall_thickness<=0||soil_conductivity<=0||soil_density<=0||soil_heat_capacity<=0||
       soil_shell_radii.empty()||far_field_radius<=0)
        throw std::invalid_argument("invalid soil-shell configuration");
    double prev=wall_outer_radius();
    for(double radius:soil_shell_radii) { if(!(radius>prev)) throw std::invalid_argument("soil shell radii must increase beyond wall radius"); prev=radius; }
    if(!(far_field_radius>prev)) throw std::invalid_argument("far-field radius must exceed outermost shell radius");
}
double swamee_jain_friction_factor(const Configuration& c,double mass_flow) { c.validate(); return friction(c,mass_flow); }
State analytical_approximation(const Configuration& c)
{
    c.validate(); const std::size_t n=c.cells; State s; s.pressure.resize(n+1); s.gas_temperature.resize(n+1);
    s.wall_temperature.resize(n+1); s.flow.assign(n,c.outlet_mass_flow); s.soil_temperature.assign(c.soil_shell_radii.size(),std::vector<double>(n+1));
    s.pressure[0]=c.inlet_pressure; s.gas_temperature[0]=c.inlet_temperature; const double u=effective_u(c);
    for(std::size_t i=1;i<=n;++i) { const double x=i*c.dx(); s.gas_temperature[i]=c.far_field_temperature+
        (c.inlet_temperature-c.far_field_temperature)*std::exp(-u*x/(c.outlet_mass_flow*cp(c))); }
    for(std::size_t i=0;i<=n;++i) set_radial_profile(c,s,i);
    const auto z=c.elevation_nodes(); for(std::size_t i=0;i<n;++i) s.pressure[i+1]=pressure_step(c,s,i,z); return s;
}
State steady_spatial_marching(const Configuration& c)
{
    c.validate(); const std::size_t n=c.cells; State s; s.pressure.resize(n+1); s.gas_temperature.resize(n+1);
    s.wall_temperature.resize(n+1); s.flow.assign(n,c.outlet_mass_flow); s.soil_temperature.assign(c.soil_shell_radii.size(),std::vector<double>(n+1));
    s.pressure[0]=c.inlet_pressure; s.gas_temperature[0]=c.inlet_temperature; const double u=effective_u(c),mc=c.outlet_mass_flow*cp(c);
    set_radial_profile(c,s,0);
    for(std::size_t i=1;i<=n;++i) { const double l=ell(c,i); s.gas_temperature[i]=(mc*s.gas_temperature[i-1]+u*l*c.far_field_temperature)/(mc+u*l); set_radial_profile(c,s,i); }
    const auto z=c.elevation_nodes(); for(std::size_t i=0;i<n;++i) s.pressure[i+1]=discrete_pressure_step(c,s,i,z); return s;
}
SolveInfo coupled_steady_newton(const Configuration& c,const State& initial,double tol,std::size_t maxit)
{ c.validate(); auto y=encode(c,initial); const auto a=newton(y,[&](const auto& x){return steady_r(c,x);},[&](const auto& x){return valid(c,x);},tol,maxit); return {decode(c,y),a.first,a.second}; }
double scaled_steady_residual(const Configuration& c,const State& s) { c.validate(); return norm(steady_r(c,encode(c,s))); }
RelaxationInfo implicit_theta_relaxation(const Configuration& c,const RelaxationOptions& o)
{
    c.validate(); if(o.theta<0.5||o.theta>1.0||o.pseudo_courant<=0||o.tolerance<=0||o.max_steps==0) throw std::invalid_argument("invalid soil-shell theta options");
    const double alpha=c.soil_conductivity/(c.soil_density*c.soil_heat_capacity);
    const double dt=o.pseudo_courant*std::max(c.dx()/std::sqrt(c.gamma*c.gas_constant*c.inlet_temperature),c.far_field_radius*c.far_field_radius/alpha);
    State s=analytical_approximation(c);
    for(auto& shell:s.soil_temperature) for(double& t:shell)
        t=c.far_field_temperature+0.95*(t-c.far_field_temperature);
    std::size_t total=0,good=0; double time=0.0,res=scaled_steady_residual(c,s);
    for(std::size_t step=1;step<=o.max_steps;++step) { const State old=s; auto y=encode(c,old);
        const auto a=newton(y,[&](const auto& x){return theta_r(c,old,x,dt,o.theta);},[&](const auto& x){return valid(c,x);},1e-9,40);
        total+=a.first; s=decode(c,y); time+=dt; res=scaled_steady_residual(c,s); good=res<=o.tolerance?good+1:0;
        if(good>=o.consecutive_converged_steps) return {s,step,total,time,dt,res};
        if(time>=o.max_pseudo_time) break;
    }
    throw std::runtime_error("soil-shell theta relaxation did not reach steady state");
}
} // namespace pipe_init::soil_shell
