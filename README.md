# Isothermal Steady Initializers

Small C++20 reference implementation of four steady initializers for the isothermal single-pipe model:

1. Analytical pressure-squared profile with Swamee–Jain Darcy friction.
2. Physical pseudo-time relaxation of a staggered low-Mach pipe model under fixed boundary conditions.
3. Damped Newton solve of the same semi-discrete equations with all time derivatives set to zero.
4. Cell-by-cell discrete steady marching using the same semi-discrete momentum equation.

The current workspace did not contain the original `isothermal-pipe-sim` checkout, so this is a standalone implementation. Its spatial state and boundary storage now follow the simulator notes: pressure/density at $N+1$ nodes, $N$ internal mass flows at staggered faces, and half-cell control volumes at both pipe ends. The boundary conditions are fixed inlet pressure and prescribed outlet mass flow. The default demonstration uses the standard 100 km, 100-cell, 0.6 m pipe at 8 MPa inlet pressure and 100 kg/s outlet flow. The reference case is flat (`hill_height = 0`); nonzero hill elevations remain configurable.

## Model implemented

For constant pipe area $A=\pi D^2/4$, fixed uniform temperature $T_0$, and constant compressibility factor $Z$,

$$
C=Z R_s T_0,\qquad \rho=\frac{p}{C},\qquad u=\frac{\dot m}{A\rho}.
$$

The code uses the low-Mach isothermal pipeline equations

$$
\frac{\partial(A\rho)}{\partial t}
+\frac{\partial\dot m}{\partial z}=0,
$$

$$
\frac{\partial\dot m}{\partial t}
+A\frac{\partial p}{\partial z}
=-A\rho g\frac{dh}{dz}
-A\frac{f_D}{2D}\rho u|u|.
$$

The convective momentum-flux term is omitted, consistent with the pressure-squared analytical initializer. This is the conventional low-Mach reduction; it is an explicit model choice and should be matched to the target solver before integration.

At steady state, continuity makes the flow constant on every internal face. The inlet half-cell balance gives the inlet flow as the first internal flow, while the outlet half-cell balance enforces the prescribed outlet flow. The discrete linepack uses trapezoidal endpoint weights, $A\,dx\,[\rho_0/2+\sum_{i=1}^{N-1}\rho_i+\rho_N/2]$, matching the simulator's nodal storage convention.

The analytical method follows from the steady equations and constant mass flow. With $y=p^2$,

$$
\frac{dy}{dz}+\frac{2g}{C}\frac{dh}{dz}y
=-\frac{f_D C\dot m|\dot m|}{D A^2}.
$$

On flat pipe this reduces to the familiar linear pressure-squared relation. For the hill, the code advances the exact integrating-factor solution over each linearly interpolated elevation segment.

## Build and run

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
ctest --test-dir build --output-on-failure
./build/pipe_init_demo
./build/pipe_init_benchmark 30 3 > initializer_timing_data.csv
python3 scripts/report_initializer_timings.py initializer_timing_data.csv
```

The benchmark warms each method, then records per-run wall time in CSV. The
analytical, discrete marching, and direct Newton methods use 30 repetitions by default; pseudo-time
relaxation uses 3 because it is substantially more expensive. Newton starts
from uniform inlet pressure and the prescribed flow, without an analytical
profile seed. Arguments override the two repetition counts, respectively.
The reporting script writes a Markdown summary with medians, ranges,
variability, work counts, residuals, and timing ratios relative to the fastest
median.

## Method behavior

- **Analytical:** fast continuum approximation; includes gravity and Swamee–Jain friction. It generally has a small mesh-dependent residual when checked against the discrete equations.
- **Discrete marching:** starts at the prescribed inlet pressure and solves the discrete momentum residual for each next pressure node, with the prescribed steady flow on every staggered face. It takes O(N) work, needs no global nonlinear solve, and reaches the same discrete steady profile as Newton for this model. It reports an error if a cell has no admissible positive-pressure root.
- **Relaxation:** initializes a uniform pressure and the prescribed outlet flow, then advances the staggered system with the simulator's implicit theta method. Each pseudo-time step solves the nonlinear theta-weighted continuity and momentum residuals with a damped Newton iteration. The default is $\theta=0.5$ (trapezoidal); set `RelaxationOptions::theta` to 1.0 for backward Euler. The pseudo-time step is `pseudo_courant * dx / sqrt(Z R_s T)`; this Courant-like factor controls convergence speed, not explicit-method stability.
- **Newton:** solves the finite-dimensional steady residual for the unknown nodal pressures $p_1,\ldots,p_N$ and staggered flows. The inlet pressure is fixed; a damped line search preserves positive pressures.

Discrete marching and steady Newton are equivalent for the current boundary conditions and constant friction factor closure; Newton remains useful as a general nonlinear-system reference. Theta relaxation should converge to the same equilibrium within its stopping tolerance. The analytical method provides a fast continuum approximation and an independent flat-pipe/hydrostatic reference.

## Adiabatic ideal-gas initializer experiment

The separate `pipe_init::adiabatic` module implements three initializer strategies for a calorically perfect ideal gas. It uses constant area, constant viscosity, Swamee–Jain Darcy friction, elevation, forward flow, prescribed inlet pressure and temperature, and prescribed outlet mass flow. The gas has constant $R_s$ and $\gamma$; there is no wall heat transfer, soil state, or real-gas EOS.

The steady cell momentum residual includes pressure, convective momentum flux, gravity, and Darcy drag. The steady energy balance conserves $h+u^2/2+gz$ between nodes. Spatial marching solves the local momentum and total-enthalpy equations cell by cell. Coupled Newton solves all continuity, momentum, and energy residuals together.

The adiabatic theta method advances the semidiscrete mass, momentum, and total-energy equations in pseudo-time under the steady boundary conditions. Pseudo-time is used only to reach the equilibrium; it is not presented as a physical transient run. The default demo uses $\theta=0.9$ for added damping of the pseudo-time iterations.

Build and run the adiabatic comparison:

```sh
cmake --build build --target adiabatic_pipe_init_demo -j
./build/adiabatic_pipe_init_demo
```

Export and plot all three adiabatic initializer profiles, with pressure, mass
flow, and temperature in the first three panels and pipe elevation in the
bottom panel:

```bash
cmake --build build --target adiabatic_pipe_init_profile_export -j
./build/adiabatic_pipe_init_profile_export > adiabatic_profiles.csv
python3 scripts/plot_adiabatic_profiles.py adiabatic_profiles.csv \\
    -o adiabatic_initializer_profiles.png
```

## Current scope

The implementation uses a constant-$Z$, constant-viscosity isothermal closure and Swamee–Jain friction. It implements forward flow for the demonstration configuration. The friction term itself is signed for reversed flow, but reversed-flow boundary handling is not yet generalized. The grid positions, endpoint half-cell storage, and boundary convention are aligned to the simulator notes. The momentum closure remains the reduced low-Mach form documented above; the original simulator's momentum flux, EOS, friction routine, and software interfaces still need direct comparison in its source checkout before claiming full physics or code-level equivalence.

This is a reference implementation for the three initializer algorithms, not yet a patch to the original simulator. The pressure/flow locations and boundary storage follow the simulator notes available here. Before merging into that repository, compare the exact momentum residual, EOS, friction routine and derivatives, and solver interfaces against its source checkout.
