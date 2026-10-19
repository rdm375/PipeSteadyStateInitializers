# Isothermal Steady Initializers

Small C++20 reference implementation of three steady initializers for the isothermal single-pipe model:

1. Analytical pressure-squared profile with Swamee–Jain Darcy friction.
2. Physical pseudo-time relaxation of a staggered low-Mach pipe model under fixed boundary conditions.
3. Damped Newton solve of the same semi-discrete equations with all time derivatives set to zero.

The current workspace did not contain the original `isothermal-pipe-sim` checkout, so this is a standalone implementation. It follows the established boundary convention: fixed inlet pressure and outlet mass flow. The elevation is a smooth hill with equal endpoint heights.

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
```

## Method behavior

- **Analytical:** fast continuum approximation; includes gravity and Swamee–Jain friction. It generally has a small mesh-dependent residual when checked against the discrete equations.
- **Relaxation:** initializes a uniform pressure and the prescribed outlet flow, then advances the staggered semi-discrete transient system with RK4 until the scaled right-hand side remains below tolerance.
- **Newton:** solves the finite-dimensional steady residual for cell pressures and staggered flows. A damped line search preserves positive pressures.

The Newton method is the discrete equilibrium reference. Relaxation should converge to the same equilibrium within its stopping tolerance. The analytical method provides a useful initial guess and a flat-pipe/hydrostatic reference.

## Current scope

The implementation uses a constant-$Z$, constant-viscosity isothermal closure and Swamee–Jain friction. It implements forward flow for the demonstration configuration. The friction term itself is signed for reversed flow, but reversed-flow boundary handling is not yet generalized. The pressure variable is cell-centered and mass flow is staggered at cell-left faces; the prescribed outlet flow closes the last cell continuity equation. Outlet pressure is recovered through a final half-cell pressure-squared step.

This is a reference implementation for the three initializer algorithms, not yet a patch to the original simulator. Before merging into that repository, align its exact spatial residual, boundary half-cell convention, EOS, friction derivatives, and solver interfaces.
