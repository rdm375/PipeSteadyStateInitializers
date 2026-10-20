# Fully thermal soil-shell initializer performance report

Source data: `soil_shell_timing_data.csv`. Values summarize repeated runs; timing is measured on the machine that ran the benchmark.

Default case: 100 km pipe, 16 axial cells, 0.6 m diameter, 8 MPa inlet pressure, 288.15 K inlet gas, 278.15 K far-field soil, 100 kg/s outlet flow, four soil shells, flat elevation. Newton uses the analytical profile as its initial guess. θ relaxation uses backward Euler (θ=1) and the configured pseudo-time step.

| Method | Runs | Median | Mean | Std. dev. | Min–max | Work | Max residual | Outlet p / Tg | vs. fastest |
|---|---:|---:|---:|---:|---:|---|---:|---|---:|
| Analytical / semi-analytical approximation | 30 | 4.309 μs | 4.334 μs | 97.901 ns | 4.280 μs–4.789 μs | direct evaluation | 3.630e-04 | 5.368596 MPa / 279.8341 K | 1.00× |
| Steady spatial marching | 30 | 8.885 μs | 9.194 μs | 1.404 μs | 8.786 μs–16.599 μs | direct evaluation | 1.746e-15 | 5.366695 MPa / 280.0944 K | 2.06× |
| Coupled steady Newton | 30 | 2.729 ms | 2.753 ms | 106.682 μs | 2.646 ms–3.121 ms | 2 Newton iterations | 2.026e-15 | 5.366695 MPa / 280.0944 K | 633.43× |
| Implicit θ pseudo-time relaxation | 5 | 99.727 ms | 99.566 ms | 429.319 μs | 98.873 ms–99.992 ms | 16 steps, 54 inner Newton iterations | 1.132e-09 | 5.366695 MPa / 280.0944 K | 23143.80× |

Timing ratios compare medians against the fastest method. The analytical approximation is a continuum estimate and is not expected to satisfy the discrete steady residual as tightly as the discrete solvers. Runtime is environment-dependent; do not compare figures across different hardware or compiler settings as if they were controlled benchmarks.
