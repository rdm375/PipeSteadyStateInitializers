# Initializer Timing Report

Source data: `initializer_timing_data.csv`

Standard case: 100 km pipe, 100 cells, 0.6 m diameter, 8 MPa inlet, 100 kg/s outlet; theta relaxation uses θ = 0.5.
Newton starts from uniform inlet pressure and prescribed flow. Relaxation runs until its scaled derivative tolerance is met.

| Method | Runs | Median | Mean | Std. dev. | Min–max | Iterations / steps | Max scaled steady residual | vs. fastest |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Analytical pressure-squared | 30 | 3.258 μs | 3.319 μs | 322.646 ns | 3.239 μs–5.025 μs | 0 Newton iterations | 1.141e-16 | 1.00× |
| Newton from flat initial state | 30 | 4.708 ms | 4.708 ms | 15.342 μs | 4.681 ms–4.734 ms | 3 Newton iterations | 3.049e-13 | 1445.12× |
| Theta pseudo-time relaxation | 3 | 3.014 s | 3.015 s | 6.663 ms | 3.008 s–3.021 s | 1801 pseudo-time steps; inner Newton count unavailable | 7.530e-08 | 925167.86× |

Timing ratios use the fastest measured median as the reference. The analytical initializer, direct nonlinear solve, and transient relaxation have different work and stopping criteria, so timing is a performance comparison rather than an accuracy ranking.
