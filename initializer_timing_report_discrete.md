# Initializer Timing Report

Source data: `initializer_timing_data_discrete.csv`

Standard case: 100 km pipe, 100 cells, 0.6 m diameter, 8 MPa inlet, 100 kg/s outlet; theta relaxation uses θ = 0.5.
Newton starts from uniform inlet pressure and prescribed flow. Relaxation runs until its scaled derivative tolerance is met.

| Method | Runs | Median | Mean | Std. dev. | Min–max | Iterations / steps | Max scaled steady residual | vs. fastest |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Analytical pressure-squared | 30 | 3.264 μs | 3.277 μs | 35.260 ns | 3.251 μs–3.400 μs | 0 Newton iterations | 1.141e-16 | 1.00× |
| Discrete steady marching | 30 | 4.047 μs | 4.307 μs | 587.440 ns | 4.026 μs–6.194 μs | 0 Newton iterations | 2.252e-16 | 1.24× |
| Newton from flat initial state | 30 | 4.726 ms | 4.850 ms | 263.501 μs | 4.663 ms–5.580 ms | 3 Newton iterations | 3.049e-13 | 1448.10× |
| Theta pseudo-time relaxation | 3 | 2.991 s | 2.992 s | 11.513 ms | 2.981 s–3.004 s | 2306 Newton iterations / 1801 pseudo-time steps | 7.530e-08 | 916385.27× |

Timing ratios use the fastest measured median as the reference. The analytical initializer, direct nonlinear solve, and transient relaxation have different work and stopping criteria, so timing is a performance comparison rather than an accuracy ranking.
