# Correctness checks and benchmarks

Run from the repository root.

## Correctness

| Script | Checks |
|---|---|
| `test_families.R` | Densities R against C++, sampler paths R against C++, and densities against `VineCopula`, for all six copula families. This is the one to run after touching any density code. |
| `test_kernel.R` | The C++ density kernel against the R implementation |
| `test_fast.R` | The full C++ sampler against the reference R sampler, same seed |

`test_families.R` is the more demanding check. Both samplers consume the random
number stream in the same order, so with `proposal = "component"` their entire
parameter paths should agree, not merely their posterior summaries. Observed
agreement is about 1e-7 over 2000 iterations. The residual is the C++ version
being more accurate in the tails, where the R version's `qnorm(pnorm(z))` round
trip for the normal margin loses precision.

## Performance

| Script | Measures |
|---|---|
| `profile_sampler.R` | `Rprof` output for the R sampler |
| `bench_ess.R` | Effective sample size per second, componentwise against joint blocks |
| `bench_crossblock.R` | Adds the joint `(alpha, xi)` block; this is the benchmark quoted in the top-level README |

A caution about `profile_sampler.R`: its output is misleading on this code.
`Rprof` attributes about 19% of self time to `pmax`, inside the clamping of
pseudo-observations, which suggests the clamping is worth optimising. Timing the
primitives directly shows otherwise. At n = 3948, 300 evaluations cost

| `pgamma` | `dgamma` | `qnorm` | `pnorm` | clamping |
|---|---|---|---|---|
| 0.250s | 0.098s | 0.025s | 0.031s | 0.009s |

so the gamma margin is roughly 70% of the cost and the clamping under 4%. The
profiler's sampling overhead inflates cheap vectorised operations called many
times. Trust the direct timings.

The benchmarks report effective sample size per second rather than seconds per
iteration, since a proposal that is cheaper per iteration is not an improvement
if the chain mixes correspondingly worse. Effective sample sizes in this model
are small enough that the estimates are themselves noisy; differences under about
1.5x should not be read as real.
