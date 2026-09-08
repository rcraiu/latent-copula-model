# Model code

Shared source for the latent variable copula model. Nothing here is run
directly; the scripts in the other directories source these files.

| File | Purpose |
|---|---|
| `copula_densities.R` | Vectorised log copula densities and the activation functions `theta(w, xi)` |
| `lvc_gamma_gauss.R` | Reference Gibbs sampler for the gamma-Gaussian model, in pure R |
| `lvc_sampler.cpp` | The same sampler in C++, with the caching and blocking described in the top-level README |
| `lvc_fast.R` | R driver for `lvc_sampler.cpp`; initialises exactly as the reference does, then hands off |
| `lvc_kernel.cpp` | Standalone C++ density kernel, used where only the density is needed |

## Copula densities

The copula parameter varies across observations, since it depends on the latent
variable, which is why these are written out rather than taken from a package.
`benchmarks/test_kernel.R` validates them against `VineCopula` at fixed parameter
values, and validates the C++ implementations against the R ones.

The activation functions map the linear predictor `gscale * (xi_0 + xi_1 w)` onto
each family's parameter space:

| Family | `theta(w, xi)` | value at `xi = 0` |
|---|---|---|
| Gaussian | `tanh(lin)` | 0, independence |
| Frank | `lin` | 0, independence |
| Clayton | `exp(lin)` | 1, Kendall's tau = 1/3 |
| Gumbel | `1 + exp(lin)` | 2, Kendall's tau = 1/2 |

Note the asymmetry in the last column: for Clayton and Gumbel the prior on `xi`
is *not* centred at independence, because those families do not contain it in
their interior. This is a modelling choice rather than a defect, but it should be
stated when priors on `xi` are compared across families.

## The two samplers

`run_lvc()` in `lvc_gamma_gauss.R` is the readable reference. `run_lvc_fast()` in
`lvc_fast.R` is what the analysis uses. They consume the random number stream in
the same order, so with the same seed and `proposal = "component"` they agree to
within floating point: `benchmarks/test_fast.R` reports agreement of about 1e-7
on the full parameter path over 2000 iterations for all five families.

The residual discrepancy is not error in the C++ code. For the Gaussian copula
the R version computes `qnorm(pnorm(z))` for the normal margin, which loses
precision in the tails, while the C++ version uses the analytic identity. The C++
result is the more accurate of the two.

Both take a `family` argument in `Gaussian`, `Clayton`, `Frank`, `Gumbel`,
`Independence`. `run_lvc_fast()` additionally takes `blockwise` and
`cross_block`, described in the top-level README.
