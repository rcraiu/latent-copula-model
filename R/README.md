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

| Family | `theta(w, xi)` | value at `xi = 0` | tail dependence |
|---|---|---|---|
| Gaussian | `tanh(lin)` | 0, independence | none |
| Frank | `lin` | 0, independence | none |
| Clayton | `exp(lin)` | 1, Kendall's tau = 1/3 | lower |
| Gumbel | `1 + exp(lin)` | 2, Kendall's tau = 1/2 | upper |
| Joe | `1 + exp(lin)` | 2, Kendall's tau = 0.36 | upper, stronger than Gumbel |
| survival Clayton | `exp(lin)` | 1, Kendall's tau = 1/3 | upper |

Note the asymmetry in the third column: for Clayton, Gumbel, Joe and survival
Clayton the prior on `xi` is *not* centred at independence, because those
families do not contain it in their interior. This is a modelling choice rather
than a defect, but it should be stated when priors on `xi` are compared across
families.

Joe and survival Clayton were added because the NHANES fits pointed clearly at
upper tail dependence. Survival Clayton is the Clayton density evaluated at
`1 - u`, which isolates the effect of tail orientation from the choice of family:
Clayton and survival Clayton differ only in which tail they concentrate in.

Both depend on the pseudo-observations only through `1 - u`. In the C++ sampler
the cache slot that holds `log u` for the other families holds `log(1 - u)` for
these two, since a run uses a single family throughout; this is what the
`need_ubar` flag selects.

## The two samplers

`run_lvc()` in `lvc_gamma_gauss.R` is the readable reference. `run_lvc_fast()` in
`lvc_fast.R` is what the analysis uses. They consume the random number stream in
the same order, so with the same seed and `proposal = "component"` they agree to
within floating point on the whole parameter path, not merely in distribution.
`benchmarks/test_families.R` checks this for every family.

The residual discrepancy for the Gaussian copula is not error in the C++ code.
There the R version computes `qnorm(pnorm(z))` for the normal margin, which loses
precision in the tails, while the C++ version uses the analytic identity. The C++
result is the more accurate of the two, so the test allows a looser tolerance for
that family alone.

Both take a `family` argument in `Gaussian`, `Clayton`, `Frank`, `Gumbel`, `Joe`,
`SurvClayton`, `Independence`. `run_lvc_fast()` additionally takes `blockwise` and
`cross_block`, described in the top-level README. `FAMILY_CODE` in `lvc_fast.R`
must stay in step with `enum Family` in `lvc_sampler.cpp`.

## Why the density test is the sensitive one

`test_families.R` checks the densities before the sampler paths, and that order
matters. These samplers are chaotic: a one-bit difference in a density can flip a
Metropolis decision, after which the two paths are unrelated. So a path
comparison that fails while the density comparison passes points at a numerical
discrepancy, whereas a density comparison that fails points at a wrong formula.

This is how a real inconsistency was caught. The R reference parameterised the
gamma margin by *rate* while the C++ used *scale*. These are reciprocals, but
`1/(alpha/b)` and `b/alpha` differ in the last bit, and for a large shape
parameter that is enough for the tail of `pgamma` to disagree by around `1e-3` in
the log density at extreme parameter values. The R code now uses scale, and all
families agree to machine precision. The lesson generalises: when two
implementations must match bit for bit, they have to use identical
parameterisations, not merely mathematically equivalent ones.
