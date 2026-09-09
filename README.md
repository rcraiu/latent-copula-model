# Latent variable copula models for bivariate surrogate responses

Code and reproduction materials for the paper *Latent Variable Copula Models
for Bivariate Surrogate Responses* (Craiu et al.).

The paper studies a hierarchical model in which a scientifically interesting
quantity `W` cannot be observed, and inference proceeds instead from a bivariate
surrogate pair `(Q, Y)` whose conditional distribution depends on `W` through
both of its margins and through the copula linking them:

```
Q_i | W_i        ~  F_q( . | phi_q(mu_q, lam_q, W_i), psi_q )
Y_i | W_i        ~  F_y( . | phi_y(mu_y, lam_y, W_i), psi_y )
(Q_i, Y_i) | W_i ~  C( F_q, F_y | theta(W_i, xi) )
W_i              ~  N( beta' x_i, 1 )
```

This repository covers two specialisations. The **mixed-response** case, with a
Poisson count margin and a normal margin, requires data augmentation and is
applied to cardiotocography data. The **continuous-continuous** case, with a
gamma margin and a normal margin, requires no augmentation and is applied to
NHANES kidney-function data, where the reciprocals of serum cystatin C and
serum beta-2 microglobulin are two surrogates for the unobservable glomerular
filtration rate. Serum creatinine, the marker used clinically, is deliberately *not* one of
the surrogates: it depends on muscle mass and so on the same age and sex that
serve as covariates for the latent variable, and it is reported on a grid coarse
enough that a continuous margin cannot represent it. See
`application-nhanes/README.md`.

## Repository layout

| Directory | Contents |
|---|---|
| `R/` | Model code: copula densities, the reference R sampler, and the Rcpp sampler |
| `application-nhanes/` | The gamma-Gaussian application: data preparation, EDA, fitting, diagnostics, model comparison |
| `simulation/` | Parameter-recovery study and the model-comparison demonstration |
| `benchmarks/` | Profiling, correctness checks of the Rcpp code against the R reference, and mixing benchmarks |
| `exploratory/` | The screening that led to the choice of surrogate pair and margins |
| `mixed-response/` | Code for the Poisson-normal case of the paper |

**All scripts are written to be run from the repository root**, so that paths
such as `R/lvc_fast.R` and `output/fits/` resolve correctly:

```sh
Rscript application-nhanes/build_data.R      # correct
cd application-nhanes && Rscript build_data.R  # will not find R/
```

Everything written to `output/` is reproducible from the scripts and is not
tracked by git.

## Requirements

R (>= 4.1) with:

```r
install.packages(c("Rcpp", "coda", "ggplot2", "haven", "dplyr", "statmod"))
```

`Rcpp` requires a working C++ compiler. On macOS this means the Xcode command
line tools (`xcode-select --install`); on Linux, `r-base-dev` or equivalent.
`VineCopula` is optional and used only to validate the copula densities in
`benchmarks/test_kernel.R`.

## Reproducing the NHANES application

```sh
# 1. Fetch the public NHANES 1999-2002 files (~40 MB) and build the analysis set.
#    Step 1 may be skipped: application-nhanes/nhanes_kidney.rds is included.
Rscript application-nhanes/download_data.R
Rscript application-nhanes/build_data.R

# 2. Descriptive figures.
Rscript application-nhanes/eda.R

# 3. Fit seven copula models, two chains each. See the runtime note below.
#    Clayton and Gaussian need M=2000000; see "What remains" below.
M=600000 ./application-nhanes/run_batch.sh \
    Gaussian:1 Gaussian:2 Clayton:1 Clayton:2 Frank:1 Frank:2 \
    Gumbel:1 Gumbel:2 Joe:1 Joe:2 SurvClayton:1 SurvClayton:2 \
    Independence:1 Independence:2

# 4. Convergence diagnostics, model comparison, figures and LaTeX tables.
Rscript application-nhanes/results.R
```

A single chain can be run directly:

```sh
Rscript application-nhanes/fit_nhanes.R <family> <chain> [M] [engine] [proposal]
Rscript application-nhanes/fit_nhanes.R Gaussian 1 600000 cpp block
```

`engine` is `cpp` (default) or `r` for the reference implementation; `proposal`
is `block` (default), `blocknox`, or `component`, as described below.

## MCMC: what is difficult here, and what was done about it

This posterior mixes badly, and the reasons are structural rather than
incidental. They are worth stating plainly, because the remedies below are what
make the reported runs feasible, and because a reader reproducing the work with
a naive sampler will not get the same answers.

**The latent variables and the marginal dispersions compete.** The model carries
one latent variable per observation. Residual dependence between the two markers
can be explained either by a tightly fitting latent variable with small marginal
dispersion, or by looser margins tied together by a strong copula. The gamma
shape `alpha` and the copula parameters `xi` are therefore strongly dependent in
the posterior: across families the posterior mean of `alpha` ran from 66 under
the Gumbel copula to 269 under independence. Because `alpha` sat in the gamma
margin's block and `xi` in the copula block, no proposal ever moved along that
ridge, and both parameters mixed at an effective sample size of about 10 per
10,000 draws.

**Four changes were made, in order of how much they helped.**

1. *A joint `(alpha, xi)` update* was appended to the block scheme, with its own
   adapted covariance, so that the chain can move along the ridge directly. This
   roughly doubled the effective sample size for `alpha`: from 10 to 19 under the
   Gaussian copula, and from 18 to 49 under Gumbel. It costs one extra gamma
   evaluation per iteration. Enabled by `proposal=block`.

2. *Joint block proposals* replaced componentwise ones. Updating the three
   parameters of a block one at a time requires three passes over the gamma
   margin per iteration; proposing them jointly from the adapted covariance
   requires one. The proposal uses the Cholesky factor of the adapted covariance,
   not merely its diagonal, and targets the acceptance rate appropriate to the
   block dimension (0.316 for a three-parameter block) rather than the scalar
   0.44. Disable with `proposal=blocknox` to keep joint blocks without the
   `(alpha, xi)` move, or `proposal=component` for the original scheme.

3. *The sampler was moved to C++* (`R/lvc_sampler.cpp`). Profiling was initially
   misleading: `Rprof` attributed 19% of self time to `pmax` inside the clamping
   of pseudo-observations, but timing the primitives directly showed the true
   cost was the gamma margin. At n = 3948, 300 evaluations cost 0.250s for
   `pgamma` and 0.098s for `dgamma`, against 0.025s for `qnorm` and 0.009s for
   the clamping. Since `pgamma` is already compiled code, a fused kernel alone
   bought only 1.3x; the gain came from caching, described next.

4. *Caching inside the loop.* The gamma pieces depend only on
   `(mu_q, lam_q, alpha, W)` and are not recomputed for proposals that leave them
   unchanged. The normal scores `z_1 = Phi^{-1}(u_1)` and the logs of the
   pseudo-observations are cached rather than recomputed inside every copula
   evaluation; the R version called `qnorm` twice per copula density and
   evaluated the copula fourteen times per iteration. For the Gaussian copula the
   normal margin makes `Phi^{-1}(F_y(y|w))` the identity, so that margin's
   `pnorm`/`qnorm` round trip is removed entirely.

Measured effect on minimum effective sample size per second, relative to the
original componentwise R sampler (`benchmarks/bench_crossblock.R`, M = 20000):

| Copula | joint blocks | joint blocks + `(alpha, xi)` |
|---|---|---|
| Gaussian | 2.73x | 3.58x |
| Gumbel | 1.95x | 2.73x |
| Frank | 3.97x | 1.81x |

The Frank ordering is within the noise of an effective sample size estimated at
these magnitudes and should not be read as a real reversal.

**What remains, and how long to run.** Even with all four changes the minimum
effective sample size is only a few tens per 10,000 retained draws, an integrated
autocorrelation time in the hundreds. This is why the analysis uses 600,000
iterations per chain with the first half discarded, rather than the tens of
thousands that the acceptance rates alone might suggest are adequate. Acceptance
rates near their targets are *not* evidence of adequate mixing in this model, and
effective sample size and the Gelman-Rubin statistic across chains should be
checked before interpreting any fit.

How long is enough turns out to depend on the copula, in a way that is itself
informative. At 600,000 iterations:

| Model | Gelman-Rubin | min ESS | posterior mean of `alpha` |
|---|---|---|---|
| Gumbel | 1.007 | 728 | 70 |
| Frank | 1.008 | 453 | 230 |
| Independence | 1.015 | 331 | 362 |
| Gaussian | 1.098 | 95 | 133 |
| Clayton | 1.261 | 122 | 265 |

The two that failed are the two that fit worst, and this is not a coincidence. A
copula that cannot represent the dependence in the data leaves that dependence to
be absorbed by the gamma shape `alpha`, which is then only weakly identified: the
95% credible interval for `alpha` under Clayton runs from 149 to 771, against 55
to 91 under Gumbel. A weakly identified `alpha` is exactly what makes the ridge
flat and the chain slow. So poor mixing here is a symptom of misspecification, and
worth reporting as such rather than only as a computational nuisance. Clayton and
Gaussian were rerun at 2,000,000 iterations, after which all models satisfy
Gelman-Rubin below 1.05.

A further improvement, not implemented here, would be a non-centered
parameterisation of `W`: `beta` is drawn conjugately given `W` while each `W_i`
moves by single-site random walk, which is the classic slow configuration for a
hierarchical model with a large number of latent variables.

## A note on the model comparison criteria

The criteria defined in the paper score a model by the predictive density
evaluated at the sampled latent variables, `p(q_i, y_i | eta, W_i)`. When both
margins are continuous *and* carry a free dispersion parameter, this quantity is
unbounded above: letting `alpha` grow while holding the conditional mean fixed
makes the gamma density at its own mean diverge, and the variance penalty does
not arrest it because that penalty is dominated by the posterior variability of
`W_i`, which is common to every copula model and so contributes about equally to
each. The criteria can then prefer a model with weak dependence and tightly
concentrated margins over one that fits the dependence well.

This does not affect the mixed-response case, where the discrete margin
contributes a probability bounded by one and the count margin has no free
dispersion parameter, nor the paper's continuous-continuous simulations, in which
the margins are held fixed. For the gamma-Gaussian case, `results.R` reports both
the criteria as defined in the paper and marginalised versions in which `W` is
integrated out by quadrature; `simulation/criteria_demo.R` demonstrates the point
on data simulated from a known model.

## Runtime

Timings are for a 2023 MacBook Pro (Apple M2 Max, 8 performance cores), R 4.4,
with each chain single-threaded.

| Task | Time |
|---|---|
| One chain, 600,000 iterations, Gaussian copula | ~45 min |
| Ten chains in parallel on 8 cores | ~65-85 min |
| `results.R` (quadrature for the marginalised criteria) | ~15 min |
| `simulation/test_recovery.R` | ~20 min |
| `benchmarks/bench_crossblock.R` | ~17 min |

Set `M` to a smaller value for a quick check; `M=20000` finishes in about a
minute per chain and is enough to confirm the code runs, though not to reproduce
the reported posterior summaries.

## License

MIT, see `LICENSE`. The NHANES data are public and distributed by the US Centers
for Disease Control and Prevention; they are downloaded by
`application-nhanes/download_data.R` rather than redistributed here.
