# Mixed response: Poisson count and normal margin

Code for the mixed-response case of the paper, in which one margin is discrete.
This is the case that requires data augmentation: because the count margin's
distribution function is a step function, the copula is not unique, and the
sampler augments the discrete observation with a continuous latent variable
through a many-to-one map.

Run everything from the repository root.

## The simulation studies of Section 5

| File | Purpose |
|---|---|
| `mixed_model.R` | Model functions for the Gaussian, Clayton, Frank, Gumbel and independence copulas: conditional distribution functions `h(u1 \| u2)` and their inverses, the observed-data likelihood, the joint density with `W` integrated out by quadrature, the WAIC/CWAIC accumulators, and a simulator |
| `test_mixed.R` | Checks `mixed_model.R` against `VineCopula` (h-functions and inverses), a finer quadrature grid, and the closed-form margin of `Y` |
| `mixed_sampler.R` | The partially collapsed Gibbs sampler, `run_mixed()`, with Andrieu-Thoms adaptation and the scale group move described below |
| `make_data.R` | Writes the two data sets: the paper's Gaussian-copula data set (read from the original `Gauss_4chain.Rdata`, path in `LVC_GAUSS_RDATA`) and a Gumbel-copula data set with the same covariates and parameters |
| `run_simulation.R` | Fits five copula models to each data set, four chains each, in parallel |
| `criteria.R` | Posterior means, Gelman-Rubin diagnostics, and the model comparison criteria in both forms; writes the LaTeX tables of Section 5 |

```sh
Rscript mixed-response/test_mixed.R
LVC_GAUSS_RDATA=/path/to/Gauss_4chain.Rdata Rscript mixed-response/make_data.R
M=300000 CORES=8 Rscript mixed-response/run_simulation.R   # ~2 h on 8 cores
CORES=8 Rscript mixed-response/criteria.R
```

Fits go to `output/mixed/fits/`, per-chain logs to `output/mixed/logs/`, tables
to `output/mixed/tables/` and a summary to `output/mixed/summary.rds`.
`run_simulation.R` skips chains whose fit file already exists, so an interrupted
batch can be resumed; `DATASETS`, `FAMILIES` and `CHAINS` restrict the batch.

### The sampler

One sweep of `run_mixed()` updates, in this order, `(mu_d, lam_d)` and
`(mu_c, lam_c, log sigma^2)` by random walk Metropolis on the observed-data
likelihood; the augmented variable `U_i` behind each count, drawn from the
conditional copula restricted to `(F_d(q_i - 1), F_d(q_i)]` by inverting the
h-function; `xi` by random walk Metropolis on the augmented likelihood; `beta`
conjugately given `W`; each `W_i` by random walk Metropolis; and finally a scale
group move

```
(W, beta, lam_d, lam_c, xi_1)  ->  (c W, c beta, lam_d / c, lam_c / c, xi_1 / c)
```

with `log c` proposed symmetrically. The observed-data likelihood is invariant
under this map, so the acceptance ratio involves only the `N(beta' x_i, 1)` prior
of the latent variables, the priors of the rescaled parameters, and the Jacobian
`c^(n + p - 3)`. The move is the mixed-model form of the parameter-expanded move
used in the continuous-continuous sampler (`R/lvc_sampler.cpp`), where the scale
is carried by `sigma_W^2` instead of `beta`. It matters for chains started away
from the scale of `W` that the likelihood wants: without it, such a chain sits
for the whole run on the ridge that trades the scale of `W` against the loadings,
with the marginal blocks never accepting a proposal.

The sampler evaluates the h-function of each family in closed form
(Gaussian, Clayton, Frank) and inverts it in closed form where possible; the
Gumbel inverse is by bisection. The Frank h-function is written in a form that
is stable for large `|theta|`.

### The criteria

The criteria as originally defined in the paper score `p(q_i, y_i | eta, W_i)`
at the sampled latent variable and are accumulated online during the run
(`fit$criteria_condW`). The criteria the paper now reports integrate `W_i` out,

```
p(q_i, y_i | eta) = int p(q_i, y_i | eta, w) phi(w - beta' x_i) dw,
```

with the conditional versions `p(q_i, y_i | eta) / p(y_i | eta)` and
`p(q_i, y_i | eta) / p(q_i | eta)`; `p(y_i | eta)` is normal in closed form and
`p(q_i | eta)` is a quadrature. The integral is computed by Simpson's rule on
`w = beta' x_i + u`, `u` in `[-8, 8]`, step 0.02, which is accurate to better
than `1e-5` in the log density; Gauss-Hermite quadrature against the normal
weight is much less accurate here because the integrand is sharply peaked for
large counts. Because this needs one pass over the grid per draw and observation,
`criteria.R` uses 500 equally spaced draws per chain (`DRAWS`).

On the simulated data the two forms disagree. The conditional WAIC favours the
independence copula for both data sets, with an effective number of parameters
of order `n/3` to `n`. The marginal criteria select the Gumbel copula on the
Gumbel data on every criterion, and on the Gaussian data place Frank and
Gaussian within 2.6 units of each other, 11 to 21 units ahead of independence;
the Poisson margin with mean `exp(-1 + W)` produces many small counts, which say
little about the shape of the copula. The mechanism behind the conditional
failure is the same as in the continuous-continuous case, in a milder form:
`W_i` is a per-observation quantity through which every model can absorb part
of the association (the independence fit to the Gumbel data has
`sigma^2 = 0.085`), and the variance penalty is dominated by the posterior
variability of `W_i`, which is common to all models.

The conditional block of the Gaussian-data table reproduces the values obtained
with the original programs (Table 4 of the earlier draft) to within 0.5.

## The original programs

| File | Purpose |
|---|---|
| `LVC_functions.R` | Densities and helper functions for the Gaussian-copula mixed model |
| `MCMC_Gauss.R` | The original Gibbs sampler with the augmentation step |
| `Generate_data.R` | Simulation of mixed-response data |
| `Copula_Criteria_Test.R` | Model comparison across copula families |

These scripts predate the rest of the repository and follow their own path
conventions; run them from within this directory. Three discrepancies relative
to `R/copula_densities.R` are worth knowing about when reading them:

- `log_dGaussian` omits the `+(z1^2 + z2^2)/2` term, returning `log phi_2`
  rather than the copula log density. The omitted term cancels in a Metropolis
  ratio taken at fixed margins, but not in the latent variable update and not in
  the WAIC or CVML computations.
- `log_dIndep` returns `log(u1) + log(u2)`; the independence copula has density
  one, so its log density is zero.
- The clamping of pseudo-observations uses `.Machine$double.xmin`, and
  `1 - 2.2e-308` is exactly 1 in double precision, so the upper clamp has no
  effect.

`mixed_model.R` and `mixed_sampler.R` supersede these for the results in the
paper.
