# Mixed response: Poisson count and normal margin

Code for the mixed-response case of the paper, in which one margin is discrete.
This is the case that requires data augmentation: because the count margin's
distribution function is a step function, the copula is not unique, and the
sampler augments the discrete observation with a continuous latent variable
through a many-to-one map.

| File | Purpose |
|---|---|
| `LVC_functions.R` | Densities and helper functions for the mixed case |
| `MCMC_Gauss.R` | Gibbs sampler with the augmentation step |
| `Generate_data.R` | Simulation of mixed-response data |
| `Copula_Criteria_Test.R` | Model comparison across copula families |

These scripts predate the continuous-continuous work and follow their own path
conventions rather than the root-relative convention used elsewhere in this
repository; run them from within this directory.

Two notes for anyone reworking this code against `R/copula_densities.R`, which
fixed three discrepancies found while writing the continuous-continuous case:

- `log_dGaussian` here omits the `+(z1^2 + z2^2)/2` term, returning `log phi_2`
  rather than the copula log density. The omitted term cancels in a Metropolis
  ratio taken at fixed margins, but not in the latent variable update and not in
  the WAIC or CVML computations.
- `log_dIndep` returns `log(u1) + log(u2)`; the independence copula has density
  one, so its log density is zero.
- The clamping of pseudo-observations uses `.Machine$double.xmin`, and
  `1 - 2.2e-308` is exactly 1 in double precision, so the upper clamp has no
  effect.

The conditional model comparison criteria are sound in this case, unlike the
continuous-continuous case: the discrete margin contributes a probability bounded
by one, and the Poisson margin has no free dispersion parameter, so the
conditional predictive density cannot be made arbitrarily large. See the
top-level README.
