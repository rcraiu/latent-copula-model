# Simulation studies

Run from the repository root.

## `test_recovery.R` — parameter recovery

Simulates from the canonical gamma-Gaussian model with a Gaussian copula at
known parameter values and checks that the sampler recovers them, in the sense
that the 95% credible interval covers the truth for each of the 14 parameters.

The copula parameters are set to `xi = (2, 4)`, which keeps
`tanh(0.1 * (xi_0 + xi_1 w))` away from saturation at plus or minus one over the
range of `w` that the simulated data actually visit. Values that saturate make
the copula parameter effectively constant and the recovery test vacuous.

Recovery is checked at n = 2000 with 150,000 iterations. Shorter runs do not
recover all parameters, which is a property of the posterior rather than a bug;
see the discussion of mixing in the top-level README. Set `LVC_M` to override the
iteration count.

## `contcont_selection.R` — copula selection, continuous case (Section 5)

The continuous-continuous counterpart of the mixed-response simulation in
`mixed-response/`. Two data sets of n = 500 are generated with the covariates
and parameters of the mixed study (`beta = (1, 0.5)`, `xi = (-1, 3)`,
`theta = kappa(0.1 (xi_0 + xi_1 W))`), one from a Gaussian and one from a
Gumbel copula, with margins `Q | W ~ Gamma(alpha = 5, mean exp(-1 + W))` and
`Y | W ~ N(W, 1)`. The five copula models are fitted with `run_lvc_fast()`
(unit-norm `beta`, free `sigma_W^2`, PX-DA group move) and scored both by the
conditional criteria accumulated by the sampler and by the criteria with `W`
integrated out by Simpson quadrature. Posterior means are compared with the
truth in the identified parametrisation, `beta / ||beta||`,
`sigma_W^2 = 1 / ||beta||^2`, `lambda ||beta||`, `xi_1 ||beta||`.

```sh
M=300000 CHAINS=2 DRAWS=500 CORES=4 Rscript simulation/contcont_selection.R
```

Output goes to `output/contcont/` (fits, LaTeX tables, `summary.rds`); the data
sets are written to `simulation/data/`.

## Checks of the identified parametrisation

`test_beta_norm1.R`, `test_wvar_free.R`, `test_beta_direction.R`,
`test_px_da.R` and `test_pxda.R` are the experiments behind Section 2.1 of the
paper: the scale ridge of the free-`beta` model, the reparametrisation with
`||beta|| = 1` and free `sigma_W^2`, the recovery of the direction of `beta`
when the generating vector has norm 2 or 1/2, and the parameter-expanded group
move that repairs the mixing of the loadings. Each script states its question
and its run length in its header.

## `criteria_demo.R` — why the conditional criteria fail here

Simulates from a model with a copula whose parameter genuinely depends on `W`, so
the correct answer is known in advance, then fits both the Gaussian copula model
and the independence model and scores them two ways: by the criteria as defined
in the paper, which condition on the sampled latent variables, and by
marginalised versions with `W` integrated out by quadrature.

The conditional criteria prefer independence even though the data were generated
with a strong copula. The marginalised criteria recover the truth. The mechanism
is that the conditional predictive density is unbounded above when both margins
carry a free dispersion parameter, so the independence model can win by shrinking
`alpha` and `sigma` and letting the latent variables track the data; the variance
penalty does not prevent this because it is dominated by the posterior
variability of the latent variables, which every model shares.

Set `DEMO_M` to override the iteration count.
