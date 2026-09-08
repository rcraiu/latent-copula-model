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
