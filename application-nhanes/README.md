# Application: kidney function

The glomerular filtration rate is the standard summary of kidney function and is
effectively unobservable in practice, since measuring it directly requires
infusing an exogenous marker and taking timed samples. Clinical practice instead
uses endogenous filtration markers, whose serum concentration rises as filtration
falls.

The marker of first resort clinically is serum creatinine, but it is unsuitable
here on two counts, and we use two others instead. Creatinine is generated in
proportion to muscle mass, so it depends on age, sex and body composition as much
as on filtration -- and those are exactly the quantities we want as *covariates*
for `W`, so using creatinine as a surrogate would confound the two roles. It is
also too coarsely rounded in these cycles to be treated as continuous at all (see
below). Serum cystatin C and serum beta-2 microglobulin avoid both problems: both
are low-molecular-weight proteins filtered freely at the glomerulus and
catabolised by the proximal tubule, neither depends on muscle mass, and both
appear in validated creatinine-free GFR equations.

What makes them a genuine surrogate *pair* rather than two copies of one
measurement is that their extra-renal determinants differ. Cystatin C responds to
thyroid status, corticosteroid use and adiposity; beta-2 microglobulin, the light
chain of the class I MHC, responds to inflammation and immune activation. Neither
is a clean measurement of filtration on its own, which is what makes their
residual dependence given `W` scientifically meaningful.

| | |
|---|---|
| `Q` | 1 / cystatin C, L/mg, gamma margin |
| `Y` | 1 / beta-2 microglobulin, L/mg, normal margin |
| `W` | latent glomerular filtration rate |
| covariates | age, sex, BMI, systolic blood pressure, standardised |

Note the orientation. Both surrogates are on the **reciprocal (clearance)
scale**, on which a renally cleared marker is linear in the filtration rate, so
both loadings are positive and `W` increases with *filtration*. Its posterior
mean correlates +0.747 with CKD-EPI eGFR and -0.898 with age.

## Data

NHANES 1999-2002. Both markers come from `SSCARD_A`, a surplus-sera file covering
1999-2004 that carries recalibrated cystatin C and beta-2 microglobulin together;
we use the two cycles for which the covariate files are downloaded. Restricting
to adults 20+ with complete data and excluding the 153 participants in CKD stages
G3b-G5 (eGFR < 45 mL/min/1.73m^2) leaves **n = 7468**, median age 46, 52.2%
female. Kendall's tau between the markers is **0.595**.

The exclusion is made because both markers stop being clean filtration
surrogates in that range -- beta-2 microglobulin especially, since it is raised by
the inflammation accompanying advanced disease, and stage G5 patients are
typically on dialysis. Staging uses creatinine-based CKD-EPI, i.e. a criterion
*external to both surrogates*. That matters for a conditional model: creatinine
proxies `W` itself rather than the residual dependence, so restricting on it
truncates the range of `W` but leaves the conditional laws of `Q` and `Y` given
`W`, and the copula between them, unaffected. Selecting on a function of the two
surrogates would distort exactly what the model estimates.

### Choice of scale: use the reciprocal, not the log

Decided before the final fits, in `eda.R`. The serum concentration of a renally
cleared marker is inversely proportional to the filtration rate, so it is the
*reciprocal* of the concentration that is linear in filtration -- the same
reasoning behind the clinical use of reciprocal-creatinine slopes to track CKD
progression. Writing every specification as a density for the observed pair of
concentrations, so the totals are comparable across scales:

| specification | total AIC | KS gamma | KS normal |
|---|---|---|---|
| Q = cystatin C, Y = log B2M | 609.1 | 0.055 | 0.054 |
| Q = B2M, Y = log cystatin C | 947.0 | 0.074 | 0.038 |
| Q = cystatin C, Y = 1/B2M | -145.2 | 0.055 | 0.021 |
| **Q = 1/cystatin C, Y = 1/B2M** (used) | **-525.6** | **0.029** | **0.021** |
| Q = 1/B2M, Y = 1/cystatin C | -30.4 | 0.035 | 0.030 |

The Box-Cox power that symmetrises the covariate residual is -0.80 for B2M and
-0.43 for cystatin C, either side of the -1 of the reciprocal; we use the
reciprocal because it is interpretable and estimates nothing from the data.

**This was not cosmetic, and the reason is worth knowing.** The normal margin is
the rigid half of the model: `Y | W` is normal and `W` is normal, so the model
implies a *symmetric* marginal `Y` and cannot do otherwise. With `Y` on the log
scale (residual skewness 0.85) `ppcheck.R` showed the model failing to reproduce
the skewness of either margin, the upper quantiles of both, or the tail
asymmetry -- and the misfit did **not** stay in the margins. Refitting on a
symmetrised scale moved the estimated dependence trend under every family and
*reversed its sign* under Gumbel. It also fixed a mixing problem: on the log
scale the Gaussian copula would not converge in 200k iterations (PSRF 1.085,
min ESS 55) and needed 600k, while on the reciprocal scale it converges as
readily as the rest. Bad mixing was a symptom of misspecification, not a
property of the family.

### Why not creatinine: measurement granularity

This is not a minor caveat but the reason creatinine is excluded as a surrogate.
In these cycles it is reported to the nearest 0.1 mg/dL:

| marker | distinct values | largest tied group | spacing / sd |
|---|---|---|---|
| 1 / cystatin C | 641 | 0.70% | 0.008 |
| 1 / beta-2 microglobulin | 384 | 1.25% | 0.011 |
| log creatinine (**not used**) | 16 | **17.6%** | **0.35** |

A continuous margin cannot represent atoms holding a sixth of the sample, and
Sklar's theorem identifies the copula uniquely only for continuous margins, so
with ties on that scale the estimand itself is not well determined.
`figs/NHANES_granularity.png` shows creatinine as a staircase beside two smooth
curves. Creatinine is retained in the data frame, but only to stage CKD for the
exclusion above.

### Covariates

Age and sex enter as covariates for `W` rather than being built into the
response, a deliberate contrast with the CKD-EPI equations that bake them into
the reported quantity. This is cleaner here than it would be with creatinine,
whose dependence on muscle mass makes age and sex determinants of the marker as
well as of kidney function.

### Tail dependence, and a caution about reading it

Nonparametric estimates of `P(U2 <= t | U1 <= t)` and its upper counterpart are
0.65 vs 0.39 at t = 0.05, and 0.54 vs 0.28 at t = 0.02: strong dependence,
markedly asymmetric, concentrated where both markers are *small*, which on the
clearance scale is where filtration is poorest. All seven families are fitted so
the comparison is not rigged toward that asymmetry.

**The caution:** Clayton has its tail dependence in exactly that tail, and still
finishes last among the copulas (620 WAIC behind Frank). The tail behaviour of
the *observed* pair is a poor guide to the family of the *conditional* copula,
because most of the marginal asymmetry is generated by `W` itself. Use the
diagnostic to choose which families to include, not to predict the winner.

## Scripts

Run from the repository root, in this order.

| Script | Purpose |
|---|---|
| `download_data.R` | Fetch the public NHANES files (~40 MB) into `application-nhanes/data/` |
| `build_data.R` | Merge, apply exclusions, and write `nhanes_kidney.rds` |
| `eda.R` | Granularity, choice of scale and margin, histograms, scatterplots, Kendall's tau by age quintile |
| `fit_nhanes.R` | Fit one copula family and one chain |
| `run_batch.sh` | Launch several chains in parallel and wait for them |
| `results.R` | Convergence diagnostics, both sets of criteria, posterior summaries, figures, LaTeX tables |
| `marginal_criteria.R` | The marginalised criteria on their own, with a quadrature accuracy check |
| `diagnostics.R` | PSRF, MPSRF and effective sample size per family |
| `ppcheck.R` | Posterior predictive check of the margins, association and tail behaviour |

`nhanes_kidney.rds` is committed, so the analysis can be reproduced without the
download step. The download-then-build path has been checked from a clean state
and reproduces the committed file exactly.

Note on the download URLs: CDC serves these from
`.../Nchs/Data/Nhanes/Public/<year>/DataFiles/<name>.xpt`, *not* the
`.../Nchs/Nhanes/<year range>/<name>.XPT` form that older documentation cites --
that form now returns an HTML "page not found" body with a 200 status, so
`download_data.R` verifies the SAS transport signature of each file rather than
trusting the status code.

## Two points worth knowing before interpreting output

**Stratify on covariates, not on the markers.** `eda.R` examines whether the
association between the markers varies, by computing Kendall's tau within age
quintiles; it rises monotonically from 0.413 to 0.645. An earlier version
stratified on `(u1 + u2)/2`, a function of both margins; conditioning on such a
quantity induces a selection effect that attenuates the within-stratum
association and produced a spurious negative trend. Stratifying on eGFR would be
circular for the same reason, since eGFR is a deterministic function of
creatinine.

**Set the quadrature window wide enough.** The marginal criteria integrate `W`
out by Simpson quadrature on a window around the integrand's mode. The error is
almost entirely *truncation*, not discretisation, and it is strongly
family-dependent: at a half-width of 1.6 the total log density was off by 2 for
the independence copula but 261 for Frank and 524 for Gaussian, which is more
than enough to reorder the families outright. Widening to 4.0 costs nothing --
the node count is irrelevant, 81 suffices -- and drops the worst error below
0.02. `marginal_criteria.R` now checks every family, not just one; checking a
single family is what hid the problem, since independence is the most accurate
case.

**Read both sets of criteria.** `results.R` reports the criteria as defined in
the paper alongside marginalised versions with `W` integrated out. In this model
they disagree, and the disagreement is systematic rather than numerical; see the
top-level README and `simulation/criteria_demo.R`. The marginalised versions have
a well-behaved penalty term, close to the number of parameters, whereas the
conditional penalty is close to n/2 for every family and so does no
discriminating between them.

## Output

Written to `output/` at the repository root: `fits/` (one `.rds` per family and
chain), `figs/`, `tables/` (LaTeX fragments for the paper), `logs/`.
