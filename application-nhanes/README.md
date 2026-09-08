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
| `Q` | serum cystatin C, mg/L, gamma margin |
| `Y` | log serum beta-2 microglobulin, mg/L, normal margin |
| `W` | latent renal impairment |
| covariates | age, sex, BMI, systolic blood pressure, standardised |

Note the orientation: both markers rise as filtration falls, and both loadings
are positive, so `W` increases with *impairment*. Its posterior mean correlates
negatively with CKD-EPI eGFR.

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

### Which marker gets which margin

Decided before fitting, in `eda.R`. Writing both assignments as densities for the
observed pair, so the totals are on a common scale:

| assignment | total AIC | KS, gamma margin | KS, normal margin |
|---|---|---|---|
| **Q = cystatin C, Y = log B2M** (used) | **609** | 0.055 | 0.054 |
| Q = B2M, Y = log cystatin C | 947 | 0.074 | 0.038 |

The chosen assignment wins on total AIC and balances the two margins instead of
making one much worse. The reason: the gamma margin has a shape parameter to
absorb positive skew and the normal margin has none, so the more skewed marker
(B2M) gains most from the log scale.

### Why not creatinine: measurement granularity

This is not a minor caveat but the reason creatinine is excluded as a surrogate.
In these cycles it is reported to the nearest 0.1 mg/dL:

| marker | distinct values | largest tied group | spacing / sd |
|---|---|---|---|
| cystatin C | 641 | 0.70% | 0.005 |
| log beta-2 microglobulin | 384 | 1.25% | 0.015 |
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

### Tail dependence

Nonparametric estimates of `P(U2 > t | U1 > t)` and its lower counterpart are
0.65 vs 0.39 at t = 0.05, and 0.54 vs 0.28 at t = 0.02: strong dependence,
markedly asymmetric, concentrated where both markers are high and filtration is
poorest. The Gaussian and Frank copulas cannot produce this and Clayton produces
it in the wrong tail, which is why Gumbel, Joe and survival Clayton are included.

## Scripts

Run from the repository root, in this order.

| Script | Purpose |
|---|---|
| `download_data.R` | Fetch the public NHANES files (~40 MB) into `application-nhanes/data/` |
| `build_data.R` | Merge, apply exclusions, and write `nhanes_kidney.rds` |
| `eda.R` | Histograms, scatterplots on the original and pseudo-observation scales, Kendall's tau by age quintile |
| `fit_nhanes.R` | Fit one copula family and one chain |
| `run_batch.sh` | Launch several chains in parallel and wait for them |
| `results.R` | Convergence diagnostics, both sets of criteria, posterior summaries, figures, LaTeX tables |
| `marginal_criteria.R` | The marginalised criteria on their own, with a quadrature accuracy check |
| `diagnostics.R` | PSRF, MPSRF and effective sample size per family |

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
