# Application: kidney function

The glomerular filtration rate is the standard summary of kidney function and is
effectively unobservable in practice, since measuring it directly requires
infusing an exogenous marker and taking timed samples. Clinical practice instead
uses endogenous filtration markers, whose serum concentration rises as filtration
falls. Two are in routine use, and they are affected by different extra-renal
factors: serum creatinine depends on muscle mass, and hence on age, sex and body
composition, while serum cystatin C is more sensitive to thyroid status,
corticosteroid use and adiposity. Neither is a clean measurement of filtration on
its own, which makes them a natural surrogate pair for a common latent quantity.

| | |
|---|---|
| `Q` | serum cystatin C, mg/L, gamma margin |
| `Y` | log serum creatinine, mg/dL, normal margin |
| `W` | latent renal impairment |
| covariates | age, sex, BMI, systolic blood pressure, standardised |

Note the orientation: both markers rise as filtration falls, and both loadings
are positive, so `W` increases with *impairment*. Its posterior mean correlates
negatively with CKD-EPI eGFR.

## Data

NHANES 1999-2002, the cycles in which cystatin C was assayed on a surplus-sera
subsample alongside the standard biochemistry panel. After restricting to adults
with complete data and excluding participants with eGFR below 15 mL/min/1.73m^2
(CKD stage G5, where dialysis decouples marker concentrations from residual
filtration), n = 3948.

The margins are chosen to match the shapes of the two distributions: cystatin C
is strictly positive with sample skewness 2.58, and log creatinine is nearly
symmetric with sample skewness 0.13. Kendall's tau between them is 0.387.

Age and sex enter as covariates for `W` rather than being built into the
response, which is a deliberate contrast with the CKD-EPI estimating equations
that bake them into the reported quantity. The surrogates stay raw measurements.

**A known limitation.** Serum creatinine in these cycles is reported on a coarse
0.1 mg/dL grid, giving only 31 distinct values among 3948 observations. It is
treated as continuous, which is a genuine approximation. Later NHANES cycles
report creatinine more finely but do not include cystatin C. The paper's
cardiotocography example has comparable granularity in `MSTV`.

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
download step.

## Two points worth knowing before interpreting output

**Stratify on covariates, not on the markers.** `eda.R` examines whether the
association between the markers varies, by computing Kendall's tau within age
quintiles. An earlier version stratified on `(u1 + u2)/2`, a function of both
margins; conditioning on such a quantity induces a selection effect that
attenuates the within-stratum association and produced a spurious negative trend.
Stratifying on eGFR would be circular for the same reason, since eGFR is a
deterministic function of creatinine.

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
