# Exploratory screening

These scripts record how the surrogate pair and the marginal families were
chosen. They are not part of the reported analysis and are kept for
transparency about the specification search.

| Script | Question |
|---|---|
| `screen_pairs.R` | Which candidate NHANES pairs behave like two surrogates for one latent quantity? |
| `screen_hom.R` | Alternative pairings in the 2017-2018 cycle |
| `screen_final.R` | Final comparison of the surviving candidates |
| `margin_choice.R`, `margin_choice2.R` | Gamma against lognormal for the positive margin |
| `check_acr.R` | Urine albumin-creatinine ratio as an alternative second surrogate |

Candidates considered were kidney function, insulin resistance, dyslipidaemia,
inflammation, iron status and liver injury. The kidney pairing was chosen for the
strength of its scientific rationale — two markers in routine clinical use for
one unobservable quantity, with different and well-documented extra-renal
determinants — and because its dependence varies with age, which is what the
`theta(W, xi)` component of the model is for.

Two pitfalls encountered here are worth recording, since both produced
misleading numbers before being caught.

**Comparing likelihoods across transformations.** An early comparison of gamma
against lognormal for cystatin C omitted the Jacobian when the lognormal was fit
on the log scale, making the gamma look much worse than it is. With the
correction, and conditioning on a latent proxy, the gap is small.

**Stratifying on a function of both margins.** Assessing whether dependence
varies by splitting on `(u1 + u2)/2` induces a selection effect that attenuates
the within-stratum association, producing a spurious negative trend. The reported
analysis stratifies on age instead.
