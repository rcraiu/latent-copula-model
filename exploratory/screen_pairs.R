## Screen candidate pairs of CONTINUOUS surrogates for a latent variable,
## scoring each on the things the gamma-Gaussian model actually needs:
##   (a) both margins genuinely continuous (few ties),
##   (b) one margin well described by a gamma with log link,
##   (c) the other margin well described by a normal after transformation,
##   (d) strong dependence, and dependence that varies with a covariate
##       (which is what the W-varying copula is for).
## NHANES 2017-2018, adults 20+.

suppressMessages({ library(haven); library(dplyr) })
rx <- function(f) as.data.frame(read_xpt(file.path("data", paste0(f, ".xpt"))))

bio  <- rx("BIOPRO_J")  %>% select(SEQN, alt = LBXSATSI, ast = LBXSASSI,
                                   ggt = LBXSGTSI, alb = LBXSAL,
                                   uric = LBXSUA, tp = LBXSTP)
ins  <- rx("INS_J")     %>% select(SEQN, insulin = LBXIN)
glu  <- rx("GLU_J")     %>% select(SEQN, glucose = LBXGLU)
tri  <- rx("TRIGLY_J")  %>% select(SEQN, trig = LBXTR)
hdl  <- rx("HDL_J")     %>% select(SEQN, hdl = LBDHDD)
crp  <- rx("HSCRP_J")   %>% select(SEQN, crp = LBXHSCRP)
fer  <- rx("FERTIN_J")  %>% select(SEQN, ferritin = LBXFER)
demo <- rx("DEMO_J")    %>% select(SEQN, age = RIDAGEYR, sex = RIAGENDR)
bmx  <- rx("BMX_J")     %>% select(SEQN, bmi = BMXBMI, waist = BMXWAIST)

D <- demo %>% filter(age >= 20) %>%
  left_join(bio, "SEQN") %>% left_join(ins, "SEQN") %>% left_join(glu, "SEQN") %>%
  left_join(tri, "SEQN") %>% left_join(hdl, "SEQN") %>% left_join(crp, "SEQN") %>%
  left_join(fer, "SEQN") %>% left_join(bmx, "SEQN") %>%
  mutate(female = as.numeric(sex == 2))

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3

granularity <- function(v) {
  u <- sort(unique(v))
  c(distinct = length(u), maxtie = 100 * max(table(v)) / length(v))
}

## gamma vs lognormal on a common scale (raw-marker density), given covariates
gamma_vs_ln <- function(v, Z) {
  out <- try({
    gl <- glm(v ~ Z, family = Gamma(link = "log"), control = list(maxit = 200))
    mu <- pmax(fitted(gl), 1e-8)
    nll <- function(la) -sum(dgamma(v, shape = exp(la), rate = exp(la) / mu, log = TRUE))
    la <- optimize(nll, c(-4, 8))$minimum
    k  <- ncol(as.matrix(Z)) + 2
    gA <- 2 * nll(la) + 2 * k
    l1 <- lm(log(v) ~ Z); s <- sqrt(mean(resid(l1)^2))
    lA <- -2 * (sum(dnorm(log(v), fitted(l1), s, log = TRUE)) - sum(log(v))) + 2 * k
    c(alpha = exp(la), gamma_minus_ln_per_obs = (gA - lA) / length(v))
  }, silent = TRUE)
  if (inherits(out, "try-error")) c(alpha = NA, gamma_minus_ln_per_obs = NA) else out
}

## how normal is the candidate normal margin, after covariates?
normality <- function(v, Z) {
  r <- resid(lm(v ~ Z))
  c(res.skew = skew(r), res.exc.kurt = kurt(r))
}

pairs <- list(
  list(name = "liver injury",        q = "alt",      y = "ast",      ylog = TRUE),
  list(name = "liver injury (GGT)",  q = "ggt",      y = "alt",      ylog = TRUE),
  list(name = "insulin resistance",  q = "insulin",  y = "glucose",  ylog = TRUE),
  list(name = "dyslipidaemia",       q = "trig",     y = "hdl",      ylog = TRUE),
  list(name = "inflammation",        q = "crp",      y = "alb",      ylog = FALSE),
  list(name = "iron status",         q = "ferritin", y = "tp",       ylog = FALSE),
  list(name = "inflammation (CRP/fer)", q = "crp",   y = "ferritin", ylog = TRUE)
)

cat(sprintf("%-22s %6s | %-9s %7s %6s | %-9s %7s %7s | %6s %7s\n",
            "latent variable", "n", "Q margin", "distinct", "maxtie",
            "Y margin", "res.skew", "kurt", "tau", "gam-ln"))
cat(strrep("-", 118), "\n")

res <- list()
for (P in pairs) {
  d <- D %>% select(all_of(c(P$q, P$y)), age, female, bmi, waist) %>%
    rename(qq = !!P$q, yy = !!P$y) %>%
    filter(complete.cases(.), qq > 0, yy > 0)
  if (nrow(d) < 500) next
  Z <- scale(as.matrix(d[, c("age", "female", "bmi", "waist")]))
  Z <- matrix(as.numeric(Z), nrow(d))
  yv <- if (P$ylog) log(d$yy) else d$yy

  gq <- granularity(d$qq); gy <- granularity(d$yy)
  gl <- gamma_vs_ln(d$qq, Z)
  nn <- normality(yv, Z)
  tau <- cor(d$qq, d$yy, method = "kendall")

  cat(sprintf("%-22s %6d | %-9s %7d %5.1f%% | %-9s %+7.3f %+7.3f | %+6.3f %+7.3f\n",
              P$name, nrow(d), P$q, gq["distinct"], gq["maxtie"],
              paste0(if (P$ylog) "log " else "", P$y),
              nn["res.skew"], nn["res.exc.kurt"], tau, gl["gamma_minus_ln_per_obs"]))

  ## does dependence vary with age?
  g <- cut(d$age, quantile(d$age, seq(0, 1, length.out = 6)),
           include.lowest = TRUE, labels = paste0("A", 1:5))
  taus <- sapply(levels(g), function(l) {
    ii <- which(g == l); cor(d$qq[ii], d$yy[ii], method = "kendall")
  })
  res[[P$name]] <- list(n = nrow(d), taus = taus, gl = gl, nn = nn,
                        gq = gq, gy = gy, tau = tau)
}

cat("\n--- Kendall's tau by age quintile (dependence varying with a covariate) ---\n")
for (nm in names(res)) {
  cat(sprintf("%-22s %s   (spread %.3f)\n", nm,
              paste(sprintf("%+.3f", res[[nm]]$taus), collapse = "  "),
              diff(range(res[[nm]]$taus))))
}

saveRDS(res, "output/screen_pairs.rds")
