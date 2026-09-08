## Final head-to-head between the surviving candidate pairs, all in
## NHANES 2017-2018 where the biomarkers are finely reported.
##
## A) latent CKD severity: ACR (gamma) + log serum creatinine (normal)
## B) latent CKD severity: serum creatinine (gamma) + log ACR (normal)
## C) latent atherogenic dyslipidaemia: triglycerides (gamma) + log HDL (normal)

suppressMessages({ library(haven); library(dplyr) })
rx <- function(f) as.data.frame(read_xpt(file.path("data", paste0(f, ".xpt"))))

skew <- function(x) mean((x - mean(x))^3) / sd(x)^3
kurt <- function(x) mean((x - mean(x))^4) / sd(x)^4 - 3
gran <- function(v) {
  u <- sort(unique(v))
  c(distinct = length(u), maxtie = 100 * max(table(v)) / length(v))
}
## gamma vs lognormal for the raw marker, on a common scale
gvl <- function(v, Z) {
  gl <- glm(v ~ Z, family = Gamma(link = "log"), control = list(maxit = 500))
  mu <- pmax(fitted(gl), 1e-10)
  nll <- function(la) -sum(dgamma(v, shape = exp(la), rate = exp(la) / mu, log = TRUE))
  la <- optimize(nll, c(-6, 8))$minimum
  k  <- ncol(as.matrix(Z)) + 2
  gA <- 2 * nll(la) + 2 * k
  l1 <- lm(log(v) ~ Z); s <- sqrt(mean(resid(l1)^2))
  lA <- -2 * (sum(dnorm(log(v), fitted(l1), s, log = TRUE)) - sum(log(v))) + 2 * k
  c(alpha = exp(la), per_obs = (gA - lA) / length(v))
}

demo <- rx("DEMO_J") %>% select(SEQN, age = RIDAGEYR, sex = RIAGENDR)
bmx  <- rx("BMX_J")  %>% select(SEQN, bmi = BMXBMI, waist = BMXWAIST)
bpx  <- rx("BPX_J")  %>% select(SEQN, BPXSY1, BPXSY2, BPXSY3, BPXSY4)
bio  <- rx("BIOPRO_J") %>% select(SEQN, screat = LBXSCR)
alb  <- rx("ALB_CR_J") %>% select(SEQN, uma = URXUMA, ucr = URXUCR)
tri  <- rx("TRIGLY_J") %>% select(SEQN, trig = LBXTR)
hdl  <- rx("HDL_J")    %>% select(SEQN, hdl = LBDHDD)

bpx$sbp <- rowMeans(bpx[, c("BPXSY1","BPXSY2","BPXSY3","BPXSY4")], na.rm = TRUE)

D <- demo %>% filter(age >= 20) %>%
  left_join(bmx, "SEQN") %>% left_join(bpx[, c("SEQN","sbp")], "SEQN") %>%
  left_join(bio, "SEQN") %>% left_join(alb, "SEQN") %>%
  left_join(tri, "SEQN") %>% left_join(hdl, "SEQN") %>%
  mutate(female = as.numeric(sex == 2),
         acr = ifelse(!is.na(uma) & !is.na(ucr) & ucr > 0, uma / ucr * 100, NA))

report <- function(lab, qv, yv, qname, yname, dsub) {
  Z <- scale(as.matrix(dsub[, c("age","female","bmi","sbp")]))
  Z <- matrix(as.numeric(Z), nrow(dsub))
  gq <- gran(qv); gy <- gran(yv)
  g  <- gvl(qv, Z)
  r  <- resid(lm(yv ~ Z))
  tau <- cor(qv, yv, method = "kendall")
  cat(sprintf("\n%s   (n = %d)\n", lab, nrow(dsub)))
  cat(sprintf("  Q = %-22s distinct %5d  maxtie %5.2f%%   gamma alpha %6.2f  gamma-ln per obs %+.4f\n",
              qname, gq["distinct"], gq["maxtie"], g["alpha"], g["per_obs"]))
  cat(sprintf("  Y = %-22s distinct %5d  maxtie %5.2f%%   resid skew %+.3f  exc.kurt %+.3f\n",
              yname, gy["distinct"], gy["maxtie"], skew(r), kurt(r)))
  cat(sprintf("  Kendall tau = %+.3f\n", tau))
  gg <- cut(dsub$age, quantile(dsub$age, seq(0, 1, length.out = 6)),
            include.lowest = TRUE)
  taus <- sapply(levels(gg), function(l) {
    ii <- which(gg == l); cor(qv[ii], yv[ii], method = "kendall") })
  cat(sprintf("  tau by age quintile: %s  (spread %.3f)\n",
              paste(sprintf("%+.3f", taus), collapse = "  "), diff(range(taus))))
  invisible(NULL)
}

## --- A and B: kidney -----------------------------------------------------
dk <- D %>% filter(!is.na(screat), !is.na(acr), screat > 0, acr > 0,
                   !is.na(bmi), !is.na(sbp), is.finite(sbp))
report("A) latent CKD severity: ACR as gamma", dk$acr, log(dk$screat),
       "urine ACR (mg/g)", "log serum creatinine", dk)
report("B) latent CKD severity: creatinine as gamma", dk$screat, log(dk$acr),
       "serum creatinine (mg/dL)", "log urine ACR", dk)

## --- C: lipids -----------------------------------------------------------
dl <- D %>% filter(!is.na(trig), !is.na(hdl), trig > 0, hdl > 0,
                   !is.na(bmi), !is.na(sbp), is.finite(sbp))
report("C) latent atherogenic dyslipidaemia", dl$trig, log(dl$hdl),
       "triglycerides (mg/dL)", "log HDL cholesterol", dl)
