## Vectorised log copula densities and conditional copulas, with the copula
## parameter allowed to vary across observations (which is what the latent
## variable model requires, and what the packaged implementations do not allow).
##
## NOTE, relative to LVC_functions.R:
##   * log_dGaussian there omits the +(z1^2+z2^2)/2 term, so it returns
##     log phi_2 up to a constant rather than the copula density. That
##     cancels in a Metropolis ratio at fixed margins, but not in the W step
##     or in WAIC/CVML.
##   * log_dIndep there returns log(u1)+log(u2); the independence copula
##     density is 1, so the log density is 0.
##   * the clamping uses eps = .Machine$double.xmin, and 1 - 2.2e-308 == 1 in
##     double precision, so the upper clamp has no effect.

CLAMP <- 1e-12

clamp01 <- function(u) pmin(pmax(u, CLAMP), 1 - CLAMP)

## --- activation functions: copula parameter as a function of w ----------
## gscale plays the role of the scaling constant gamma in the paper.
theta_of_w <- function(w, xi, family, gscale = 0.1) {
  lin <- gscale * (xi[1] + xi[2] * w)
  switch(family,
    Gaussian    = pmax(pmin(tanh(lin), 1 - 1e-8), -1 + 1e-8),
    Frank       = ifelse(abs(lin) < 1e-6, 1e-6, lin),
    Gumbel      = 1 + exp(lin),
    Clayton     = exp(lin),
    ## Joe shares Gumbel's parameter space, survival Clayton shares Clayton's
    Joe         = 1 + exp(lin),
    SurvClayton = exp(lin),
    Independence = rep(0, length(w)),
    stop("unknown family: ", family)
  )
}

## --- log copula densities ----------------------------------------------
log_c_Gaussian <- function(u1, u2, rho) {
  z1 <- qnorm(clamp01(u1)); z2 <- qnorm(clamp01(u2))
  rho <- pmax(pmin(rho, 1 - 1e-8), -1 + 1e-8)
  om <- 1 - rho^2
  d <- -0.5 * log(om) - (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * om) + (z1^2 + z2^2) / 2
  d[!is.finite(d)] <- -Inf
  d
}

log_c_Clayton <- function(u1, u2, th) {
  u1 <- clamp01(u1); u2 <- clamp01(u2)
  th <- pmax(th, 1e-10)
  A <- u1^(-th) + u2^(-th) - 1
  d <- log1p(th) - (1 + th) * (log(u1) + log(u2)) - (2 + 1 / th) * log(A)
  d[!is.finite(d)] <- -Inf
  d
}

log_c_Gumbel <- function(u1, u2, th) {
  u1 <- clamp01(u1); u2 <- clamp01(u2)
  th <- pmax(th, 1 + 1e-10)
  x <- -log(u1); y <- -log(u2)
  S <- x^th + y^th
  A <- S^(1 / th)
  d <- -A - log(u1) - log(u2) + (th - 1) * (log(x) + log(y)) +
    (1 - 2 * th) * log(A) + log(A + th - 1)
  d[!is.finite(d)] <- -Inf
  d
}

log_c_Frank <- function(u1, u2, th) {
  u1 <- clamp01(u1); u2 <- clamp01(u2)
  th <- ifelse(abs(th) < 1e-8, 1e-8, th)
  ## c = th*(1-e^-th)*e^{-th(u+v)} / [ (1-e^-th) - (1-e^{-th u})(1-e^{-th v}) ]^2
  g  <- function(x) expm1(-th * x)          # e^{-th x} - 1
  num <- -th * g(1) * (1 + g(u1 + u2))      # = th (1-e^-th) e^{-th(u+v)}
  den <- g(u1) * g(u2) + g(1)
  d <- log(num) - 2 * log(abs(den))
  d[!is.finite(d)] <- -Inf
  d
}

## Joe: upper tail dependence 2 - 2^(1/th), no lower tail dependence.
## c = A^(1/th - 2) ub^(th-1) vb^(th-1) (th - 1 + A),
##   A = ub^th + vb^th - ub^th vb^th,  ub = 1 - u1, vb = 1 - u2.
## A is accumulated as ub^th (1 - vb^th) + vb^th, a sum of positive terms, to
## avoid the cancellation that the written form suffers when th is large.
log_c_Joe <- function(u1, u2, th) {
  a <- log1p(-clamp01(u1)); b <- log1p(-clamp01(u2))   # log ub, log vb
  th <- pmax(th, 1 + 1e-10)
  ea <- exp(th * a); eb <- exp(th * b)
  A <- ea * (-expm1(th * b)) + eb
  d <- (1 / th - 2) * log(A) + (th - 1) * (a + b) + log(th - 1 + A)
  d[!is.finite(A) | A <= 0] <- -Inf
  d[!is.finite(d)] <- -Inf
  d
}

## Survival (180 degree rotated) Clayton: the Clayton density in 1 - u.
## Included because Clayton's lower tail dependence is the wrong orientation
## for these data, so this isolates orientation from family.
log_c_SurvClayton <- function(u1, u2, th) {
  a <- log1p(-clamp01(u1)); b <- log1p(-clamp01(u2))
  th <- pmax(th, 1e-10)
  A <- exp(-th * a) + exp(-th * b) - 1
  d <- log1p(th) - (1 + th) * (a + b) - (2 + 1 / th) * log(A)
  d[!is.finite(A) | A <= 0] <- -Inf
  d[!is.finite(d)] <- -Inf
  d
}

log_c_Independence <- function(u1, u2, th) rep(0, length(u1))

log_copula_density <- function(u1, u2, th, family) {
  switch(family,
    Gaussian     = log_c_Gaussian(u1, u2, th),
    Clayton      = log_c_Clayton(u1, u2, th),
    Gumbel       = log_c_Gumbel(u1, u2, th),
    Frank        = log_c_Frank(u1, u2, th),
    Joe          = log_c_Joe(u1, u2, th),
    SurvClayton  = log_c_SurvClayton(u1, u2, th),
    Independence = log_c_Independence(u1, u2, th),
    stop("unknown family: ", family)
  )
}

## --- validation against VineCopula at fixed parameter values ------------
if (sys.nframe() == 0) {
  suppressMessages(library(VineCopula))
  set.seed(1)
  u1 <- runif(200, 0.01, 0.99); u2 <- runif(200, 0.01, 0.99)
  chk <- function(family, vfam, par) {
    mine <- log_copula_density(u1, u2, rep(par, length(u1)), family)
    theirs <- log(VineCopula::BiCopPDF(u1, u2, family = vfam, par = par))
    cat(sprintf("%-13s par=%5.2f  max abs diff = %.3e\n",
                family, par, max(abs(mine - theirs))))
  }
  chk("Gaussian", 1,  0.6);  chk("Gaussian", 1, -0.4)
  chk("Clayton",  3,  2.0);  chk("Clayton",  3,  0.5)
  chk("Gumbel",   4,  1.8);  chk("Gumbel",   4,  3.0)
  chk("Frank",    5,  4.0);  chk("Frank",    5, -3.0)
  ## VineCopula: 6 = Joe, 13 = survival (180 degree) Clayton
  chk("Joe",         6,  2.0); chk("Joe",         6,  5.0)
  chk("SurvClayton", 13, 2.0); chk("SurvClayton", 13, 0.5)
  cat("Independence max abs =",
      max(abs(log_copula_density(u1, u2, rep(0, 200), "Independence"))), "\n")
}
