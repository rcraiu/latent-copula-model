library(tidyverse)
library(mvtnorm)
library(bayesplot)
library(truncnorm)
library(zipfR)
library(MHN)
library(doParallel)
library(foreach)
library(doRNG)
library(stableGR)
library(scales)
library(ggpubr)
library(copula)

rm(list = ls())
gc()
set.seed(1729)

source("LVC_functions.R")

##### FIXED PARAMETERS

nu <- 1
gam <- 0.1

##### PRIOR PARAMETERS

# for beta
betadim <- 3
mubeta <- rep(0, times=betadim)
Sbeta <- diag(rep(1000, times=betadim))

# for xi
muxi <- c(-1,3)
Sxi0 <- 1
Sxi1 <- 1

# for lambda
mulamc <- 0
mulamd <- 0
ssqlamc <- 100
ssqlamd <- 100

# for mu
ssqmuc <- 100
ssqmud <- 100

# for ssq
a <- 1000
b <- 1000

M <- 100000 # number of iterations

BB <- 1 # number of data replicates

load("Gauss-FOE.Rdata")


# cl <- makePSOCKcluster(BB) # set the number of processor cores
# registerDoParallel(cl)
# registerDoRNG(1729, once=FALSE)

# MCMC <- foreach(bb = 1:BB, .errorhandling = "pass", .options.RNG = 1729, .packages = c("truncnorm", "zipfR", "MHN", "mvtnorm", "invgamma")) %dopar% {

MCMC <- vector("list", BB)
for (bb in 1:BB) {
  
  DAT <- DAT.Ga[[bb]]
  
  y <- DAT$y
  q <- DAT$q
  X <- DAT$X
  n <- length(y)
  w <- DAT$w
  
  betapost_cov <- solve(t(X)%*%X) # covariance for beta posterior
  
  
  ##### PREPARE THE GIBBS SAMPLER
  
  parnames <- c("beta0", "beta1", "beta2", "mud", "lamd", "muc", "lamc", "ssq", "xi0", "xi1")
  PAR <- matrix(0L, nrow=length(parnames), ncol=M)
  rownames(PAR) <- parnames
  
  mean_pM <- rep(0, n)        # Initialize mean of p_M(q_i, y_i | ...)
  mean_logpM <- rep(0, n)     # Initialize mean of log(p_M(q_i, y_i | ...))
  var_logpM <- rep(0, n)      # Initialize variance of log(p_M(q_i, y_i | ...))
  
  mean_pM_cond <- rep(0, n)   # Initialize mean of p_M(q_i | y_i, ...)
  mean_logpM_cond <- rep(0, n) # Initialize mean of log(p_M(q_i | y_i, ...))
  var_logpM_cond <- rep(0, n)  # Initialize variance of log(p_M(q_i | y_i, ...))
  
  B <- floor(4*M/5)       # Start updating after burn-in (4/5 of M)
  
  
  
  ##### INITIALIZE THINGS
  
  
  ## Initialize beta
  
  beta_OLS <- lm(y - mean(y) ~ 0 + ., data=as.data.frame(X))$coefficients
  PAR[c("beta0", "beta1", "beta2"), 1] <- beta_OLS
  beta_cur <- PAR[c("beta0", "beta1", "beta2"), 1]
  

  ## Initialize W
  
  W <- as.vector(X %*% beta_cur)
  W <- w
  
  
  ## Initialize lambda_d and mu_d
  
  ll_mud_lamd <- function(mmud, llamd) {
    sum(dpois(q, lambda = exp(mmud + llamd*W), log=TRUE))
  }
  
  PAR[c("mud", "lamd"), 1] <- optim(par=c(1,1), fn=function(pr) ll_mud_lamd(pr[1], pr[2]), control=list(fnscale=-1))$par
  mud_cur <- PAR["mud", 1]
  lamd_cur <- PAR["lamd", 1]
  
  
  ## Initialize lambda_c, mu_c, and sigma^2
  
  ll_muc_lamc_ssq <- function(mmuc, llamc, logsd) {
    sum(dnorm(y, mean = mmuc + llamc*W, sd = exp(logsd), log=TRUE))
  }
  
  PAR[c("muc", "lamc", "ssq"), 1] <- optim(par=c(0,0,0), fn=function(pr) ll_muc_lamc_ssq(pr[1], pr[2], pr[3]), control=list(fnscale=-1))$par
  PAR["ssq", 1] <- exp(2*PAR["ssq", 1])
  muc_cur <- PAR["muc", 1]
  lamc_cur <- PAR["lamc", 1]
  ssq_cur <- PAR["ssq", 1]
  
  
  ## Initialize xi_0 and xi_1
  
  ll_xi <- function(xxi0, xxi1) {
    u1 <- ppois(q, lambda = exp(mud_cur + lamd_cur*W))
    u1m1 <- ppois(q-1, lambda = exp(mud_cur + lamd_cur*W))
    u1m1[which(u1m1 == 1)] <- 1 - .Machine$double.eps
    u2 <- pnorm(y, mean=muc_cur + lamc_cur*W)
    thet <- thetac_Gauss(W, c(xxi0, xxi1))
    
    lp <- sum(log(P_Gauss(u1, u2, thet) - P_Gauss(u1m1, u2, thet)))
    
    return(lp)
  }
  
  PAR[c("xi0", "xi1"), 1] <- optim(par=c(1,1), fn=function(pr) ll_xi(pr[1], pr[2]), control=list(fnscale=-1))$par
  
  
  
  ##### PREPARE ADAPTATION
  
  gam_AM <- function(t) {t^(-0.5)} 

  mu_xi <- rep(0, 2)  # Mean vector for adaptation
  Sig_xi <- diag(2)   # Covariance matrix for adaptation
  lam_xi <- rep(2.38^2 / 2, 2)  # Scaling factors
  acc_xi <- 0  # Initialize acceptance counter

  # mu_Z <- rep(0, times=n)
  # sig_Z <- rep(1, times=n)
  # lam_Z <- rep(2.38, times=n)
  # 
  mu_W <- rep(0, times=n)  # Mean vector for adaptation
  sig_W <- rep(1, times=n) # Variance for adaptation
  lam_W <- rep(2.38, times=n) # Scaling factor
  acc_W <- 0*1:n # Initialize acceptance counter
  
  acc_d <- 0

  
  #### SAMPLE SOME STUFF
  
  for (t in 2:M) {
    
    beta_cur <- PAR[c("beta0", "beta1", "beta2"), t-1]
    mud_cur <- PAR["mud", t-1]
    muc_cur <- PAR["muc", t-1]
    lamd_cur <- PAR["lamd", t-1]
    lamc_cur <- PAR["lamc", t-1]
    xi0_cur <- PAR["xi0", t-1]
    xi1_cur <- PAR["xi1", t-1]
    xi_cur <- PAR[c("xi0", "xi1"), t-1]
    ssq_cur <- PAR["ssq", t-1]
    theta_cur <- thetac_Gauss(W, xi_cur)
    
    
    
    
    ## Step 1: Sample discrete marginal parameters
    
    lp_mud_lamd <- function(mmud, llamd) {
      u1 <- ppois(q, lambda = exp(mmud + llamd*W))
      u1m1 <- ppois(q-1, lambda = exp(mmud + llamd*W))
      u2 <- pnorm(y, mean=muc_cur + lamc_cur*W, sd=sqrt(ssq_cur))
      thet <- thetac_Gauss(W, xi_cur)
      
      lp <- sum(log(P_Gauss(u1, u2, thet) - P_Gauss(u1m1, u2, thet))) +
        dnorm(mmud, sd=sqrt(ssqmud), log=TRUE) + 
        dnorm(llamd, mean=mulamc, sd=sqrt(ssqlamc), log=TRUE)
      
      return(lp)
    }
    
    d_prop <- optim(par=c(mud_cur, lamd_cur), fn=function(pr) lp_mud_lamd(pr[1], pr[2]), control=list(fnscale=-1), hessian=TRUE)
    d_prop_means <- d_prop$par
    d_prop_ihess <- diag(solve(d_prop$hessian))
    
    mud_prop <- extraDistr::rlst(n=1, df=7, mu=d_prop_means[1], sigma=-d_prop_ihess[1])
    lamd_prop <- extraDistr::rlst(n=1, df=7, mu=d_prop_means[2], sigma=-d_prop_ihess[2])
    
    # Current log-posterior
    log_posterior_cur <- lp_mud_lamd(mud_cur, lamd_cur)
    
    # Compute the log-posterior for the proposed values
    log_posterior_prop <- lp_mud_lamd(mud_prop, lamd_prop)
    
    # Compute proposal densities for the current and proposed values
    log_prop_density_cur_to_prop <- extraDistr::dlst(mud_prop, df=7, mu=d_prop_means[1], sigma=-d_prop_ihess[1], log=TRUE) + 
      extraDistr::dlst(lamd_prop, df=7, mu=d_prop_means[2], sigma=-d_prop_ihess[2], log=TRUE)
    
    log_prop_density_prop_to_cur <- extraDistr::dlst(mud_cur, df=7, mu=d_prop_means[1], sigma=-d_prop_ihess[1], log=TRUE) + 
      extraDistr::dlst(lamd_cur, df=7, mu=d_prop_means[2], sigma=-d_prop_ihess[2], log=TRUE)
    
    # Compute the full Metropolis-Hastings acceptance ratio
    log_accept_ratio_d <- (log_posterior_prop - log_posterior_cur) + 
      (log_prop_density_prop_to_cur - log_prop_density_cur_to_prop)
    
    # Accept or reject the proposal
    if (log(runif(1)) < log_accept_ratio_d) {
      PAR[c("mud", "lamd"), t] <- c(mud_prop, lamd_prop)
      acc_d <- acc_d + 1
      accrate_d <- acc_d / t
    } else {
      PAR[c("mud", "lamd"), t] <- PAR[c("mud", "lamd"), t-1]
    }
    
    mud_cur <- PAR["mud", t]
    lamd_cur <- PAR["lamd", t]
    
    
    
    
    ## Step 2: Sample Z
    
    Z <- rtruncnorm(n=n, 
                    a = qnorm(ppois(q-1, lambda=exp(mud_cur + lamd_cur*W))),
                    b = qnorm(ppois(q, lambda=exp(mud_cur + lamd_cur*W))),
                    mean=theta_cur*(y - muc_cur - lamc_cur*W)/sqrt(ssq_cur),
                    sd=sqrt(1 - theta_cur^2))
    
    
    
    
    ## Step 3: Sample continuous marginal parameters
    
    lamcvar <- 1/(sum(W^2/(1 - theta_cur^2))/ssq_cur + 1/ssqlamc)
    lamcmean <- (sum(W*(1/(1 - theta_cur^2))*(y - muc_cur - sqrt(ssq_cur)*theta_cur*Z))/ssq_cur + mulamc/ssqlamc)*lamcvar
    
    PAR["lamc", t] <- rnorm(n=1, mean=lamcmean, sd=sqrt(lamcvar))
    lamc_cur <- PAR["lamc", t]
    
    
    mucvar <- 1/(sum(1/(1 - theta_cur^2))/ssq_cur + 1/ssqmuc)
    mucmean <- (sum( (1/(1 - theta_cur^2))*(y - lamc_cur*W - sqrt(ssq_cur)*theta_cur*Z))/ssq_cur)*mucvar
    
    PAR["muc", t] <- rnorm(n=1, mean=mucmean, sd=sqrt(mucvar))
    muc_cur <- PAR["muc", t]
    
    
    PAR["ssq", t] <- 1/rMHN(n=1, 
                            alpha=n+a,
                            beta=(sum( (y - muc_cur - lamc_cur*W)^2/(1 - theta_cur^2) ) + b)/2,
                            gamma=sum( (y - muc_cur - lamc_cur*W)*(1/(1 - theta_cur^2))*theta_cur*Z))^2
    ssq_cur <- PAR["ssq", t]
    
    
    
    
    ## Step 4: Sample copula parameters
    
    # Log-posterior function for xi_0 and xi_1
    lp_xi <- function(x0, x1) {
      lp <- sum(log_dGauss(cbind(pnorm(Z), pnorm((y - muc_cur - lamc_cur*W) / sqrt(ssq_cur))), 
                       thetac_Gauss(W, c(x0, x1)))) +
        dnorm(x0, mean = muxi[1], sd = sqrt(Sxi0), log = TRUE) + 
        dnorm(x1, mean = muxi[2], sd = sqrt(Sxi1), log = TRUE)
      return(lp)
    }
    
    # Propose new values for xi_0 and xi_1 from multivariate normal distribution
    xi_prop <- xi_cur + rmvnorm(n = 1, sigma = diag(sqrt(lam_xi)) %*% Sig_xi %*% diag(sqrt(lam_xi)))

    # Current log-posterior
    log_posterior_xi_cur <- lp_xi(xi0_cur, xi1_cur)
    
    # Compute the log-posterior for the proposed values
    log_posterior_xi_prop <- lp_xi(xi_prop[1], xi_prop[2])
    
    
    # Compute the full Metropolis-Hastings acceptance ratio
    log_accept_ratio_xi <- log_posterior_xi_prop - log_posterior_xi_cur
    
    # Accept or reject the proposal
    if (is.infinite(log_posterior_xi_cur)) {
      # If the log acceptance difference is infinite, retain current values
      PAR[c("xi0", "xi1"), t] <- xi_cur
    } else if (log(runif(1)) < log_accept_ratio_xi) {
      PAR[c("xi0", "xi1"), t] <- xi_prop
      acc_xi <- acc_xi + 1
    } else {
      PAR[c("xi0", "xi1"), t] <- xi_cur
    }
    
    xi_cur <- PAR[c("xi0", "xi1"), t]
    xi0_cur <- PAR[c("xi0"), t]
    xi1_cur <- PAR[c("xi1"), t]
    theta_cur <- thetac_Gauss(W, c(xi0_cur, xi1_cur))
    
    # Calculate the acceptance rate
    accrate_xi <- acc_xi / t
    
    # Adaptation of lam_xi for both xi_0 and xi_1
    lam_xi[1] <- lam_xi[1] * exp(gam_AM(t) * (exp(lp_xi(xi_prop[1], xi1_cur)) - 0.44))
    lam_xi[2] <- lam_xi[2] * exp(gam_AM(t) * (exp(lp_xi(xi0_cur, xi_prop[2])) - 0.44))
    
    # Update mu_xi and Sig_xi based on the new xi_cur
    mu_xi_old <- mu_xi
    Sig_xi_old <- Sig_xi
    
    mu_xi <- mu_xi_old + gam_AM(t) * (xi_cur - mu_xi_old)
    Sig_xi <- Sig_xi_old + gam_AM(t) * ((xi_cur - mu_xi_old) %*% t(xi_cur - mu_xi_old) - Sig_xi_old)
    

    
    
    ## Step 5: Sample beta
    
    betapost_mean <- betapost_cov %*% (t(X)%*%W)
    PAR[c("beta0", "beta1", "beta2"), t] <- rmvnorm(n=1, mean=betapost_mean, sigma=betapost_cov)
    beta_cur <- PAR[c("beta0", "beta1", "beta2"), t]
    
    
    

    ## Step 6: Sample W

    # Log-posterior for all W_i's (vectorized)
    lp_W <- function(WW) {
      u1 <- ppois(q, lambda = exp(mud_cur + lamd_cur*WW))
      u1m1 <- ppois(q-1, lambda = exp(mud_cur + lamd_cur*WW))
      u2 <- pnorm(y, mean=muc_cur + lamc_cur*WW)
      thet <- thetac_Gauss(WW, xi_cur)

      lp <- log(P_Gauss(u1, u2, thet) - P_Gauss(u1m1, u2, thet)) +
        dnorm(y, mean=muc_cur + lamc_cur*WW, log=TRUE) +
        dnorm(WW, mean=as.vector(X%*%beta_cur), sd=1, log=TRUE)

      return(lp)
    }

    # Propose new values for W using random walk Metropolis
    W_prop <- W + rnorm(n, mean = 0, sd = sqrt(lam_W*sig_W))  # Univariate proposals
    W_prop <- W + rnorm(n, mean = 0, sd = 1.5)  # Univariate proposals

    # Compute the log-posterior for the current and proposed W values
    log_posterior_W_cur <- lp_W(W)
    log_posterior_W_prop <- lp_W(W_prop)

    # Compute the log acceptance ratio (vectorized for all W_i's)
    log_accept_ratio_W <- pmin(0, log_posterior_W_prop - log_posterior_W_cur)
    log_accept_ratio_W[which(is.nan(log_accept_ratio_W))] <- -Inf

    # Accept or reject the proposals
    MH_bool <- log(runif(n)) < log_accept_ratio_W
    W <- ifelse(MH_bool, W_prop, W)

    # Calculate the acceptance rate
    accrate_W <- acc_W / t
    acc_W <- acc_W + as.numeric(MH_bool)

    # Adapt the proposal scale (lam_W) for each W_i
    lam_W <- lam_W * exp(gam_AM(t) * (exp(log_accept_ratio_W) - 0.44))
    sig_W_old <- sig_W

    # Update adaptation parameters (mu_W, Sig_W) for each W_i (univariate updates)
    mu_W_old <- mu_W
    mu_W <- mu_W_old + gam_AM(t) * (W - mu_W_old)  # Update mean
    sig_W <- sig_W_old + gam_AM(t) * ((W - mu_W_old)^2 - sig_W_old)  # Update variance (scalar for univariate case)# # Log-posterior for all W_i's (vectorized)
    # lp_W <- function(WW) {
    #   u1 <- ppois(q, lambda = exp(mud_cur + lamd_cur*WW))
    #   u1m1 <- ppois(q-1, lambda = exp(mud_cur + lamd_cur*WW))
    #   u2 <- pnorm(y, mean=muc_cur + lamc_cur*WW)
    #   thet <- thetac_Gauss(WW, xi_cur)
    # 
    #   lp <- log(P_Gauss(u1, u2, thet) - P_Gauss(u1m1, u2, thet)) +
    #     dnorm(y, mean=muc_cur + lamc_cur*WW, log=TRUE) +
    #     dnorm(WW, mean=as.vector(X%*%beta_cur), sd=1, log=TRUE)
    # 
    #   return(lp)
    # }
    # 
    # # Propose new values for W using random walk Metropolis
    # W_prop <- W + rnorm(n, mean = 0, sd = sqrt(lam_W*sig_W))  # Univariate proposals
    # W_prop <- W + rnorm(n, mean = 0, sd = 1.5)  # Univariate proposals
    # 
    # # Compute the log-posterior for the current and proposed W values
    # log_posterior_W_cur <- lp_W(W)
    # log_posterior_W_prop <- lp_W(W_prop)
    # 
    # # Compute the log acceptance ratio (vectorized for all W_i's)
    # log_accept_ratio_W <- pmin(0, log_posterior_W_prop - log_posterior_W_cur)
    # log_accept_ratio_W[which(is.nan(log_accept_ratio_W))] <- -Inf
    # 
    # # Accept or reject the proposals
    # MH_bool <- log(runif(n)) < log_accept_ratio_W
    # W <- ifelse(MH_bool, W_prop, W)
    # 
    # # Calculate the acceptance rate
    # accrate_W <- acc_W / t
    # acc_W <- acc_W + as.numeric(MH_bool)
    # 
    # # Adapt the proposal scale (lam_W) for each W_i
    # lam_W <- lam_W * exp(gam_AM(t) * (exp(log_accept_ratio_W) - 0.44))
    # sig_W_old <- sig_W
    # 
    # # Update adaptation parameters (mu_W, Sig_W) for each W_i (univariate updates)
    # mu_W_old <- mu_W
    # mu_W <- mu_W_old + gam_AM(t) * (W - mu_W_old)  # Update mean
    # sig_W <- sig_W_old + gam_AM(t) * ((W - mu_W_old)^2 - sig_W_old)  # Update variance (scalar for univariate case)

    
    
    # # Update stuff for the WAIC and CWAIC
    # 
    # if (t == (B + 1)) {
    #   # Initialize joint density updates
    #   p_t <- p_Gauss(yi = y, Wi = W, Zi = Z, muc = muc_cur, lamc = lamc_cur, ssq = ssq_cur, xi0 = xi0_cur, xi1 = xi1_cur)
    #   mean_pM <- p_t
    #   mean_logpM <- log(p_t)
    #   var_logpM <- 0
    #   
    #   # Initialize conditional density updates
    #   p_t_cond <- pcond_Gauss(yi = y, Wi = W, Zi = Z, muc = muc_cur, lamc = lamc_cur, ssq = ssq_cur, xi0 = xi0_cur, xi1 = xi1_cur)
    #   mean_pM_cond <- p_t_cond
    #   mean_logpM_cond <- log(p_t_cond)
    #   var_logpM_cond <- 0
    #   
    #   # Error check to stop if NaN occurs
    #   if (any(is.nan(mean_logpM)) || any(is.nan(mean_logpM_cond))) stop("something is NaN")
    #   
    # } else if (t > B) {
    #   # Update joint density after burn-in
    #   p_t <- p_Gauss(yi = y, Wi = W, Zi = Z, muc = muc_cur, lamc = lamc_cur, ssq = ssq_cur, xi0 = xi0_cur, xi1 = xi1_cur)
    #   
    #   # Iterative update of joint mean_pM
    #   mean_pM <- mean_pM + (p_t - mean_pM) / (t - B)  # Updating the mean
    #   
    #   # Iterative update of joint mean_logpM
    #   log_p_t <- log(p_t)
    #   delta_logp <- log_p_t - mean_logpM
    #   mean_logpM <- mean_logpM + delta_logp / (t - B)  # Updating the log-mean
    #   
    #   # Iterative update of joint var_logpM (Welford's method)
    #   var_logpM <- var_logpM + delta_logp * (log_p_t - mean_logpM)  # Updating the variance
    #   
    #   # Update conditional density after burn-in
    #   p_t_cond <- pcond_Gauss(yi = y, Wi = W, Zi = Z, muc = muc_cur, lamc = lamc_cur, ssq = ssq_cur, xi0 = xi0_cur, xi1 = xi1_cur)
    #   
    #   # Iterative update of conditional mean_pM_cond
    #   mean_pM_cond <- mean_pM_cond + (p_t_cond - mean_pM_cond) / (t - B)  # Updating the mean
    #   
    #   # Iterative update of conditional mean_logpM_cond
    #   log_p_t_cond <- log(p_t_cond)
    #   delta_logp_cond <- log_p_t_cond - mean_logpM_cond
    #   mean_logpM_cond <- mean_logpM_cond + delta_logp_cond / (t - B)  # Updating the log-mean
    #   
    #   # Iterative update of conditional var_logpM_cond (Welford's method)
    #   var_logpM_cond <- var_logpM_cond + delta_logp_cond * (log_p_t_cond - mean_logpM_cond)  # Updating the variance
    #   
    #   # Error check to stop if NaN occurs
    #   if (any(is.nan(mean_logpM)) || any(is.nan(var_logpM)) || 
    #       any(is.nan(mean_logpM_cond)) || any(is.nan(var_logpM_cond))) stop("something is NaN")
    # }
    
    print(t)
    print(PAR[,t])
    
  }
  
  # # Compute WAIC
  # WAIC <- -2 * (sum(log(mean_pM)) - sum(var_logpM) / (M - B))
  # 
  # # Compute conditional WAIC
  # WAIC_cond <- -2 * (sum(log(mean_pM_cond)) - sum(var_logpM_cond) / (M - B))
  #   
  # list(PAR=PAR, WAIC=WAIC, WAIC_cond=WAIC_cond)
  
}

# stopImplicitCluster()
# registerDoSEQ()
# bayesplot::mcmc_trace(t(PAR))
