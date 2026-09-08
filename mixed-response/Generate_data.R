library(tidyverse)
library(mvtnorm)
library(bayesplot)
library(MHN)
library(readxl)
library(ggpubr)
library(Cairo)
library(viridis)
library(truncnorm)

rm(list = ls())
gc()
source("LVC_functions.R")
set.seed(1729)


##### FIXED PARAMETERS

gam <- 0.1
alpha <- 1
nu <- 1

BB <- 1 # number of independent datasets to generate for each model (with fixed covariates)

n <- 500 # number of observations per dataset


##### FIXED (TRUE) PARAMETERS

alpha <- 1
beta <- c(1, 1, 0.5)
betadim <- length(beta)
xi0 <- -1
xi1 <- 3
xi <- c(xi0, xi1)
lamd <- 1
mud <- -1
lamc <- 1
muc <- 0
ssq <- 1


# Generate covariates

X <- matrix(0L, nrow=n, ncol=betadim)
for (i in 1:n) {
  X[i,] <- c(1, runif(n=1, min=-3, max=1), runif(n=1, min=-1, max=1))
}

#### GAUSSIAN COPULA


DAT.Ga <- vector("list", length=BB)

# Generate dataset replicates

for (j in 1:BB) {
  
  y <- 0*1:n
  q <- 0*1:n
  z <- 0*1:n
  w <- 0*1:n
  tau.Ga <- 0*1:n
  
  
  for (i in 1:n) {
    
    w[i] <- rnorm(n=1, mean=t(beta)%*%X[i,], sd=1)
    
    ui <- copula::rCopula(n=1, copula=copula::normalCopula(param=thetac_Gauss(w[i], xi), dim = 2))
    tau.Ga[i] <- copula::tau(copula::normalCopula(thetac_Gauss(w[i], xi)))
    
    # p_pois <- ppois((-1):100, lambda=exp(mud + lamd*w[i]), lower.tail=F)
    # p_pois <- pmin(pmax(p_pois, .Machine$double.eps), 1 - .Machine$double.eps) # Clamping
    # alphF_i <- qnorm(p_pois, sd=nu, lower.tail=F)
    # z[i] <- qnorm(ui[1], sd=nu)
    # q[i] <- findInterval(z[i], alphF_i) - 1
    q[i] <- qpois(ui[1], lambda=exp(mud + lamd*w[i]))
    
    y[i] <- qnorm(ui[2], mean=mud + lamc*w[i], sd=sqrt(ssq))
  }
  
  DAT.Ga[[j]] <- list(y=y, q=q, X=X, z=z, w=w)
}

save(DAT.Ga, file="Gauss-FOE.Rdata")




color_scheme_set("viridis")
Gauss.dat <- data.frame(y, q, tau=tau.Ga)

Gaussian.orig <- ggplot(data = Gauss.dat, mapping = aes(x = y, y = q)) +
  # Layer 1: Points filled and colored by `tau`
  geom_point(aes(color = tau, fill = tau), size = 3, stroke = 0) +
  # Layer 2: Red outline for "Testing" points with `tau` fill
  scale_color_viridis(discrete = FALSE, end = 1, name = "Kendall's tau",limits=c(-1,1)) +
  scale_fill_viridis(discrete = FALSE, end = 1, name = "Kendall's tau",limits=c(-1,1)) +
  scale_shape_manual(values = c("Training" = 24, "Testing" = 21)) +
  labs(x = "y", y = "q", shape = "Observation Type", fill = "Observation Type") +
  theme(legend.position = c(0.15, 0.7), 
        axis.title = element_text(size = 24, family = "CMU Serif", face = 'italic'), 
        title = element_text(size = 24, family = "CMU Serif"), 
        legend.text = element_text(family = "CMU Serif", size=14),
        legend.title = element_text(size = 16)) +
  ggtitle("Gaussian copula")

Gaussian.orig

ggsave("sim_Gaussian_dat.pdf", plot=Gaussian.orig, units="in", device=cairo_pdf, width=10, height=7)


### See what the MLE is

logpostunorm <- function(ALL) {
  bbeta <- ALL[1:3]
  mmuc <- ALL[4]
  llamc <- ALL[5]
  ssq <- ALL[6]
  xxi <- ALL[7:8]
  
  if( FALSE )  {
    return(-.Machine$double.xmax)
  } else
  {
    l <- sum(dnorm(z, log=TRUE)) +
      sum(dnorm(y, mean=mmuc + w*llamc, sd=sqrt(ssq), log=TRUE)) +
      sum(log_dGauss(cbind(pnorm(z), pnorm(y, mean=mmuc + w*llamc, sd=sqrt(ssq))), thetac_Gauss(w, xxi))) +
      -0.5*sum( (w - as.vector(X %*% bbeta))^2 ) +
      dmvnorm(bbeta, mean=c(0,0,0), sigma=diag(rep(1000, times=3)), log=TRUE) +
      dnorm(llamc, mean=0, sd=sqrt(100), log=TRUE) +
      # dnorm(llamd, mean=0, sd=sqrt(100), log=TRUE) +
      invgamma::dinvgamma(ssq, shape=1000, rate=1000, log=TRUE) +
      dmvnorm(xxi, mean=c(0,0), sigma=diag(rep(500, times=2)), log=TRUE) +
      dnorm(mmuc, mean=0, sd=sqrt(100), log=TRUE)
    1
  }
  return(l)
}

init <- optim(par=c(0,0,0,
                    # 0,
                    # lamd,
                    0,
                    0,
                    0.5,
                    -2, 2), 
              fn=logpostunorm, control = list(fnscale=-1))

data.frame(name=c("beta0", "beta1", "beta2", "mu_c", "lam_c", "ssq", "xi0", "xi1"),
           MLE=init$par, true=c(beta, muc, lamc, ssq, xi))


logpostunorm_xi <- function(xxi) {
  bbeta <- beta
  mmuc <- muc
  llamc <- lamc
  ssq <- ssq
  
  if( FALSE )  {
    return(-.Machine$double.xmax)
  } else
  {
    l <- sum(dnorm(z, log=TRUE)) +
      sum(dnorm(y, mean=mmuc + w*llamc, sd=sqrt(ssq), log=TRUE)) +
      sum(log_dGauss(cbind(pnorm(z), pnorm(y, mean=mmuc + w*llamc, sd=sqrt(ssq))), thetac_Gauss(w, xxi))) +
      -0.5*sum( (w - as.vector(X %*% bbeta))^2 ) +
      dmvnorm(bbeta, mean=c(0,0,0), sigma=diag(rep(1000, times=3)), log=TRUE) +
      dnorm(llamc, mean=0, sd=sqrt(100), log=TRUE) +
      # dnorm(llamd, mean=0, sd=sqrt(100), log=TRUE) +
      invgamma::dinvgamma(ssq, shape=1000, rate=1000, log=TRUE) +
      # dmvnorm(xxi, mean=c(0,0), sigma=diag(rep(500, times=2)), log=TRUE) +
      dnorm(mmuc, mean=0, sd=sqrt(100), log=TRUE)
    1
  }
  return(l)
}

#### Log-likelihood surface
{
  # Create a sequence of values for xi0 and xi1
  xi0_vals <- seq(xi0-1, xi0+1, length.out = 100)
  xi1_vals <- seq(xi1-1, xi1+1, length.out = 100)
  
  # Create a grid of xi0 and xi1 values
  grid <- expand.grid(xi0 = xi0_vals, xi1 = xi1_vals)
  
  # Assuming logpostunorm_xi is your function that takes a vector xi = c(xi0, xi1)
  log_likelihood <- apply(grid, 1, function(x) logpostunorm_xi(c(x[1], x[2])))
  
  # Add the log-likelihood values to the grid
  grid$log_likelihood <- log_likelihood
  
  # Reshape the data for persp
  log_likelihood_matrix <- matrix(log_likelihood, nrow = length(xi0_vals), ncol = length(xi1_vals))
  
  # Create the 3D perspective plot and store the transformation matrix
  persp_matrix <- persp(xi0_vals, xi1_vals, log_likelihood_matrix,
                        theta = 30, phi = 30, expand = 0.5, col = "lightblue",
                        xlab = expression(xi[0]), ylab = expression(xi[1]), zlab = "Log-Likelihood",
                        xlim = range(xi0_vals), ylim = range(xi1_vals), zlim = range(log_likelihood_matrix),
                        ticktype = "detailed",  # Adds numbers to the axes
                        main = "Log-Likelihood Surface - Gauss")
  
  # First line at xi0 = -1.5, varying xi1
  xi0_line1 <- rep(xi0, length(xi1_vals))  # xi0 is constant at -1.5
  xi1_line1 <- xi1_vals  # xi1 ranges over its values
  z_line1 <- apply(cbind(xi0_line1, xi1_line1), 1, function(x) logpostunorm_xi(c(x[1], x[2])))
  
  # Transform the 3D points into 2D for the first line
  line1_points <- trans3d(xi0_line1, xi1_line1, z_line1, persp_matrix)
  
  # Add the first solid red line (xi0 = -1.5)
  lines(line1_points, col = "red", lwd = 2)
  
  # Second line at xi1 = 2, varying xi0
  xi0_line2 <- xi0_vals  # xi0 ranges over its values
  xi1_line2 <- rep(xi1, length(xi0_vals))  # xi1 is constant at 2
  z_line2 <- apply(cbind(xi0_line2, xi1_line2), 1, function(x) logpostunorm_xi(c(x[1], x[2])))
  
  # Transform the 3D points into 2D for the second line
  line2_points <- trans3d(xi0_line2, xi1_line2, z_line2, persp_matrix)
  
  # Add the second solid red line (xi1 = 2)
  lines(line2_points, col = "red", lwd = 2)
}
####

