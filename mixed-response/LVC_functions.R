eps <- .Machine$double.xmin

thetac_Gaussian <- function(w, x) {
  
  r <- tanh(gam * (x[1] + x[2] * w)) # This is already vectorized
  max_r <- 1 - eps
  r <- pmax(pmin(r, max_r), -max_r) # Vectorized clamping
  return(r)
}

thetac_Frank <- function(w, x) {
  gam * (x[1] + x[2] * w)
}

thetac_Gumbel <- function(w, x) {
  exp(gam * (x[1] + x[2] * w)) + 1
}

thetac_Clayton <- function(w, x) {
  exp(gam * (x[1] + x[2] * w))
}


# regularized incomplete gamma function
Q <- function(j, y) {zipfR::Rgamma(j, y, lower=FALSE)}
Qinv <- function(j, y) {zipfR::Rgamma.inv(j, y, lower=FALSE)}





##### LOGGED JOINT DENSITY FUNCTIONS


log_dGaussian <- function(u, rho) {

  # Ensure u is strictly between 0 and 1
  u[,1] <- pmin(u[,1], 1 - eps)
  u[,2] <- pmin(u[,2], 1 - eps)
  
  # Transform to z-scores
  z1 <- qnorm(u[,1])
  z2 <- qnorm(u[,2])
  
  # Clamp rho for numerical stability
  rho <- pmax(pmin(rho, 1 - eps), -1 + eps)
  
  # Compute log-density
  d <- -0.5 * log(1 - rho^2) - (z1^2 - 2 * rho * z1 * z2 + z2^2) / (2 * (1 - rho^2))
  
  # Handle NaN and Inf values
  d[!is.finite(d)] <- -Inf
  
  return(d)
}


log_dIndep <- function(u) {
  
  # Ensure u is strictly between 0 and 1
  u[,1] <- pmin(u[,1], 1 - eps)
  u[,2] <- pmin(u[,2], 1 - eps)
  
  # Compute log-density
  d <- log(u[,1]) + log(u[,2])
  
  # Handle NaN and Inf values
  d[!is.finite(d)] <- -Inf
  
  return(d)
}


log_dFrank <- function(u, theta) {
  
  # Ensure u is strictly between 0 and 1
  u[,1] <- pmin(u[,1], 1 - eps)
  u[,2] <- pmin(u[,2], 1 - eps)
  
  # Helper function g
  g <- function(x) {expm1(-theta*x)}
  
  # Compute the numerator and denominator
  numerator <- -theta * g(1) * (1 + g(u[,1] + u[,2]))
  denominator <- abs( g(u[,1]) * g(u[,2]) + g(1) )
  
  # Compute log-density
  d <- log(numerator) - 2*log(denominator)
  
  # Handle NaN and Inf values
  d[!is.finite(d)] <- -Inf
  
  return(d)
}


log_dGumbel <- function(u, theta) {
  
  # Ensure u is strictly between 0 and 1
  u[,1] <- pmin(u[,1], 1 - eps)
  u[,2] <- pmin(u[,2], 1 - eps)
  
  # Compute the negative log of the inputs
  mlog_u1 <- -log(u[,1])
  mlog_u2 <- -log(u[,2])
  
  # Compute terms
  s <- cbind((-log(u[, 1]))^theta, (-log(u[, 2]))^theta)
  sum_s <- s[, 1] + s[, 2]
  C <- -(sum_s)^(1/theta)
  
  # Compute log-density
  d <- C - log(u[, 1]) - log(u[, 2]) +
    (-2 + 2/theta) * log(sum_s) +
    (theta - 1) * log(log(u[, 1]) * log(u[, 2])) +
    log(1 + (theta - 1) * sum_s^(-1/theta))
  
  # Handle NaN and Inf values
  d[!is.finite(d)] <- -Inf
  
  return(d)
}


log_dClayton <- function(u, theta) {
  
  # Ensure u is strictly between 0 and 1
  u[,1] <- pmin(u[,1], 1 - eps)
  u[,2] <- pmin(u[,2], 1 - eps)
  
  # Compute the negative log of the inputs
  mlog_u1 <- -log(u[,1])
  mlog_u2 <- -log(u[,2])
  
  # Compute terms
  A <- u[,1]^(-theta) + u[,2]^(-theta) - 1

  # Compute log-density
  d <- log(1 + theta) + (-1 - theta) * log(u[,1] * u[,2]) - (1/theta + 2) * log(A)
  
  # Handle NaN and Inf values
  d[!is.finite(d)] <- -Inf
  
  return(d)
}







##### CONDITIONAL DISTRIBUTION FUNCTIONS


P_cGaussian <- function(u1, u2, theta) {
  
  # Ensure inputs are within valid range
  if (any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1)) {
    stop("u1 and u2 must be in the interval [0, 1]")
  }
  
  # Initialize result vector with zeros (same length as u1)
  result <- numeric(length(u1))
  
  # Identify indices where either u1 or u2 is zero
  zero_mask <- (u1 == 0 | u2 == 0)
  
  # Identify indices where u1 or u2 is 1
  u1_one_mask <- (u1 == 1)
  u2_one_mask <- (u2 == 1)
  
  # Assign values based on boundary conditions
  result[zero_mask] <- 0
  result[u1_one_mask & !u2_one_mask] <- 1
  result[u2_one_mask & !u1_one_mask] <- u1[u2_one_mask & !u1_one_mask]
  
  # For remaining cases, compute the conditional Gaussian copula
  regular_mask <- !(zero_mask | u1_one_mask | u2_one_mask)
  if (any(regular_mask)) {
    z2 <- qnorm(u2[regular_mask])
    z1 <- qnorm(u1[regular_mask])
    
    # Compute conditional mean and variance
    mu_cond <- theta[regular_mask] * z2
    sigma_cond <- sqrt(1 - theta[regular_mask]^2)
    
    # Compute the conditional CDF
    result[regular_mask] <- pnorm(z1, mean = mu_cond, sd = sigma_cond)
  }
  
  # Return the result vector
  return(result)
}



P_cIndep <- function(u1, u2) {
  
  # Ensure inputs are within valid range
  if (any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1)) {
    stop("u1 and u2 must be in the interval [0, 1]")
  }
  
  # Initialize result vector with zeros (same length as u1)
  result <- numeric(length(u1))
  
  # Identify indices where either u1 or u2 is zero
  zero_mask <- (u1 == 0 | u2 == 0)
  
  # Identify indices where u1 or u2 is 1
  u1_one_mask <- (u1 == 1)
  u2_one_mask <- (u2 == 1)
  
  # Assign values based on boundary conditions
  result[zero_mask] <- 0
  result[u1_one_mask & !u2_one_mask] <- 1
  result[u2_one_mask & !u1_one_mask] <- u1[u2_one_mask & !u1_one_mask]
  
  # For remaining cases, compute the conditional Gaussian copula
  regular_mask <- !(zero_mask | u1_one_mask | u2_one_mask)
  if (any(regular_mask)) {
    
    # Compute the conditional CDF
    result[regular_mask] <- u1[regular_mask]
  }
  
  # Return the result vector
  return(result)
}



P_cFrank <- function(u1, u2, theta) {
  
  # Ensure inputs are within valid range
  if (any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1)) {
    stop("u1 and u2 must be in the interval [0, 1]")
  }
  
  # Initialize result vector with zeros (same length as u1)
  result <- numeric(length(u1))
  
  # Identify indices where either u1 or u2 is zero
  zero_mask <- (u1 == 0 | u2 == 0)
  
  # Identify indices where u1 or u2 is 1
  u1_one_mask <- (u1 == 1)
  u2_one_mask <- (u2 == 1)
  
  # Assign values based on boundary conditions
  result[zero_mask] <- 0
  result[u1_one_mask & !u2_one_mask] <- 1
  result[u2_one_mask & !u1_one_mask] <- u1[u2_one_mask & !u1_one_mask]
  
  # Helper function g
  g <- function(x) {expm1(-theta*x)}
  
  # For remaining cases, compute the conditional Gaussian copula
  regular_mask <- !(zero_mask | u1_one_mask | u2_one_mask)
  if (any(regular_mask)) {
    num <- g(u1) * g(u2) + g(u1)
    denom <- g(u1) * g(u2) + g(1)
    
    # Compute the conditional CDF
    result[regular_mask] <- num[regular_mask]/denom[regular_mask]
  }
  
  # Return the result vector
  return(result)
}



P_cGumbel <- function(u1, u2, theta) {
  
  # Ensure inputs are within valid range
  if (any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1)) {
    stop("u1 and u2 must be in the interval [0, 1]")
  }
  
  # Initialize result vector with zeros (same length as u1)
  result <- numeric(length(u1))
  
  # Identify indices where either u1 or u2 is zero
  zero_mask <- (u1 == 0 | u2 == 0)
  
  # Identify indices where u1 or u2 is 1
  u1_one_mask <- (u1 == 1)
  u2_one_mask <- (u2 == 1)
  
  # Assign values based on boundary conditions
  result[zero_mask] <- 0
  result[u1_one_mask & !u2_one_mask] <- 1
  result[u2_one_mask & !u1_one_mask] <- u1[u2_one_mask & !u1_one_mask]
  
  # Helper function g
  g <- function(x) {expm1(-theta*x)}
  
  # For remaining cases, compute the conditional Gaussian copula
  regular_mask <- !(zero_mask | u1_one_mask | u2_one_mask)
  if (any(regular_mask)) {
    
    # Negative log of u1 and u2
    mlogu1 <- -log(u1)
    mlogu2 <- -log(u2)
    
    C <- exp( -( mlogu1^theta + mlogu2^theta )^(1/theta)  )
    num <- (mlogu1^theta + mlogu2^theta)^(1/theta - 1)
    denom <- theta*mlogu2^(1-theta)
    
    # Compute the conditional CDF
    result[regular_mask] <- C[regular_mask]*num[regular_mask]/denom[regular_mask]
  }
  
  # Return the result vector
  return(result)
}


P_cClayton <- function(u1, u2, theta) {
  
  # Ensure inputs are within valid range
  if (any(u1 < 0 | u1 > 1) || any(u2 < 0 | u2 > 1)) {
    stop("u1 and u2 must be in the interval [0, 1]")
  }
  
  # Initialize result vector with zeros (same length as u1)
  result <- numeric(length(u1))
  
  # Identify indices where either u1 or u2 is zero
  zero_mask <- (u1 == 0 | u2 == 0)
  
  # Identify indices where u1 or u2 is 1
  u1_one_mask <- (u1 == 1)
  u2_one_mask <- (u2 == 1)
  
  # Assign values based on boundary conditions
  result[zero_mask] <- 0
  result[u1_one_mask & !u2_one_mask] <- 1
  result[u2_one_mask & !u1_one_mask] <- u1[u2_one_mask & !u1_one_mask]
  
  # Helper function g
  g <- function(x) {expm1(-theta*x)}
  
  # For remaining cases, compute the conditional Gaussian copula
  regular_mask <- !(zero_mask | u1_one_mask | u2_one_mask)
  if (any(regular_mask)) {
    
    # Compute the conditional CDF
    cC <- u2^(-theta - 1) * (u1^(-theta) + u2^(-theta) - 1)^(-(theta + 1)/theta)
    
    result[regular_mask] <- cC[regular_mask]
  }
  
  # Return the result vector
  return(result)
}
