# R implementation for tracking WAIC, CVML, and their conditional versions
# during an MCMC algorithm for the bivariate copula model with normal margins

# Load external functions from LVC_functions.R
source("LVC_functions.R")

# Load required libraries
library(VineCopula)
library(ggplot2)
library(gridExtra)

# Function to generate bivariate data from a specified copula family
# Args:
# - n: Number of observations
# - copula_family: Copula family ("Gumbel", "Frank", "Gaussian", "Clayton", or "Independence")
# - tau: Desired Kendall's tau
# Returns: Matrix of bivariate samples
generate_copula_data <- function(n, copula_family, tau) {
  if (copula_family %in% c("Gumbel", "Clayton") && (tau <= 0 || tau >= 1)) {
    stop("For Gumbel and Clayton copulas, Kendall's tau must be in (0, 1).")
  }
  if (copula_family == "Gaussian" && (tau <= -1 || tau >= 1)) {
    stop("For Gaussian copula, Kendall's tau must be in (-1, 1).")
  }
  if (copula_family == "Frank" && (tau <= -1 || tau >= 1)) {
    stop("For Frank copula, Kendall's tau must be in (-1, 1).")
  }
  
  theta <- solve_copula_param(copula_family, tau)
  
  copula <- switch(copula_family,
                   "Gumbel" = VineCopula::BiCop(family = 4, par = theta),
                   "Frank" = VineCopula::BiCop(family = 5, par = theta),
                   "Gaussian" = VineCopula::BiCop(family = 1, par = theta),
                   "Clayton" = VineCopula::BiCop(family = 3, par = theta),
                   "Independence" = VineCopula::BiCop(family = 0)
  )
  
  if (is.null(copula)) {
    stop("Unsupported copula family: ", copula_family)
  }
  
  u <- VineCopula::BiCopSim(n, copula)
  x1 <- qnorm(u[, 1])
  x2 <- qnorm(u[, 2])
  samples <- cbind(x1, x2)
  colnames(samples) <- c("X1", "X2")
  return(samples)
}

# Initialization: Setup for online updates
initialize_updates <- function(n, B) {
  list(
    mean_p = rep(0, n),     # Mean of P(y_i | omega), initialized with small positive values
    mean_log_p = rep(log(1e-1), n),   # rep(log(1e-15), n),  Mean of log(P(y_i | omega)), initialized with small positive log values
    var_log_p = rep(0, n),  # Variance of log(P(y_i | omega))
    count = B               # Iteration count for Welford updates
  )
}

# Welford's online update function
update_online <- function(state, p_i, log_p_i) {
  t <- state$count + 1
  delta_m <- p_i - state$mean_p
  state$mean_p <- state$mean_p + delta_m / t
  delta_l <- log_p_i - state$mean_log_p
  state$mean_log_p <- state$mean_log_p + delta_l / t
  state$var_log_p <- state$var_log_p + delta_l * (log_p_i - state$mean_log_p)
  state$count <- t
  return(state)
}

# Final computation of WAIC, CVML, and their conditional versions
compute_measures <- function(state, M, B) {
  n <- length(state$mean_p)
  waic <- -2 * sum(log(pmax(state$mean_p, 1e-15))) + 2 * sum(state$var_log_p / (M - B + 1))
  return(waic)
}

# Function to solve for copula parameters given Kendall's tau using VineCopula
solve_copula_param <- function(copula_family, tau) {
  if (copula_family == "Gumbel") {
    return(1 / (1 - tau))  # Gumbel: theta = 1 / (1 - tau)
  } else if (copula_family == "Frank") {
    return(VineCopula::BiCopTau2Par(5, tau))  # Frank: use VineCopula for calculation
  } else if (copula_family == "Gaussian") {
    return(sin(tau * pi / 2))  # Gaussian: theta = sin(tau * pi / 2)
  } else if (copula_family == "Clayton") {
    return(2 * tau / (1 - tau))  # Clayton: theta = 2 * tau / (1 - tau)
  } else {
    return(0)  # Independence or unsupported family
  }
}

# Vectorized density computation functions using external functions
compute_density <- function(data, params, copula_family) {
  u1 <- pmax(pnorm(data[, 1], mean = params$mu1, sd = params$sigma1), 1e-15)
  u2 <- pmax(pnorm(data[, 2], mean = params$mu2, sd = params$sigma2), 1e-15)
  if (copula_family == "Gumbel") {
    return(exp(log_dGumbel(cbind(u1, u2), params$theta)))
  } else if (copula_family == "Frank") {
    return(exp(log_dFrank(cbind(u1, u2), params$theta)))
  } else if (copula_family == "Gaussian") {
    return(exp(log_dGaussian(cbind(u1, u2), params$rho)))
  } else if (copula_family == "Clayton") {
    return(exp(log_dClayton(cbind(u1, u2), params$theta)))
  } else if (copula_family == "Independence") {
    return(rep(1, length(u1)))
  } else {
    stop("Unsupported copula family")
  }
}

compute_conditional_density <- function(data, params, copula_family, condition_on_y2 = TRUE) {
  u1 <- pmax(pnorm(data[, 1], mean = params$mu1, sd = params$sigma1), 1e-15)
  u2 <- pmax(pnorm(data[, 2], mean = params$mu2, sd = params$sigma2), 1e-15)
  if (copula_family == "Gumbel") {
    log_density <- log_dGumbel(cbind(u1, u2), params$theta)
  } else if (copula_family == "Frank") {
    log_density <- log_dFrank(cbind(u1, u2), params$theta)
  } else if (copula_family == "Gaussian") {
    log_density <- log_dGaussian(cbind(u1, u2), params$rho)
  } else if (copula_family == "Clayton") {
    log_density <- log_dClayton(cbind(u1, u2), params$theta)
  } else if (copula_family == "Independence") {
    log_density <- rep(0, length(u1))  # log(1) = 0
  } else {
    stop("Unsupported copula family")
  }
  if (condition_on_y2) {
    return(exp(log_density) * dnorm(data[, 2], mean = params$mu2, sd = params$sigma2))
  } else {
    return(exp(log_density) * dnorm(data[, 1], mean = params$mu1, sd = params$sigma1))
  }
}

# Main MCMC function with Random Walk Metropolis for sampling
run_mcmc <- function(data, n_iter, burn_in, copula_family) {
  n <- nrow(data)  # Number of data points
  
  # Initialize updates for all measures
  waic_state <- initialize_updates(n, burn_in)
  cvml_state <- initialize_updates(n, burn_in)
  cwaic1_state <- initialize_updates(n, burn_in)
  cwaic2_state <- initialize_updates(n, burn_in)
  ccvml1_state <- initialize_updates(n, burn_in)
  ccvml2_state <- initialize_updates(n, burn_in)
  
  # Initialize parameter sampling
  current_tau <- 0.5  # Initial Kendall's tau
  current_theta <- solve_copula_param(copula_family, current_tau)
  params <- list(theta = current_theta, mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1)
  if (copula_family == "Gaussian") {
    params$rho <- current_theta
  }
  
  # Placeholder for parameter samples
  samples <- numeric(n_iter)
  acceptances <- 0
  
  for (iter in 1:n_iter) {
    # Random Walk Proposal for Kendall's Tau
    proposed_tau <- rnorm(1, mean = current_tau, sd = 0.1)  # Proposal with small random step
    proposed_tau <- if (copula_family %in% c("Gumbel", "Clayton")) {
      pmin(pmax(proposed_tau, 0.01), 0.99)  # Restrict tau to [0, 1) for Gumbel and Clayton
    } else {
      pmin(pmax(proposed_tau, -0.99), 0.99)  # Allow negative dependence for other families
    }
    
    proposed_theta <- solve_copula_param(copula_family, proposed_tau)
    proposed_params <- list(theta = proposed_theta, mu1 = 0, sigma1 = 1, mu2 = 0, sigma2 = 1)
    if (copula_family == "Gaussian") {
      proposed_params$rho <- proposed_theta
    }
    
    # Compute log-likelihood for current and proposed states
    current_log_lik <- sum(log(compute_density(data, params, copula_family)))
    proposed_log_lik <- sum(log(compute_density(data, proposed_params, copula_family)))
    
    # Metropolis acceptance step
    log_acceptance_ratio <- proposed_log_lik - current_log_lik
    if (log(runif(1)) < log_acceptance_ratio) {
      current_tau <- proposed_tau
      current_theta <- proposed_theta
      params <- proposed_params
      acceptances <- acceptances + 1
    }
    
    # Store the current tau value
    samples[iter] <- current_tau
    
    # Step 2: Compute log-densities
    p_i <- compute_density(data, params, copula_family)
    log_p_i <- log(p_i)
    
    p_y1_given_y2 <- compute_conditional_density(data, params, copula_family, condition_on_y2 = TRUE)
    log_p_y1_given_y2 <- log(p_y1_given_y2)
    
    p_y2_given_y1 <- compute_conditional_density(data, params, copula_family, condition_on_y2 = FALSE)
    log_p_y2_given_y1 <- log(p_y2_given_y1)
    
    # Update WAIC, CVML, and conditional measures
    waic_state <- update_online(waic_state, exp(log_p_i), log_p_i)
    cvml_state <- update_online(cvml_state, 1 / p_i, -log_p_i)
    
    cwaic1_state <- update_online(cwaic1_state, exp(log_p_y1_given_y2), log_p_y1_given_y2)
    cwaic2_state <- update_online(cwaic2_state, exp(log_p_y2_given_y1), log_p_y2_given_y1)
    ccvml1_state <- update_online(ccvml1_state, 1 / p_y1_given_y2, -log_p_y1_given_y2)
    ccvml2_state <- update_online(ccvml2_state, 1 / p_y2_given_y1, -log_p_y2_given_y1)
  }
  
  # Final computation of measures after MCMC
  waic <- compute_measures(waic_state, n_iter, burn_in)
  cvml <- compute_measures(cvml_state, n_iter, burn_in)
  cwaic <- list(
    cwaic1 = compute_measures(cwaic1_state, n_iter, burn_in),
    cwaic2 = compute_measures(cwaic2_state, n_iter, burn_in)
  )
  ccvml <- list(
    ccvml1 = compute_measures(ccvml1_state, n_iter, burn_in),
    ccvml2 = compute_measures(ccvml2_state, n_iter, burn_in)
  )
  
  # Calculate acceptance rate
  acceptance_rate <- acceptances / n_iter
  
  # Traceplot of samples after burn-in
  traceplot <- data.frame(Iteration = (burn_in + 1):n_iter, Tau = samples[(burn_in + 1):n_iter], Family = copula_family)
  
  return(list(waic = waic, cvml = cvml, cwaic = cwaic, ccvml = ccvml, samples = samples, 
              acceptance_rate = acceptance_rate, traceplot = traceplot))
}






# Example usage of generate_copula_data
set.seed(23)
data <- generate_copula_data(n = 1000, copula_family = "Gumbel", tau = 0.45)
print(head(data))

n_iter <- 10000

results_gumbel <- run_mcmc(data, n_iter = n_iter, burn_in = n_iter/4, copula_family = "Gumbel")
results_frank <- run_mcmc(data, n_iter = n_iter, burn_in = n_iter/4, copula_family = "Frank")
results_independence <- run_mcmc(data, n_iter = n_iter, burn_in = n_iter/4, copula_family = "Independence")
results_normal <- run_mcmc(data, n_iter = n_iter, burn_in = n_iter/4, copula_family = "Gaussian")
results_clayton <- run_mcmc(data, n_iter = n_iter, burn_in = n_iter/4, copula_family = "Clayton")

# Combine traceplots for all families into a grid
traceplots_list <- list(
  ggplot(results_gumbel$traceplot, aes(x = Iteration, y = Tau)) +
    geom_line(color = "blue") + theme_minimal() + labs(title = "Gumbel", x = "Iteration", y = "Kendall's Tau"),
  ggplot(results_frank$traceplot, aes(x = Iteration, y = Tau)) +
    geom_line(color = "red") + theme_minimal() + labs(title = "Frank", x = "Iteration", y = "Kendall's Tau"),
  ggplot(results_independence$traceplot, aes(x = Iteration, y = Tau)) +
    geom_line(color = "green") + theme_minimal() + labs(title = "Independence", x = "Iteration", y = "Kendall's Tau"),
  ggplot(results_normal$traceplot, aes(x = Iteration, y = Tau)) +
    geom_line(color = "purple") + theme_minimal() + labs(title = "Gaussian", x = "Iteration", y = "Kendall's Tau"),
  ggplot(results_clayton$traceplot, aes(x = Iteration, y = Tau)) +
    geom_line(color = "orange") + theme_minimal() + labs(title = "Clayton", x = "Iteration", y = "Kendall's Tau")
)

combined_traceplot <- gridExtra::grid.arrange(grobs = traceplots_list, ncol = 2)

# Display combined traceplot grid
print(combined_traceplot)

# Compile model results into a table
model_results <- list(
  Gumbel = results_gumbel,
  Frank = results_frank,
  Independence = results_independence,
  Gaussian = results_normal,
  Clayton = results_clayton
)

formatted_results <- lapply(names(model_results), function(name) {
  res <- model_results[[name]]
  c(
    Model = name,
    WAIC = round(res$waic, 1),
    CVML = round(res$cvml, 1),
 #   CWAIC1 = round(res$cwaic$cwaic1, 1),
 #   CWAIC2 = round(res$cwaic$cwaic2, 1),
 #   CCVML1 = round(res$ccvml$ccvml1, 1),
 #   CCVML2 = round(res$ccvml$ccvml2, 1),
    AcceptanceRate = round(res$acceptance_rate, 2)
  )
})

results_table <- do.call(rbind, formatted_results)
print(results_table)
