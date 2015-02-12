# univariate.R

mcmc.univariate <- function(iterations, vx, a, b)
{
  # Initialise
  n = length(vx)
  zero.set = which(vx == 0)
  nonzero.set = which(vx != 0)
  vlambda = rep(NA, iterations+1)
  vrho = rep(NA, iterations+1)
  vlambda[1] = mean(vx)
  vrho[1] = length(zero.set)/length(vx)
  veta = rep(NA, n)
  vr = rep(NA, n)
  # Build up a list of the zero observations. We'll only generate
  # r[i]'s for those. The non-zero observations will always have r[i] = 1.
  vr[nonzero.set] = 1
  # Iterate ----
  for (i in 1:iterations) {
    varg = exp(-vlambda[i] + logit(vrho[i]))
    veta = varg/(1 + varg)
    
    # The following code is correct, but 100 times slower as R is interpreted.
    # So we use the vectorised code instead.
    #for (j in 1:length(zero_observation_idx)) {
    #  r[zero_observation_idx[j]] = rbinom(1, 1, eta)
    #}
    
    vr[zero.set] = rbinom(length(zero.set), 1, veta) 
    vlambda[i+1] = rgamma(1, a + sum(vx), b + sum(vr))
    vrho[i+1] = rbeta(1, sum(vr) + 1, n - sum(vr) + 1)
  }
  return(list(vlambda=vlambda, vrho=vrho))
}

calculate_lower_bound.univariate <- function(univariate)
{
  vx = univariate$vx
  vp = univariate$vp
  a_lambda = univariate$a_lambda
  b_lambda = univariate$b_lambda
  a_rho = univariate$a_rho
  b_rho = univariate$b_rho
  
  zero.set <- which(vx==0)
  
  E_lambda = a_lambda/b_lambda
  E_log_lambda = digamma(a_lambda) - log(b_lambda)  
  
  E_r = ifelse(vx == 0, vp, 1)
  E_xi_log_lambda_r = ifelse(vx == 0, 0, vx*E_log_lambda)
  
  E_log_rho = digamma(a_rho) - digamma(a_rho + b_rho)
  E_log_one_minus_rho = digamma(b_rho) - digamma(a_rho + b_rho)
  
  E_log_q_r = (vp[zero.set]*log(vp[zero.set]) + (1-vp[zero.set])*log(1-vp[zero.set]))
  E_log_q_lambda = -gamma_entropy(a_lambda, b_lambda)
  E_log_q_rho = -beta_entropy(a_rho, b_rho) # FIXME: Why is this entropy positive?
  
  result = a_lambda * log(b_lambda) + (a_lambda-1) * E_log_lambda - b_lambda * E_lambda - lgamma(a_lambda)
  result = result - E_lambda * sum(E_r)
  result = result + sum(E_xi_log_lambda_r) - sum(lgamma(vx+1))
  result = result + sum(E_r) * E_log_rho + sum(1 - E_r) * E_log_one_minus_rho
  result = result - sum(E_log_q_r) - E_log_q_lambda - E_log_q_rho
  
  return(result)
}

create_univariate <- function(vx, a, b)
{
  # Initialise
  n = length(vx)
  vp = rep(1, n)
  a_lambda = a + sum(vx)
  b_lambda = b
  univariate = list(vx=vx, vp=vp, a_lambda=a_lambda, b_lambda=b_lambda)
  class(univariate) = "univariate"
  return(univariate)
}

zero_infl_var.univariate <- function(univariate, verbose=FALSE, plot_lower_bound=FALSE)
{
  vx = univariate$vx
  a = univariate$a_lambda
  b = univariate$b_lambda
  n = length(univariate$vx)
  zero.set = which(univariate$vx == 0)
  nonzero.set = which(univariate$vx != 0)
  vlower_bound <- c()
  
  i = 0
  # Iterate ----
  while (i <= 2 || vlower_bound[i] > vlower_bound[i-1]) {
    i = i+1
    
    # Update parameter for q_lambda
    univariate$b_lambda = b + sum(univariate$vp)
    
    # Update parameters for q_rho
    univariate$a_rho = 1 + sum(univariate$vp)
    univariate$b_rho = n - sum(univariate$vp) + 1
    
    # Update parameters for q_vr
    univariate$vp[zero.set] = expit(-expected_lambda(univariate) + digamma(univariate$a_rho) - digamma(univariate$b_rho))
    
    vlower_bound[i] <- calculate_lower_bound(univariate)
    
    if (verbose && i > 1)
      cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
          vlower_bound[i] - vlower_bound[i-1], " parameters ", "a_lambda", univariate$a_lambda,
          "b_lambda", univariate$b_lambda, "a_rho", univariate$a_rho, "b_rho", univariate$b_rho, "\n")
  }
  
  if (plot_lower_bound)
    plot(lower_bound_vector,type="l")
  
  params = list(a_lambda=univariate$a_lambda, b_lambda=univariate$b_lambda, a_rho=univariate$a_rho, b_rho=univariate$b_rho)
  return(params)
}

