# variational_approximation_to_zero_inflated_model.R

generate_test_data = function (n, rho, lambda)
{
  x = rep(NA, n)
  for (i in 1:n) {
    if (runif(1) >= rho) {
      x[i] = rpois(1, lambda)
    } else {
      x[i] = 0
    }
  }
  return(x)
}

logit = function(p) log(p/(1-p))
expit = function(eta) 1/(1+exp(-eta))

zero_infl_mcmc = function(iterations, x, a, b)
{
  # Initialise
  n = length(x)
  lambda = rep(NA, iterations+1)
  rho = rep(NA, iterations+1)
  lambda[1] = mean(x)
  rho[1] = sum(x == 0)/length(x)
  eta = rep(NA, n)
  r = rep(NA, n)
  # Build up a list of the zero observations. We'll only generate
  # r[i]'s for those. The non-zero observations will always have r[i] = 1.
  r[x != 0] = 1
  zero_observation_idx = which(x == 0)
  # Iterate ----
  for (i in 1:iterations) {
    arg = exp(-lambda[i] + logit(rho[i]))
    eta = arg/(1 + arg)
    
    # The following code is correct, but 100 times slower as R is interpreted.
    # So we use the vectorised code instead.
    #for (j in 1:length(zero_observation_idx)) {
    #  r[zero_observation_idx[j]] = rbinom(1, 1, eta)
    #}
    
    r[zero_observation_idx] = rbinom(length(zero_observation_idx), 1, eta) 
    lambda[i+1] = rgamma(1, a + sum(x), b + sum(r))
    rho[i+1] = rbeta(1, sum(r) + 1, n - sum(r) + 1)
  }
  return(list(lambda=lambda, rho=rho))
}

lower_bound = function(x, p, a, b, a_lambda, b_lambda, a_rho, b_rho)
{
  gamma_entropy = function(a, b)
  {
    a - log(b) + lgamma(a) + (1-a) * digamma(a)
  }

  E_log_lambda = -gamma_entropy(a, b)
  E_lambda = a/b
  E_r = ifelse(x == 0, p, 1)
  E_xi_log_lambda_r = ifelse(x == 0, E_log_lambda, 0)

  beta_entropy = function(a, b)
  {
    lbeta(a, b) - (a-1)*digamma(a) - (b-1)*digamma(b)+ (a+b-2)*digamma(a+b)
  }

  E_log_rho = 0 # rho is Unif([0, 1]), so the answer is digamma(1) - digamma(1) = 0
  E_log_one_minus_rho = E_log_rho
  
  log_q_r = ifelse(x == 0, log(p), 0)
  log_q_lambda = -gamma_entropy(a_lambda, b_lambda)
  log_q_rho = -beta_entropy(a_rho, b_rho)
  
  # FIXME: Sign error?
  result = a * log(b) + (a-1) * E_log_lambda - b * E_lambda - lgamma(a)
  result = result + - E_lambda * sum(E_r)
  result = result + sum(E_xi_log_lambda_r) - sum(log(factorial(x)))
  result = result + sum(E_r) * E_log_rho + sum((1 - E_r)) * E_log_one_minus_rho
  result = result - sum(log_q_r) - log_q_lambda - log_q_rho
  
  browser()
  return(result)
}

zero_infl_var = function(iterations, x, a, b)
{
  n = length(x)
  p = rep(sum(x == 0)/length(x), n)
  eta = rep(NA, n)
  a_lambda = a + sum(x)
  # Iterate ----
  for (i in 1:iterations) {
    b_lambda = b + sum(p)
    a_rho = 1 + sum(p)
    b_rho = n - sum(p) + 1
    # TODO: You could vectorise the loop below.
    for (i in 1:n) {
      if (x[i] != 0)
        p[i] = 1.0
      else {
        eta[i] = - a_lambda/b_lambda + digamma(a_rho) - digamma(b_rho)
        p[i] = expit(eta[i])
      }
    }
    # TODO: Lower bound? ----
    lower_bd = lower_bound(x, p, a, b, a_lambda, b_lambda, a_rho, b_rho)
    print(lower_bd)
    params = list(a_lambda=a_lambda, b_lambda=b_lambda, a_rho=a_rho, b_rho=b_rho)
    #print(params)
  }
  return(params)
}

# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
calculate_accuracy = function(result_mcmc, result_var)
{
  density_mcmc_rho = density(result_mcmc$rho)
  integrand = function(x)
  {
    fn = splinefun(density_mcmc_rho$x, density_mcmc_rho$y)
    return(abs(fn(x) - dbeta(x, result_var$a_rho, result_var$b_rho)))
  }
  integrate(integrand, min(density_mcmc_rho$x), max(density_mcmc_rho$x), subdivisions = length(density_mcmc_rho$x))
  
  density_mcmc_lambda = density(result_mcmc$lambda)
  integrand = function(x)
  {
    fn = splinefun(density_mcmc_lambda$x, density_mcmc_lambda$y)
    return(abs(fn(x) - dgamma(x, result_var$a_lambda, result_var$b_lambda)))
  }
  result = integrate(integrand, min(density_mcmc_lambda$x), max(density_mcmc_lambda$x), subdivisions = length(density_mcmc_lambda$x))
  accuracy = 1 - .5 * result$value
  return(accuracy)
}

check_accuracy = function(n, rho, lambda)
{
  x = generate_test_data(n, rho, lambda)
  
  a = 10
  b = 10
  
  # MCMC ----
  iterations = 1e6
  burnin = 1e3
  
  start = Sys.time()
  result_mcmc = zero_infl_mcmc(iterations+burnin, x, a, b)
  mcmc_runtime = Sys.time() - start
  # Throw away burn-in samples.
  # Brackets turn out to be incredibly important here!!!
  result_mcmc$rho = result_mcmc$rho[(burnin+1):(burnin+iterations+1)]
  result_mcmc$lambda = result_mcmc$lambda[(burnin+1):(burnin+iterations+1)]
  # 1000 iterations in 0.07461715 seconds.
  
  # Variational approximation ----
  start = Sys.time()
  result_var = zero_infl_var(10, x, a, b)
  var_runtime = Sys.time() - start
  # Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
  # or 200 iterations a second.
  var_lambda = result_var$a_lambda / result_var$b_lambda
  var_rho = result_var$a_rho / (result_var$a_rho + result_var$b_lambda)
  
  accuracy = calculate_accuracy(result_mcmc, result_var)
  return(list(n=n, rho=rho, lambda=lambda, accuracy=accuracy))
}

# Check accuracy of the approximation for a range of parameter values ----
main_check_accuracy = function()
{
  n = 1000
  rho = .5
  lambda = 100
  
  for (rho in .1*1:9)
    for (lambda in seq(5, 100, 5))
      print(check_accuracy(n, rho, lambda))
}
#main_check_accuracy()

# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC

x = generate_test_data(n, rho, lambda)

a = 10
b = 10

# MCMC ----
iterations = 1e6
burnin = 1e3

start = Sys.time()
result_mcmc = zero_infl_mcmc(iterations+burnin, x, a, b)
mcmc_runtime = Sys.time() - start
# Throw away burn-in samples.
# Brackets turn out to be incredibly important here!!!
result_mcmc$rho = result_mcmc$rho[(burnin+1):(burnin+iterations+1)]
result_mcmc$lambda = result_mcmc$lambda[(burnin+1):(burnin+iterations+1)]
# 1000 iterations in 0.07461715 seconds.

# Variational approximation ----
start = Sys.time()
result_var = zero_infl_var(10, x, a, b)
var_runtime = Sys.time() - start
# Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
# or 200 iterations a second.
var_lambda = result_var$a_lambda / result_var$b_lambda
var_rho = result_var$a_rho / (result_var$a_rho + result_var$b_lambda)
