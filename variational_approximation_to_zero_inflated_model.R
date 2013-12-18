# variational_approximation_to_zero_inflated_model.R
require(parallel)

# I really have to come up with a better way of naming these ----

# Simulate the data ----
rho = .5
lambda = 100
n = 1000
x = rep(NA, n)
for (i in 1:n) {
  if (runif(1) >= rho) {
    x[i] = rpois(1, lambda)
  } else {
    x[i] = 0
  }
}
a = 10
b = 10
# MCMC ----
logit = function(p) log(p/(1-p))
expit = function(eta) 1/(1+exp(-eta))
zero_infl_mcmc = function(iterations)
{
  # Initialise
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
    
    #for (j in 1:length(zero_observation_idx)) {
    #  r[zero_observation_idx[j]] = rbinom(1, 1, eta)
    #}
    r[zero_observation_idx] = rbinom(length(zero_observation_idx), 1, eta) 
    lambda[i+1] = rgamma(1, a + sum(x), b + sum(r))
    rho[i+1] = rbeta(1, sum(r) + 1, n - sum(r) + 1)
  }
  return(list(lambda=lambda, rho=rho))
}

start = Sys.time()
result_mcmc = zero_infl_mcmc(1e3)
Sys.time() - start
# 1000 iterations in 0.07461715 seconds.

# Variational approximation ----
# Initialise ----
zero_infl_var = function(iterations)
{
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
    params = list(a_lambda=a_lambda, b_lambda=b_lambda, a_rho=a_rho, b_rho=b_rho)
    #print(params)
  }
  return(params)
}
start = Sys.time()
result_var = zero_infl_var(10)
Sys.time() - start
# Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
# or 200 iterations a second.
result_var$a_lambda / result_var$b_lambda
result_var$a_rho / (result_var$a_rho + result_var$b_lambda)

# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
