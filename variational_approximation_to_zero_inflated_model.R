# variational_approximation_to_zero_inflated_model.R
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
# Initialise
lambda = mean(x)
rho = sum(x == 0)/length(x)
eta = rep(NA, n)
r = rep(NA, n)
# Iterate ----
logit = function(p) log(p/(1-p))
iterations = 10000
for (j in 1:iterations) {
  # TODO: You could vectorise the loop below.
  for (i in 1:n) {
    zero_ind = as.numeric(x[i]==0)
    nonzero_ind = as.numeric(x[i]!=0)
    arg = exp(-lambda + logit(rho))
    eta[i] = arg/(zero_ind + arg)
    if (x[i] == 0)
      r[i] = rbinom(1, 1, eta[i])
    else
      r[i] = 1
  }
  # head(cbind(x, eta, r), 100)
  lambda = rgamma(1, a + sum(x), b + sum(r))
  rho = rbeta(1, sum(r) + 1, n - sum(r) + 1)
}

# Variational approximation ----
# Initialise ----
p = rep(.5, n)
# Iterate ----
for (i in 1:iterations) {
  a_lambda = a + sum(x)
  b_lambda = b + sum(p)
  a_rho = 1 + sum(p)
  b_rho = n - sum(p) + 1
  # TODO: You could vectorise the loop below.
  for (i in 1:n) {
    p[i] = - a_lambda/b_lambda + digamma(a_rho) - digamma(b_rho)
  }
# TODO: Lower bound? ----
}
# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
