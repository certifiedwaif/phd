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
iterations = 10000
for (j in 1:iterations) {
  # TODO: You could vectorise the loop below.
  for (i in 1:n) {
    zero_ind = as.numeric(x[i]==0)
    nonzero_ind = as.numeric(x[i]!=0)
    eta[i] = ((exp(-lambda) * rho))/(zero_ind + ((exp(-lambda) * rho)))
    r[i] = rbinom(1, 1, eta[i])
  }
  # head(cbind(x, eta, r), 100)
  lambda = rgamma(1, a + sum(x), b + sum(r))
  rho = rbeta(1, sum(r) + 1, n - sum(r) + 1)
}

# Variational approximation ----
# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
