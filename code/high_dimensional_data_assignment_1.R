# High dimensional data Assignment 1 Problem 2.R
# Simulate data from a linear regression model with p_alpha_f = 5 and
# n = 16 and n = 64.

# Generate the test data
generate_test_data = function(variables=1:5, n=16)
{
  # True model:
  beta = c(12, -8, 9, -7, 12)
  X = matrix(NA, n, 5)
  for (j in 1:5) {
    X[1:n, j] = rnorm(n, 0, 1)
  }
  y = X[variables] %*% beta[variables]
  return(list(X=X, y=y))
}
# Consider data generating models with one, three and four non-zero parameters.
# What exactly does this mean?
# Pick a subset of the five variables. Generate test data using that subset.
# Then enumerate all possible sub-models and generate them. See which
# is selected by the AIC etc.

variables = sort(sample(1:5, 4))
test_data = generate_test_data(variables=variables, n=16)

sigma_a = function(y, yhat, a)
{
  n = length(y)
  rss = sum((y - yhat)^2)
  return(sqrt(rss/(n - a)))
}

log_lik = function(X, y, beta, sigma)
{
  residuals = y - X %*% beta
  return(sum(log(pnorm(residuals, 0, sigma^2))))
}

# AIC
# BIC
# AIC_c
# AIC_c^*

# Calculate these model selection measures and use them to rank models
# Pick the best one
