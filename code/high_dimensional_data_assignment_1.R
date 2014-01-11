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
  y = X[,variables] %*% beta[variables]
  return(list(X=X, y=y))
}

sigma_a = function(y, yhat, a)
{
  n = length(y)
  rss = sum((y - yhat)^2)
  return(sqrt(rss/(n - a)))
}

# Probably a better idea to use the log likelihood functions built into
# R.
log_lik = function(X, y, beta, sigma)
{
  residuals = y - X %*% beta
  return(sum(log(pnorm(residuals, 0, sigma^2))))
}

l_n = function(y, X, theta_a)
# theta_a is implicit in the X matrix which is passed.  
{
  # Fit a linear model to y using the model in theta_a. Return the
  # log likelihood of the model.
  data_tmp = as.data.frame(cbind(1, X))
  fit = lm(y~., data=data_tmp)
  return(logLik(fit))
}

# leaps will find the best subset using an exhaustive search
# Use the regsubsets function

# AIC
aic = function(loglik, n, k)
{
  -2 * loglik + 2 * k
}

# BIC
bic = function(loglik, n, k)
{
  -2 * loglik + log(n) * k  
}

# AIC_c
aic_c = function(loglik, n, k)
{
  # AIC with a finite sample correction
  aic(loglik, n, k) + (2*k*(k+1))/(n - k - 1)
}

# AIC_c^*
aic_cstar = function(loglik, n, k)
{
  # ???
}

# Consider data generating models with one, three and four non-zero parameters.
# What exactly does this mean?
# Pick a subset of the five variables. Generate test data using that subset.
# Then enumerate all possible sub-models and generate them. See which
# is selected by the AIC etc.
variables = sort(sample(1:5, 4))
test_data = generate_test_data(variables=variables, n=16)

# Calculate these model selection measures and use them to rank models
# Pick the best one
X = test_data$X
y = test_data$y
subsets = regsubsets(x=X, y=y, nvmax=4, nbest=2^4, method="exhaustive")
names(subsets)
summ_subsets = summary(subsets)
summ_subsets$which
dimnames(summ_subsets$which)
# Extract the models of interest.
models = summ_subsets$which[dimnames(summ_subsets$which)[[1]] == "4",]
models = models[,-1] # Drop the intercept, which is always in
# Calculate the measures.
aic_v = NULL
bic_v = NULL
aic_c_v = NULL
for (i in 1:dim(models)[1]) {
  fit = lm(y~X[,models[i,]])
  theta = coef(fit)
  data_tmp = as.data.frame(cbind(1, X[,models[i,]]))
  fit = lm(y~., data=data_tmp)
  loglik = logLik(fit)
  n = length(y)
  k = 4

  aic_v[i] = aic(loglik, n, k)
  bic_v[i] = bic(loglik, n, k)
  aic_c_v[i] = aic_c(loglik, n, k)
}
# Pick the best one according to our measure.

# Idea: Use regsubsets to calculate many of the models. Then calculate the
# AIC/BIC etc. for all of them, and select the best one yourself.
