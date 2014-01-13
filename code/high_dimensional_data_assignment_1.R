# High dimensional data Assignment 1 Problem 2.R
# Simulate data from a linear regression model with p_alpha_f = 5 and
# n = 16 and n = 64.
require(leaps)

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

sigma = function(RSS, n, a)
{
  sqrt(RSS/(n-a))
}

l_n = function(y, X, beta, sigma)
{
  #print(beta)
  k = length(beta)
  - .5*k * log(2*pi) - k * log(sigma) - .5 * sum((y - X %*% beta)^2)
}

pick_best_model = function(X, y, vars=4, true_model)
{
  # leaps will find the best subset using an exhaustive search
  # Use the regsubsets function
  subsets = regsubsets(x=X, y=y, nvmax=vars, nbest=2^vars, method="exhaustive")
  names(subsets)
  summ_subsets = summary(subsets)
  summ_subsets$which
  dimnames(summ_subsets$which)
  # TODO: Return true model
  # true_model = some subset, which one?
  
  # Extract the models of interest.
  models = summ_subsets$which[dimnames(summ_subsets$which)[[1]] == as.character(vars),]
  models = models[,-1] # Drop the intercept, which is always in
  # Calculate the measures.
  aic_v = NULL
  bic_v = NULL
  aic_c_v = NULL
  aic_cstar_v = NULL
  for (i in 1:dim(models)[1]) {
    # We don't care about models that don't include the intercept
    if (models[i,1] == FALSE)
      next
    X_subset = cbind(1, X[,models[i,]])
    fit = lm(y~., data=as.data.frame(X_subset))
    #beta_ML = coef(fit)[models[i,]]
    beta_ML = coef(fit)#[models[i,]]
    print(beta_ML)
    beta_ML = beta_ML[!is.na(beta_ML)]
    print(beta_ML)
    print(models[i,])
    RSS = sum(residuals(fit)^2)
    #loglik = logLik(fit)
    # TODO: Calculate the scale estimator
    # TODO: Calculate the log likelihood
    n = length(y)
    k = vars
  
    # AIC: -2 * loglik + 2 * k
    # BIC: -2 * loglik + log(n) * k
    # AIC_c: -2 l_n(beta_ML, sigma_0) + 2(p_alpha + 1) (n)/(n - p_alpha - 2)
    # AIC_c*: -2 l_n(beta_ML, sigma_{p_alpha + 2}) + 2 (p_alpha + 1)
    #print(beta_ML)
    print(RSS)
    aic_v[i] = -2 * l_n(y, X_subset, beta_ML, sigma(RSS, n, k+1)) + 2 * k
    bic_v[i] = -2 * l_n(y, X_subset, beta_ML, sigma(RSS, n, k+1)) + log(n) * k
    aic_c_v[i] = -2 * l_n(y, X_subset, beta_ML, sigma(RSS, n, 0)) + (2*(k + 1) * n)/(n - k - 2)
    aic_cstar_v[i] = -2 * l_n(y, X_subset, beta_ML, sigma(RSS, n, k + 2)) + 2 * (k + 1)
  }
  print(aic_v)
  print(bic_v)
  print(aic_c_v)
  print(aic_cstar_v)
  # Pick the best one according to our measure.
  best_aic_model = which.min(aic_v)
  best_bic_model = which.min(bic_v)
  best_aic_c_model = which.min(aic_c_v)
  best_aic_cstar_model = which.min(aic_cstar_v)
  
  # Did you pick the best one?
  # In my test case, we always picked the same one. Is something wrong?
  
  # TODO: What should I return here?
  return(list(best_aic_model=best_aic_model,
              best_bic_model=best_bic_model,
              best_aic_c_model=best_aic_model,
              best_aic_cstar_model=best_aic_cstar_model))
}

# Consider data generating models with one, three and four non-zero parameters.
# What exactly does this mean?
# Pick a subset of the five variables. Generate test data using that subset.
# Then enumerate all possible sub-models and generate them. See which
# is selected by the AIC etc.

# Idea: Use regsubsets to calculate many of the models. Then calculate the
# AIC/BIC etc. for all of them, and select the best one yourself.

# Calculate these model selection measures and use them to rank models
# Pick the best one

# TODO: One, three and four non-zero regression parameters.
# Bootstrap of 1000 runs.
true_model = sort(sample(1:5, 4))
aic_models = NULL
bic_models = NULL
aic_c_models = NULL
aic_cstar_models = NULL
for (B in 1:1000) {
  # Resample
  test_data = generate_test_data(variables=true_model, n=16)
  X = test_data$X
  y = test_data$y
  results = pick_best_model(X, y, 4, true_model)
  aic_models[B] = results$best_aic_model
  bic_models[B] = results$best_bic_model
  aic_c_models[B] = results$best_aic_c_model
  aic_cstar_models[B] = results$best_aic_cstar_model
}
