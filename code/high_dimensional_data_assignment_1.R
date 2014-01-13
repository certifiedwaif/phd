# High dimensional data Assignment 1 Problem 2.R
# Simulate data from a linear regression model with p_alpha_f = 5 and
# n = 16 and n = 64.
require(leaps)

# Generate the test data
generate_test_data = function(variables=1:5, n=16)
{
  # True model:
  beta = c(0.1, 5, 9, 7, 12)
  X = matrix(NA, n, 5)
  for (j in 1:5) {
    X[1:n, j] = runif(n, 0, 1)
  }
  
  if (length(variables) == 1) {
    y = X[,variables] * beta[variables] + rnorm(n, 0, 10)    
  } else {
    y = X[,variables] %*% beta[variables] + rnorm(n, 0, 10)
  }
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

# Slide 12, lecture 6

# ASM : function to select all possible submodels
# input : k.s = number of slope parameters
# intcpt = TRUE, hence all models include an intercept
# output: indicator matrix where rows represent models and columns variables
ASM = function(k.s,intcpt = TRUE){
  if (intcpt == FALSE){k.s = k.s+1}
  A = NULL
  for(i in 1:k.s){
    A0 = cbind(0L,A)
    A1 = cbind(1L,A)
    A = rbind(A0,A1)
  } # end i-for-loop
  if (intcpt == TRUE){A = cbind(1L,A)}
  return(A[order(apply(A,1,sum)),])
} # end ASM function

#A = ASM(dim(dat2)[2]-1) # with ASM function from Lecture 3
# input : A from ASM; y; X with first column a column of 1’s
lm.loglik = function(A,y,X) {
  RES = NULL;
  nr.a = dim(A)[1];
  k =dim(A)[2];
  n = length(y);
  A.c = A %*% diag(1:k)
  for(a in 1:nr.a){
    a.c = A.c[a,]; a.c = a.c[a.c>0];
    X2 = as.matrix(X[,a.c])
    lm.fit.a = lm(y ~ -1 + X2)
    RES = rbind(RES,logLik(lm.fit.a))
  } # a-for end
  return(RES)
}

lm.rss = function(A,y,X) {
  RES = NULL;
  nr.a = dim(A)[1];
  k =dim(A)[2];
  n = length(y);
  A.c = A %*% diag(1:k)
  for(a in 1:nr.a){
    a.c = A.c[a,]; a.c = a.c[a.c>0];
    X2 = as.matrix(X[,a.c])
    lm.fit.a = lm(y ~ -1 + X2)
    RES = rbind(RES,sum(residuals(lm.fit.a)^2))
  } # a-for end
  return(RES)
}

lm.df = function(A,y,X) {
  RES = NULL;
  nr.a = dim(A)[1];
  k =dim(A)[2];
  n = length(y);
  A.c = A %*% diag(1:k)
  for(a in 1:nr.a){
    a.c = A.c[a,]; a.c = a.c[a.c>0];
    X2 = as.matrix(X[,a.c])
    lm.fit.a = lm(y ~ -1 + X2)
    RES = rbind(RES,summary(lm.fit.a)$df[1])
  } # a-for end
  return(RES)
}


#y = dat2$Bodyfat
#X = cbind(1L,dat2[,-1])
#Loss = -2*lm.loglik(A,y,X)
#Size = apply(A,1,sum)

pick_best_model = function(X, y, vars=4, true_model)
{
  # leaps will find the best subset using an exhaustive search
  # Use the regsubsets function
  subsets = regsubsets(x=X, y=y, nvmax=vars, nbest=choose(5, vars), method="exhaustive")
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

pick_best_model2 = function(X, y, k=4, true_model)
{
  n = length(y)
  A = ASM(5)
  rss = lm.rss(A, y, X)
  sigma = function(a) sqrt(rss/(n - a))
  df = lm.df(A, y, X)
  models_we_want = apply(A, 1, sum) %in% 1 + c(1, 3, 4)
  rss = rss[models_we_want]
  df = df[models_we_want]
  print(l_n)
  print(rss)
  print(df)
  print(A[models_we_want,])
  
  # Calculate the measures.
  p = df - 1
  aic_v = 2 * n * log(sigma(0)) + .5 * n  + 2 * df  
  bic_v = 2 * n * log(sigma(0)) + .5 * n  + log(n) * df  
  aic_c_v = 2 * n * log(sigma(0)) + .5 * n + (2*(p + 1) * n)/(n - p - 2)
  aic_cstar_v = 2 * n * log(sigma(p + 2)) + .5 * (n - (p + 2)) + 2 * (p + 1)
  
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
run = function(n=16, vars=4)
{
  true_model = sort(sample(1:5, vars))
  aic_models = NULL
  bic_models = NULL
  aic_c_models = NULL
  aic_cstar_models = NULL
  for (B in 1:100) {
    # Resample
    test_data = generate_test_data(variables=true_model, n=16)
    X = test_data$X
    y = test_data$y
    results = pick_best_model2(X, y, vars, true_model)
    aic_models[B] = results$best_aic_model
    bic_models[B] = results$best_bic_model
    aic_c_models[B] = results$best_aic_c_model
    aic_cstar_models[B] = results$best_aic_cstar_model
  }
  print(true_model)
  return(list(aic_models=aic_models,
              bic_models=bic_models,
              aic_c_models=aic_c_models,
              aic_cstar_models=aic_cstar_models))
}
run(n=16, vars=1)
run(n=64, vars=1)
run(n=16, vars=3)
run(n=64, vars=3)
run(n=16, vars=4)
run(n=64, vars=4)
