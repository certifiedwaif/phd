# High dimensional data Assignment 1 Problem 2.R
# Simulate data from a linear regression model with p_alpha_f = 5 and
# n = 16 and n = 64.
require(leaps)

# Generate the test data
generate_test_data = function(variables=1:5, n=16, sigma=1)
{
  # True model:
  beta = c(0.1, 5, 9, 7, 12)
  X = matrix(NA, n, 5)
  for (j in 1:5) {
    X[1:n, j] = runif(n, 0, 1)
  }
  
  if (length(variables) == 1) {
    y = X[,variables] * beta[variables] + rnorm(n, 0, sigma)    
  } else {
    y = X[,variables] %*% beta[variables] + rnorm(n, 0, sigma)
  }
  return(list(X=X, y=y))
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

pick_best_model2 = function(X, y, k=4, true_model)
{
  n = length(y)
  A = ASM(5)
  X = cbind(1, X)
  rss = lm.rss(A, y, X)
  sigma = function(a) sqrt(rss/(n - a))
  df = lm.df(A, y, X)
  models_we_want = apply(A, 1, sum) %in% (1 + c(1, 3, 4))
  #print(A[models_we_want])
  rss = rss[models_we_want]
  df = df[models_we_want]
  #print(l_n)
  #print(rss)
  #print(df)
  #print(A[models_we_want,])
  
  # Calculate the measures.
  p = df - 1
  aic_v = 2 * n * log(sigma(0)) + .5 * n  + 2 * df  
  bic_v = 2 * n * log(sigma(0)) + .5 * n  + log(n) * df  
  aic_c_v = 2 * n * log(sigma(0)) + .5 * n + (2*(p + 1) * n)/(n - p - 2)
  aic_cstar_v = 2 * n * log(sigma(p + 2)) + .5 * (n - (p + 2)) + 2 * (p + 1)
  
  #print(aic_v)
  #print(bic_v)
  #print(aic_c_v)
  #print(aic_cstar_v)

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
              best_aic_cstar_model=best_aic_cstar_model,
              A=A[models_we_want,]))
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
run = function(MAX_B=1000, n=16, vars=4, sigma=1)
{
  true_model = sort(sample(1:5, vars))
  aic_models = NULL
  bic_models = NULL
  aic_c_models = NULL
  aic_cstar_models = NULL
  for (B in 1:MAX_B) {
    # Resample
    test_data = generate_test_data(variables=true_model, n=n, sigma=sigma)
    X = test_data$X
    y = test_data$y
    results = pick_best_model2(X, y, vars, true_model)
    aic_models[B] = results$best_aic_model
    bic_models[B] = results$best_bic_model
    aic_c_models[B] = results$best_aic_c_model
    aic_cstar_models[B] = results$best_aic_cstar_model
    A = results$A
  }
  
  true_model_row = c(1, 0, 0, 0, 0, 0)
  true_model_row[true_model+1] = 1
  print("true_model_row")
  print(true_model_row)
  true_model_number = NA
  for (row in 1:dim(A)[1]) {
    if (all(A[row,] == true_model_row))
      true_model_number = row
  }
  
  prop_true_model_aic = sum(as.numeric(aic_models==true_model_number))/MAX_B
  prop_true_model_bic = sum(as.numeric(bic_models==true_model_number))/MAX_B
  prop_true_model_aic_c = sum(as.numeric(aic_c_models==true_model_number))/MAX_B
  prop_true_model_aic_cstar = sum(as.numeric(aic_cstar_models==true_model_number))/MAX_B
    
  return(list(aic_models=aic_models,
              bic_models=bic_models,
              aic_c_models=aic_c_models,
              aic_cstar_models=aic_cstar_models,
              A=A,
              true_model_row=true_model_row,
              true_model_number=true_model_number,
              prop_true_model_aic=prop_true_model_aic,
              prop_true_model_bic=prop_true_model_bic,
              prop_true_model_aic_c=prop_true_model_aic_c,
              prop_true_model_aic_cstar=prop_true_model_aic_cstar))
}

display = function(result, pdf_filename)
{
  pdf(pdf_filename)
  par(mfrow=c(2, 2))
  hist(result$aic_models, main="AIC models selected")
  hist(result$bic_models, main="BIC models selected")
  hist(result$aic_c_models, main="AIC C models selected")
  hist(result$aic_cstar_models, main="AIC C* models selected")
  dev.off()
}

result = run(n=16, vars=1, sigma=1.5)
display(result, "sixteen_samples_one_var_hist.pdf")
result = run(n=64, vars=1, sigma=1.5)
display(result, "sixty_four_samples_one_var_hist.pdf")
result = run(n=16, vars=3, sigma=1.5)
display(result, "sixteen_samples_three_vars_hist.pdf")
result = run(n=64, vars=3, sigma=1.5)
display(result, "sixty_four_samples_three_vars_hist.pdf")
result = run(n=16, vars=4, sigma=1.5)
display(result, "sixteen_samples_four_vars_hist.pdf")
result = run(n=64, vars=4, sigma=1.5)
display(result, "sixty_four_samples_four_vars_hist.pdf")
