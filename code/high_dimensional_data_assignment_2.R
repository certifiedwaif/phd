# high_dimensional_data_assignment2.R
require(MASS)

data_set = function()
{
  set.seed(2014)
  n = 64
  x1 = sort(rnorm(n))
  x2 = runif(n, -5, 5)
  beta = c(10, 0, 1)
  error = rnorm(n)
  X = cbind(1L, x1, x2)
  y = X %*% beta + error
  y[1:2] = 100
  dat = data.frame(y, x1, x2)
  dat
}

MAD = function(x)
{
  1.483 * median(abs(x - median(x)))
}

rho = function(x)
{
  min(x^2, 2^2)
}

# ASM : function to select all possible submodels
# input : k.s = number of slope parameters
# intcpt = TRUE, hence all models include an intercept
# output: indicator matrix where rows represent models and columns variables
ASM = function(k.s,intcpt = TRUE)
{
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

fit_model = function(y, X, regression_method="lm")
{
  if (regression_method == "lm") {
    fit = lm(y~X)
  }
  if (regression_method == "MM") {
    fit = rlm(X, y, method="MM")
  }
  resid = residuals(fit)
  return(list(resid=resid, sigma=MAD(resid), coef=coef(fit)))
}

M_1 = function(result)
{
  resid = result$resid
  sigma = result$sigma
  n = result$n
  p = result$p
  (sum(rho(resid^2/sigma)) + log(n) * p)/n
}

M_2 = function(y, X, m=16, regression_method="lm")
{
  n = length(y)
  M_2 = NULL
  # Bootstrap for M_2
  for (B in 1:100) {
    bootstrap_rows = sample(1:64, m, replace=TRUE)
    bootstrap_y = y[bootstrap_rows]
    bootstrap_X = X[bootstrap_rows,]
    
    model_result = fit_model(bootstrap_y, as.matrix(bootstrap_X),
                                        regression_method=regression_method)

    beta_star_vec = model_result$coef
    beta_star_vec = beta_star_vec[!is.na(beta_star_vec)]
    beta_star = matrix(beta_star_vec, nrow=length(beta_star_vec), ncol=1)

    #resid = model_result$resid
    #sigma = model_result$sigma
    #print(dim(y))
    #print(dim(X))
    #print(dim(beta_star))
    #print(beta_star)
    #browser()
    
    resid = y - X %*% beta_star
    sigma = MAD(resid)
    M_2[B] = sum(rho(resid^2/sigma)) # Bootstrapped version
  }
  return(sum(M_2)/(n*B))
  # Not quite. Get the coefficient. Then apply that to the entire data set.
}

main = function()
{
  m = c(16, 24, 32, 64)
  dat = data_set()
  y = dat[,1]
  X = cbind(1, dat[,2:3])

  A = ASM(2)
  A = A[2:4,] # Exclude the intercept model
  
  M_1_lm = NULL
  M_1_MM = NULL
  M_2_lm = matrix(NA, 3, 4)
  M_2_MM = matrix(NA, 3, 4)  
  M_lm = matrix(NA, 3, 1)
  M_MM = matrix(NA, 3, 1)
  
  for (model_num in 1:dim(A)[1]) {
    #print("model_num")
    #print(model_num)
    #print(A[model_num,])
    columns = A[model_num,]*1:3
    columns = sort(columns[columns != 0])
    #print(columns)
    X_alpha = as.matrix(X[,columns])
    n = length(y)
    p = sum(A[model_num,])    
    
    lm_result = fit_model(y, X_alpha, regression_method="lm")
    lm_result = c(lm_result, n=n, p=p)
    M_1_lm[model_num] = M_1(lm_result)
    
    MM_result = fit_model(y, X_alpha, regression_method="MM")
    MM_result = c(MM_result, n=n, p=p)
    M_1_MM[model_num] = M_1(MM_result)

    # For M_2 we need to bootstrap
    for (i in 1:4) {
      M_2_lm[model_num, i] = M_2(y, X_alpha, m=m[i], regression_method="lm")
      M_2_MM[model_num, i] = M_2(y, X_alpha, m=m[i], regression_method="MM")
    }
    # As the bootstrap sample size increases, M2 goes down. Why?
    
    #for (i in 1:4) {
      M_lm[model_num, 1] = MAD(M_1_lm[model_num] + M_2_lm[model_num,])    
      M_MM[model_num, 1] = MAD(M_1_MM[model_num] + M_2_MM[model_num,])    
    #}
  }
  return(list(M_1_lm=M_1_lm,
              M_1_MM=M_1_MM,
              M_2_lm_16=M_2_lm[,1],
              M_2_MM_16=M_2_MM[,1],
              M_2_lm_24=M_2_lm[,2],
              M_2_MM_24=M_2_MM[,2],
              M_2_lm_32=M_2_lm[,3],
              M_2_MM_32=M_2_MM[,3],
              M_2_lm_64=M_2_lm[,4],
              M_2_MM_64=M_2_MM[,4],
              M_lm=M_lm,
              M_MM=M_MM))
}
main()

# Theory: The reason that the RSS is so much bigger for the full data set
# is that the whole regression line fits the whole data set very badly, whereas
# fitting to a subsample is much less egregious.

# Basic checks ----
ds = data_set()
y = ds[,1]
x = ds[,2]
plot(x, y)
lm_fit = lm(y~x)
summary(lm_fit)
# Large, spurious x coefficient
plot(lm_fit)
# Outlying responses exert very high leverage
mm_fit = rlm(y~x, method="MM")
summary(mm_fit)

x = ds[,3]
plot(x, y)
lm_fit = lm(y~x)
summary(lm_fit)
# Overestimate of the true co-efficients.
mm_fit = rlm(y~x, method="MM")
summary(mm_fit)
# Much closer to the true co-efficient values.
