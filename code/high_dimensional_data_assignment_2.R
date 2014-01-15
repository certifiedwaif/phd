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
  sum(abs(x - median(x)))/length(x)
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

find_best_model = function(A, y, X, method="lm")
{
  n = length(y)
  M = NULL
  #print(X)
  #print(A)
  for (model_num in 1:dim(A)[1]) {
    print("model_num")
    print(model_num)
    #print(A[model_num,])
    print(A[model_num,]*1:3)
    columns = A[model_num,]*1:3
    columns = columns[columns != 0]
    X_alpha = X[,columns]
    #print(X_alpha)
    #print(dim(X_alpha))
    #print(X_alpha)
    if (method == "lm") {
      fit = lm(y~X_alpha)
    }
    if (method == "MM") {
      fit = rlm(X_alpha, y, method="MM")
    }
    resid = residuals(fit)
    print(sum(resid^2))
    sigma = MAD(resid)
    p = sum(A[model_num,]) - 1
    M[model_num] = sum(rho(resid))/sigma + log(n) * p
  }
  print(M)
  return(which.min(M))
}

main = function(method="lm")
{
  dat = data_set()
  y = dat[,1]
  X = cbind(1, dat[,2:3])
  
  A = ASM(2)
  A = A[2:4,] # Exclude the intercept model
  model_selected = NULL
  for (B in 1:100) {
    bootstrap_rows = sample(1:64, 64, replace=TRUE)
    bootstrap_y = y[bootstrap_rows]
    bootstrap_X = X[bootstrap_rows,]
    
    model_selected[B] = find_best_model(A, bootstrap_y, as.matrix(bootstrap_X),
                                        method=method)
  }
  return(model_selected)
}
model_selected_lm = main(method="lm")
model_selected_MM = main(method="MM")
par(mfrow=c(2, 1))
hist(model_selected_lm)
hist(model_selected_MM)
