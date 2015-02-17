# generate.R

generate_univariate_test_data <- function (n, rho, lambda)
{
  vx = rep(NA, n)
  for (i in 1:n) {
    if (runif(1) <= rho) {
      vx[i] = rpois(1, lambda)
    } else {
      vx[i] = 0
    }
  }
  return(vx)
}

gen_mult_data <- function (mX, mZ, m, n, rho, vbeta, sigma2_u, verbose=FALSE)
{
  if (is.null(mZ)) {
    mC = mX
  } else {
    mC = cbind(mX, mZ)
  }
  vy = rep(NA, sum(n))
  vu = rep(NA, length(n))
  if (verbose)
    cat("vy", vy, "\n")
  idx = 0
  for (i in 1:length(n)) {
    if (sigma2_u == 0)
      vu[i] = 0
    else
      vu[i] = rnorm(1, 0, sqrt(sigma2_u))
    
    for (j in 1:n[i]) {
      idx = idx + 1
      if (verbose)
        cat("idx", idx, "\n")
      
      if (runif(1) <= rho) {
        veta = mX[idx,] %*% vbeta + vu[i]
        if (verbose)
          cat("veta", veta, "\n")
        vy[idx] = rpois(1, exp(veta))
        if (verbose)
          cat("Generated vy[idx]", vy[idx], "\n")
      } else {
        vy[idx] = 0
      }
    }
  }
  if (verbose)
    cat("vy", vy, "\n")
  if(NA %in% vy)
    stop("NAs in vy")
  result = list(vy=vy, vu=vu)
  return(result)
}

gen_mult_data_inv_wish <- function (mX, mZ, m, n, rho, vbeta, v, psi, verbose=FALSE)
{
  if (is.null(mZ)) {
    mC = mX
  } else {
    mC = cbind(mX, mZ)
  }
  vy = rep(NA, sum(n))
  vu = rep(NA, length(n))
  if (verbose)
    cat("vy", vy, "\n")
  idx = 0
  for (i in 1:length(n)) {
    sigma2_u = solve(rWishart(1, v, solve(psi)))
    vu = rmvnorm(1, rep(0, ncol(mZ)), sigma2_u)
    #if (sigma2_u == 0)
    #  vu[i] = 0
    #else
    #  vu[i] = rnorm(1, 0, sqrt(sigma2_u))
    
    for (j in 1:n[i]) {
      idx = idx + 1
      if (verbose)
        cat("idx", idx, "\n")
      
      if (runif(1) <= rho) {
        veta = mX[idx,] %*% vbeta + vu[i]
        if (verbose)
          cat("veta", veta, "\n")
        vy[idx] = rpois(1, exp(veta))
        if (verbose)
          cat("Generated vy[idx]", vy[idx], "\n")
      } else {
        vy[idx] = 0
      }
    }
  }
  if (verbose)
    cat("vy", vy, "\n")
  if(NA %in% vy)
    stop("NAs in vy")
  result = list(vy=vy, vu=vu)
  return(result)
}

