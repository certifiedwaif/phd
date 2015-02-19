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

# Create mZ matrix for random slopes
# There's probably a better way to do this
makeZ <- function(mX, m, ni, p=1)
{
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  mZ <- kronecker(diag(1,m),rep(1,ni))
  # Create mZ matrix for random slopes
  # There's probably a better way to do this
  mZ2 = matrix(0, nrow=m*ni, ncol=2*m)
  for (i in 1:m) {
    row_idx = ni*(i-1)+1:ni
    mZ2[row_idx,((i-1)*p+1):(i*p)] = mX[row_idx,]
  }
  mZ2
}

generate_test_data = function(m, ni)
{
  m = m
  ni = ni
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  #print("mX=")
  #print(mX)
  #cat("dim(mX)", dim(mX), "\n")
  
  #v = c(rep(1, g), rep(0, g))
  # Indicator variables for groups
  
  #mZ = matrix(cbind(v, 1-v), sum(n), 2)
  #mZ <- matrix(0,sum(n),m)
  #count <- 0
  #for (i in 1:m) {
  #  mZ[count + (1:n[i]),i] <- 1
  #  count <- count + n[i]
  #}
  
  mZ <- kronecker(diag(1,m),rep(1,ni))
  
  #print("mZ=")
  #print(mZ)
  #cat("dim(mZ)", dim(mZ), "\n")
  
  expected_rho = 0.5
  expected_beta = c(2, 1)
  expected_sigma2_u = .5^2
  a_sigma = 1e-2
  b_sigma = 1e-2
  
  tau = 1.0E2
  
  sigma2.beta <- 1.0E3
  
  test_data = gen_mult_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u)
  vy = test_data$vy
  
  # Test accuracy
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau, m=m, blocksize=1, spline_dim=0)
  
  return(mult)
}

generate_slope_test_data = function() {
  m = 10
  ni = 5
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  mZ = makeZ(mX, m, ni, p=2)
  
  expected_rho = 0.5
  expected_beta = c(2, 1)
  expected_sigma2_u = .5^2
  a_sigma = 1e-2
  b_sigma = 1e-2
  
  tau = 1.0E2
  
  sigma2.beta <- 1.0E5
  
  test_data = gen_mult_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=TRUE)
  vy = test_data$vy
  
  # Test model fitting
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau, m=m, blocksize=2, spline_dim=0)
  return(mult)
}

generate_spline_test_data = function()
{
  n = 5000
  vx = matrix(sort(runif(n, -1, 1))) 
  
  mX = cbind(1,vx)
  
  expected_rho = 1
  #expected_mu = c(0, 1)
  expected_sigma2_u = 0
  sigma2.beta = 1e5
  a_sigma = 1e5
  b_sigma = 1e5
  tau = 1.0E-5
  
  sigma2.true = 0.01
  expected_beta = c(0, 1)
  vf = 5+2*sin(pi*vx)
  vy = rpois(n,exp(vf))
  
  source("ZOsull.r")
  numIntKnots <- 35
  intKnots <- quantile(unique(vx),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
  
  mZ = ZOSull(vx,range.x=c(-1.1,1.1),intKnots=intKnots,drv=0)
  #vy = 2+mX[,1]^3+rnorm(m)*.1
  #result = fit_spline(vx, vy)
  #result = fit_spline(mX[,1], vy)
  #mZ = result$Z
  
  #mZ <- mZ/max(mZ)
  
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau, m=0, blocksize=1, spline_dim=37)
  
  return(mult)
}

