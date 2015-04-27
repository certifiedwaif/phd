# generate.R
source("splines.R")

generate_univariate_test_data <- function (n, rho, lambda)
{
  vx <- rep(NA, n)
  for (i in 1:n) {
    if (runif(1) <= rho) {
      vx[i] <- rpois(1, lambda)
    } else {
      vx[i] <- 0
    }
  }
  return(vx)
}

gen_mult_data <- function (mX, mZ, m, n, rho, vbeta, sigma2_u, verbose=FALSE)
{
  if (is.null(mZ)) {
    mC <- mX
  } else {
    mC <- cbind(mX, mZ)
  }
  vy <- rep(NA, sum(n))
  vu <- rep(NA, length(n))
  if (verbose)
    cat("vy", vy, "\n")
  idx <- 0
  for (i in 1:length(n)) {
    if (sigma2_u == 0)
      vu[i] <- 0
    else
      vu[i] <- rnorm(1, 0, sqrt(sigma2_u))
    
    for (j in 1:n[i]) {
      idx <- idx + 1
      if (verbose)
        cat("idx", idx, "\n")
      
      if (runif(1) <= rho) {
        veta <- mX[idx,] %*% vbeta + vu[i]
        if (verbose)
          cat("veta", veta, "\n")
        vy[idx] <- rpois(1, exp(veta))
        if (verbose)
          cat("Generated vy[idx]", vy[idx], "\n")
      } else {
        vy[idx] <- 0
      }
    }
  }
  if (verbose)
    cat("vy", vy, "\n")
  if(NA %in% vy)
    stop("NAs in vy")
  result <- list(vy=vy, vu=vu)
  return(result)
}

gen_slope_data <- function(vx, vbeta, vu, rho, m, ni)
{
  eta <- matrix(NA, m, ni)
  zeta <- matrix(NA, m, ni)
  r <- matrix(NA, m, ni)
  y <- matrix(NA, m, ni)
  
  for (i in 1:m) {
    for (j in 1:ni) {
      eta[i, j] <- vbeta[1] + vbeta[2]*vx[i, j] + vu[i, 1] + vu[i, 2]*vx[i, j]
      zeta[i, j] <- rpois(1, exp(eta[i, j]))
      r[i, j] <- rbinom(1, 1, rho)
      y[i, j] <- r[i, j] * zeta[i, j]
    }
  }
  
  return(y)
}

generate_test_data <- function(m, ni)
{
  m <- m
  ni <- ni
  n <- rep(ni,m)
  mX <- matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  #print("mX=")
  #print(mX)
  #cat("dim(mX)", dim(mX), "\n")
  
  #v <- c(rep(1, g), rep(0, g))
  # Indicator variables for groups
  
  #mZ <- matrix(cbind(v, 1-v), sum(n), 2)
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
  
  expected_rho <- 0.5
  expected_beta <- c(2, 1)
  expected_sigma2_u <- .5^2
  a_sigma <- 1e-2
  b_sigma <- 1e-2
  
  tau <- 1.0E2
  
  sigma2.beta <- 1.0E3
  
  test_data <- gen_mult_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u)
  vy <- test_data$vy
  
  # Test accuracy
  mult <- create_mult(vy, mX, mZ, sigma2.beta, m=m, blocksize=1, spline_dim=0)
  
  return(mult)
}

generate_slope_test_data <- function(m=5, ni=20)
{
  n <- m*ni
  x <- rnorm(n, 0, 1)
  vx <- t(matrix(x, ni, m))
  groups <- gl(m, ni)
  mC <- model.matrix(~1+x+groups*x)
  p <- 2
  mX <- mC[,1:p]
  mZ <- mC[,p+(1:((m-1)*p))]
  # TODO: Re-order columns of z so that columns for same groups are adjacent
  # This will ensure banded structure of t(mC) %*% mC
  # Take 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  # to   1, 6, 2, 7, 3, 8, 4, 9, 5, 10
  # This works because R stores its matrices using the Fortran convention of
  # column-major ordering. Yes, I was surprised too!
  ordering <- rbind(1:(m-1), (m-1)+1:(m-1))
  mZ_reordered <- mZ[, as.vector(ordering)]
  
  # Centre slope term?
  #mX <- cbind(mX[,1], scale(mX[,2]))
  vbeta <- c(1.5, 0.5)
  mSigma_0 <- matrix(c( 1.0, -0.3,
                      -0.3,  1.0), 2, 2)
  vu <- rmvnorm(m, sigma <- mSigma_0)
  rho <- 1.0
  vy <- as.vector(gen_slope_data(vx, vbeta, vu, rho, m, ni))
  
  # Create mult object
  sigma2.beta <- 1.0E5
  mult <- create_mult(vy, mX, mZ_reordered, sigma2.beta, m=m, blocksize=2, spline_dim=0)
  return(mult)
}

generate_spline_test_data <- function(n=100)
{
  vx <- matrix(sort(runif(n, -1, 1))) 
  
  mX <- cbind(1,vx)
  
  expected_rho <- 1
  expected_sigma2_u <- 0
  sigma2_beta <- 1e5
  a_sigma <- 1e5
  b_sigma <- 1e5
  tau <- 1.0E-5
  
  sigma2.true <- 0.01
  expected_beta <- c(0, 1)
  vf <- 2 + sin(pi * vx)
  vy <- rpois(n, exp(vf))
  
  numIntKnots <- 10
  intKnots <- quantile(unique(vx), seq(0, 1, length=(numIntKnots+2))[-c(1, (numIntKnots+2))])
  
  mZ <- ZOSull(vx, range.x=c(-1.1, 1.1), intKnots=intKnots, drv=0)
  #vy <- 2+mX[,1]^3+rnorm(m)*.1
  #result <- fit_spline(vx, vy)
  #result <- fit_spline(mX[,1], vy)
  #mZ <- result$Z
  
  #mZ <- mZ/max(mZ)
  
  mult <- create_mult(vy, mX, mZ, sigma2_beta, m=1, blocksize=0, spline_dim=12)
  
  # Check whether we've accidentally created a data matrix with repeated
  # columns. This can happen when, for instance, the basis vectors [1 x]
  # are inadvertenty repeated.
  if (abs(det(tcrossprod(mult$mC))) < 1e-7) {
    stop("The cross-product of mC^T is singular. Perhaps you repeated some basis vectors?")
  }
  
  return(mult)
}

