# generate.R
library(mvtnorm)
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

gen_int_data <- function (vx, vbeta, vu, rho, m, ni)
{
  eta <- matrix(NA, m, ni)
  zeta <- matrix(NA, m, ni)
  r <- matrix(NA, m, ni)
  y <- matrix(NA, m, ni)
  
  for (i in 1:m) {
    for (j in 1:ni) {
      eta[i, j] <- vbeta[1] + vbeta[2] * vx[i, j] + vu[i, 1]
      zeta[i, j] <- rpois(1, exp(eta[i, j]))
      r[i, j] <- rbinom(1, 1, rho)
      y[i, j] <- r[i, j] * zeta[i, j]
    }
  }
  
  return(as.vector(y))
}

generate_int_test_data <- function(m, ni, expected_beta = c(2, 1), expected_rho = 1.0)
{
  n <- m * ni
  x <- rnorm(n, 0, 1)
  vx <- t(matrix(x, ni, m))
  mX <- cbind(rep(1, n), x)
  
  mZ <- kronecker(diag(1,m), rep(1,ni))
  mZ <- mZ[, 2:m]

  mSigma_0 <- matrix(c(1.0), 1, 1)
  vu <- rmvnorm(m, sigma=mSigma_0)
  
  expected_sigma2_u <- .5^2
  
  vy <- gen_int_data(vx, expected_beta, vu, expected_rho, m, ni)
  
  sigma2.beta <- 1.0E5
  # Test accuracy
  mult <- create_mult(vy, mX, mZ, sigma2.beta, m=m, blocksize=1, spline_dim=0, v=4)
  
  return(mult)
}

gen_slope_data <- function(mX, vbeta, vu, rho, m, ni)
{
  eta <- matrix(NA, m, ni)
  zeta <- matrix(NA, m, ni)
  r <- matrix(NA, m, ni)
  y <- matrix(NA, m, ni)
  
  for (i in 1:m) {
    for (j in 1:ni) {
      eta[i, j] <- vbeta[1] + vbeta[2]*mX[i, j] + vu[i, 1] + vu[i, 2] * mX[i, j]
      zeta[i, j] <- rpois(1, exp(eta[i, j]))
      r[i, j] <- rbinom(1, 1, rho)
      y[i, j] <- r[i, j] * zeta[i, j]
    }
  }
  
  return(as.vector(y))
}

generate_slope_test_data <- function(m=10, ni=20, expected_beta=c(2, 1), expected_rho=0.5)
{
  n <- m * ni
  x <- rnorm(n, 0, 1)
  vx <- t(matrix(x, ni, m))
  # TODO: Rewrite this
  groups <- gl(m, ni)
  mC <- model.matrix(~1+x+groups*x)
  p <- 2
  mX <- mC[,1:p]
  mZ <- mC[,p+(1:((m-1)*p))]
  # Re-order columns of z so that columns for same groups are adjacent
  # This will ensure banded structure of t(mC) %*% mC
  # Take 1, 2, 3, 4, 5, 6, 7, 8, 9, 10
  # to   1, 6, 2, 7, 3, 8, 4, 9, 5, 10
  # This works because R stores its matrices using the Fortran convention of
  # column-major ordering. Yes, I was surprised too!
  # TODO: This is ghastly. Rewrite to use the Kronecker product.
  ordering <- rbind(1:(m-1), (m-1)+1:(m-1))
  mZ_reordered <- mZ[, as.vector(ordering)]
  
  # Centre slope term?
  #mX <- cbind(mX[,1], scale(mX[,2]))
  vbeta <- c(2, 1)
  mSigma_0 <- matrix(c(1.0, -0.3,
                      -0.3,  1.0), 2, 2)
  vu <- rmvnorm(m, sigma=mSigma_0)
  rho <- 0.5
  vy <- gen_slope_data(vx, vbeta, vu, rho, m, ni)
  
  # Create mult object
  sigma2.beta <- 1.0E5
  mult <- create_mult(vy, mX, mZ_reordered, sigma2.beta, m=m, blocksize=2, spline_dim=0, v=2+2)
  return(mult)
}

generate_spline_test_data <- function(n=100)
{
  a <- -1
  b <- 1
  vx <- matrix(sort(runif(n, a, b))) 
  
  mX <- cbind(rep(1, n), vx)
  # mX <- NULL
  sigma2_beta <- 1.0E5
  
  vf <- 4 + sin(pi * vx)
  vy <- rpois(n, exp(vf))
  
  numIntKnots <- 10
  intKnots <- quantile(unique(vx), seq(0, 1, length=(numIntKnots + 2))[-c(1, (numIntKnots + 2))])
  
  result <- formOmega(a, b, intKnots)
  mZ <- ZOSull2(vx, range.x=range(vx), intKnots=intKnots, allKnots=result$allKnots,
                Omega=result$Omega, drv=0)
  #mZ <- mZ/max(mZ)
  
  # mult <- create_mult(vy, mX, mZ, sigma2_beta, m=1, blocksize=0, spline_dim=12)
  sigma2_beta <- 1.0E5
  mult <- create_mult(vy, mX, mZ, sigma2_beta, m=1, blocksize=0, spline_dim=(numIntKnots + 2),
                      v=(numIntKnots + 2) + 1)
  # mult$vmu <- glm.fit(mult$mC, mult$vy, family=poisson())$coefficients
  # TODO: You need to set mPsi <- mOmega
  mult$prior$mPsi <- result$Omega
  mult$mPsi <- result$Omega
  
  # Check whether we've accidentally created a data matrix with repeated
  # columns. This can happen when, for instance, the basis vectors [1 x]
  # are inadvertenty repeated.
  # if (abs(det(crossprod(mult$mC))) < 1e-99) {
  #   stop("The cross-product of mC is singular. Perhaps you repeated some basis vectors?")
  # }
  
  return(list(mult=mult, allKnots=result$allKnots))
}

