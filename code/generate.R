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

generate_slope_test_data <- function(m=10, ni=20)
{
  n <- m*ni
  x <- rnorm(n, 0, 1)
  vx <- t(matrix(x, ni, m))
  # TODO: Rewrite this
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
  mSigma_0 <- matrix(c(1.0, -0.3,
                      -0.3,  1.0), 2, 2)
  vu <- rmvnorm(m, sigma <- mSigma_0)
  rho <- 1.0
  vy <- as.vector(gen_slope_data(vx, vbeta, vu, rho, m, ni))
  
  # Create mult object
  sigma2.beta <- 1.0E5
  mult <- create_mult(vy, mX, mZ_reordered, sigma2.beta, m=m, blocksize=2, spline_dim=0, v=2+2)
  mult$vmu <- lm.fit(mult$mC, log(mult$vy + 1))$coefficients
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
  mult$vmu <- glm.fit(mult$mC, mult$vy, family=poisson())$coefficients
  # TODO: You need to set mPsi <- mOmega
  mult$prior$mPsi <- result$Omega
  mult$mPsi <- result$Omega
  
  # Check whether we've accidentally created a data matrix with repeated
  # columns. This can happen when, for instance, the basis vectors [1 x]
  # are inadvertenty repeated.
  # if (abs(det(crossprod(mult$mC))) < 1e-99) {
  #   stop("The cross-product of mC is singular. Perhaps you repeated some basis vectors?")
  # }
  
  return(mult)
}

