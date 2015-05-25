# zero_inflated_model.R
library(limma)
source("common.R")
source("gaussian.R")


calculate_lower_bound <- function(mult, verbose=FALSE)
{
	vy <- mult$vy
  vp <- mult$vp
  vmu <- mult$vmu
  mX <- mult$mX
  mC <- mult$mC
  m <- mult$m
  p <- mult$p
  v <- mult$v
  d <- mult$d
  mPsi <- mult$mPsi
  a_rho <- mult$a_rho
  b_rho <- mult$b_rho
  mLambda <- mult$mLambda
  prior <- mult$prior

	if (!is.null(mX)) {
  	beta_idx <- 1:p
  }
  
	zero.set <- which(vy == 0)

	# Terms for (beta, u)
	T1 <- t(vy * vp) %*% mC %*% vmu
	T1 <- T1 - t(vp) %*% exp(mC %*% vmu + .5 * diag(mC %*% mLambda %*% t(mC)))
	T1 <- T1 - sum(lgamma(vy + 1))
  if (p > 0) {
	 T1 <- T1 - .5 * p * log(prior$sigma2.beta) - .5 * (sum(vmu[beta_idx] ^ 2) + tr(mLambda[beta_idx, beta_idx])) / prior$sigma2.beta
  }
	T1 <- T1 + .5 * (p + m) + .5 * log(det(mLambda))
	eps <- 1e-10

	E_log_det_mSigma_u <- sum(digamma(.5 * (v + 1 - 1:d))) + d * log(2) - log(det(mPsi))
	E_mSigma_u_inv <- solve(mPsi) / (v + d + 1)
	T2 <- .5 * (prior$v * log(det(prior$mPsi)) - v * log(det(mPsi)))
	T2 <- T2 - lgammap(.5 * prior$v, d) + lgammap(.5 * v, d)
	T2 <- T2 + .5 * E_log_det_mSigma_u - .5 * tr(E_mSigma_u_inv * (prior$mPsi - mPsi))

	#if (verbose) cat("calculate_lower_bound: ", vp[zero.set], "\n")
	# 0 log 0 returns NaN. We sidestep this by adding epsilon to vp[zero.set]
	T3 <- sum(-vp[zero.set] * log(vp[zero.set] + eps) - (1 - vp[zero.set]) * log(1 - vp[zero.set] + eps))
	T3 <- T3 - lbeta(prior$a_rho, prior$b_rho) + lbeta(a_rho, b_rho)

	if (verbose) {
	  	cat("calculate_lower_bound: T1", T1, "T2", T2, "T3", T3, "\n")
	}
  
	return(T1 + T2 + T3)
}

create_mult <- function(vy, mX, mZ, sigma2.beta, m=ncol(mZ), blocksize=1, spline_dim=0, v=0)
{
  # Initialise
  n <- length(vy)
  vp <- rep(1, n)
  if (is.null(mZ)) {
    mC <- mX
  } else if (is.null(mX)) {
    mC <- mZ
  }
  else {
    mC <- cbind(mX, mZ)
  }
  vmu <- rep(0, ncol(mC))
  
  if (!is.null(ncol(mX))) {
    p <- ncol(mX)
    mSigma.beta.inv <- diag(1 / sigma2.beta, ncol(mX))
  } else {
    p <- 0
    mSigma.beta.inv <- NULL
  }
  d <- blocksize + spline_dim
  mLambda <- diag(rep(1, ncol(mC)))
  a_rho <- 1 + sum(vp)
  b_rho <- n - sum(vp) + 1
  
  # Set prior parameters for Inverse Wishart distribution
  #mPsi=diag(1, rep(blocksize))
  mPsi <- 1e-5 * diag(1, d)
  
  prior <- list(v=v, mPsi=mPsi,
                a_rho=1, b_rho=1,
                sigma2.beta=sigma2.beta)
  v <- prior$v + m
  
  if (!is.null(ncol(mZ))) {
    # TODO: We don't yet handle the case where there are splines and random intecepts/
    # slopes.
    if (m == 1) {
      mSigma.u.inv <- diag(1, spline_dim)
    } else {
      mSigma.u.inv <- kronecker(diag(1, (m - 1)), mPsi)
    }
  } else {
    mSigma.u.inv <- NULL
  }

  if (is.null(mSigma.u.inv)) {
    mSigma.inv <- mSigma.beta.inv
  } else if (is.null(mSigma.beta.inv)) {
    mSigma.inv <- mSigma.u.inv
  } else {
    mSigma.inv <- blockDiag(mSigma.beta.inv, mSigma.u.inv)
  }
  
  mult <- list(vy=vy, vp=vp, vmu=vmu,
                mX=mX, mZ=mZ, mC=mC,
                m=m, blocksize=blocksize, spline_dim=spline_dim,
                p=p, d=d, v=v, mPsi=mPsi,
                a_rho=a_rho, b_rho=b_rho,
                mLambda=mLambda,
                prior=prior,
                mSigma.beta.inv=mSigma.beta.inv,
                mSigma.u.inv=mSigma.u.inv)
  class(mult) <- "multivariate"
  return(mult)
}

zero_infl_var <- function(mult, method="gva", verbose=FALSE, plot_lower_bound=FALSE)
{
  MAXITER <- 5
  
  # Initialise variables from mult
  N <- length(mult$vy)
  p <- mult$p
  vy <- mult$vy
  vp <- mult$vp
  mZ <- mult$mZ
  mC <- mult$mC
  mSigma.beta.inv <- mult$mSigma.beta.inv
  mSigma.u.inv <- mult$mSigma.u.inv
  mSigma.inv <- mult$mSigma.inv
  vmu <- mult$vmu
  mLambda <- mult$mLambda
  a_rho <- mult$a_rho
  b_rho <- mult$b_rho
  prior <- mult$prior
  mPsi <- mult$mPsi
  v <- mult$v
  d <- mult$d
  m <- mult$m
  spline_dim <- mult$spline_dim
  blocksize <- mult$blocksize

  if (verbose) {
    cat("N", N, "\n")
    cat("p", p, "\n")
  }

  if (!is.null(mZ)) {
    if (verbose) {
      cat("m", m, "\n")
      cat("blocksize", blocksize, "\n")
      cat("spline_dim", spline_dim, "\n")
    }
  }
  
  zero.set <- which(vy == 0)
  vlower_bound <- c()
  
  i <- 0
  # Iterate ----
  # TODO: Add check on whether parameters are still changing
  while ((i <= 1) || is.nan(vlower_bound[i] - vlower_bound[i - 1]) ||
        (vlower_bound[i] > vlower_bound[i - 1])) {
  #for (i in 1:MAXITER) {
    if (i >= MAXITER) {
      cat("Iteration limit reached, breaking ...")
      break
    }
    
    i <- i + 1
    
    if (is.null(mSigma.u.inv)) {
      mSigma.inv <- mSigma.beta.inv
    } else if (is.null(mSigma.beta.inv)) {
      mSigma.inv <- mSigma.u.inv
    } else {
      mSigma.inv <- blockDiag(mSigma.beta.inv, mSigma.u.inv)
    }
    
    # Update parameter for q_vnu by maximising using the Gaussian Variational Approximation from Dr Ormerod's Poisson mixed model code
    if (method == "laplace") {
      fit1 <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
    } else if (method == "gva") {	
      fit_lap <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
      fit1 <- fit.GVA(fit_lap$vmu, fit_lap$mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B")
      # fit1 <- fit.GVA(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B")
    } else if (method == "gva2") {
      fit_lap <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
      fit1 <- fit.GVA_new(fit_lap$vmu, fit_lap$mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B", p=p, m=m,
                         blocksize=blocksize, spline_dim=spline_dim)
    } else if (method == "gva_nr") {
      fit1 <- fit.GVA_nr(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B", p=p, m=m, 
                        blocksize=blocksize, spline_dim=spline_dim)
    } else {
      stop("method must be either laplace, gva, gva2 or gva_nr")
    }
    
    vmu <- fit1$vmu
    mLambda <- fit1$mLambda
    f <- fit1$res$value
    
    # Update parameters for q_vr
    if (length(zero.set) != 0) {
      vp[zero.set] <- expit((vy[zero.set] * mC[zero.set,]) %*% vmu - exp(mC[zero.set,] %*% vmu + 
                            0.5 * diag((mC[zero.set,]) %*% mLambda %*% t(mC[zero.set,])) + 
                            digamma(a_rho) - digamma(b_rho)))
    }
    
    # Update parameters for q_rho
    a_rho <- prior$a_rho + sum(vp)
    b_rho <- prior$b_rho + N - sum(vp)
    
    # Update parameters for q_sigma_u^2 if we need to
    if (!is.null(mZ)) {
      u_dim <- (m - 1) * blocksize + spline_dim
      u_idx <- p + 1:u_dim  # (ncol(mX)+1):ncol(mC)

      # a_sigma is fixed
      #a_sigma <- prior$a_sigma + u_dim/2
      #b_sigma <- prior$b_sigma + sum(vu^2)/2 + (tr_mSigma)/2
      #tr_mSigma <- ncol(mZ) * prior$a_sigma/prior$b_sigma

      # Extract right elements of mLambda
      #b_sigma <- prior$b_sigma + sum(vmu[u_idx]^2)/2 + tr(mLambda[u_idx, u_idx])/2    
      
      if (spline_dim == 0) {
        acc <- matrix(0, blocksize, blocksize)
        for (j in 1:(m - 1)) {
          j_idx <- p + (j - 1) * blocksize + (1:blocksize)
          acc <- acc + vmu[j_idx] %*% t(vmu[j_idx]) + mLambda[j_idx, j_idx]
        }
      } else {
        acc <- matrix(0, spline_dim, spline_dim)
        spline_idx <- (p + 1):length(vmu)
        acc <- vmu[spline_idx] %*% t(vmu[spline_idx]) + mLambda[spline_idx, spline_idx]
      }
      mPsi <- prior$mPsi + acc
      
      #tau_sigma <- a_sigma/b_sigma
      #mSigma.u.inv <- diag(tau_sigma, u_dim)
      # TODO: This doesn't handle the case where you have splines and random intercepts
      # or slopes.
      if (m == 1) {
        mSigma.u.inv <- solve(mPsi)/(v - d - 1)
      } else {
        mSigma.u.inv <- kronecker(diag(1, m - 1), solve(mPsi)/(v - d - 1))
      }
      #mSigma.u.inv <- solve(mPsi) # What multiplicative factor for psi?
    }
    
    # Update all variables in mult so that we can calculate the lower
    # bound.
    mult$mSigma.beta.inv <- mSigma.beta.inv
    mult$mSigma.u.inv <- mSigma.u.inv
    mult$mSigma.inv <- mSigma.inv
    mult$vmu <- vmu
    mult$mLambda <- mLambda
    mult$a_rho <- a_rho
    mult$b_rho <- b_rho
    mult$mPsi <- mPsi
    mult$f <- f

    vlower_bound[i] <- calculate_lower_bound(mult, verbose=verbose)
    
    if (verbose) {
      if (i > 1)
        cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
            vlower_bound[i] - vlower_bound[i - 1], " parameters ", "vmu", round(vmu, 2),
            "diag(mLambda)", round(diag(mLambda), 2), "a_rho", round(a_rho, 2), "b_rho",
            round(b_rho, 2))
      if (!is.null(mZ)) {
        # cat(" mSigma.u.inv ", mSigma.u.inv)
      }
      cat("\n")
    }
  }
  
  if (plot_lower_bound)
    plot(vlower_bound, type="l")
  
  params <- list(vmu=vmu, mLambda=mLambda, a_rho=a_rho, b_rho=b_rho,
                 mSigma=solve(mSigma.inv), vlower_bound=vlower_bound)
  return(params)
}
