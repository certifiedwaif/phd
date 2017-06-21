# zero_inflated_model.R

blockDiag <- function(mA, mB)
{
  stopifnot(nrow(mA) == ncol(mA))
  stopifnot(nrow(mB) == ncol(mB))
  mResult <- matrix(0, nrow(mA) + nrow(mB), ncol(mA) + ncol(mB))
  mResult[1:nrow(mA), 1:ncol(mA)] <- mA
  mResult[(nrow(mA) + 1):(nrow(mA) + nrow(mB)), (ncol(mA) + 1):(ncol(mA) + ncol(mB))] <- mB
  return(mResult)
}

#' Create a mult object which can be fit using the zipvb function.
#'
#' @param vy The vector of responses.
#' @param mX The matrix of fixed effects covariates.
#' @param mZ The matrix of random effects covariates.
#' @param sigma2.beta The variance of the fixed effects, beta.
#' @param m The number of groups of random effects.
#' @param blocksize The block size of the random effects e.g. 1 for random intercepts, 2 for random slopes.
#' @param spline_dim The dimension of the splines
#' @param v The degrees of freedom of the Inverse WIshart distribution
#' @return A mult object, initialised to contain all of the parameters.
#' @export
create_mult <- function(vy, mX, mZ, sigma2.beta, mPsi=NULL, m=ncol(mZ), blocksize=1, spline_dim=0, v=0)
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
  if (is.null(mPsi)) {
    mPsi <- 1e-5 * diag(1, d)
  }

  prior <- list(v=v, mPsi=mPsi,
                a_rho=1, b_rho=1,
                sigma2.beta=sigma2.beta)
  v <- prior$v

  if (!is.null(ncol(mZ))) {
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
  # cat("T1", T1, "\n")
	T1 <- T1 - t(vp) %*% exp(mC %*% vmu + .5 * diag(mC %*% mLambda %*% t(mC)))
  # cat("T1", T1, "\n")
  # cat("diag(mC %*% mLambda %*% t(mC)))", diag(mC %*% mLambda %*% t(mC)), "\n")
	T1 <- T1 - sum(lgamma(vy + 1))
  # cat("T1", T1, "\n")
  if (p > 0) {
	 T1 <- T1 - .5 * p * log(prior$sigma2.beta) - .5 * (sum(vmu[beta_idx] ^ 2) + tr(mLambda[beta_idx, beta_idx])) / prior$sigma2.beta
    # cat("T1", T1, "\n")
  }
	T1 <- T1 + .5 * (p + m) + .5 * log(det(mLambda))
  # cat("T1", T1, "\n")
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

#' Fit a zero-inflated Poisson model using Variational Bayes/Gaussian Variational Approximation.
#'
#' @param mult A mult object, created by create_mult.
#' @param method The algorithm used for fitting the model. One of "laplace", "gva", "gva2",
#'               "gva_nr".
#' @param verbose A logical flag, for controlling the amount of debugging output.
#' @param plot_lower_bound A logical flag, for controlling whether the lower bound is plotted
#'                        as the fitting algorithm is being run.
#' @param glm_init A logical flag, for controlling whether the model fitting algoritym is
#'                 initialised using co-efficients from a Poisson model fit using glm.
#' @return A list detailing the model fit.
#' @export
zipvb <- function(mult, method="gva", verbose=FALSE, plot_lower_bound=FALSE, glm_init=TRUE)
{
  # browser()
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
  d <- mult$d
  m <- mult$m
  spline_dim <- mult$spline_dim
  blocksize <- mult$blocksize
  v <- prior$v + m + p

  zero.set <- which(vy == 0)
  vlower_bound <- c()

  # Make an initial guess for vmu
  if (glm_init) {
    fit <- glm.fit(mult$mC, mult$vy, family=poisson())
    vmu <- fit$coefficients
    #mLambda <- vcov(fit)
  }
  # vmu[is.na(vmu)] <- 0
  i <- 0
  # Iterate ----
  # TODO: Add check on whether parameters are still changing
  while ((i <= 2) || is.nan(vlower_bound[i] - vlower_bound[i - 1]) ||
        (vlower_bound[i] > vlower_bound[i - 1]) ||
        sqrt(sum((mult$vmu - old_vmu)^2)) > 1e-8) {
  # for (i in 1:MAXITER) {
    if (i > 0 && verbose) {
      cat("vmu", mult$vmu, "\n")
      cat("old_vmu", old_vmu, "\n")
      cat("Diff in norm of vmu", sqrt(sum((mult$vmu - old_vmu)^2)), "\n")
    }
    if (i >= MAXITER) {
      if (verbose) {
        cat("Iteration limit reached, breaking ...")
      }
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
    old_vmu <- vmu
    if (method == "laplace") {
      fit1 <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
    } else if (method == "gva") {
      if (i <= 1) {
        fit1 <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
        eig <- eigen(fit1$mLambda)
        if (any(is.complex(eig$values)) || any(eig$values < 1e-8)) {
          fit1 <- fit.GVA(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B")
        }
      } else {
        # Idea: If L-BFGS-B fails due to solve a computationally singular linear system,
        # try another Laplace update instead.
        # fit_lap <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
        # fit1 <- fit.GVA(fit_lap$vmu, fit_lap$mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B")
        fit1 <- fit.GVA(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B")
      }
      # fit1 <- fit.GVA(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B")
    } else if (method == "gva2") {
      if (i <= 1) {
        # fit1 <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
        fit1 <- fit.GVA_nr(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B", p=p, m=m,
                          blocksize=blocksize, spline_dim=spline_dim)
        eig <- eigen(fit1$mLambda)
        # Very occasionally, fit.Lap will return a non-invertible mLambda. In that case,
        # fall back to L-BFGS.
        if (any(is.complex(eig$values)) || any(eig$values < 1e-8)) {
          fit1 <- fit.GVA_new(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B", p=p, m=m,
                        blocksize=blocksize, spline_dim=spline_dim)
        }
      } else {
        fit1 <- fit.GVA_new(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B", p=p, m=m,
                        blocksize=blocksize, spline_dim=spline_dim)
        # If you use L-BFGS-B, this code will overflow frequently - about one time in ten.
        # If you use BFGS or CG, you'll get strange errors about one time in a thousand.
        # CG also isn't very accurate.
        # fit1 <- fit.GVA_new(vmu, mLambda, vy, vp, mC, mSigma.inv, "BFGS")
      }
    } else if (method == "gva_nr") {
      if (i <= 1) {
        fit1 <- fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
      } else {
        fit1 <- fit.GVA_nr(vmu, mLambda, vy, vp, mC, mSigma.inv, "L-BFGS-B", p=p, m=m,
                          blocksize=blocksize, spline_dim=spline_dim)
      }
    } else {
      stop("method must be either laplace, gva, gva2 or gva_nr")
    }

    vmu <- fit1$vmu
    mLambda <- fit1$mLambda
    f <- fit1$res$value

    # Update parameters for q_vr
    if (length(zero.set) != 0) {
      #cat("zero.set", zero.set, "\n")
      vp[zero.set] <- expit((vy[zero.set] * mC[zero.set,]) %*% vmu - exp(mC[zero.set,] %*% vmu +
                            0.5 * diag(mC[zero.set,] %*% mLambda %*% t(mC[zero.set,])) +
                            digamma(a_rho) - digamma(b_rho)))
      #cat("vp[zero.set]", vp[zero.set], "\n")
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
        mSigma.u.inv <- solve(mPsi/(v - d - 1))
      } else {
        mSigma.u.inv <- kronecker(diag(1, m - 1), solve(mPsi/(v - d - 1)))
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

  mult$mSigma <- solve(mSigma.inv + diag(1e-8, ncol(mSigma.inv)))
  mult$vlower_bound <- vlower_bound
  mult$iterations <- i

  return(mult)
}
