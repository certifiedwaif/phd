# variational_approximation_to_zero_inflated_model.R
library(limma)
source("common.R")
source("gaussian.R")

mcmc <- function(object)
{
  UseMethod("zero_infl_mcmc", object)
}

gamma_entropy <- function(alpha, beta)
{
  alpha - log(beta) + lgamma(alpha) - (alpha-1) * digamma(alpha)
}

beta_entropy <- function(alpha, beta)
{
  lbeta(alpha, beta) - (alpha-1)*digamma(alpha) - (beta-1)*digamma(beta) + (alpha+beta-2)*digamma(alpha+beta)
}	

calculate_lower_bound.multivariate <- function(multivariate, verbose=FALSE)
{
  result = with(multivariate, {
    p = ncol(mX)
    beta_idx = 1:p
    if (!is.null(mZ)) {
      m = ncol(mZ) 
      u_idx = p + (1:m)  # (ncol(mult$mX)+1):ncol(mult$mC)
      vu = vmu[u_idx]
    } else {
      m = 0
    }
    
    zero.set <- which(vy==0)
    
    # Terms for (beta, u)
    #result = result + (vy*vp) %*% mC %*% vmu
    #mat <- matrix(0.5 * (t(vmu) %*% mLambda %*% vmu), nrow = nrow(mC), ncol = 1)
    #result = result - vp %*% exp(mC %*% vmu + mat) - sum(lgamma(vy + 1))
    #result = result + 0.5 * (det(2*pi*mLambda) + t(vmu) %*% solve(mLambda) %*% vmu)
    
    # Terms for (beta, u)
    T1 = t(vy*vp) %*% mC %*% vmu
    T1 = T1 - t(vp) %*% exp(mC %*% vmu + .5 * diag(mC%*%mLambda%*%t(mC)))
    T1 = T1 - sum(lgamma(vy + 1))
    T1 = T1 - .5*p*log(prior$sigma2.beta) - .5*(sum(vmu[beta_idx]^2) + tr(mLambda[beta_idx, beta_idx]))/prior$sigma2.beta
    if (!is.null(mZ)) {
      #T1 = T1 + prior$a_sigma * log (prior$b_sigma) - lgamma(prior$a_sigma)
      #T1 = T1 - a_sigma * log(b_sigma) + lgamma(a_sigma)
    }
    T1 = T1 + .5*(p+m) + .5*log(det(mLambda))
    eps = 1e-10
    
    lgammap = function(a, p)
    {
      .25*p*(p-1)*log(pi) + sum(lgamma(a + .5*(1 - 1:p)))
    }
    
    E_log_det_mSigma_u = sum(digamma(.5*(v + 1 - 1:p))) + p*log(2) - log(det(mPsi))
    E_mSigma_u_inv = solve(mPsi)/(v+p+1)
    T2 = .5*(prior$v*log(det(prior$mPsi)) - v*log(det(mPsi)))
    T2 = T2 - lgammap(.5*prior$v, p) + lgammap(.5*v, p)
    T2 = T2 + .5*E_log_det_mSigma_u - .5*tr(E_mSigma_u_inv*(prior$mPsi - mPsi))
    
    # Something is wrong in T3. It sometimes goes backwards as we're optimised.
    # This should be unchanged from the univariate lower bound
    #if (verbose) cat("calculate_lower_bound: ", vp[zero.set], "\n")
    # 0 log 0 returns NaN. We sidestep this by adding epsilon to vp[zero.set]
    vp[zero.set] = vp[zero.set] + eps
    T3 = sum(-vp[zero.set]*log(vp[zero.set]) - (1-vp[zero.set])*log(1-vp[zero.set]))
    T3 = T3 - lbeta(prior$a_rho, prior$b_rho) + lbeta(a_rho, b_rho)
    
	if (verbose)
    	cat("calculate_lower_bound: T1", T1, "T2", T2, "T3", T3, "\n")

    result = T1 + T2 + T3
    result
  })
  
  return(result)
}

calculate_lower_bound <- function(object, verbose=verbose)
{
  UseMethod("calculate_lower_bound", object)
}

create_multivariate <- function(vy, mX, mZ, sigma2.beta, m=ncol(mZ), blocksize=1, spline_dim=NA, v=blocksize+1)
{
  # Initialise
  n = length(vy)
  vp = rep(1, n)
  if (is.null(mZ)) {
    mC = mX
  } else {
    mC = cbind(mX, mZ)
  }
  vmu = rep(0, ncol(mC))
  
  p = ncol(mX)
  mSigma.beta.inv = diag(1/sigma2.beta, ncol(mX))
  mLambda = diag(rep(1, ncol(mC)))
  a_rho = 1 + sum(vp)
  b_rho = n - sum(vp) + 1
  
  # Set prior parameters for Inverse Wishart distribution
  #mPsi=diag(1, rep(blocksize))
  mPsi = 1e-5*matrix(c(1, 0, 0, 1), blocksize, blocksize)
  
  prior = list(v=v, mPsi=mPsi,
               a_rho=1, b_rho=1,
               sigma2.beta=sigma2.beta)
  v=prior$v + m
  
  if (!is.null(ncol(mZ))) {
    mSigma.u.inv = kronecker(diag(1, (m-1)), mPsi)
  } else {
    mSigma.u.inv = NULL
  }
  
  multivariate = list(vy=vy, vp=vp, vmu=vmu,
                      mX=mX, mZ=mZ, mC=mC,
                      m=m, blocksize=blocksize, spline_dim=spline_dim,
                      p=p, v=v, mPsi=mPsi,
                      a_rho=a_rho, b_rho=b_rho,
                      mLambda=mLambda,
                      prior=prior,
                      mSigma.beta.inv=mSigma.beta.inv,
                      mSigma.u.inv=mSigma.u.inv,
                      mSigma.beta=solve(mSigma.beta.inv))
  class(multivariate) = "multivariate"
  return(multivariate)
}

zero_infl_var.multivariate <- function(mult, method="gva", verbose=FALSE, plot_lower_bound=FALSE)
{
  MAXITER <- 30
  
  # Initialise
  N = length(mult$vy)

  if (verbose) cat("N", N, "\n")
  if (!is.null(mult$mX)) {
    p = ncol(mult$mX) 
    if (verbose) cat("p", p, "\n")
  }	else {
    p = 0
  }
  if (!is.null(mult$mZ)) {
    m = mult$m
    blocksize = mult$blocksize
    spline_dim = mult$spline_dim
    if (verbose) {
      cat("m", m, "\n")
      cat("blocksize", blocksize, "\n")
      cat("spline_dim", spline_dim, "\n")
    }
  } else {
    m = 0
    blocksize = 0
  }
  
  #mC_sp = Matrix(mult$mC, sparse=TRUE)

  zero.set = which(mult$vy == 0)
  nonzero.set = which(mult$vy != 0)
  vlower_bound <- c()
  
  i = 0
  # Iterate ----
  while ( (i <= 1) || is.nan(vlower_bound[i] - vlower_bound[i - 1]) || (vlower_bound[i] > vlower_bound[i-1])  ) {
  #while ( (i <= MAXITER)  ) {	
    if (i >= MAXITER) {
      cat("Iteration limit reached, breaking ...")
      break
    }
    
    i = i+1
    
    if (!is.null(mult$mSigma.u.inv)) {
      mult$mSigma.inv <- blockDiag(mult$mSigma.beta.inv,mult$mSigma.u.inv)
    } else {
      mult$mSigma.inv <- mult$mSigma.beta.inv
    }
    
    # Update parameter for q_vnu
    # Maximise the Gaussian Variational Approximation using Dr Ormerod's Poisson mixed model code
    if (method == "laplacian") {
      fit1 = fit.Lap(mult$vmu, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, mult$mLambda)
    } else if (method == "gva") {	
      #fit2 = fit.Lap(mult$vmu, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, mult$mLambda)
      fit1 = fit.GVA(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B")
    } else if (method == "gva2") {
      fit2 = fit.Lap(mult$vmu, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, mult$mLambda)
      fit1 = fit.GVA_new(fit2$vmu, fit2$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B", p=p, m=m, blocksize=mult$blocksize)
    } else if (method == "gva2new") {
      #fit2 = fit.Lap(mult$vmu, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, mult$mLambda)
      fit1 = fit.GVA_new2(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B", p=p, m=m, blocksize=mult$blocksize, spline_dim=spline_dim, mC_sp=mC_sp)
      #fit1 = fit.GVA_new2(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "BFGS", p=p, m=m)
    } else if (method == "gva_nr") {
      fit1 = fit.GVA_nr(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B", p=p, m=m, blocksize=mult$blocksize, spline_dim=spline_dim)
    } else {
      stop("method must be either laplacian, gva, gva2, gva2new or gva_nr")
    }
    
    mult$vmu = fit1$vmu
    mult$mLambda = fit1$mLambda
    mult$f = fit1$res$value
    
    # Update parameters for q_vr
    if (length(zero.set) != 0) {
      mult$vp[zero.set] = expit((mult$vy[zero.set]*mult$mC[zero.set,])%*%mult$vmu-exp(mult$mC[zero.set,]%*%mult$vmu + 0.5*diag((matrix(mult$mC[zero.set,], length(zero.set), p+(m-1)*blocksize+spline_dim))%*%mult$mLambda%*%t(matrix(mult$mC[zero.set,], length(zero.set), p+(m-1)*blocksize+spline_dim))) + digamma(mult$a_rho) - digamma(mult$b_rho)))
    }
    
    # Update parameters for q_rho
    mult$a_rho = mult$prior$a_rho + sum(mult$vp)
    mult$b_rho = mult$prior$b_rho + N - sum(mult$vp)
    
    # Update parameters for q_sigma_u^2 if we need to
    if (!is.null(mult$mZ)) {
      # a_sigma is fixed
      u_dim = (m-1)*blocksize+spline_dim
      mult$a_sigma = mult$prior$a_sigma + u_dim/2
      u_idx = p + 1:u_dim  # (ncol(mult$mX)+1):ncol(mult$mC)
      # u_idx = p + 1:(m*blocksize) + spline_dim
      #tr_mSigma = ncol(mult$mZ) * mult$prior$a_sigma/mult$prior$b_sigma
      #mult$b_sigma = mult$prior$b_sigma + sum(vu^2)/2 + (tr_mSigma)/2
      mult$b_sigma = mult$prior$b_sigma + sum(mult$vmu[u_idx]^2)/2 + tr(mult$mLambda[u_idx, u_idx])/2    # Extract right elements of mLambda
      
      # FIXME: This code is hard to read
      acc = matrix(0, blocksize, blocksize)
      for (j in 1:(m-1)) {
        j_idx = p + (j-1)*blocksize+(1:blocksize)
        acc = acc + with(mult, vmu[j_idx] %*% t(vmu[j_idx]) + mLambda[j_idx, j_idx])
      }
      mult$mPsi = with(mult, prior$mPsi + acc)
      
      #tau_sigma = mult$a_sigma/mult$b_sigma
      
      #mult$mSigma.u.inv = diag(tau_sigma, u_dim)
  
      mult$mSigma.u.inv = with(mult, kronecker(diag(1, m-1), solve(mPsi/(v - blocksize - 1))))
      
      #mult$mSigma.u.inv = with(mult, solve(mPsi)) # What multiplicative factor for psi?
    }
    
    vlower_bound[i] <- 0 # calculate_lower_bound(mult)
    vlower_bound[i] <- calculate_lower_bound(mult, verbose=verbose)
    
    if (verbose && i > 1)
      cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
          vlower_bound[i] - vlower_bound[i-1], " parameters ", "vmu", mult$vmu,
          "diag(mLambda)", diag(mult$mLambda), "a_rho", mult$a_rho, "b_rho", mult$b_rho)
    if (verbose && !is.null(mult$mZ)) {
      cat(" a_sigma", mult$a_sigma, "b_sigma", mult$b_sigma)
    }
    if (verbose) cat("\n")
  }
  
  if (plot_lower_bound)
    plot(vlower_bound,type="l")
  
  params = list(vmu=mult$vmu, mLambda=mult$mLambda, a_rho=mult$a_rho, b_rho=mult$b_rho,
                a_sigma=mult$a_sigma, b_sigma=mult$b_sigma, vlower_bound=vlower_bound)
  return(params)
}

zero_infl_var <- function(object, method="laplacian", verbose=FALSE, plot_lower_bound=FALSE)
{
  UseMethod("zero_infl_var", object)
}
