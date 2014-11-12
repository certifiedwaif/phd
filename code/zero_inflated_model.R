# variational_approximation_to_zero_inflated_model.R
source("common.R")
source("zero_inflated_poisson_model.R")

mcmc.univariate <- function(iterations, vx, a, b)
{
  # Initialise
  n = length(vx)
  zero.set = which(vx == 0)
  nonzero.set = which(vx != 0)
  vlambda = rep(NA, iterations+1)
  vrho = rep(NA, iterations+1)
  vlambda[1] = mean(vx)
  vrho[1] = length(zero.set)/length(vx)
  veta = rep(NA, n)
  vr = rep(NA, n)
  # Build up a list of the zero observations. We'll only generate
  # r[i]'s for those. The non-zero observations will always have r[i] = 1.
  vr[nonzero.set] = 1
  # Iterate ----
  for (i in 1:iterations) {
    varg = exp(-vlambda[i] + logit(vrho[i]))
    veta = varg/(1 + varg)
    
    # The following code is correct, but 100 times slower as R is interpreted.
    # So we use the vectorised code instead.
    #for (j in 1:length(zero_observation_idx)) {
    #	r[zero_observation_idx[j]] = rbinom(1, 1, eta)
    #}
    
    vr[zero.set] = rbinom(length(zero.set), 1, veta) 
    vlambda[i+1] = rgamma(1, a + sum(vx), b + sum(vr))
    vrho[i+1] = rbeta(1, sum(vr) + 1, n - sum(vr) + 1)
  }
  return(list(vlambda=vlambda, vrho=vrho))
}

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

calculate_lower_bound.univariate <- function(univariate)
{
  vx = univariate$vx
  vp = univariate$vp
  a_lambda = univariate$a_lambda
  b_lambda = univariate$b_lambda
  a_rho = univariate$a_rho
  b_rho = univariate$b_rho
  
  zero.set <- which(vx==0)
  
  E_lambda = a_lambda/b_lambda
  E_log_lambda = digamma(a_lambda) - log(b_lambda)	
  
  E_r = ifelse(vx == 0, vp, 1)
  E_xi_log_lambda_r = ifelse(vx == 0, 0, vx*E_log_lambda)
  
  E_log_rho = digamma(a_rho) - digamma(a_rho + b_rho)
  E_log_one_minus_rho = digamma(b_rho) - digamma(a_rho + b_rho)
  
  E_log_q_r = (vp[zero.set]*log(vp[zero.set]) + (1-vp[zero.set])*log(1-vp[zero.set]))
  E_log_q_lambda = -gamma_entropy(a_lambda, b_lambda)
  E_log_q_rho = -beta_entropy(a_rho, b_rho) # FIXME: Why is this entropy positive?
  
  result = a_lambda * log(b_lambda) + (a_lambda-1) * E_log_lambda - b_lambda * E_lambda - lgamma(a_lambda)
  result = result - E_lambda * sum(E_r)
  result = result + sum(E_xi_log_lambda_r) - sum(lgamma(vx+1))
  result = result + sum(E_r) * E_log_rho + sum(1 - E_r) * E_log_one_minus_rho
  result = result - sum(E_log_q_r) - E_log_q_lambda - E_log_q_rho
  
  return(result)
}

calculate_lower_bound.multivariate <- function(multivariate)
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
      T1 = T1 + prior$a_sigma * log (prior$b_sigma) - lgamma(prior$a_sigma)
      T1 = T1 - a_sigma * log(b_sigma) + lgamma(a_sigma)
    }
    T1 = T1 + .5*(p+m) + .5*log(det(mLambda))
    
    # Something is wrong in T2. It sometimes goes backwards as we're optimised.
    # This should be unchanged from the univariate lower bound
    #cat("calculate_lower_bound: ", vp[zero.set], "\n")
    # 0 log 0 returns NaN. We sidestep this by adding epsilon to vp[zero.set]
    eps = 1e-10
    vp[zero.set] = vp[zero.set] + eps
    T2 = sum(-vp[zero.set]*log(vp[zero.set]) - (1-vp[zero.set])*log(1-vp[zero.set]))
    T2 = T2 - lbeta(prior$a_rho, prior$b_rho) + lbeta(a_rho, b_rho)
    
    cat("calculate_lower_bound: T1", T1, "T2", T2, "\n")
    result = T1 + T2
    result
  })
  
  return(result)
}

calculate_lower_bound <- function(object)
{
  UseMethod("calculate_lower_bound", object)
}

create_univariate <- function(vx, a, b)
{
  # Initialise
  n = length(vx)
  vp = rep(1, n)
  a_lambda = a + sum(vx)
  b_lambda = b
  univariate = list(vx=vx, vp=vp, a_lambda=a_lambda, b_lambda=b_lambda)
  class(univariate) = "univariate"
  return(univariate)
}

zero_infl_var.univariate <- function(univariate, verbose=FALSE, plot_lower_bound=FALSE)
{
  vx = univariate$vx
  a = univariate$a_lambda
  b = univariate$b_lambda
  n = length(univariate$vx)
  zero.set = which(univariate$vx == 0)
  nonzero.set = which(univariate$vx != 0)
  vlower_bound <- c()
  
  i = 0
  # Iterate ----
  while (i <= 2 || vlower_bound[i] > vlower_bound[i-1]) {
    i = i+1
    
    # Update parameter for q_lambda
    univariate$b_lambda = b + sum(univariate$vp)
    
    # Update parameters for q_rho
    univariate$a_rho = 1 + sum(univariate$vp)
    univariate$b_rho = n - sum(univariate$vp) + 1
    
    # Update parameters for q_vr
    univariate$vp[zero.set] = expit(-expected_lambda(univariate) + digamma(univariate$a_rho) - digamma(univariate$b_rho))
    
    vlower_bound[i] <- calculate_lower_bound(univariate)
    
    if (verbose && i > 1)
      cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
          vlower_bound[i] - vlower_bound[i-1], " parameters ", "a_lambda", univariate$a_lambda,
          "b_lambda", univariate$b_lambda, "a_rho", univariate$a_rho, "b_rho", univariate$b_rho, "\n")
  }
  
  if (plot_lower_bound)
    plot(lower_bound_vector,type="l")
  
  params = list(a_lambda=univariate$a_lambda, b_lambda=univariate$b_lambda, a_rho=univariate$a_rho, b_rho=univariate$b_rho)
  return(params)
}

create_multivariate <- function(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau)
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
  
  mSigma.beta.inv = diag(1/sigma2.beta, ncol(mX))
  if (!is.null(ncol(mZ))) {
    mSigma.u.inv = diag(tau, ncol(mZ))
  } else {
    mSigma.u.inv = NULL
  }
  mLambda = diag(rep(1, ncol(mC)))
  a_rho = 1 + sum(vp)
  b_rho = n - sum(vp) + 1
  
  prior = list(a_sigma=a_sigma, b_sigma=b_sigma, a_rho=1, b_rho=1, sigma2.beta=sigma2.beta)
  multivariate = list(vy=vy, mX=mX, mZ=mZ, mC=mC, vp=vp, vmu=vmu,
                      a_sigma=a_sigma, b_sigma=b_sigma,
                      a_rho=a_rho, b_rho=b_rho,
                      mLambda=mLambda,
                      prior=prior,
                      mSigma.beta.inv=mSigma.beta.inv,
                      mSigma.u.inv=mSigma.u.inv,
                      mSigma.beta=solve(mSigma.beta.inv),
                      mSigma.vu=mSigma.u.inv)
  class(multivariate) = "multivariate"
  return(multivariate)
}

library(limma)

zero_infl_var.multivariate <- function(mult, method="gva", verbose=FALSE, plot_lower_bound=FALSE)
{
  MAXITER <- 20
  
  # Initialise
  N = length(mult$vy)
  if (verbose) cat("N", N)
  if (!is.null(mult$mX)) {
    p = ncol(mult$mX) 
    if (verbose) cat("p", p)
  }	else {
    p = 0
  }
  if (!is.null(mult$mZ)) {
    m = ncol(mult$mZ) 
    if (verbose) cat("m", m)
  } else {
    m = 0
  }
  cat("\n")
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
      fit1 = fit.GVA(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B")
    } else if (method == "gva2")
    {
      fit1 = fit.GVA_new(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B")
    } else if (method == "gva_nr") {
      fit1 = fit.GVA_nr(mult$vmu, mult$mLambda, mult$vy, mult$vp, mult$mC, mult$mSigma.inv, "L-BFGS-B")
    } else {
      stop("method must be either laplacian, gva, gva2 or gva_nr")
    }
    
    mult$vmu = fit1$vmu
    mult$mLambda = fit1$mLambda
    mult$f = fit1$res$value
    
    # Update parameters for q_vr
    if (length(zero.set) != 0) {
      mult$vp[zero.set] = expit((mult$vy[zero.set]*mult$mC[zero.set,])%*%mult$vmu-exp(mult$mC[zero.set,]%*%mult$vmu + 0.5*diag((matrix(mult$mC[zero.set,], length(zero.set), p+m))%*%mult$mLambda%*%t(matrix(mult$mC[zero.set,], length(zero.set), p+m))) + digamma(mult$a_rho) - digamma(mult$b_rho)))
    }
    
    # Update parameters for q_rho
    mult$a_rho = mult$prior$a_rho + sum(mult$vp)
    mult$b_rho = mult$prior$b_rho + N - sum(mult$vp)
    
    # Update parameters for q_sigma_u^2 if we need to
    if (!is.null(mult$mZ)) {
      # a_sigma is fixed
      mult$a_sigma = mult$prior$a_sigma + m/2
      u_idx = p + (1:m)  # (ncol(mult$mX)+1):ncol(mult$mC)
      #tr_mSigma = ncol(mult$mZ) * mult$prior$a_sigma/mult$prior$b_sigma
      #mult$b_sigma = mult$prior$b_sigma + sum(vu^2)/2 + (tr_mSigma)/2
      mult$b_sigma = mult$prior$b_sigma + sum(mult$vmu[u_idx]^2)/2 + tr(mult$mLambda[u_idx, u_idx])/2    # Extract right elements of mLambda
      
      tau_sigma = mult$a_sigma/mult$b_sigma
      
      mult$mSigma.u.inv = diag(tau_sigma, m)	
    }
    
    vlower_bound[i] <- 0 # calculate_lower_bound(mult)
    vlower_bound[i] <- calculate_lower_bound(mult)
    
    if (verbose && i > 1)
      cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
          vlower_bound[i] - vlower_bound[i-1], " parameters ", "vmu", mult$vmu,
          "a_rho", mult$a_rho, "b_rho", mult$b_rho)
    if (!is.null(mult$mZ)) {
      cat(" a_sigma", mult$a_sigma, "b_sigma", mult$b_sigma)
    }
    cat("\n")
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