# mcmc.R
library(rstan)
library(parallel)
# Andrew Gelman says that this magical line of code automatically makes Stan
# run in parallel and cache compiled models.
# Two minutes later: Hey, it actually works!
#source("http://mc-stan.org/rstan/stan.R")

mcmc_approximation <- function(mult, seed=1, iterations=NA, warmup=NA, mc.cores=1,
															 stan_file="multivariate_zip.stan")
{
  # Use Stan to create MCMC samples, because Stan deals much better with highly
  # correlated posteriors.
  
  m <- mult$m
  blocksize <- mult$blocksize
  spline_dim <- mult$spline_dim
  mX <- mult$mX
  mZ <- mult$mZ
  vy <- mult$vy
  mPsi <- mult$prior$mPsi
  mSigma.beta <- mult$mSigma.beta
  
  zip_data <- list(N=length(vy), P=2, M=m, B=blocksize, spline_dim=spline_dim,
  								 y=vy, X=mX, Z=mZ,
                   psi=mPsi, BetaPrior=mSigma.beta) #, v=prior$v))

	fit <- stan(stan_file, seed=seed, data=zip_data, iter=iterations, warmup=warmup, chains=1)
  mcmc_samples <- extract(fit)

  return(mcmc_samples)
}
