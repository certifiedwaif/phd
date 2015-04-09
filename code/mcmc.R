# mcmc.R
library(rstan)
library(parallel)
# Andrew Gelman says that this magical line of code automatically makes Stan
# run in parallel and cache compiled models.
# Two minutes later: Hey, it actually works!
#source("http://mc-stan.org/rstan/stan.R")

mcmc_approximation <- function(mult, seed=1, iterations=1e3, warmup=floor(iterations/2), mc.cores=1)
{
  # Use Stan to create MCMC samples, because Stan deals much better with highly
  # correlated posteriors.
  
  u_dim <- with(mult, m*blocksize+spline_dim)
  zip_data <- with(mult, list(N=length(vy), P=2, M=m, B=blocksize, y=vy, X=mX, Z=mZ,
                              psi=prior$mPsi, BetaPrior=mSigma.beta)) #, v=prior$v))
	fit <- stan("multivariate_zip.stan", seed=seed, data=zip_data, iter=iterations, warmup=warmup, chains = 1)
  mcmc_samples <- extract(fit)

  return(mcmc_samples)
}
