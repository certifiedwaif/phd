# mcmc.R
mcmc_approximation <- function(mult, iterations=1e3, warmup=floor(iterations/2), mc.cores = 1)
{
  #mcmc_result = mcmc(mult, iterations=1e5+2000, burnin=2000, thinning=1)
  # Use Stan to create MCMC samples, because Stan deals much better with highly
  # correlated posteriors.
  require(rstan)
  require(parallel)
  #source("multivariate_stan.R")
  
  #browser()
  mS = matrix(c(1, 0, 0, 1), 2, 2)
  u_dim = with(mult, m*blocksize+spline_dim)
  zip_data <- with(mult, list(N=length(vy), P=2, M=m, B=blocksize, y=vy, X=mX, Z=mZ,
                              S=mS, BetaPrior=mSigma.beta))
  #print(str(zip_data))
  rng_seed <- 5;
  fit <- stan("multivariate_zip.stan", data=zip_data, iter=iterations, warmup=warmup, chains = 1)
  #sflist <- 
  #  mclapply(1:4, mc.cores = mc.cores, 
  #           function(i) stan(fit=foo, data=zip_data, seed = rng_seed, 
  #                            chains = 1, chain_id = i, refresh = -1,
  #                            iter=iterations))
  #fit <- sflist2stanfit(sflist)
  
  #fit <- stan(model_code = zip_code, data = zip_dat, 
  #            iter = 1e5, chains = 4)  
  
  mcmc_samples = extract(fit)
  return(mcmc_samples)
}
