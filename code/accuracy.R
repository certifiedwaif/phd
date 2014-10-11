# accuracy.R
setwd("~/phd/code")
source("zero_inflated_model.R")
source("test_zero_inflated_model.R")
source("rwmh.R")

generate_test_data = function(m, ni)
{
  m = m
  ni = ni
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  #print("mX=")
  #print(mX)
  #cat("dim(mX)", dim(mX), "\n")
  
  #v = c(rep(1, g), rep(0, g))
  # Indicator variables for groups
  
  #mZ = matrix(cbind(v, 1-v), sum(n), 2)
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
  
  expected_rho = 0.5
  expected_beta = c(2, 1)
  expected_sigma2_u = .5^2
  a_sigma = 1e-2
  b_sigma = 1e-2
  
  tau = 1.0E2
  
  sigma2.beta <- 1.0E3
  
  test_data = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u)
  vy = test_data$vy
  
  # Test accuracy
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau)
  mult$n = n
  mult$m = m
  mult$ni = ni
  
  return(mult)
}

# This code runs too slow.
test_multivariate_accuracy <- function()
{
  m = 20
  ni = 5
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  #print("mX=")
  #print(mX)
  #cat("dim(mX)", dim(mX), "\n")
  
  #v = c(rep(1, g), rep(0, g))
  # Indicator variables for groups
  
  #mZ = matrix(cbind(v, 1-v), sum(n), 2)
  #mZ <- matrix(0,sum(n),m)
  #count <- 0
  #for (i in 1:m) {
  #	mZ[count + (1:n[i]),i] <- 1
  #	count <- count + n[i]
  #}
  
  mZ <- kronecker(diag(1,m),rep(1,ni))
  
  #print("mZ=")
  #print(mZ)
  #cat("dim(mZ)", dim(mZ), "\n")
  
  expected_rho = 0.5
  expected_beta = c(2, 1)
  expected_sigma2_u = .5^2
  a_sigma = 1e-2
  b_sigma = 1e-2
  
  tau = 1.0E2
  
  sigma2.beta <- 1.0E3
  
  test_data = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u)
  vy = test_data$vy
  
  # Test accuracy
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau)
  mcmc_result = mcmc(mult, iterations=1e5+2000, burnin=2000, thinning=1)
  
  par(mfrow=c(2, 1))
  hist(mcmc_result$vnu[1,])
  hist(mcmc_result$vnu[2,])
  par(mfrow=c(1, 1))
  print(summary(mcmc_result$vnu[1,]))
  print(summary(mcmc_result$vnu[2,]))
  var_result = zero_infl_var(mult, method="gva", verbose=TRUE)
  
  # Compare MCMC distribution with variational approximation for each parameter
  # vnu[i] ~ Normal, dnorm
  # sigma2_u ~ IG, dgamma(1/x)
  # rho ~ Beta, dbeta
  # vr[i] ~ Bernoulli, dbinom
  # For each parameter of interest,
  # * estimate density of MCMC
  # * compare with q distribution using L_1 norm
  
  # Kernel density estimates of MCMC-estimated posteriors
  # Use L_1 distance to compare against variational approximations of posteriors
  par(mfrow=c(2,1))
  # Pseudocode:
  # Estimate density
  # Define integrand, using distribution and parameters.
  # Calculate accuracy integral
  # Plots
  calculate_accuracy = function(mcmc_samples, dist_fn, param1, param2)
  {
    #cat("dist_fn", deparse(substitute(dist_fn)), "\n")
    #cat("param1", param1, "\n")
    #cat("param2", param2, "\n")
    mcmc_density = density(mcmc_samples)
    mcmc_fn = splinefun(mcmc_density$x, mcmc_density$y)
    
    integrand <- function(x)
    {
      return(abs(mcmc_fn(x) - dist_fn(x, param1, param2)))
    }
    #cat("min(mcmc_density$x)", min(mcmc_density$x), "\n")
    #cat("max(mcmc_density$x)", max(mcmc_density$x), "\n")
    result = integrate(integrand, min(mcmc_density$x), max(mcmc_density$x),
                       subdivisions = length(mcmc_density$x))
    accuracy = 1 - .5 * result$value
    return(accuracy)
  }
  # vnu accuracy
  for (t in 1:nrow(mcmc_result$vnu)) {
    accuracy = calculate_accuracy(mcmc_result$vnu[t,], dnorm,
                                  var_result$vmu[t], sqrt(var_result$mLambda[t,t]))
    cat("vnu[", t, "] accuracy: ", accuracy, "\n")
  }
  # sigma2_u accuracy
  accuracy = calculate_accuracy(1/mcmc_result$sigma2_u, dgamma,
                                var_result$a_sigma, var_result$b_sigma)
  cat("sigma2_u accuracy: ", accuracy, "\n")
  
  # rho accuracy
  accuracy = calculate_accuracy(mcmc_result$rho, dbeta,
                                var_result$a_rho, var_result$b_rho)
  cat("rho accuracy: ", accuracy, "\n")
  
  accuracy_plot = function(mcmc_samples, dist_fn, param1, param2)
  {
    mcmc_density = density(mcmc_samples)
    plot(mcmc_density)
    curve(dist_fn(x, param1, param2),
          from=min(mcmc_density$x), to=max(mcmc_density$x),
          add=TRUE, lty=2, col="blue")
  }
  par(mfrow=c(2,1))
  accuracy_plot(mcmc_result$vnu[t,], dnorm, 
                var_result$vmu[t], sqrt(var_result$mLambda[t,t]))
  plot(mcmc_result$vnu[t,], type="l")
  par(mfrow=c(1,1))
  t = 6
  plot(mcmc_result$vnu[t,1:(1e4-1)], mcmc_result$vnu[t,2:1e4])
  acf(mcmc_result$vnu[t,])
  pacf(mcmc_result$vnu[t,])
  
  par(mfrow=c(2,1))
  accuracy_plot(1/mcmc_result$sigma2_u, dgamma,
                var_result$a_sigma, var_result$b_sigma)
  plot(1/mcmc_result$sigma2_u, type="l")
  par(mfrow=c(1,1))
  plot(mcmc_result$sigma2_u[1:(1e4-1)], mcmc_result$sigma2_u[2:1e4])
  
  par(mfrow=c(2,1))
  accuracy_plot(mcmc_result$rho, dbeta,
                var_result$a_rho, var_result$b_rho)
  plot(mcmc_result$rho, type="l")
  plot(mcmc_result$rho[1:(1e4-1)], mcmc_result$rho[2:1e4])
  par(mfrow=c(1,1))
  
  par(mfrow=c(2,1))
  density_mcmc_sigma2_u_inv = density(1/mcmc_result$sigma2_u)
  plot(density_mcmc_sigma2_u_inv)
  curve(dgamma(x, result_var$a_sigma, result_var$b_sigma),
        from=min(density_mcmc_sigma2_u_inv$x), to=max(density_mcmc_sigma2_u_inv$x),
        add=TRUE, lty=2, col="blue")
  par(mfrow=c(1,1))
  
  par(mfrow=c(2,1))
  density_mcmc_rho = density(mcmc_result$rho)
  plot(density_mcmc_rho)
  curve(dbeta(x, result_var$a_rho, result_var$b_rho),
        from=min(density_mcmc_rho$x), to=max(density_mcmc_rho$x),
        add=TRUE, lty=2, col="blue")
  par(mfrow=c(1,1))
}

mcmc_approximation <- function(mult, iterations=1e3)
{
  #mcmc_result = mcmc(mult, iterations=1e5+2000, burnin=2000, thinning=1)
  # Use Stan to create MCMC samples, because Stan deals much better with highly
  # correlated posteriors.
  require(rstan)
  require(parallel)
  #source("multivariate_stan.R")
  
  zip_data <- with(mult, list(N=sum(n), P=2, M=m, y=vy, X=mX, Z=mZ))
  #print(str(zip_data))
  rng_seed <- 5;
  foo <- stan("multivariate_zip.stan", data=zip_data, chains = 0)
  sflist <- 
    mclapply(1:4, mc.cores = 1, 
             function(i) stan(fit=foo, data=zip_data, seed = rng_seed, 
                              chains = 1, chain_id = i, refresh = -1,
                              iter=iterations))
  fit <- sflist2stanfit(sflist)
  
  #fit <- stan(model_code = zip_code, data = zip_dat, 
  #            iter = 1e5, chains = 4)  
  
  mcmc_samples = extract(fit)
  return(mcmc_samples)
}

test_accuracy = function(mult, mcmc_samples, approximation)
{
  cat("approximation", approximation, "\n")
  pdf("accuracy_plots.pdf")
  var_result = zero_infl_var(mult, method=approximation)
  # vbeta accuracy
  calculate_accuracy3 = function(mcmc_samples, dist_fn, param1, param2)
  {
    mcmc_density = density(mcmc_samples)
    mcmc_fn = splinefun(mcmc_density$x, mcmc_density$y)
    
    integrand <- function(x)
    {
      return(abs(mcmc_fn(x) - dist_fn(x, param1, param2)))
    }
    #cat("min(mcmc_density$x)", min(mcmc_density$x), "\n")
    #cat("max(mcmc_density$x)", max(mcmc_density$x), "\n")
    result = integrate(integrand, min(mcmc_density$x), max(mcmc_density$x),
                       subdivisions = length(mcmc_density$x))
    accuracy = 1 - .5 * result$value
    return(accuracy)
  }
  
  # Compare MCMC distribution with variational approximation for each parameter
  # vnu[i] ~ Normal, dnorm
  # sigma2_u ~ IG, dgamma(1/x)
  # rho ~ Beta, dbeta
  # vr[i] ~ Bernoulli, dbinom
  # For each parameter of interest,
  # * estimate density of MCMC
  # * compare with q distribution using L_1 norm
  
  # Kernel density estimates of MCMC-estimated posteriors
  # Use L_1 distance to compare against variational approximations of posteriors
  par(mfrow=c(1,1))
  # Pseudocode:
  # Estimate density
  # Define integrand, using distribution and parameters.
  # Calculate accuracy integral
  # Plots
  
  accuracy_plot = function(mcmc_samples, dist_fn, param1, param2)
  {
    mcmc_density = density(mcmc_samples)
    plot(mcmc_density)
    curve(dist_fn(x, param1, param2),
          from=min(mcmc_density$x), to=max(mcmc_density$x),
          add=TRUE, lty=2, col="blue")
  }
  
  for (i in 1:ncol(mult$mX)) {
    accuracy = calculate_accuracy3(mcmc_samples$vbeta[,i], dnorm,
                                  var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
    cat("vbeta[", i, "]", approximation, "accuracy:", accuracy, "\n")
    
    par(mfrow=c(1,1))
    param_name = sprintf("vbeta[%d]", 1)
    pdf(sprintf("~/phd/code/vbeta_accuracy_%s.pdf", i))
    accuracy_plot(mcmc_samples$vbeta[,i], dnorm,
                  var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
    dev.off()
    #plot(mcmc_samples$vbeta[,i], type="l")
    par(mfrow=c(1,1))
  }
  
  # vu accuracy
  for (i in 1:ncol(mult$mZ)) {
    accuracy = calculate_accuracy3(mcmc_samples$u[,i], dnorm,
                                  var_result$vmu[i+2], sqrt(var_result$mLambda[i+2,i+2]))
    cat("vu[", i, "]", approximation, "accuracy:", accuracy, "\n")
    par(mfrow=c(1,1))
    pdf(sprintf("~/phd/code/vu_accuracy_%s.pdf", i))
    accuracy_plot(mcmc_samples$u[,i], dnorm,
                  var_result$vmu[i+2], sqrt(var_result$mLambda[i+2,i+2]))
    dev.off()
    #plot(mcmc_samples$u[,i], type="l")
    par(mfrow=c(1,1))
  }
  
  # sigma2_u accuracy
  accuracy = calculate_accuracy3(1/mcmc_samples$sigma_u^2, dgamma,
                                 var_result$a_sigma, var_result$b_sigma)
  cat("sigma2_u", approximation, "accuracy:", accuracy, "\n")
  par(mfrow=c(1,1))
  pdf("~/phd/code/sigma2_u_accuracy.pdf")
  accuracy_plot(1/mcmc_samples$sigma_u^2, dgamma,
                var_result$a_sigma, var_result$b_sigma)
  dev.off()
  #plot(mcmc_samples$sigma_u, type="l")
  par(mfrow=c(1,1))
  
  # rho accuracy
  accuracy = calculate_accuracy3(mcmc_samples$rho, dbeta,
                                 var_result$a_rho, var_result$b_rho)
  cat("rho", approximation, "accuracy: ", accuracy, "\n")
  par(mfrow=c(1,1))
  pdf("~/phd/code/rho_accuracy.pdf")
  accuracy_plot(mcmc_samples$rho, dbeta,
                var_result$a_rho, var_result$b_rho)
  dev.off()
  #plot(mcmc_samples$rho, type="l")
  par(mfrow=c(1,1))
  return(var_result)
}

# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
calculate_accuracy2 <- function(result_mcmc, result_var)
{
  density_mcmc_rho = density(result_mcmc$vrho)
  integrand <- function(x)
  {
    fn = splinefun(density_mcmc_rho$x, density_mcmc_rho$y)
    return(abs(fn(x) - dbeta(x, result_var$a_rho, result_var$b_rho)))
  }
  result = integrate(integrand, min(density_mcmc_rho$x), max(density_mcmc_rho$x), subdivisions = length(density_mcmc_rho$x))
  rho_accuracy = 1 - .5 * result$value
  
  density_mcmc_lambda = density(result_mcmc$vlambda)
  integrand <- function(x)
  {
    fn = splinefun(density_mcmc_lambda$x, density_mcmc_lambda$y)
    return(abs(fn(x) - dgamma(x, result_var$a_lambda, result_var$b_lambda)))
  }
  result = integrate(integrand, min(density_mcmc_lambda$x), max(density_mcmc_lambda$x), subdivisions = length(density_mcmc_lambda$x))
  lambda_accuracy = 1 - .5 * result$value
  return(list(rho_accuracy=rho_accuracy, lambda_accuracy=lambda_accuracy))
}

check_accuracy <- function(n, rho, lambda)
{
  x = generate_univariate_test_data(n, rho, lambda)
  
  a = 10
  b = 10
  
  # MCMC ----
  iterations = 1e6
  burnin = 1e3
  
  start = Sys.time()
  result_mcmc = mcmc.univariate(iterations+burnin, x, a, b)
  mcmc_runtime = Sys.time() - start
  # Throw away burn-in samples.
  # Brackets turn out to be incredibly important here!!!
  result_mcmc$rho = result_mcmc$rho[(burnin+1):(burnin+iterations+1)]
  result_mcmc$lambda = result_mcmc$lambda[(burnin+1):(burnin+iterations+1)]
  # 1000 iterations in 0.07461715 seconds.
  
  # Variational approximation ----
  start = Sys.time()
  univariate = create_univariate(x, a, b)
  result_var = zero_infl_var(univariate)
  var_runtime = Sys.time() - start
  # Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
  # or 200 iterations a second.
  var_lambda = result_var$a_lambda / result_var$b_lambda
  var_rho = result_var$a_rho / (result_var$a_rho + result_var$b_lambda)
  
  accuracy = calculate_accuracy(result_mcmc, result_var)
  return(list(n=n, rho=rho, lambda=lambda, accuracy=accuracy))
}

# Check accuracy of the approximation for a range of parameter values ----
main_check_accuracy <- function()
{
  n = 1000
  rho = .5
  lambda = 100
  
  sink("~/phd/code/univariate_accuracy_results.txt")
  for (rho in seq(.1,.9, 0.05))
    for (lambda in seq(0.1, 10, 0.05))
      print(check_accuracy(n, rho, lambda))
  sink()
}

#main_check_accuracy()

# Need to be able to compare the solution paths of each approximation

# Generate data
for (i in 1:100) {
  set.seed(i)
  mult = generate_test_data(10, 100)
  # Monte Carlo Markov Chains approximation
  mcmc_samples = mcmc_approximation(mult, iterations=1e3)
  # Save the results, because this takes such a long time to run.
}
save(mult, mcmc_samples, file="accuracy.RData")
# Test all other approximations against it

# Test multivariate approximation's accuracy
test_accuracy(mult, mcmc_samples, "laplacian")
test_accuracy(mult, mcmc_samples, "gva")
test_accuracy(mult, mcmc_samples, "gva2")
test_accuracy(mult, mcmc_samples, "gva_nr")
