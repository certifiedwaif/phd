# accuracy.R
setwd("~/phd/code")
source("zero_inflated_model.R")
source("test_zero_inflated_model.R")
source("rwmh.R")

generate_int_test_data = function(m, ni)
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
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau, m=m, blocksize=1, spline_dim=0)
  
  return(mult)
}

generate_slope_test_data = function() {
  m = 20
  ni = 10
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  mZ = makeZ(mX, m, ni, p=2)
  
  expected_rho = 0.5
  expected_beta = c(2, 1)
  expected_sigma2_u = .5^2
  a_sigma = 1e-2
  b_sigma = 1e-2
  
  tau = 1.0E2
  
  sigma2.beta <- 1.0E3
  
  test_data = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=TRUE)
  vy = test_data$vy
  
  # Test model fitting
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau, m=m, blocksize=2, spline_dim=0)
  return(mult)
}

generate_spline_test_data = function()
{
  n = 5000
  vx = matrix(sort(runif(n, -1, 1))) 
  
  mX = cbind(1,vx)
  
  expected_rho = 1
  #expected_mu = c(0, 1)
  expected_sigma2_u = 0
  sigma2.beta = 1e5
  a_sigma = 1e5
  b_sigma = 1e5
  tau = 1.0E-5
  
  sigma2.true = 0.01
  expected_beta = c(0, 1)
  vf = 5+2*sin(pi*vx)
  vy = rpois(n,exp(vf))
  
  source("ZOsull.r")
  numIntKnots <- 35
  intKnots <- quantile(unique(vx),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
  
  mZ = ZOSull(vx,range.x=c(-1.1,1.1),intKnots=intKnots,drv=0)
  #vy = 2+mX[,1]^3+rnorm(m)*.1
  #result = fit_spline(vx, vy)
  #result = fit_spline(mX[,1], vy)
  #mZ = result$Z
  
  #mZ <- mZ/max(mZ)
  
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau, m=0, blocksize=1, spline_dim=37)
  
  return(mult)
}

calculate_accuracy = function(mcmc_samples, var_result, print_flag=FALSE, plot_flag=FALSE)
{
  # TODO: Add support for checking the accuracy over multiple dimensions
  # cubature$adaptIntegrate
  
  if (plot_flag) pdf(paste0("accuracy_plots_", approximation, ".pdf"))
  #return(var_result)
  # vbeta accuracy
  calculate_accuracy3 = function(mcmc_samples, dist_fn, param1, param2)
  {
    mcmc_density = density(mcmc_samples)
    mcmc_fn = splinefun(mcmc_density$x, mcmc_density$y)
    
    integrand <- function(x)
    {
      return(abs(mcmc_fn(x) - dist_fn(x, param1, param2)))
    }
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
  
  accuracy_plot = function(mcmc_samples, dist_fn, param1, param2)
  {
    mcmc_density = density(mcmc_samples)
    plot(mcmc_density)
    curve(dist_fn(x, param1, param2),
          from=min(mcmc_density$x), to=max(mcmc_density$x),
          add=TRUE, lty=2, col="blue")
  }
  
  vbeta_accuracy = rep(NA, ncol(mult$mX))
  for (i in 1:ncol(mult$mX)) {
    vbeta_accuracy[i] = calculate_accuracy3(mcmc_samples$vbeta[,i], dnorm,
                                            var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
    if (print_flag) cat("vbeta[", i, "]", approximation, "accuracy:", vbeta_accuracy[i], "\n")
    
    param_name = sprintf("vbeta[%d]", 1)
    if (plot_flag) accuracy_plot(mcmc_samples$vbeta[,i], dnorm,
                            var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
  }
  
  # vu accuracy
  # FIXME: To check for random slopes accuracy, this section will have
  # to get more complex. Or not, if John's right.
  vu_accuracy = rep(NA, ncol(mult$mZ))
  for (i in 1:ncol(mult$mZ)) {
    vu_accuracy[i] = calculate_accuracy3(mcmc_samples$u[,i], dnorm,
                                         var_result$vmu[i+2], sqrt(var_result$mLambda[i+2,i+2]))
    if (print_flag) cat("vu[", i, "]", approximation, "accuracy:", vu_accuracy[i], "\n")
    if (plot_flag) accuracy_plot(mcmc_samples$u[,i], dnorm,
                            var_result$vmu[i+2], sqrt(var_result$mLambda[i+2,i+2]))
  }
  
  # sigma2_u accuracy
  sigma2_u_accuracy = calculate_accuracy3(1/mcmc_samples$sigma_u^2, dgamma,
                                          var_result$a_sigma, var_result$b_sigma)
  if (print_flag) cat("sigma2_u", approximation, "accuracy:", sigma2_u_accuracy, "\n")
  if (plot_flag) accuracy_plot(1/mcmc_samples$sigma_u^2, dgamma,
                          var_result$a_sigma, var_result$b_sigma)
  
  # rho accuracy
  rho_accuracy = calculate_accuracy3(mcmc_samples$rho, dbeta,
                                     var_result$a_rho, var_result$b_rho)
  if (print_flag) cat("rho", approximation, "accuracy: ", rho_accuracy, "\n")
  if (plot_flag) accuracy_plot(mcmc_samples$rho, dbeta,
                          var_result$a_rho, var_result$b_rho)
  if (plot_flag) dev.off()
  return(list(var_result=var_result,
              vbeta_accuracy=vbeta_accuracy,
              vu_accuracy=vu_accuracy,
              sigma2_u_accuracy=sigma2_u_accuracy,
              rho_accuracy=rho_accuracy))
}

test_accuracy = function(mult, mcmc_samples, approximation, plot=FALSE)
{
  var_result = zero_infl_var(mult, method=approximation, verbose=TRUE)
  return(calculate_accuracy(mcmc_samples, var_result))
}

# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
calculate_univariate_accuracy <- function(result_mcmc, var_result)
{
  density_mcmc_rho = density(result_mcmc$vrho)
  integrand <- function(x)
  {
    fn = splinefun(density_mcmc_rho$x, density_mcmc_rho$y)
    return(abs(fn(x) - dbeta(x, var_result$a_rho, var_result$b_rho)))
  }
  result = integrate(integrand, min(density_mcmc_rho$x), max(density_mcmc_rho$x), subdivisions = length(density_mcmc_rho$x))
  rho_accuracy = 1 - .5 * result$value
  
  density_mcmc_lambda = density(result_mcmc$vlambda)
  integrand <- function(x)
  {
    fn = splinefun(density_mcmc_lambda$x, density_mcmc_lambda$y)
    return(abs(fn(x) - dgamma(x, var_result$a_lambda, var_result$b_lambda)))
  }
  result = integrate(integrand, min(density_mcmc_lambda$x), max(density_mcmc_lambda$x), subdivisions = length(density_mcmc_lambda$x))
  lambda_accuracy = 1 - .5 * result$value
  return(list(rho_accuracy=rho_accuracy, lambda_accuracy=lambda_accuracy))
}

test_accuracies = function()
{
  # Need to be able to compare the solution paths of each approximation
  
  # Generate data
  #for (i in 1:100) {
  #   set.seed(i)
  #   mult = generate_int_test_data(20, 100)
  #   # Monte Carlo Markov Chains approximation
  #   mcmc_samples = mcmc_approximation(mult, iterations=1e6)
  #   # Save the results, because this takes such a long time to run.
  # }
  # save(mult, mcmc_samples, file="accuracy_good.RData")
  set.seed(1)
#   mult = generate_test_data(10, 100)
#   # Monte Carlo Markov Chains approximation
#   mcmc_samples = mcmc_approximation(mult, iterations=1e4)
#   # Save the results, because this takes such a long time to run.
#   #save(mult, mcmc_samples, file="accuracy.RData")
#   save(mult, mcmc_samples, file="accuracy_int.RData")
  load(file="accuracy_int.RData")
  #load(file="accuracy.RData")
  # Test all other approximations against it
  #load(file="accuracy.RData")
  
  # Test multivariate approximation's accuracy
}
