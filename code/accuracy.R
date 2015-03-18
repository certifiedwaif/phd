# accuracy.R
source("zero_inflated_model.R")
source("generate.R")
source("mcmc.R")
source("rwmh.R")

calculate_accuracy <- function(mcmc_samples, dist_fn, param1, param2)
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

accuracy_plot <- function(mcmc_samples, dist_fn, param1, param2)
{
  mcmc_density = density(mcmc_samples)
  plot(mcmc_density)
  curve(dist_fn(x, param1, param2),
        from=min(mcmc_density$x), to=max(mcmc_density$x),
        add=TRUE, lty=2, col="blue")
}

calculate_accuracies <- function(mult, mcmc_samples, var_result, approximation, print_flag=FALSE, plot_flag=FALSE)
{
  # TODO: Add support for checking the accuracy over multiple dimensions
  # cubature$adaptIntegrate
  
  if (plot_flag) pdf(paste0("accuracy_plots_", approximation, ".pdf"))
  #return(var_result)
  # vbeta accuracy
  
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
 
  vbeta_accuracy = rep(NA, ncol(mult$mX))
  for (i in 1:ncol(mult$mX)) {
    vbeta_accuracy[i] = calculate_accuracy(mcmc_samples$vbeta[,i], dnorm,
                                            var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
    if (print_flag) cat("vbeta[", i, "]", approximation, "accuracy:", vbeta_accuracy[i], "\n")
    
    param_name = sprintf("vbeta[%d]", i)
    if (plot_flag) accuracy_plot(mcmc_samples$vbeta[,i], dnorm,
                            var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
  }
  
  # vu accuracy
  # FIXME: To check for random slopes accuracy, this section will have
  # to get more complex.
  print(dim(mult$mZ))
  print(dim(mcmc_samples$u))
  vu_accuracy = rep(NA, ncol(mult$mZ))
  B = mult$blocksize
  b_idx = 1
  for (i in 1:ncol(mult$mZ)) {
    m_idx = ceiling(i/B)
    vu_accuracy[i] = calculate_accuracy(mcmc_samples$vu[,m_idx,b_idx], dnorm,
                                         var_result$vmu[i+mult$p], sqrt(var_result$mLambda[i+mult$p,i+mult$p]))
    if (print_flag) cat("vu[", i, "]", approximation, "accuracy:", vu_accuracy[i], "\n")
    if (plot_flag) accuracy_plot(mcmc_samples$vu[,m_idx,b_idx], dnorm,
                            var_result$vmu[i+mult$p], sqrt(var_result$mLambda[i+mult$p,i+mult$p]))

    b_idx = b_idx + 1
    if (b_idx > B)
      b_idx=1
  }
  
  # sigma2_u accuracy
  # FIXME - this may be wrong for blocksize != 1?
  # This is totally wrong for the Inverse Wishart model?
  #sigma2_u_accuracy = calculate_accuracy(1/mcmc_samples$sigma_u^2, dgamma,
  #                                        var_result$a_sigma, var_result$b_sigma)
  #if (print_flag) cat("sigma2_u", approximation, "accuracy:", sigma2_u_accuracy, "\n")
  #if (plot_flag) accuracy_plot(1/mcmc_samples$sigma_u^2, dgamma,
  #                        var_result$a_sigma, var_result$b_sigma)
  
  # rho accuracy
  rho_accuracy = calculate_accuracy(mcmc_samples$rho, dbeta,
                                     var_result$a_rho, var_result$b_rho)
  if (print_flag) cat("rho", approximation, "accuracy: ", rho_accuracy, "\n")
  if (plot_flag) accuracy_plot(mcmc_samples$rho, dbeta,
                          var_result$a_rho, var_result$b_rho)
  if (plot_flag) dev.off()
  return(list(var_result=var_result,
              vbeta_accuracy=vbeta_accuracy,
              vu_accuracy=vu_accuracy,
              #sigma2_u_accuracy=sigma2_u_accuracy,
              rho_accuracy=rho_accuracy))
}

test_accuracy <- function(mult, mcmc_samples, approximation, plot=FALSE)
{
  var_result = zero_infl_var(mult, method=approximation, verbose=TRUE)
  return(calculate_accuracies(mult, mcmc_samples, var_result, approximation, plot_flag=plot))
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

test_accuracies <- function()
{
  # Need to be able to compare the solution paths of each approximation
  
  # Generate data
  # for (i in 1:100) {
  #   set.seed(i)
  #   mult = generate_test_data(20, 100)
  #   # Monte Carlo Markov Chains approximation
  #   mcmc_samples = mcmc_approximation(mult, iterations=1e6)
  #   # Save the results, because this takes such a long time to run.
  # }
  # save(mult, mcmc_samples, file="accuracy_good.RData")
  #set.seed(1)
  #mult = generate_test_data(10, 100)
  # Monte Carlo Markov Chains approximation
  #mcmc_samples = mcmc_approximation(mult, iterations=1e6, warmup = 1e4)
#   # Save the results, because this takes such a long time to run.
#   #save(mult, mcmc_samples, file="accuracy.RData")
  #save(mult, mcmc_samples, file="data/accuracy_int.RData")
  load(file="data/accuracy_int_2015_02_17.RData")
  #mult$spline_dim = 0
  #load(file="accuracy.RData")
  # Test all other approximations against it
  #load(file="accuracy.RData")
  
  # Test multivariate approximation's accuracy
  now = Sys.time()
  var1 = test_accuracy(mult, mcmc_samples, "laplacian")
  print(Sys.time() - now)
  #print(image(Matrix(var1$var_result$mLambda)))
  print(var1)
  
  now = Sys.time()
  var2 = test_accuracy(mult, mcmc_samples, "gva")
  print(Sys.time() - now)
  #print(image(Matrix(var2$var_result$mLambda)))
  print(var2)
  
  now = Sys.time()
  var3 = test_accuracy(mult, mcmc_samples, "gva2", plot=TRUE)
  print(Sys.time() - now)
  #print(image(Matrix(var3$mLambda)))
  print(var3)
  
  #Rprof()
  #now = Sys.time()
  #var4 = test_accuracy(mult, mcmc_samples, "gva2new")
  #print(Sys.time() - now)
  #print(image(Matrix(var4$var_result$mLambda)))
  #print(var4)
  
  #Rprof(NULL)
  #print(summaryRprof())
  #print(image(Matrix(var3_new$mLambda)))
  
  now = Sys.time()
  var5 = test_accuracy(mult, mcmc_samples, "gva_nr")
  print(Sys.time() - now)
  #print(image(Matrix(var4$result_var$mLambda)))
  print(var5)
  
  #save(var1, var2, var3, var4, var5, file="accuracy_results_int.RData")
  #for (i in 1:100) {
  #  set.seed(i)
  #  mult = generate_test_data(20, 100)
  #  mcmc_samples = mcmc_approximation(mult, iterations=1e4)
  #  
  #  var1 = test_accuracy(mult, mcmc_samples, "laplacian")
  #  var2 = test_accuracy(mult, mcmc_samples, "gva")
  #  var3 = test_accuracy(mult, mcmc_samples, "gva2")
  #  var4 = test_accuracy(mult, mcmc_samples, "gva_nr")
  #}
  
}
#test_accuracies()

test_accuracies_slope <- function()
{
  # Monte Carlo Markov Chains approximation
  #seed = 1
  #set.seed(seed)
  #mult = generate_slope_test_data()
  #mcmc_samples = mcmc_approximation(mult, seed=seed, iterations=1e5, warmup=1e4)
  #save(mult, mcmc_samples, file="data/accuracy_slope_2015_03_17.RData")  
  load(file="data/accuracy_slope_2015_03_03.RData")
  
  now = Sys.time()
  var1 = test_accuracy(mult, mcmc_samples, "laplacian", plot=TRUE)
  print(Sys.time() - now)
  print(var1)
  
  now = Sys.time()
  var2 = test_accuracy(mult, mcmc_samples, "gva", plot=TRUE)
  print(Sys.time() - now)
  #print(image(Matrix(var2$var_result$mLambda)))
  print(var2)
  
  #browser()
  now = Sys.time()
  #mult$mLambda = diag(1, ncol(mult$mC))
  var3 = test_accuracy(mult, mcmc_samples, "gva2", plot=TRUE)
  #print(image(Matrix(var3$var_result$mLambda)))
  print(Sys.time() - now)
  print(var3)

  #now = Sys.time()
  #Rprof(line.profiling=TRUE)
  #var4 = test_accuracy(mult, mcmc_samples, "gva2new", plot=TRUE)
  #Rprof(NULL)
  #print(summaryRprof(lines = "both"))
  #print(image(Matrix(var4$var_result$mLambda)))  
  #var4 = test_accuracy(mult, mcmc_samples, "gva2new")
  #print(Sys.time() - now)
  #print(var4)
  
  now = Sys.time()
  # GVA NR is unstable, and sometimes fails with an error
  #var5 = test_accuracy(mult, mcmc_samples, "gva_nr")
  var5 = test_accuracy(mult, mcmc_samples, "gva_nr", plot=TRUE)
  print(Sys.time() - now)
  #print(image(Matrix(var5$var_result$mLambda)))  
  print(var5)
  
  browser()
  
  #save(var1, var2, var3, var4, var5, file="accuracy_results_slope.RData")
  
  #print(image(Matrix(var3_new$mLambda)))
  
  #now = Sys.time()
  #var4 = test_accuracy(mult, mcmc_samples, "gva_nr")
  #Sys.time() - now
  #print(image(Matrix(var4$mLambda)))

#   mult = generate_spline_test_data()
#   mcmc_samples = mcmc_approximation(mult, iterations=1e4)
#   save(mult, mcmc_samples, file="accuracy_spline.RData")    

  #now = Sys.time()
  #var1 = test_accuracy(mult, mcmc_samples, "laplacian")
  #Sys.time() - now
  #print(image(Matrix(var1$var_result$mLambda)))
  #var1
  
  #now = Sys.time()
  #var2 = test_accuracy(mult, mcmc_samples, "gva")
  #Sys.time() - now
  #print(image(Matrix(var2$var_result$mLambda)))
  #var2
  
  #now = Sys.time()
  #var3 = test_accuracy(mult, mcmc_samples, "gva2")
  #Sys.time() - now
  #print(image(Matrix(var3$var_result$mLambda)))
  #var3
  
  #now = Sys.time()
  #var4 = test_accuracy(mult, mcmc_samples, "gva2new")
  #Sys.time() - now
  #print(image(Matrix(var4$var_result$mLambda)))  
  #var4
  
  #Rprof()
  #now = Sys.time()
  #var5 = test_accuracy(mult, mcmc_samples, "gva_nr")
  #print(Sys.time() - now)
  #print(image(Matrix(var5$var_result$mLambda)))  
  #var5
}
test_accuracies_slope()
