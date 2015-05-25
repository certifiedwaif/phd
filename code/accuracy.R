# accuracy.R
source("zero_inflated_model.R")
source("generate.R")
library(rstan)

# Andrew Gelman says that this magical line of code automatically makes Stan
# run in parallel and cache compiled models.
# Two minutes later: Hey, it actually works!
# Some time later: Sometimes it works.
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
  v <- mult$prior$v
  mSigma.beta <- mult$mSigma.beta
  
  zip_data <- list(N=length(vy), P=2, M=m, B=blocksize, spline_dim=spline_dim,
                   y=vy, X=mX, Z=as.matrix(mZ),
                   v=v, psi=mPsi, BetaPrior=mSigma.beta)

  fit <- stan(stan_file, seed=seed, data=zip_data, iter=iterations, warmup=warmup, chains=1)
  mcmc_samples <- extract(fit)

  return(list(fit=fit,
              mcmc_samples=mcmc_samples))
}

calculate_accuracy <- function(mcmc_samples, dist_fn, param1, param2)
{
  mcmc_density <- density(mcmc_samples)
  mcmc_fn <- splinefun(mcmc_density$x, mcmc_density$y)
  
  integrand <- function(x)
  {
    return(abs(mcmc_fn(x) - dist_fn(x, param1, param2)))
  }
  result <- integrate(integrand, min(mcmc_density$x), max(mcmc_density$x),
                     subdivisions = length(mcmc_density$x))
  accuracy <- 1 - .5 * result$value
  return(accuracy)
}

accuracy_plot <- function(mcmc_samples, dist_fn, param1, param2)
{
  mcmc_density <- density(mcmc_samples)
  plot(mcmc_density)
  curve(dist_fn(x, param1, param2),
        from=min(mcmc_density$x), to=max(mcmc_density$x),
        add=TRUE, lty=2, col="blue")
}

calculate_accuracies <- function(mult, mcmc_samples, var_result, approximation, print_flag=FALSE, plot_flag=FALSE)
{
  # TODO: Add support for checking the accuracy over multiple dimensions
  # cubature$adaptIntegrate
  
  if (plot_flag) pdf(paste0("results/accuracy_plots_", approximation, ".pdf"))
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
 
  if (mult$p > 0) {
    vbeta_accuracy <- rep(NA, ncol(mult$mX))
    vbeta_means <- rep(NA, ncol(mult$mX))
    for (i in 1:ncol(mult$mX)) {
      vbeta_accuracy[i] <- calculate_accuracy(mcmc_samples$vbeta[,i], dnorm,
                                              var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
      vbeta_means[i] <- mean(mcmc_samples$vbeta[, i])
      if (print_flag) cat("vbeta[", i, "]", approximation, "accuracy:", vbeta_accuracy[i], "\n")
      
      param_name <- sprintf("vbeta[%d]", i)
      if (plot_flag) accuracy_plot(mcmc_samples$vbeta[,i], dnorm,
                              var_result$vmu[i], sqrt(var_result$mLambda[i,i]))
    }
  } else {
    vbeta_accuracy <- NULL
    vbeta_means <- NULL
  }
  
  if (mult$m > 0) {
    # vu accuracy
    vu_accuracy <- rep(NA, ncol(mult$mZ))
    vu_means <- rep(NA, ncol(mult$mZ))
    if (mult$spline_dim == 0) {
      # FIXME: Does this really have to be this complicated?
      B <- mult$blocksize
      b_idx <- 1
      for (i in 1:ncol(mult$mZ)) {
        m_idx <- ceiling(i / B)
        vu_mean <- var_result$vmu[i + mult$p]
        vu_means[i] <- mean(mcmc_samples$vu[, m_idx, b_idx])
        vu_sd <- sqrt(var_result$mLambda[i + mult$p, i + mult$p])
        vu_accuracy[i] <- calculate_accuracy(mcmc_samples$vu[, m_idx, b_idx], dnorm, vu_mean, vu_sd)
        if (print_flag) cat("vu[", i, "]", approximation, "accuracy:", vu_accuracy[i], "\n")
        if (plot_flag) accuracy_plot(mcmc_samples$vu[, m_idx, b_idx], dnorm, vu_mean, vu_sd)

        b_idx <- b_idx + 1
        if (b_idx > B)
        b_idx=1
      }
    } else {
      for (i in 1:ncol(mult$mZ)) {
        vu_mean <- var_result$vmu[i + mult$p]
        vu_sd <- sqrt(var_result$mLambda[i + mult$p, i + mult$p])
        vu_accuracy[i] <- calculate_accuracy(mcmc_samples$vu[, i], dnorm, vu_mean, vu_sd)
        vu_means[i] <- mean(mcmc_samples$vu[, i])
        if (print_flag) cat("vu[", i, "]", approximation, "accuracy:", vu_accuracy[i], "\n")
        if (plot_flag) accuracy_plot(mcmc_samples$vu[, i], dnorm, vu_mean, vu_sd)
      }
    }
  } else {
    vu_accuracy <- NULL
    vu_means <- NULL
  }

  # sigma2_u accuracy
  # FIXME - this may be wrong for blocksize != 1?
  # This is totally wrong for the Inverse Wishart model?
  #sigma2_u_accuracy <- calculate_accuracy(1/mcmc_samples$sigma_u^2, dgamma,
  #                                        var_result$a_sigma, var_result$b_sigma)
  #if (print_flag) cat("sigma2_u", approximation, "accuracy:", sigma2_u_accuracy, "\n")
  #if (plot_flag) accuracy_plot(1/mcmc_samples$sigma_u^2, dgamma,
  #                        var_result$a_sigma, var_result$b_sigma)
  # TODO: Can I just modify this to check each of the diagonals of the covariance matrices?
  # According to Wikipedia's article on Wishart distributions, the diagonal elements follow
  # chi-squared distributions
  
  # rho accuracy
  rho_accuracy <- calculate_accuracy(mcmc_samples$rho, dbeta,
                                     var_result$a_rho, var_result$b_rho)
  if (print_flag) cat("rho", approximation, "accuracy: ", rho_accuracy, "\n")
  if (plot_flag) accuracy_plot(mcmc_samples$rho, dbeta,
                          var_result$a_rho, var_result$b_rho)
  if (plot_flag) dev.off()
  return(list(var_result=var_result,
              vbeta_means=vbeta_means,
              vbeta_accuracy=vbeta_accuracy,
              vu_means=vu_means,
              vu_accuracy=vu_accuracy,
              #sigma2_u_accuracy=sigma2_u_accuracy,
              rho_accuracy=rho_accuracy))
}

# FIXME: This function should be unnecessary. Stan's fit object already contains this
# information and more.
mcmc_means <- function(mult, mcmc_samples)
{
  mX <- mult$mX
  mZ <- mult$mZ
  
  vbeta_means <- rep(NA, ncol(mX))
  for (p in 1:ncol(mX)) {
    vbeta_means[p] <- mean(mcmc_samples$vbeta[, p])
  }

  vu_means <- rep(NA, ncol(mZ))
  for (m in 1:ncol(mZ)) {
    vu_means[m] <- mean(mcmc_samples$vu[, m])
  }

  return(list(vbeta_means=vbeta_means,
              vu_means=vu_means,
              vmu=c(vbeta_means, vu_means)))
}

test_spline_accuracy <- function(mult, allKnots, mcmc_samples, approximation, plot=FALSE)
{
  var_result <- zero_infl_var(mult, method=approximation, verbose=TRUE)
  # Calculate the mean for vbeta, vu
  mcmc_means <- mcmc_means(mult, mcmc_samples)
  # Construct a BSpline matrix over the range we wish to plot
  # Plot the function using our MCMC and VB estimates
  # mCtilde %*% vmu
  xtilde <- seq(from=-1, to=1, by=1e-2)
  result <- spline.des(allKnots, xtilde, derivs=rep(0, length(xtilde)), outer.ok=TRUE)
  mC_tilde <- cbind(1, xtilde, result$design)
  f_hat_vb <- mC_tilde %*% var_result$vmu
  f_hat_mcmc <- mC_tilde %*% mcmc_means$vmu
  # TODO: Plot true function as well
  if (plot) {
    plot(mult$mX[,2], mult$vy)
    vf <- 4 + sin(pi * xtilde)
    lines(xtilde, exp(vf), type="l", col="black")
    lines(xtilde, exp(f_hat_mcmc), type="l", col="red")
    lines(xtilde, exp(f_hat_vb), type="l", col="blue")
    legend("topleft", c("True function", "MCMC", "VB"), fill=c("black", "red", "blue"))
  }
  #return(calculate_accuracies(mult, mcmc_samples, var_result, approximation, plot_flag=plot))
  return()
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
  now <- Sys.time()
  var1_result <- zero_infl_var(mult, method="laplace", verbose=TRUE)
  var1_accuracy <- calculate_accuracies(mult, mcmc_samples, var1_result, "laplace", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var1$var_result$mLambda)))
  print(var1_accuracy)
  
  now <- Sys.time()
  var2_result <- zero_infl_var(mult, method="gva", verbose=TRUE)
  var2_accuracy <- calculate_accuracies(mult, mcmc_samples, var2_result, "gva", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var2$var_result$mLambda)))
  print(var2_accuracy)

  now <- Sys.time()
  var3_result <- zero_infl_var(mult, method="gva2", verbose=TRUE)
  var3_accuracy <- calculate_accuracies(mult, mcmc_samples, var3_result, "gva2", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var3$var_result$mLambda)))
  print(var3_accuracy)

  now <- Sys.time()
  var4_result <- zero_infl_var(mult, method="gva_nr", verbose=TRUE)
  var4_accuracy <- calculate_accuracies(mult, mcmc_samples, var4_result, "gva_nr", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var4$var_result$mLambda)))
  print(var4_accuracy)
}
#test_accuracies()

test_accuracies_slope <- function()
{
  # Monte Carlo Markov Chains approximation
  # seed <- 1
  # set.seed(seed)
  # mult <- generate_slope_test_data(m=20, ni=10)
  # mcmc_samples <- mcmc_approximation(mult, seed=seed, iterations=1e5, warmup=1e4)
  # save(mult, mcmc_samples, file="data/accuracy_slope_2015_05_04.RData")  
  load(file="data/accuracy_slope_2015_05_04.RData")
  # load(file="data_macbook/accuracy_slope_2015_03_30.RData")
  
  now <- Sys.time()
  var1_result <- zero_infl_var(mult, method="laplace", verbose=TRUE)
  var1_accuracy <- calculate_accuracies(mult, mcmc_samples, var1_result, "laplace", plot_flag=TRUE)
  print(Sys.time() - now)
  #print(image(Matrix(var1$var_result$mLambda)))
  print(var1_accuracy)
  
  now <- Sys.time()
  var2_result <- zero_infl_var(mult, method="gva", verbose=TRUE)
  var2_accuracy <- calculate_accuracies(mult, mcmc_samples, var2_result, "gva", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var2$var_result$mLambda)))
  print(var2_accuracy)

  now <- Sys.time()
  var3_result <- zero_infl_var(mult, method="gva2", verbose=TRUE)
  var3_accuracy <- calculate_accuracies(mult, mcmc_samples, var3_result, "gva2", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var3$var_result$mLambda)))
  print(var3_accuracy)

  now <- Sys.time()
  var4_result <- zero_infl_var(mult, method="gva_nr", verbose=TRUE)
  var4_accuracy <- calculate_accuracies(mult, mcmc_samples, var4_result, "gva_nr", plot_flag=plot)
  print(Sys.time() - now)
  #print(image(Matrix(var4$var_result$mLambda)))
  print(var4_accuracy)
}
# test_accuracies_slope()

test_accuracies_spline <- function()
{
  # Monte Carlo Markov Chains approximation
  # seed <- 1
  # set.seed(seed)
  # result <- generate_spline_test_data()
  # mult <- result$mult
  # allKnots <- result$allKnots
  # mcmc_result <- mcmc_approximation(mult, seed=seed, iterations=1e5, warmup=1e3,
  #                                   stan_file="multivariate_zip_splines.stan")
  # mcmc_samples <- mcmc_result$mcmc_samples
  # fit <- mcmc_result$fit
  # print(fit)
  # save(mult, mcmc_samples, fit, allKnots, file="data/accuracy_spline_2015_05_19.RData")
  load(file="data/accuracy_spline_2015_05_19.RData")
  
  # now <- Sys.time()
  # var1 <- test_accuracy(mult, mcmc_samples, "laplace", plot=TRUE)
  # print(Sys.time() - now)
  # print(var1)
  
  # now <- Sys.time()
  # var2 <- test_accuracy(mult, mcmc_samples, "gva", plot=TRUE)
  # print(Sys.time() - now)
  # #print(image(Matrix(var2$var_result$mLambda)))
  # print(var2)

  # now <- Sys.time()
  # var3 <- test_accuracy(mult, mcmc_samples, "gva2", plot=TRUE)
  # print(Sys.time() - now)
  # #print(image(Matrix(var3$var_result$mLambda)))
  # print(var3)

  # now <- Sys.time()
  # var4 <- test_accuracy(mult, mcmc_samples, "gva_nr", plot=TRUE)
  # print(Sys.time() - now)
  # #print(image(Matrix(var4$var_result$mLambda)))
  # print(var4)

  # John says I should look at the functions
  test_spline_accuracy(mult, allKnots, mcmc_samples, "gva2", plot=TRUE)
  # test_spline_accuracy(mult, allKnots, mcmc_samples, "gva_nr", plot=TRUE)
}
test_accuracies_spline()