#!/usr/bin/env Rscript
# accuracy.R
library(zipvb)
library(rstan)
library(optparse)
library(mvtnorm)
library(limma)
library(latex2exp)

source("generate.R")

# Andrew Gelman says that this magical line of code automatically makes Stan
# run in parallel and cache compiled models.
# Two minutes later: Hey, it actually works!
# Some time later: Sometimes it works.
#source("http://mc-stan.org/rstan/stan.R")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

mcmc <- function(mult, seed=1, iterations=NA, warmup=NA, mc.cores=1, p=2,
                               stan_file="multivariate_zip.stan", stan_fit=NA)
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
  
  zip_data <- list(N=length(vy), P=p, M=m, B=blocksize, spline_dim=spline_dim,
                   y=vy, X=mX, Z=as.matrix(mZ),
                   v=v, psi=mPsi, BetaPrior=mSigma.beta)

  options(mc.cores = mc.cores)
  fit <- stan(stan_file, fit=stan_fit, seed=seed, data=zip_data, iter=iterations, warmup=warmup,
              chains=1)
  mcmc_samples <- extract(fit)

  return(list(fit=fit,
              mcmc_samples=mcmc_samples))
}

##### Trapezoidal Integration
##### By Eric Cai - The Chemical Statistician
##### define the function for trapezoidal integration

trapezoidal.integration = function(x, f)
{
### 3 checks to ensure that the arguments are numeric and of equal lengths
# check if the variable of integration is numeric
  if (!is.numeric(x))
  {
    stop('The variable of integration "x" is not numeric.')
  }

# check if the integrand is numeric
  if (!is.numeric(f))
  {
    stop('The integrand "f" is not numeric.')
  }

  # check if the variable of integration and the integrand have equal lengths
  if (length(x) != length(f))
  {
    stop('The lengths of the variable of integration and the integrand do not match.')
  }

  ### finish checks

  # obtain length of variable of integration and integrand
  n <- length(x)

  # integrate using the trapezoidal rule
  integral <- 0.5*sum((x[2:n] - x[1:(n-1)]) * (f[2:n] + f[1:(n-1)]))

  # print the definite integral
  return(integral)
}

integrate2 <- function(f, min_x, max_x, subdivisions = 0)
{
  # Calculate x and f for trapezoidal.integration
  vx <- seq(min_x, max_x, length.out = subdivisions)
  vf <- f(vx)
  list(value=trapezoidal.integration(vx, vf))
}

calculate_accuracy_normalised <- function(mcmc_samples, dist_fn, ...)
{
  mcmc_density <- density(mcmc_samples)
  opt <- optimize(function(x) dist_fn(x, ...), interval=c(min(mcmc_density$x), max(mcmc_density$x)), maximum=TRUE)
  mcmc_fn <- splinefun(mcmc_density$x, mcmc_density$y * opt$object / max(mcmc_density$y))
  result1 <- integrate(mcmc_fn, min(mcmc_density$x), max(mcmc_density$x),
                     subdivisions = length(mcmc_density$x))
  result2 <- integrate(function(x) dist_fn(x, ...), min(mcmc_density$x), max(mcmc_density$x),
                     subdivisions = length(mcmc_density$x))
  
  integrand <- function(x)
  {
    return(abs(mcmc_fn(x)/result1$value - dist_fn(x, ...)/result2$value))
  }
  result <- integrate2(integrand, min(mcmc_density$x), max(mcmc_density$x),
                     subdivisions = length(mcmc_density$x))
  accuracy <- 1 - .5 * result$value
  return(100 * accuracy)
}

calculate_accuracy <- function(mcmc_samples, dist_fn, ...)
{
  mcmc_density <- density(mcmc_samples)
  mcmc_fn <- splinefun(mcmc_density$x, mcmc_density$y)
  
  integrand <- function(x)
  {
    return(abs(mcmc_fn(x) - dist_fn(x, ...)))
  }
  result <- integrate2(integrand, min(mcmc_density$x), max(mcmc_density$x),
                     subdivisions = length(mcmc_density$x))
  accuracy <- 1 - .5 * result$value
  return(100 * accuracy)
}

calculate_accuracy_spline <- function(vx, mcmc_fn, vb_fn)
{
  integrand <- function(x)
  {
    return(abs(mcmc_fn(x) - vb_fn(x)))
  }
  result <- integrate2(integrand, min(vx), max(vx),
                     subdivisions = length(vx))
  accuracy <- 1 - .5 * result$value
  return(accuracy)
  
}

accuracy_plot <- function(title, mcmc_samples, dist_fn, ...)
{
  mcmc_density <- density(mcmc_samples)
  # Find the maximum of the function and use it to normalise plots.
  opt <- optimize(function(x) dist_fn(x, ...), interval=c(min(mcmc_density$x), max(mcmc_density$x)), maximum=TRUE)
  # browser()
  curve(dist_fn(x, ...),
        from=min(mcmc_density$x), to=max(mcmc_density$x),
        lty=2, col="blue", ylab="")
  lines(mcmc_density$x, mcmc_density$y * opt$object / max(mcmc_density$y))
  legend("topright", c("MCMC estimate", "VB estimate"), lty=c(1, 2))
  title(title)
  # mcmc_fn <- splinefun(mcmc_density$x, mcmc_density$y)
  # curve(mcmc_fn(x, ...),
  #       from=min(mcmc_density$x), to=max(mcmc_density$x),
  #       lty=1, col="black", add=TRUE)
}

calculate_accuracies <- function(test, mult, mcmc_samples, var_result, approximation, print_flag=FALSE, plot_flag=FALSE)
{
  # TODO: Add support for checking the accuracy over multiple dimensions
  # cubature$adaptIntegrate
  
  if (plot_flag) pdf(paste0("results/accuracy_plots_", test, "_", approximation, ".pdf"))
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
      title <- latex2exp(sprintf("%s $\\textbf{\\beta_%d}$ accuracy: %2.0f%%", approximation, i, vbeta_accuracy[i]))
      if (print_flag) print(title)
      if (plot_flag) accuracy_plot(title, mcmc_samples$vbeta[,i], dnorm,
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
        title <- sprintf("%s vu[%d] accuracy: %2.0f%%", approximation, i, vu_accuracy[i])
        if (print_flag) print(title)
        if (plot_flag) accuracy_plot(title, mcmc_samples$vu[, m_idx, b_idx], dnorm, vu_mean, vu_sd)

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
        title <- sprintf("%s vu[%d] accuracy: %2.0f%%", approximation, i, vu_accuracy[i], "\n")
        if (print_flag) print(title)
        if (plot_flag) accuracy_plot(title, mcmc_samples$vu[, i], dnorm, vu_mean, vu_sd)
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
  # TODO: Add mSigma.beta accuracy
  # \sigma_Z^2 <- solve(mSigma)[i, i]^2
  # sigma2_z * dchisq(x, df = v)
  # TODO: Add mSigma.u accuracy
  # TODO: Need to generalise this to the splines case
  # We need to iterate through the sigma_u samples, inverting each matrix
  n <- dim(mcmc_samples$sigma_u)[1]
  sigma_u_inv <- array(0, dim(mcmc_samples$sigma_u))
  # for (i in 1:n) {
  #   sigma_u_inv[i, , ] <- solve(sigma_u[i, , ])
  # }
  sigma2_vu_accuracy <- rep(0, B)
  v <- mult$v
  if (B == 2) {
    # E_mPsi_inv <- var_result$mult$mPsi / (4 - 2 - 1)
    sigma_u_inv <- array(0, c(dim(mcmc_samples$sigma_u)[1], 2, 2))
    a <- mcmc_samples$sigma_u[, 1, 1]
    b <- mcmc_samples$sigma_u[, 1, 2]
    c <- mcmc_samples$sigma_u[, 2, 1]
    d <- mcmc_samples$sigma_u[, 2, 2]
    sigma_u_inv[, 1, 1] <- d / (a * d - b * c)
    sigma_u_inv[, 1, 2] <- -b / (a * d - b * c)
    sigma_u_inv[, 2, 1] <- -c / (a * d - b * c)
    sigma_u_inv[, 2, 2] <- a / (a * d - b * c)
    # sigma_u_inv <- mcmc_samples$sigma_u
  } else {
    sigma_u_inv <- array(0, c(dim(mcmc_samples$sigma_u)[1], 1, 1))
    sigma_u_inv[, 1, 1] <- 1 / mcmc_samples$sigma_u[, 1, 1]
    # sigma_u_inv[, 1, 1] <- mcmc_samples$sigma_u[, 1, 1]
  }
  mPsi_inv <- solve(var_result$mPsi)
  for (i in 1:B) {
    # sigma2 <- E_mPsi_inv[i, i]
    # sigma2_inv <- solve(var_result$mPsi)[i, i]
    sigma_vu <- sqrt(mPsi_inv[i, i])
    # sigma_vu <- sqrt(var_result$mPsi[i, i])
    sigma2_vu_accuracy[i] <- calculate_accuracy_normalised(sigma_u_inv[, i, i],
                                               function(x, ...) dchisq(x / sigma_vu, df = v, ...))
    title <- latex2exp(sprintf("%s $\\sigma^2_{u_%d}$ accuracy: %2.0f%%",
                       approximation,
                       i,
                       sigma2_vu_accuracy[i]))
    if (print_flag) print(title)
    if (plot_flag) {
      accuracy_plot(title, sigma_u_inv[, i, i], function(x, ...) dchisq(x / sigma_vu, df = v, ...))
    }
  }
  
  # rho accuracy
  rho_accuracy <- calculate_accuracy(mcmc_samples$rho, dbeta,
                                     var_result$a_rho, var_result$b_rho)
  title <- sprintf("%s rho accuracy: %2.0f%%", approximation, rho_accuracy, "\n")
  if (print_flag) print(title)
  if (plot_flag) accuracy_plot(title, mcmc_samples$rho, dbeta,
                          var_result$a_rho, var_result$b_rho)
  if (plot_flag) dev.off()
  return(list(var_result=var_result,
              vbeta_means=vbeta_means,
              vbeta_accuracy=vbeta_accuracy,
              vu_means=vu_means,
              vu_accuracy=vu_accuracy,
              sigma2_vu_accuracy=sigma2_vu_accuracy,
              rho_accuracy=rho_accuracy))
}

test_accuracies_intercept <- function(save=FALSE)
{
  # Need to be able to compare the solution paths of each approximation
  
  # Generate data
  # for (i in 1:100) {
  #   set.seed(i)
  #   mult = generate_test_data(20, 100)
  #   # Monte Carlo Markov Chains approximation
  #   mcmc_samples = mcmc(mult, iterations=1e6)
  #   # Save the results, because this takes such a long time to run.
  # }
  # save(mult, mcmc_samples, file="accuracy_good.RData")
  if (save) {
    set.seed(1)
    m <- 20
    ni <- 10
    mult <- generate_int_test_data(m, ni, expected_beta = c(2, 1), expected_rho = 0.5)
    # Monte Carlo Markov Chains approximation
    result <- mcmc(mult, iterations=1e5, warmup = 1e4)
    fit <- result$fit
    mcmc_samples <- result$mcmc_samples
  #   # Save the results, because this takes such a long time to run.
    #save(mult, mcmc_samples, file="accuracy.RData")
    #save(mult, mcmc_samples, file="data/accuracy_int.RData")
    save(mult, mcmc_samples, fit, file="data/accuracy_int_2015_05_27.RData")
  } else {
    load(file="data/accuracy_int_2015_05_27.RData")
  }
  # load(file="data/accuracy_int_2015_02_17.RData")
  #mult$spline_dim = 0
  #load(file="accuracy.RData")
  # Test all other approximations against it
  #load(file="accuracy.RData")
  
  # Test multivariate approximation's accuracy
  now <- Sys.time()
  var1_result <- zipvb(mult, method="laplace", verbose=FALSE)
  print(Sys.time() - now)
  var1_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var1_result, "laplace", plot_flag=TRUE)
  # #print(image(Matrix(var1$var_result$mLambda)))
  # print(var1_accuracy)
  print(var1_accuracy$vbeta_accuracy)
  print(mean(var1_accuracy$vu_accuracy))
  print(var1_accuracy$sigma2_vu_accuracy)
  print(var1_accuracy$rho_accuracy)
  
  now <- Sys.time()
  var2_result <- zipvb(mult, method="gva", verbose=FALSE)
  print(Sys.time() - now)
  var2_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var2_result, "gva", plot_flag=TRUE)
  # #print(image(Matrix(var2$var_result$mLambda)))
  # print(var2_accuracy)
  print(var2_accuracy$vbeta_accuracy)
  print(mean(var2_accuracy$vu_accuracy))
  print(var2_accuracy$sigma2_vu_accuracy)
  print(var2_accuracy$rho_accuracy)

  now <- Sys.time()
  var3_result <- zipvb(mult, method="gva2", verbose=FALSE)
  print(Sys.time() - now)
  var3_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var3_result, "gva2", plot_flag=TRUE)
  #print(image(Matrix(var3$var_result$mLambda)))
  # print(var3_accuracy)
  print(var3_accuracy$vbeta_accuracy)
  print(mean(var3_accuracy$vu_accuracy))
  print(var3_accuracy$sigma2_vu_accuracy)
  print(var3_accuracy$rho_accuracy)

  now <- Sys.time()
  var4_result <- zipvb(mult, method="gva_nr", verbose=FALSE)
  print(Sys.time() - now)
  var4_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var4_result, "gva_nr", plot_flag=TRUE)
  # #print(image(Matrix(var4$var_result$mLambda)))
  # print(var4_accuracy)
  print(var4_accuracy$vbeta_accuracy)
  print(mean(var4_accuracy$vu_accuracy))
  print(var4_accuracy$sigma2_vu_accuracy)
  print(var4_accuracy$rho_accuracy)
}
# test_accuracies()

test_accuracies_slope <- function(save=FALSE)
{
  # Monte Carlo Markov Chains approximation
  if (save) {
    seed <- 3
    set.seed(seed)
    mult <- generate_slope_test_data(m=20, ni=10)
    result <-  mcmc(mult, iterations=3e5, warmup = 5e4)
    fit <- result$fit
    print(fit)
    mcmc_samples <- result$mcmc_samples
    save(mult, mcmc_samples, fit, file="data/accuracy_slope_2015_06_22.RData")  
    # load(file="data/accuracy_slope_2015_05_04.RData")
    # load(file="data_macbook/accuracy_slope_2015_03_30.RData")
  } else {
    load(file="data/accuracy_slope_2015_06_22.RData")  
  }
  # m <- 20
  # mult$vmu <- c(2, 1, rep(0, (m-1) * 2))
  
  now <- Sys.time()
  var1_result <- zipvb(mult, method="laplace", verbose=FALSE)
  cat("Laplace", Sys.time() - now, "\n")
  var1_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var1_result, "laplace", print_flag=FALSE, plot_flag=TRUE)
  print(var1_accuracy$vbeta_accuracy)
  print(mean(var1_accuracy$vu_accuracy))
  print(var1_accuracy$sigma2_vu_accuracy)
  print(var1_accuracy$rho_accuracy)
  
  now <- Sys.time()
  var2_result <- zipvb(mult, method="gva", verbose=FALSE)
  cat("GVA", Sys.time() - now, "\n")
  var2_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var2_result, "gva", print_flag=FALSE, plot_flag=TRUE)
  print(var2_accuracy$vbeta_accuracy)
  print(mean(var2_accuracy$vu_accuracy))
  print(var2_accuracy$sigma2_vu_accuracy)
  print(var2_accuracy$rho_accuracy)

  now <- Sys.time()
  var3_result <- zipvb(mult, method="gva2", verbose=TRUE)
  cat("GVA2", Sys.time() - now, "\n")
  var3_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var3_result, "gva2", print_flag=FALSE, plot_flag=TRUE)
  print(var3_accuracy$vbeta_accuracy)
  print(mean(var3_accuracy$vu_accuracy))
  print(var3_accuracy$sigma2_vu_accuracy)
  print(var3_accuracy$rho_accuracy)

  now <- Sys.time()
  var4_result <- zipvb(mult, method="gva_nr", verbose=FALSE)
  cat("GVA NR", Sys.time() - now, "\n")
  var4_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var4_result, "gva_nr", print_flag=FALSE, plot_flag=TRUE)
  print(var4_accuracy$vbeta_accuracy)
  print(mean(var4_accuracy$vu_accuracy))
  print(var4_accuracy$sigma2_vu_accuracy)
  print(var4_accuracy$rho_accuracy)
}
# test_accuracies_slope()

# Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
sourceCpp(file="spline_ci.cpp")

test_spline_accuracy <- function(mult, allKnots, fit, approximation, plot=FALSE)
{
  # Change the knots
  # Small count example
  var_result <- zipvb(mult, method=approximation, verbose=TRUE)
  # John said this isn't good enough.
  # We should do a grid search from -1 to 1, simulating from the normal that
  # VB
  vmu_vb <- var_result$vmu
  mLambda_vb <- var_result$mLambda

  vbeta_mcmc <- apply(fit$vbeta, 2, mean)
  vu_mcmc <- apply(fit$vu, 2, mean)
  vmu_mcmc <- c(vbeta_mcmc, vu_mcmc)
  vmu_mcmc_samples <- cbind(fit$vbeta, fit$vu)

  x <- seq(from=-1, to=1, by=1e-2)
  vb_lci <- rep(0, length(x))
  vb_uci <- rep(0, length(x))
  mcmc_lci <- rep(0, length(x))
  mcmc_uci <- rep(0, length(x))
  for (i in 1:length(x)) {
    NUM_REPS <- 1000
    f_hat_vb <- rep(NA, NUM_REPS)
    f_hat_mcmc <- rep(NA, NUM_REPS)
    # Create an appropriate mC
    mC_x <- cbind(1, x[i], spline.des(allKnots, x[i], derivs=c(0), outer.ok=TRUE)$design)

    result <- spline_ci(NUM_REPS, mC_x, vmu_vb, mLambda_vb, fit$vbeta, fit$vu)
    f_hat_vb <- result$f_hat_vb
    f_hat_mcmc <- result$f_hat_mcmc
    
    # for (j in 1:NUM_REPS) {

    #   # Generate samples - simulate vbeta, vu
    #   f_hat_vb[j] <- t(rmvnorm(1, mC_x %*% vmu_vb, mC_x %*% mLambda_vb %*% t(mC_x)))
    #   sample_idx <- sample(dim(fit$vbeta)[1], 1)
    #   vbeta_mcmc <- fit$vbeta[sample_idx, ]
    #   vu_mcmc <- fit$vu[sample_idx, ]
    #   vnu_mcmc <- c(vbeta_mcmc, vu_mcmc)
    #   f_hat_mcmc[j] <- t(mC_x %*% vnu_mcmc)
    # }
    
    vb_quantiles <- quantile(f_hat_vb, c(.025, .975))
    vb_lci[i] <- vb_quantiles[1]
    vb_uci[i] <- vb_quantiles[2]
    mcmc_quantiles <- quantile(f_hat_mcmc, c(.025, .975))
    mcmc_lci[i] <- mcmc_quantiles[1]
    mcmc_uci[i] <- mcmc_quantiles[2]
  }

  if (plot) {
    pdf(sprintf("results/accuracy_plots_spline_%s.pdf", approximation))
    # pdf(sprintf("results/splines_ci_%s.pdf", approximation))
    within_bounds_idx <- (mult$mX[, 2] > -1) && (mult$mX[, 2] < 1)
    plot(x[within_bounds_idx], exp(vb_lci[within_bounds_idx]), type="l", col="blue",
         xlab="x", ylab=expression(3 + 3 * sin(pi * x)),
         xlim=c(-1.0, 1.0), ylim=c(0.0, exp(6.0)), lty=2)
    lines(x[within_bounds_idx], exp(vb_uci[within_bounds_idx]), col="blue", lty=2)
    lines(x[within_bounds_idx], exp(mcmc_lci[within_bounds_idx]), col="red", lty=2)
    lines(x[within_bounds_idx], exp(mcmc_uci[within_bounds_idx]), col="red", lty=2)
    # legend("topleft", c("VB", "MCMC"), fill=c("blue", "red"))
    
    # Calculate the mean for vbeta, vu
    # Construct a BSpline matrix over the range we wish to plot
    # Plot the function using our MCMC and VB estimates
    # mCtilde %*% vmu
    xtilde <- seq(from=-1, to=1, by=1e-2)
    result <- spline.des(allKnots, xtilde, derivs=rep(0, length(xtilde)), outer.ok=TRUE)
    mC_tilde <- cbind(1, xtilde, result$design)
    f_hat_vb <- mC_tilde %*% var_result$vmu
    f_hat_mcmc <- mC_tilde %*% vmu_mcmc

    points(mult$mX[within_bounds_idx,2], mult$vy[within_bounds_idx],
           xlim=c(-1, 1))
    vf <- 3 + 3 * sin(pi * xtilde)
    lines(xtilde[1:(length(xtilde)-1)], exp(vf[1:(length(xtilde)-1)]), type="l", col="black")
    lines(xtilde[1:(length(xtilde)-1)], exp(f_hat_mcmc[1:(length(xtilde)-1)]), type="l", col="red")
    lines(xtilde[1:(length(xtilde)-1)], exp(f_hat_vb[1:(length(xtilde)-1)]), type="l", col="blue")
    legend("topleft", c("True function", "MCMC estimate", "VB estimate"),
           fill=c("black", "red", "blue"))
    dev.off()
  }
  #return(calculate_accuracies(mult, mcmc_samples, var_result, approximation, plot_flag=plot))
  return()
}

test_accuracies_spline <- function(save=FALSE)
{
  # Monte Carlo Markov Chains approximation
  if (save) {
    seed <- 1
    set.seed(seed)
    result <- generate_spline_test_data(n=100)
    mult <- result$mult
    allKnots <- result$allKnots
    mcmc_result <- mcmc(mult, seed=seed, iterations=1e4, warmup=1e3,
                                      stan_file="multivariate_zip_splines.stan")
    mcmc_samples <- mcmc_result$mcmc_samples
    fit <- mcmc_result$fit
    print(fit)
    save(mult, mcmc_samples, fit, allKnots, file="data/accuracy_spline_2015_05_19.RData")
    # save(mult, mcmc_samples, fit, allKnots, file="/tmp/accuracy_spline_2015_05_19.RData")
  } else {
    load(file="data/accuracy_spline_2015_05_19.RData")
    # load(file="/tmp/accuracy_spline_2015_05_19.RData")
  }
  # John says I should look at the functions
  test_spline_accuracy(mult, allKnots, mcmc_samples, "laplace", plot=TRUE)
  test_spline_accuracy(mult, allKnots, mcmc_samples, "gva", plot=TRUE)
  test_spline_accuracy(mult, allKnots, mcmc_samples, "gva2", plot=TRUE)
  test_spline_accuracy(mult, allKnots, mcmc_samples, "gva_nr", plot=TRUE)
}
# test_accuracies_spline()

main <- function()
{
  option_list <- list(make_option(c("-t", "--test"), default="slope"),
                      make_option(c("-s", "--save"), action="store_true", default=FALSE))
  opt <- parse_args(OptionParser(option_list=option_list))
  test <- opt$test
  if (test == "intercept") {
    test_accuracies_intercept(opt$save)
  }
  if (test == "slope") {
    test_accuracies_slope(opt$save)
  }
  if (test == "spline") {
    test_accuracies_spline(opt$save)
  }
}

main()
