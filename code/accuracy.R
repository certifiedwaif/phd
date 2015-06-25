#!/usr/bin/env Rscript
# accuracy.R
source("zero_inflated_model.R")
source("generate.R")
library(rstan)
library(optparse)

# Andrew Gelman says that this magical line of code automatically makes Stan
# run in parallel and cache compiled models.
# Two minutes later: Hey, it actually works!
# Some time later: Sometimes it works.
#source("http://mc-stan.org/rstan/stan.R")

mcmc <- function(mult, seed=1, iterations=NA, warmup=NA, mc.cores=1,
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
  
  zip_data <- list(N=length(vy), P=2, M=m, B=blocksize, spline_dim=spline_dim,
                   y=vy, X=mX, Z=as.matrix(mZ),
                   v=v, psi=mPsi, BetaPrior=mSigma.beta)

  fit <- stan(stan_file, fit=stan_fit, seed=seed, data=zip_data, iter=iterations, warmup=warmup,
              chains=1)
  mcmc_samples <- extract(fit)

  return(list(fit=fit,
              mcmc_samples=mcmc_samples))
}

calculate_accuracy <- function(mcmc_samples, dist_fn, ...)
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
  result <- integrate(integrand, min(mcmc_density$x), max(mcmc_density$x),
                     subdivisions = length(mcmc_density$x))
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
        lty=2, col="blue")
  lines(mcmc_density$x, mcmc_density$y * opt$object / max(mcmc_density$y))
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
      title <- sprintf("%s vbeta[%d] accuracy: %f", approximation, i, vbeta_accuracy[i])
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
        title <- sprintf("%s vu[%d] accuracy: %f", approximation, i, vu_accuracy[i])
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
        title <- sprintf("%s vu[%d] accuracy: %f", approximation, i, vu_accuracy[i], "\n")
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
  sigma2_vu_accuracy <- rep(0, B)
  if (B == 2) {
    v <- 4
    d <- var_result$mult$d
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
  } else {
    v <- 2
    sigma_u_inv <- array(0, c(dim(mcmc_samples$sigma_u)[1], 1, 1))
    sigma_u_inv[, 1, 1] <- 1 / mcmc_samples$sigma_u[, 1, 1]
    print(dim(sigma_u_inv))
  }
  for (i in 1:B) {
    # sigma2 <- E_mPsi_inv[i, i]
    sigma2 <- var_result$mPsi[i, i]
    sigma2_vu_accuracy[i] <- calculate_accuracy(sigma_u_inv[, i, i],
                                               function(x, ...) dgamma(x/sqrt(sigma2), ...), v/2, sigma2/2)
    title <- sprintf("%s sigma2_u[%d] accuracy: %f", approximation, i, sigma2_vu_accuracy[i])
    if (print_flag) print(title)
    if (plot_flag)
      accuracy_plot(title, sigma_u_inv[, i, i], function(x, ...) dgamma(x/sqrt(sigma2), ...), v/2, sigma2/2)
    
  }
  
  # rho accuracy
  rho_accuracy <- calculate_accuracy(mcmc_samples$rho, dbeta,
                                     var_result$a_rho, var_result$b_rho)
  title <- sprintf("%s rho accuracy: %f", approximation, rho_accuracy, "\n")
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

test_spline_accuracy <- function(mult, allKnots, fit, approximation, plot=FALSE)
{
  var_result <- zipvb(mult, method=approximation, verbose=TRUE)
  # Calculate the mean for vbeta, vu
  # Construct a BSpline matrix over the range we wish to plot
  # Plot the function using our MCMC and VB estimates
  # mCtilde %*% vmu
  xtilde <- seq(from=-1, to=1, by=1e-2)
  result <- spline.des(allKnots, xtilde, derivs=rep(0, length(xtilde)), outer.ok=TRUE)
  mC_tilde <- cbind(1, xtilde, result$design)
  f_hat_vb <- mC_tilde %*% var_result$vmu
  # FIXME: This is probably broken
  f_hat_mcmc <- mC_tilde %*% fit$vmu$mean
  if (plot) {
    pdf(sprintf("results/accuracy_plots_spline_%s.pdf", approximation))
    plot(mult$mX[,2], mult$vy)
    vf <- 4 + sin(pi * xtilde)
    lines(xtilde, exp(vf), type="l", col="black")
    lines(xtilde, exp(f_hat_mcmc), type="l", col="red")
    lines(xtilde, exp(f_hat_vb), type="l", col="blue")
    legend("topleft", c("True function", "MCMC", "VB"), fill=c("black", "red", "blue"))
    dev.off()
  }
  #return(calculate_accuracies(mult, mcmc_samples, var_result, approximation, plot_flag=plot))
  return()
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
  var1_result <- zipvb(mult, method="laplace", verbose=TRUE)
  print(Sys.time() - now)
  var1_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var1_result, "laplace", plot_flag=TRUE)
  # #print(image(Matrix(var1$var_result$mLambda)))
  # print(var1_accuracy)
  print(var1_accuracy$vbeta_accuracy)
  print(var1_accuracy$vu_accuracy)
  print(var1_accuracy$rho_accuracy)
  
  now <- Sys.time()
  var2_result <- zipvb(mult, method="gva", verbose=TRUE)
  print(Sys.time() - now)
  var2_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var2_result, "gva", plot_flag=TRUE)
  # #print(image(Matrix(var2$var_result$mLambda)))
  # print(var2_accuracy)
  print(var2_accuracy$vbeta_accuracy)
  print(var2_accuracy$vu_accuracy)
  print(var2_accuracy$rho_accuracy)

  now <- Sys.time()
  var3_result <- zipvb(mult, method="gva2", verbose=TRUE)
  print(Sys.time() - now)
  var3_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var3_result, "gva2", plot_flag=TRUE)
  #print(image(Matrix(var3$var_result$mLambda)))
  # print(var3_accuracy)
  print(var3_accuracy$vbeta_accuracy)
  print(var3_accuracy$vu_accuracy)
  print(var3_accuracy$rho_accuracy)

  now <- Sys.time()
  var4_result <- zipvb(mult, method="gva_nr", verbose=TRUE)
  print(Sys.time() - now)
  var4_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var4_result, "gva_nr", plot_flag=TRUE)
  # #print(image(Matrix(var4$var_result$mLambda)))
  # print(var4_accuracy)
  print(var4_accuracy$vbeta_accuracy)
  print(var4_accuracy$vu_accuracy)
  print(var4_accuracy$rho_accuracy)
}
# test_accuracies()

test_accuracies_slope <- function(save=FALSE)
{
  # Monte Carlo Markov Chains approximation
  # 1 good
  # 2, 3 bad - good with more samples
  # 4 good
  # 5 bad
  if (save) {
    seed <- 3
    set.seed(seed)
    mult <- generate_slope_test_data(m=20, ni=10)
    result <-  mcmc(mult, iterations=3e4, warmup = 5e3)
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
  var1_result <- zipvb(mult, method="laplace", verbose=TRUE)
  print(Sys.time() - now)
  var1_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var1_result, "laplace", print_flag=TRUE, plot_flag=TRUE)
  print(var1_accuracy$vbeta_accuracy)
  print(var1_accuracy$vu_accuracy)
  print(var1_accuracy$rho_accuracy)
  
  now <- Sys.time()
  var2_result <- zipvb(mult, method="gva", verbose=TRUE)
  print(Sys.time() - now)
  var2_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var2_result, "gva", print_flag=TRUE, plot_flag=TRUE)
  print(var2_accuracy$vbeta_accuracy)
  print(var2_accuracy$vu_accuracy)
  print(var2_accuracy$rho_accuracy)

  now <- Sys.time()
  var3_result <- zipvb(mult, method="gva2", verbose=TRUE)
  print(Sys.time() - now)
  var3_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var3_result, "gva2", print_flag=TRUE, plot_flag=TRUE)
  print(var3_accuracy$vbeta_accuracy)
  print(var3_accuracy$vu_accuracy)
  print(var3_accuracy$rho_accuracy)

  now <- Sys.time()
  var4_result <- zipvb(mult, method="gva_nr", verbose=TRUE)
  print(Sys.time() - now)
  var4_accuracy <- calculate_accuracies("slope", mult, mcmc_samples, var4_result, "gva_nr", print_flag=TRUE, plot_flag=TRUE)
  print(var4_accuracy$vbeta_accuracy)
  print(var4_accuracy$vu_accuracy)
  print(var4_accuracy$rho_accuracy)
}
# test_accuracies_slope()

test_accuracies_spline <- function(save=FALSE)
{
  # Monte Carlo Markov Chains approximation
  if (save) {
    seed <- 1
    set.seed(seed)
    result <- generate_spline_test_data()
    mult <- result$mult
    allKnots <- result$allKnots
    mcmc_result <- mcmc(mult, seed=seed, iterations=1e5, warmup=1e3,
                                      stan_file="multivariate_zip_splines.stan")
    mcmc_samples <- mcmc_result$mcmc_samples
    fit <- mcmc_result$fit
    print(fit)
    save(mult, mcmc_samples, fit, allKnots, file="data/accuracy_spline_2015_05_19.RData")
  } else {
    load(file="data/accuracy_spline_2015_05_19.RData")
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
  option_list <- list(make_option(c("-a", "--accuracy"), default="slope"),
                      make_option(c("-s", "--save"), action="store_true", default=FALSE))
  opt <- parse_args(OptionParser(option_list=option_list))
  accuracy <- opt$accuracy
  if (accuracy == "intercept") {
    test_accuracies_intercept(opt$save)
  }
  if (accuracy == "slope") {
    test_accuracies_slope(opt$save)
  }
  if (accuracy == "spline") {
    test_accuracies_spline(opt$save)
  }
}

# main()
