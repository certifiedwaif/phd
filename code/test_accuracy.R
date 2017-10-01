#!/usr/bin/env Rscript
# test_accuracy.R
source("accuracy.R")
source("generate.R")
library(optparse)

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
    save(mult, mcmc_samples, fit, file="data/accuracy_int_2017_05_30.RData")
  } else {
    load(file="data/accuracy_int_2017_05_30.RData")
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
    mult <- generate_slope_test_data(m=50, ni=10)
    result <-  mcmc(mult, iterations=1e5, warmup = 5e4)
    fit <- result$fit
    print(fit)
    mcmc_samples <- result$mcmc_samples
    save(mult, mcmc_samples, fit, file="data/accuracy_slope_2017_05_26.RData")  
    # load(file="data/accuracy_slope_2015_05_04.RData")
    # load(file="data_macbook/accuracy_slope_2015_03_30.RData")
  } else {
    # load(file="data/accuracy_slope_2017_05_26.RData")  
    load(file="data/accuracy_slope_2017_05_24.RData")  
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
  var3_result <- zipvb(mult, method="gva2", verbose=FALSE)
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

test_stability <- function(save=FALSE)
{
  load(file="data/accuracy_int_2015_05_27.RData")
  thresholds <- seq(from=0.01, to=0.3, by=0.01)
  results <- matrix(NA, length(thresholds), 24)
  for (i in 1:length(thresholds)) {
    options(threshold=thresholds[i])
    var3_result <- zipvb(mult, method="gva2", verbose=FALSE)
    var3_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var3_result, "gva2", plot_flag=TRUE)
    results[i, ] <- c(thresholds[i], var3_accuracy$vbeta_accuracy, var3_accuracy$vu_accuracy, var3_accuracy$sigma2_vu_accuracy, var3_accuracy$rho_accuracy)
  }
  write.csv(results, file="stability_intercept.csv")
  return()
  load(file="data/accuracy_slope_2015_06_22.RData")
  for (threshold in seq(from=0.1, to=5, by=0.1)) {
    options(safe_exp=TRUE, threshold=threshold)
    var3_result <- zipvb(mult, method="gva2", verbose=FALSE)
    var3_accuracy <- calculate_accuracies("intercept", mult, mcmc_samples, var3_result, "gva2", plot_flag=TRUE)
    cat(threshold, var3_accuracy$vbeta_accuracy, var3_accuracy$vu_accuracy, "\n")
  }
}



main <- function()
{
  option_list <- list(make_option(c("-t", "--test"), default="slope"),
                      make_option(c("-s", "--save"), action="store_true", default=FALSE))
  opt <- parse_args(OptionParser(option_list=option_list))
  print(str(opt))
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
  if (test == "stability") {
    test_stability(opt$save)
  }
}

main()
