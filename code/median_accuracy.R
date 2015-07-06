#!/usr/bin/env Rscript
# median_accuracy.R
library(optparse)
library(parallel)
source("generate.R")
source("zero_inflated_model.R")
source("accuracy.R")

# Repeatedly run trials and compare accuracy. Plot boxplots.
generate_test_case <- function(i, test) {
  set.seed(i)
  # Run code
  if (test == "intercept") {
    m <- 20
    ni <- 10
    mult <- generate_int_test_data(m, ni, expected_beta = c(2, 1), expected_rho = 0.7)
  } else if (test == "slope") {
    m <- 20
    ni <- 10
    mult <- generate_slope_test_data(m=20, ni=10, expected_beta=c(2, 1), expected_rho=0.7)
  }
  mult
}

save_mcmc_data <- function(test)
{
  ITER <- 100

  # Saving the fit allows us to skip the recompilation of the C++ for the model
  # Update: Well, it's supposed to. Stan seems to want to keep recompiling the
  # model anyway, regardless of what I do.
  accuracy <- mclapply(1:ITER, function(i) {
    # Idea: Save generated data so we can reproduce problem data sets?
    mult <- generate_test_case(i, test)
    stan_fit <- mcmc(mult, iterations=1e5, warmup=1e4, mc.cores = 1)
    save(stan_fit, file = sprintf("/dskh/nobackup/markg/phd/stan_fit_%s_%d.RData", test, i))
  }, mc.cores = 32)
}

load_stan_data <- function(i, test)
{
  hostname <- Sys.info()["nodename"]
  if (hostname == "verona.maths.usyd.edu.au")
    load(file = sprintf("/dskh/nobackup/markg/phd/stan_fit_%s_%d.RData", test, i))
  else if (hostname == "markg-OptiPlex-9020")
    load(file = sprintf("/home/markg/phd_data/stan_fit_%s_%d.RData", test, i))
  else
    stop("Cannot find hostname")
  return(stan_fit)
}

median_accuracy <- function(approximation="gva", test="intercept")
{
  ITER <- 100

  # Saving the fit allows us to skip the recompilation of the C++ for the model
  # Update: Well, it's supposed to. Stan seems to want to keep recompiling the
  # model anyway, regardless of what I do.
  accuracy <- lapply(1:ITER, function(i) {
    cat("Test case", i, "\n")
    mult <- generate_test_case(i, test)
    var_result <- zipvb(mult, method=approximation, verbose=TRUE)
    stan_fit <- load_stan_data(i, test)
    mcmc_samples <- stan_fit$mcmc_samples
    # Must handle the case where integrate() throws a non-finite function value error
    # If the integration fails, something must have gone fundamentally wrong, and we should
    # figure out why.
    result <- calculate_accuracies("", mult, mcmc_samples, var_result, approximation, print_flag=TRUE)
    save(result, file = sprintf("results/accuracy_result_%s_%s_%d.RData", approximation, test, i))
    result
  }) # , mc.cores = 32)
  
  save(accuracy, file = sprintf("results/accuracy_list_%s_%s.RData", approximation, test))
  # median_accuracy_graph_all(test)
}

create_accuracy_df <- function(accuracy, make_colnames=FALSE) {
  # Construct a data frame
  vbeta_accuracy <- t(sapply(accuracy, function(x) {x$vbeta_accuracy}))
  vu_accuracy <- t(sapply(accuracy, function(x) {x$vu_accuracy}))
  sigma2_vu_accuracy <- t(sapply(accuracy, function(x) {x$sigma2_vu_accuracy}))
  rho_accuracy <- sapply(accuracy, function(x) {x$rho_accuracy})
  accuracy_df <- cbind(vbeta_accuracy, vu_accuracy, sigma2_vu_accuracy, rho_accuracy)
  if (make_colnames) {
    colnames(accuracy_df) <- c("vbeta_0", "vbeta_1", 
      sapply(1:ncol(vu_accuracy), function(x) paste0("vu_", x)),
      sapply(1:ncol(sigma2_vu_accuracy), function(x) paste0("sigma2_vu_", x)), "rho")
  }
  accuracy_df
}

median_accuracy_graph <- function(approximation="gva", test="intercept") {
  load(file = sprintf("results/accuracy_list_%s_%s.RData", approximation, test))
  
  accuracy_df <- create_accuracy_df(accuracy, make_colnames=TRUE)
  browser()
  pdf(sprintf("results/median_accuracy_%s_%s.pdf", approximation, test))
  # Should be plotting mean, and comparing approximations
  if (test == "intercept") {
    mean_vu <- apply(accuracy_df[, 3:22], 1, mean)
    accuracy_df2 <- cbind(accuracy_df[, c("vbeta_0", "vbeta_1", "sigma2_vu", "rho")])
    boxplot(accuracy_df2, ylim=c(0, 1))
    axis(1, at=1:5, labels=c(expression(bold(beta)[0], bold(beta)[1], bold(u)[1], bold(sigma[u]^2), rho)))
  } else if (test == "slope") {
    mean_vu <- matrix(0, 100, 2)
    mean_vu[, 1] <- apply(accuracy_df[, 3 + 0:19 * 2], 1, mean)
    mean_vu[, 2] <- apply(accuracy_df[, 3 + 1 + 0:19 * 2], 1, mean)
    mean_vu_df <- data.frame(mean_vu_0=mean_vu[, 1], mean_vu_1=mean_vu[, 2])
    accuracy_df2 <- cbind(accuracy_df[, c("vbeta_0", "vbeta_1")],
                          mean_vu_df,
                          accuracy_df[, c("sigma2_vu_1", "sigma2_vu_2", "rho")])
    boxplot(accuracy_df2, ylim=c(0, 1))
    axis(1, at=1:7, labels=c(expression(bold(beta)[0], bold(beta)[1], bold(u)[1], bold(u)[2], bold(sigma[u1]^2), bold(sigma[u2]^2), rho)))
  }
  title(sprintf("%s %s median accuracy", approximation, test))
  dev.off()
}

median_accuracy_graph_all <- function(test) {
  # for (approximation in c("gva", "gva2")) {
  for (approximation in c("laplace", "gva", "gva2", "gva_nr")) {
    median_accuracy_graph(approximation=approximation, test=test)
    median_accuracy_csv(approximation=approximation, test=test)
  }
}

median_accuracy_csv <- function(approximation="gva", test=test) {
  load(file = sprintf("results/accuracy_list_%s.RData", approximation))
  accuracy_df <- create_accuracy_df(accuracy, make_colnames=TRUE)
  write.csv(accuracy_df, file=sprintf("results/accuracy_list_%s.csv", approximation))
}

# Graph of Var_q(theta) against Var(theta|y)
# How to get this?
# Run fits for a range of theta values?

# median_accuracy()

is.between <- function(x, a, b) x >= a && x <= b

coverage_percentage <- function(approximation="gva")
{
	# Percentage coverage of the true parameter values by approximate 95%
	# credible intervals based on variational Bayes approximate posterior
	# density functions. The percentages are based on 10,000 replications.
	
	counter <- rep(0, 3)
  p <- 2
  # If ni is substantially smaller than m, you're probably going to get a lot of
  # overflow errors.
  m <- 10
  ni <- 25
  expected_beta <- c(2, 1)
  expected_rho <- 0.5
	for (i in 1:1e4) {
    set.seed(i)
      cat("i", i, "counter", counter, "percentage", round(100*counter/i, 2), "\n")
		mult <- generate_test_data(m, ni, expected_beta = expected_beta, expected_rho = expected_rho)
    
		# var_result <- tryCatch(zipvb(mult, method=approximation, verbose=FALSE),
    #                        error <- function (E) { return(NULL) })
    var_result <- zipvb(mult, method=approximation, verbose=FALSE)
  	# If there was an error, re-try with another generated data set. Sometimes the
  	# optimiser gets passed a non-finite value.
  	if (is.null(var_result)) {
    		i <- i - 1
    		next
  	}

    # TODO: Fix intercepts
    for (j in 1:p) {
      # Check that true parameter is within 95% credible interval.
      lci <- var_result$vmu[j] - 1.96*sqrt(var_result$mLambda[j, j])
      uci <- var_result$vmu[j] + 1.96*sqrt(var_result$mLambda[j, j])
      # cat("vbeta lci[", j, "]", lci, "\n")
      # cat("vbeta uci[", j, "]", uci, "\n")
      if (is.between(expected_beta[j], lci, uci)) {
        # Increment counter
        counter[j] <- counter[j] + 1
      }
    }

    u_dim <- (p + 1):(p + m - 1)
  	sum_vmu_vu <- sum(var_result$vmu[u_dim])
  	sum_mLambda_vu <- sum(sqrt(diag(var_result$mLambda[u_dim, u_dim]))) # I think this will be an over-estimate
    lci <- sum_vmu_vu - 1.96*sum_mLambda_vu
    uci <- sum_vmu_vu + 1.96*sum_mLambda_vu
    # cat("vu lci", lci, "\n")
    # cat("vu uci", uci, "\n")

    # cat("vbeta[1] + sum_vu", var_result$vmu[1] + sum_vmu_vu, "\n")

  	if (is.between(0, sum_vmu_vu - 1.96 * sum_mLambda_vu, sum_vmu_vu + 1.96 * sum_mLambda_vu)) {
    		counter[p + 1] <- counter[p + 1] + 1
  	}
	}
	print(counter)
}
# coverage_percentage()

mean_var <- function(vbeta)
{
  result <- compare_approximations(vbeta)
  return(with(result, {
  var_approx_mean <- var_result$vmu[2]
  mcmc_approx_mean <- mean(mcmc_samples$vbeta[,2])
  var_approx_var <- var_result$mLambda[2,2]
  mcmc_approx_var <- var(mcmc_samples$vbeta[,2])
  list(var_approx_mean=var_approx_mean,
		mcmc_approx_mean=mcmc_approx_mean,
  		var_approx_var=var_approx_var,
		mcmc_approx_var=mcmc_approx_var)
  }))
}

# This is most definitely a verona job
print_mean_var <- function()
{
  for (theta in seq(1, 2, by=.01))
	  print(mean_var(c(1, theta)))
}
# Graph of E_q(theta) against E(theta|y)

# For most important parameters: beta_1, beta_2

plot_graphs <- function()
{
  # Read files, plot graphs
  require(stringr)
  fn <- function(fname)
  {
    lines <- readLines(fname)
    num_str <- str_match(lines, "[0-9]*\\.[0-9]*")
    return(na.omit(as.numeric(num_str)))
  }
  
  var_mean <- fn("var_approx_mean.txt")
  var_var <- fn("var_approx_var.txt")
  mcmc_mean <- fn("mcmc_approx_mean.txt")
  mcmc_var <- fn("mcmc_approx_var.txt")
  
  plot(mcmc_mean)
  points(var_mean, col=2)
  plot(mcmc_mean, var_mean)
  
  plot(mcmc_var)
  points(var_var, col=2)
  plot(mcmc_var, var_var)
}

main <- function()
{
  option_list <- list(make_option(c("-a", "--approximation"), default="gva"),
                      make_option(c("-t", "--test"), default="intercept"),
                      make_option(c("-c", "--coverage_percentage"), action="store_true", default=FALSE),
                      make_option(c("-m", "--median_accuracy"),  action="store_true", default=FALSE))
  opt <- parse_args(OptionParser(option_list=option_list))
  approximation <- opt$approximation
  test <- opt$test
  if (opt$median_accuracy) {
    median_accuracy(approximation=approximation, test=test)
  }

  if (opt$coverage_percentage) {
    coverage_percentage(approximation=approximation)
  }
}

main()
