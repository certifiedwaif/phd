#!/usr/bin/env Rscript
# median_accuracy.R
library(optparse)
library(parallel)
source("test_zero_inflated_model.R")
source("accuracy.R")

# Repeatedly run trials and compare accuracy. Plot boxplots.
median_accuracy <- function(approximation="gva")
{
  ITER <- 100

  # Saving the fit allows us to skip the recompilation of the C++ for the model
  # Update: Well, it's supposed to. Stan seems to want to keep recompiling the
  # model anyway, regardless of what I do.
  stan_fit <- NA
  accuracy <- mclapply(1:ITER, function(i) {
    set.seed(i)
    # Run code
    m <- 10
    ni <- 100
    done <- FALSE
    while (!done) {
      done <- TRUE
      mult <- generate_int_test_data(m, ni, expected_beta = c(2, 1), expected_rho = 0.5)
      # Make an initial guess for vmu
      var_result <- zero_infl_var(mult, method=approximation, verbose=TRUE)
      stanfit <- mcmc_approximation(mult, iterations=1e5, warmup=5e3, mc.cores = 1,
                                    stan_fit=stan_fit)
      mcmc_samples <- stanfit$mcmc_samples
      result <- tryCatch(calculate_accuracies("", mult, mcmc_samples, var_result, approximation, print_flag=TRUE),
                                error=function(e) {
                                  print(e)
                                  NULL
                                })
      # print(result)
      # result <- calculate_accuracies("", mult, mcmc_samples, var_result, approximation,
      #                                print_flag=TRUE)

      # Check whether integrate() failed. If so, generate another data set and re-try.
      if (is.null(result)) {
        done <- FALSE
      }
    }
    save(result, file = sprintf("results/accuracy_result_%d_%s.RData", i, approximation))
    result
  }, mc.cores = 32)
  
  # For name in names(accuracy)
  save(accuracy, file = sprintf("results/accuracy_list_%s.RData", approximation))
  
  # for (approximation in c("laplace", "gva", "gva2", "gva_nr")) {
  # for (approximation in c("gva", "gva2")) {
  #   load(file = sprintf("results/accuracy_list_%s.RData", approximation))
    
  #   # Construct a data frame
  #   vbeta_accuracy <- t(sapply(accuracy, function(x) {x$vbeta_accuracy}))
  #   vu_accuracy <- t(sapply(accuracy, function(x) {x$vu_accuracy}))
  #   sigma2_vu_accuracy <- sapply(accuracy, function(x) {x$sigma2_vu_accuracy})
  #   rho_accuracy <- sapply(accuracy, function(x) {x$rho_accuracy})
  #   accuracy_df <- cbind(vbeta_accuracy, vu_accuracy, sigma2_vu_accuracy, rho_accuracy)
  #   colnames(accuracy_df) <- c("vbeta_0", "vbeta_1", 
  #     sapply(1:19, function(x) paste0("vu_", x)),
  #     "sigma2_vu", "rho")
  #   pdf(sprintf("results/median_accuracy_%s.pdf", approximation))
  #   boxplot(accuracy_df, ylim=c(0, 1))
  #   axis(1, at=1:23, labels=c(expression(bold(beta)[0], bold(beta)[1], bold(u)[1], bold(u)[2], bold(u)[3], bold(u)[4], bold(u)[5], bold(u)[6], bold(u)[7], bold(u)[8], bold(u)[9], bold(u)[11], bold(u)[12], bold(u)[13], bold(u)[14], bold(u)[15], bold(u)[16], bold(u)[17], bold(u)[18], bold(u)[19], bold(u)[20], bold(sigma^2[u]), rho)))
  #   title(sprintf("%s median accuracy", approximation))
  #   dev.off()
  # }
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
    
		# var_result <- tryCatch(zero_infl_var(mult, method=approximation, verbose=FALSE),
    #                        error <- function (E) { return(NULL) })
    var_result <- zero_infl_var(mult, method=approximation, verbose=FALSE)
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
                      make_option(c("-c", "--coverage_percentage"), action="store_true", default=FALSE),
                      make_option(c("-m", "--median_accuracy"),  action="store_true", default=FALSE))
  opt <- parse_args(OptionParser(option_list=option_list))
  approximation <- opt$approximation
  if (opt$median_accuracy) {
    median_accuracy(approximation=approximation)
  }

  if (opt$coverage_percentage) {
    coverage_percentage(approximation=approximation)
  }
}

main()
