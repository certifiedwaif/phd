# median_accuracy.R
source("test_zero_inflated_model.R")
source("accuracy.R")

# Repeatedly run trials and compare accuracy. Plot boxplots.
median_accuracy <- function()
{
  ITER <- 100

  set.seed(1234)
  accuracy <- list()
  for (i in 1:ITER) {
    # Run code
    result <- compare_approximations(c(2, 1))
    accuracy[[i]] <- with(result, calculate_accuracy(mcmc_samples, var_result, print_flag=TRUE))
  }
  # For name in names(accuracy)
  
  result_df <- data.frame()
  for (i in 1:length(accuracy[[1]]$vu_accuracy)) {
    result_df[1:ITER, paste0("vu_accuracy_", i)] <- sapply(accuracy, function (x) x$vu_accuracy[i])
  }
  for (i in 1:length(accuracy[[1]]$vbeta_accuracy)) {
    result_df[1:ITER, paste0("vbeta_accuracy_", i)] <- sapply(accuracy, function (x) x$vbeta_accuracy[i])
  }
  result_df[1:ITER, "sigma2_u_accuracy"] <- sapply(accuracy, function (x) x$sigma2_u_accuracy)
  result_df[1:ITER, "rho_accuracy"] <- sapply(accuracy, function (x) x$rho_accuracy)
  pdf("median_accuracy.pdf")
  boxplot(result_df)
  dev.off()
}

# Graph of Var_q(theta) against Var(theta|y)
# How to get this?
# Run fits for a range of theta values?
compare_approximations <- function(vbeta, approximation="gva")
{
  multivariate <- generate_test_mult(expected_beta = c(2, 1))
  var_result <- zero_infl_var(multivariate, method=approximation, verbose=TRUE)
  mcmc_samples <- mcmc_approximation(multivariate, iterations=1e4, mc.cores = 4)
  return(list(multivariate=multivariate, var_result=var_result, mcmc_samples=mcmc_samples))
}

is.between <- function(x, a, b) x >= a && x <= b

coverage_percentage <- function()
{
	# Percentage coverage of the true parameter values by approximate 95%
	# credible intervals based on variational Bayes approximate posterior
	# density functions. The percentages are based on 1000 replications.
	
	set.seed(12345)
	counter <- rep(0, 3)
  p <- 2
  # If ni is substantially smaller than m, you're probably going to get a lot of
  # overflow errors.
  m <- 10
  ni <- 20
  expected_beta <- c(2, 1)
  expected_rho <- c(0.5, 1)
	for (i in 1:1000) {
    cat("i", i, "counter", counter, "percentage", round(100*counter/i, 2), "\n")
		mult <- generate_test_data(m, ni, expected_beta = expected_beta, expected_rho = 0.5)
		approximation <- "gva2"
    
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
			if (is.between(expected_beta[j], var_result$vmu[j] - 1.96*sqrt(var_result$mLambda[j, j]), var_result$vmu[j] + 1.96*sqrt(var_result$mLambda[j, j]))) {
				# Increment counter
				counter[j] <- counter[j] + 1
			}
  	}
    u_dim <- (p + 1):(p + m - 1)
  	sum_vmu_vmu <- sum(var_result$vmu[u_dim])
  	sum_mLambda_vmu <- sum(sqrt(diag(var_result$mLambda[u_dim, u_dim]))) # I think this will be an over-estimate
  	if (is.between(0, sum_vmu_vmu - 1.96 * sum_mLambda_vmu, sum_vmu_vmu + 1.96 * sum_mLambda_vmu)) {
    		counter[p + 1] <- counter[p + 1] + 1
  	}
	}
	print(counter)
}
coverage_percentage()

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
