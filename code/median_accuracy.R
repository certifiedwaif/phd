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
    result <- generate_int_test_data(m, ni, expected_beta = c(2, 1), expected_rho = 0.7)
  } else if (test == "slope") {
    m <- 20
    ni <- 10
    result <- generate_slope_test_data(m=20, ni=10, expected_beta=c(2, 1), expected_rho=0.7)
  } else if (test == "spline") {
    # result <- generate_spline_test_data()
    # mult <- result$mult
    # allKnots <- result$allKnots
    result <- generate_spline_test_data()
  }
  result
}

save_mcmc <- function(test)
{
  ITER <- 100

  # Saving the fit allows us to skip the recompilation of the C++ for the model
  # Update: Well, it's supposed to. Stan seems to want to keep recompiling the
  # model anyway, regardless of what I do.
  accuracy <- mclapply(1:ITER, function(i) {
  # accuracy <- lapply(97, function(i) {
    # Idea: Save generated data so we can reproduce problem data sets?
    result <- generate_test_case(i, test)
    if (test == "spline") {
      mult <- result$mult
      stan_fit <- mcmc(mult, seed=i, iterations=1e5, warmup=1e3,
                          stan_file="multivariate_zip_splines.stan")
    } else {
      mult <- result
      stan_fit <- mcmc(mult, iterations=1e5, warmup=1e4, mc.cores = 1)
    }
    hostname <- Sys.info()["nodename"]
    if (hostname == "verona.maths.usyd.edu.au")
      save(stan_fit, file = sprintf("/dskh/nobackup/markg/phd/stan_fit_%s_%d.RData", test, i))
    else if (hostname == "markg-OptiPlex-9020")
      save(stan_fit, file = sprintf("/home/markg/phd_data/stan_fit_%s_%d.RData", test, i))
    else
      stop("Cannot find hostname")
  }, mc.cores = 40)
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

median_accuracy <- function(approximation="gva", test="intercept", save_flag=FALSE)
{
  if (save_flag) {
    ITER <- 100

    # Saving the fit allows us to skip the recompilation of the C++ for the model
    # Update: Well, it's supposed to. Stan seems to want to keep recompiling the
    # model anyway, regardless of what I do.
    accuracy <- lapply(1:ITER, function(i) {
      cat("Test case", i, "\n")
      result <- generate_test_case(i, test)
      if (test == "spline") {
        allKnots <- result$allKnots
        mult <- result$mult
      } else {
        mult <- result
      }
      var_result <- zipvb(mult, method=approximation, verbose=TRUE)
      stan_fit <- load_stan_data(i, test)
      mcmc_samples <- stan_fit$mcmc_samples
      # Must handle the case where integrate() throws a non-finite function value error
      # If the integration fails, something must have gone fundamentally wrong, and we should
      # figure out why.
      if (test == "spline") {
        vy <- mult$vy
        vx <- mult$mX[,1]
        xtilde <- seq(from=-1, to=1, by=1e-2)
        result <- spline.des(allKnots, xtilde, derivs=rep(0, length(xtilde)), outer.ok=TRUE)
        mC_tilde <- cbind(1, xtilde, result$design)
        f_hat_vb <- mC_tilde %*% var_result$vmu
        pars <- stan_fit$fit@sim$pars_oi
        probs <- c(0.025, 0.25, 0.5, 0.75, 0.975)
        s <- summary(stan_fit$fit, pars, probs)
        mcmc_vmu <- s$summary[1:14, 1]
        f_hat_mcmc <- mC_tilde %*% mcmc_vmu
        mcmc_spline_fn <- splinefun(xtilde, f_hat_mcmc)
        vb_spline_fn <- splinefun(xtilde, f_hat_vb)
        result <- calculate_accuracy_spline(xtilde, mcmc_spline_fn, vb_spline_fn)
      } else {
        result <- calculate_accuracies("", mult, mcmc_samples, var_result, approximation, print_flag=TRUE)
      }
      # In spline case, should calculate accuracy of the fit functions
      save(result, file = sprintf("results/accuracy_result_%s_%s_%d.RData", approximation, test, i))
      result
    }) # , mc.cores = 32)
    
    save(accuracy, file = sprintf("results/accuracy_list_%s_%s.RData", approximation, test))
  }
  median_accuracy_graph_all(test)
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

# Somehow do this across all approximations.
# So for each parameter, plot the Laplace, GVA, GVA2 and GVA2 accuracy
# But label the x axis only for the parameter. Use a legend to indicate the approximation
median_accuracy_graph <- function(approximation="gva", test="intercept") {
  accuracy_df <- data.frame()
  # Create accuracy data frame
  # Merge with other accuracy data frames
  cat("median_accuracy_graph: approximation", approximation, ", test", test, "\n")
  load(file = sprintf("results/accuracy_list_%s_%s.RData", approximation, test))
  
  accuracy_df_delta <- create_accuracy_df(accuracy, make_colnames=TRUE)
  accuracy_df <- rbind(accuracy_df, cbind(approximation, accuracy_df_delta))

  # Create graph
  # Order by parameter, then approximation within that

  pdf_fname <- sprintf("results/median_accuracy_%s_%s.pdf", test, approximation)
  cat("Creating", pdf_fname, "\n")
  pdf(pdf_fname)
  # For all parameters of interest
  par(mfrow=c(1, 4))
  # Just label the graphs and move on.
  boxplot(as.numeric(as.character(vbeta_0))~approximation, main="vbeta_0", data=accuracy_df)
  boxplot(as.numeric(as.character(vbeta_1))~approximation, main="vbeta_1",  data=accuracy_df)
  boxplot(as.numeric(as.character(sigma2_vu_1))~approximation, main="sigma2_vu_1",  data=accuracy_df)
  boxplot(as.numeric(as.character(sigma2_vu_2))~approximation, main="sigma2_vu_2",  data=accuracy_df)
  par(mfrow=c(1, 1))
  if (test == "intercept") {
    mean_vu <- apply(accuracy_df[, 3:22], 1, mean)
    accuracy_df2 <- cbind(accuracy_df[, c("vbeta_0", "vbeta_1", "sigma2_vu", "rho")])
    boxplot(accuracy_df2, ylim=c(0, 1))
    axis(1, at=1:5, labels=c(expression(bold(beta)[0], bold(beta)[1], bold(u)[1], bold(sigma[u]^2), rho)))
  } else if (test == "slope") {
    mean_vu <- matrix(0, nrow(accuracy_df), 2)
    mean_vu[, 1] <- apply(accuracy_df[, 3 + 0:19 * 2], 2, function(x) mean(as.numeric(x)))
    mean_vu[, 2] <- apply(accuracy_df[, 3 + 1 + 0:19 * 2], 2,  function(x) mean(as.numeric(x)))
    mean_vu_df <- data.frame(mean_vu_0=mean_vu[, 1], mean_vu_1=mean_vu[, 2])
    accuracy_df2 <- cbind(accuracy_df[, c("vbeta_0", "vbeta_1")],
                          mean_vu_df,
                          accuracy_df[, c("sigma2_vu_1", "sigma2_vu_2", "rho")])
    colnames(accuracy_df2) <- rep("", 7)
    boxplot(accuracy_df2) # , ylim=c(0, 1))
    axis(1, at=1:7, labels=c(expression(bold(beta)[0], bold(beta)[1], bold(u)[1], bold(u)[2], bold(sigma[u1]^2), bold(sigma[u2]^2), rho)))
  } else if (test == "spline") {
    # TODO: Base on random slopes case
    stop("I have no idea what to do")
  }
  title(sprintf("%s %s median accuracy", approximation, test))
  dev.off()
}

median_accuracy_graph_all <- function(test) {
  # for (approximation in c("gva", "gva2")) {
  # for (approximation in c("laplace", "gva", "gva2", "gva_nr")) {
  #   median_accuracy_graph(approximation=approximation, test=test)
  #   median_accuracy_csv(approximation=approximation, test=test)
  # }

  # Load all of the data files
  accuracy_df <- data.frame()
  for (approximation in c("laplace", "gva", "gva2", "gva_nr")) {
    load(file = sprintf("results/accuracy_list_%s_%s.RData", approximation, test))
    
    accuracy_df_delta <- create_accuracy_df(accuracy, make_colnames=TRUE)
    accuracy_df <- rbind(accuracy_df, cbind(approximation, accuracy_df_delta))
  }

  # Take the mean of the random effects components for intercept and slope  
  mean_vu <- matrix(0, nrow(accuracy_df), 2)
  mean_vu[, 1] <- apply(accuracy_df[, 3 + 0:19 * 2], 2, function(x) mean(as.numeric(x)))
  mean_vu[, 2] <- apply(accuracy_df[, 3 + 1 + 0:19 * 2], 2,  function(x) mean(as.numeric(x)))
  mean_vu_df <- data.frame(mean_vu_0=mean_vu[, 1], mean_vu_1=mean_vu[, 2])
  accuracy_df2 <- cbind(accuracy_df[, c("vbeta_0", "vbeta_1")],
                        mean_vu_df,
                        accuracy_df[, c("sigma2_vu_1", "sigma2_vu_2", "rho")])
  parameters <- colnames(accuracy_df2)
  
  # Construct a data frame with a column for each combination of approximation and parameter
  acc <- data.frame(a=rep(NA, 100))
  for (param in parameters) {
    for (approximation in c("laplace", "gva", "gva2", "gva_nr")) {
      # Take the mean of the vu accuracies
      accuracy_numbers <- as.numeric(as.character(accuracy_df2[accuracy_df$approximation == approximation, param]))
      acc <- cbind(acc, accuracy_numbers)
      colnames(acc)[length(colnames(acc))] <- sprintf("%s_%s", approximation, param)
    }
  }
  # acc <- acc[, 2:ncol(acc)]

  # Plot the combined graph
  pdf("results/median_accuracy_combined.pdf")
  boxplot(acc,
          col=1:4,
          xaxt="n",
          ylab="Accuracy")
  axis(1,
       labels = c(expression(bold(beta)[0]), 
                  expression(bold(beta)[1]),
                  expression(bold(u)[0]),
                  expression(bold(u)[1]),
                  expression(sigma[bold(u)[0]]),
                  expression(sigma[bold(u)[1]]),
                  expression(rho)),
       at = 1:7 * 4 - 0.25,
       tick=TRUE)
  legend("bottomright", c("Laplace", "GVA", "GVA2", "GVA FP"), lty=1, col=1:4)
  title("Combined median accuracy graph")
  dev.off()
  # TODO: Label the axes better
  # TODO: Why are the accuracies for the variance components so low, although the estimation for the vbeta
  # components is better?
}

median_accuracy_csv <- function(approximation="gva", test=test) {
  load(file = sprintf("results/accuracy_list_%s_%s.RData", approximation, test))
  accuracy_df <- create_accuracy_df(accuracy, make_colnames=TRUE)
  write.csv(accuracy_df, file=sprintf("results/accuracy_list_%s.csv", approximation))
}

# Graph of Var_q(theta) against Var(theta|y)
# How to get this?
# Run fits for a range of theta values?

# median_accuracy()

is.between <- function(x, a, b) x >= a && x <= b

coverage_percentage <- function(approximation="gva", test="intercept")
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
  expected_rho <- 0.7
  result <- matrix(0, 1e4, 42)
	for (i in 1:1e4) {
    set.seed(i)
      cat("i", i, "counter", counter, "percentage", round(100*counter/i, 2), "\n")
		# mult <- generate_test_data(m, ni, expected_beta = expected_beta, expected_rho = expected_rho)
    mult <- generate_test_case(i, test)
    
		# var_result <- tryCatch(zipvb(mult, method=approximation, verbose=FALSE),
    #                        error <- function (E) { return(NULL) })
    var_result <- zipvb(mult, method=approximation, verbose=FALSE)
    cat("vmu", var_result$vmu, "\n")
    result[i, ] <- var_result$vmu
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
  write.csv(result, file="results/coverage_results.csv")
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
                      make_option(c("-m", "--median_accuracy"),  action="store_true", default=FALSE),
                      make_option(c("-s", "--save_flag"),  action="store_true", default=FALSE),
                      make_option(c("-d", "--save_mcmc"),  action="store_true", default=FALSE))
  opt <- parse_args(OptionParser(option_list=option_list))
  test <- opt$test
  if (opt$save_mcmc) {
    save_mcmc(test)
  } else if (opt$median_accuracy) {
    approximation <- opt$approximation
    median_accuracy(approximation=approximation, test=test, save_flag=opt$save_flag)
  } else if (opt$coverage_percentage) {
    approximation <- opt$approximation
    coverage_percentage(approximation=approximation)
  }
}
# main()
