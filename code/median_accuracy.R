# median_accuracy.R
source("test_zero_inflated_model.R")
source("accuracy.R")

# Repeatedly run trials and compare accuracy. Plot boxplots.
median_accuracy = function()
{
  ITER = 100

  set.seed(1234)
  accuracy = list()
  for (i in 1:ITER) {
    # Run code
    result = compare_approximations(c(1, 2))
    accuracy[[i]] = with(result, calculate_accuracy(mcmc_samples, var_result))
  }
  # For name in names(accuracy)
  boxplot(accuracy)
}

# Graph of Var_q(theta) against Var(theta|y)
# How to get this?
# Run fits for a range of theta values?
compare_approximations = function(vbeta)
{
  m = 20
  ni = 10
  n = rep(ni,m)
  mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  mZ <- kronecker(diag(1,m),rep(1,ni))
  
  expected_rho = 0.5
  expected_beta = vbeta
  expected_sigma2_u = .5^2
  a_sigma = 1e-2
  b_sigma = 1e-2
  
  tau = 1.0E2
  
  sigma2.beta <- 1.0E3
  
  test_data = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=FALSE)
  vy = test_data$vy
  
  multivariate = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau)
  approximation = "gva2new"
  result_var = zero_infl_var(multivariate, method=approximation, verbose=FALSE)
  mcmc_samples = mcmc_approximation(multivariate, iterations=1e4, mc.cores = 32)
  return(list(multivariate=multivariate, result_var=result_var, mcmc_samples=mcmc_samples))
}

mean_var = function(vbeta)
{
  result = compare_approximations(vbeta)
  return(with(result, {
  var_approx_mean = result_var$vmu[2]
  mcmc_approx_mean = mean(mcmc_samples$vbeta[,2])
  var_approx_var = result_var$mLambda[2,2]
  mcmc_approx_var = var(mcmc_samples$vbeta[,2])
  list(var_approx_mean=var_approx_mean,
		mcmc_approx_mean=mcmc_approx_mean,
  		var_approx_var=var_approx_var,
		mcmc_approx_var=mcmc_approx_var)
  }))
}

# This is most definitely a verona job
for (theta in seq(1, 2, by=.01))
	print(mean_var(c(1, theta)))

# Graph of E_q(theta) against E(theta|y)

# For most important parameters: beta_1, beta_2

plot_graphs = function()
{
  # Read files, plot graphs
  require(stringr)
  fn = function(fname)
  {
    lines = readLines(fname)
    num_str = str_match(lines, "[0-9]*\\.[0-9]*")
    return(na.omit(as.numeric(num_str)))
  }
  
  var_mean = fn("var_approx_mean.txt")
  var_var = fn("var_approx_var.txt")
  mcmc_mean = fn("mcmc_approx_mean.txt")
  mcmc_var = fn("mcmc_approx_var.txt")
  
  plot(mcmc_mean)
  points(var_mean, col=2)
  plot(mcmc_mean, var_mean)
  
  plot(mcmc_var)
  points(var_var, col=2)
  plot(mcmc_var, var_var)
}