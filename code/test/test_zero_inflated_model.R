# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC
library(zipvb)
# source("rwmh.R")
require(testthat)
#require(Matrix)
source("generate.R")

test_univariate_zip <- function()
{
	# Simulate data
	m = 10000
	expected_rho = .5
	expected_lambda = 1

	vx = generate_univariate_test_data(m, expected_rho, expected_lambda)

	a_lambda = 0.01
	b_lambda = 0.01

	# Test model fitting
	univariate = create_univariate(vx, a_lambda, b_lambda)
	var_result = zero_infl_var(univariate)
  return(var_result)
}

test_multivariate_zip_no_zeros <- function(approximation="gva")
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	# FIXME: You have serious overflow issues
	m = 50
	n = rep(1, m)
	mX = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
	cat("mX", mX, "\n")
	mZ = NULL
	expected_rho = 1
	expected_mu = c(2, 1)
	expected_sigma2_u = 0
  sigma2.beta = 1e5
	a_sigma = 1e5
	b_sigma = 1e5
  test_data = gen_mult_test_data(mX, NULL, m, n, expected_rho, expected_mu, expected_sigma2_u)
	vy = test_data$vy

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta)
	var_result = zero_infl_var(multivariate, verbose=TRUE, method=approximation)
  return(var_result)
}

test_multivariate_zip_half_zeros <- function(approximation="gva")
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	m = 100
	n = rep(1, m)
	mX = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
	mZ = NULL
	expected_rho = .5
	expected_mu = c(2, 1)
	expected_sigma2_u = 0
	sigma2.beta = 1e5  
	a_sigma = 1e5
	b_sigma = 1e5
	test_data = gen_mult_test_data(mX, NULL, m, n, expected_rho, expected_mu, expected_sigma2_u)
	vy = test_data$vy

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta)
	var_result = zero_infl_var(multivariate, verbose=TRUE, method=approximation)
  return(var_result)
}

test_multivariate_zip_no_zeros_random_intercept <- function(approximation="gva")
{
	# Simulate data
	m = 20
	ni = 10
	n = rep(ni,m)
	mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
	
	mZ <- kronecker(diag(1,m),rep(1,ni))
	
	expected_rho = 1
	expected_beta = c(2, 1)
	expected_sigma2_u = .5^2
	a_sigma = 1e-2
	b_sigma = 1e-2
	
	sigma2.beta <- 1.0E3
	
	tau = 1.0E2
	
	test_data = gen_mult_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=TRUE)
	vy = test_data$vy
	
	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta)
	var_result = zero_infl_var(multivariate, method=approximation, verbose=TRUE)
  return(var_result)
}

test_multivariate_zip_half_zeros_random_slope <- function(approximation="gva")
{
	m = 20
	ni = 10
	n = rep(ni,m)
	mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
  mZ = makeZ(mX, m, ni, p=2)
  
	expected_rho = 0.5
	expected_beta = c(2, 1)
	expected_sigma2_u = .5^2
	a_sigma = 1e-2
	b_sigma = 1e-2
	
	tau = 1.0E2
	
	sigma2.beta <- 1.0E3
	
	test_data = gen_mult_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=TRUE)
	vy = test_data$vy
	
	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, blocksize=2, sigma2.beta)
	var_result = zero_infl_var(multivariate, method=approximation, verbose=TRUE)
  return(var_result)
}

# Idea: Run each of the tests for convergence repeatedly.

test_spline = function(approximation="gva")
{
  n = 5000
  vx = matrix(sort(runif(n, -1, 1))) 
  
  mX = cbind(1,vx)
 
  expected_rho = 1
  #expected_mu = c(0, 1)
  expected_sigma2_u = 0
  sigma2.beta = 1e5
  a_sigma = 1e5
  b_sigma = 1e5
  tau = 1.0E-5
  
  sigma2.true = 0.01
  expected_beta = c(0, 1)
  vf = 5+2*sin(pi*vx)
  vy = rpois(n,exp(vf))
 
  source("ZOsull.r")
  numIntKnots <- 35
  intKnots <- quantile(unique(vx),seq(0,1,length=(numIntKnots+2))[-c(1,(numIntKnots+2))])
  
  mZ = ZOSull(vx,range.x=c(-1.1,1.1),intKnots=intKnots,drv=0)
  #vy = 2+mX[,1]^3+rnorm(m)*.1
  #result = fit_spline(vx, vy)
  #result = fit_spline(mX[,1], vy)
  #mZ = result$Z
  
  #mZ <- mZ/max(mZ)
 
  mult = create_multivariate(vy, mX, mZ, sigma2.beta, m=0, blocksize=1, spline_dim=37)
  
  var_result = zero_infl_var(mult, method=approximation, verbose=TRUE)
}

main <- function()
{
	set.seed(5)
	options(recover = traceback)
  
	test_univariate_zip()
  
	# Tests with multivariate fixed effects
	test_multivariate_zip_no_zeros("gva")
	test_multivariate_zip_no_zeros("gva2")
	test_multivariate_zip_half_zeros("gva")
	test_multivariate_zip_half_zeros("gva2")
	
	# Tests with multivariate fixed effects and random intercepts
	test_multivariate_zip_no_zeros_random_intercept("gva")
	test_multivariate_zip_no_zeros_random_intercept("gva2")
	test_multivariate_zip_half_zeros_random_intercept("gva")
	test_multivariate_zip_half_zeros_random_intercept("gva2")
	test_multivariate_zip_half_zeros_random_intercept("gva_nr")
  test_spline("gva")
}

#main()
