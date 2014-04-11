# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC
source("zero_inflated_model.R")
require(testthat)

generate_univariate_test_data <- function (n, rho, lambda)
{
	vx = rep(NA, n)
	for (i in 1:n) {
		if (runif(1) <= rho) {
			vx[i] = rpois(1, lambda)
		} else {
			vx[i] = 0
		}
	}
	return(vx)
}

generate_multivariate_test_data <- function (mX, rho, vbeta, sigma2_u)
{
	n = nrow(mX)
	vy = rep(NA, n)
	for (i in 1:n) {
		if (runif(1) <= rho) {
			vy[i] = rpois(1, exp(mX[i,] %*% vbeta) + rnorm(1, 0, sigma2_u))
		} else {
			vy[i] = 0
		}
	}
	return(vy)
}

test_univariate_zip <- function()
{
	# Simulate data
	n = 10000
	expected_rho = .5
	expected_lambda = 1

	vx = generate_univariate_test_data(n, expected_rho, expected_lambda)

	a_lambda = 0.01
	b_lambda = 0.01

	# Test model fitting
	univariate = create_univariate(vx, a_lambda, b_lambda)
	result_var = zero_infl_var(univariate)

	expect_equal(result_var$a_lambda / result_var$b_lambda, expected_lambda, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

test_multivariate_zip_no_zeros <- function()
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	n = 100
	mX = matrix(as.vector(cbind(rep(1, n), rnorm(n))), n, 2)
	expected_rho = 1
	expected_beta = c(1, 2)
	expected_sigma2_u = 0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mX, expected_rho, expected_beta, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mX, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate)

	expect_equal(as.vector(result_var$vbeta), expected_beta, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

test_multivariate_zip_half_zeros <- function()
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	n = 100
	mX = matrix(as.vector(cbind(rep(1, n), rnorm(n))), n, 2)
	expected_rho = .5
	expected_beta = c(1, 2)
	expected_sigma2_u = 0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mX, expected_rho, expected_beta, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mX, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate, trace=TRUE)

	expect_equal(as.vector(result_var$vbeta), expected_beta, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

test_multivariate_zip_no_zeros_random_intercept <- function()
{
	stop("Not implemented yet")
}

test_multivariate_zip_half_zeros_random_intercept <- function()
{
	stop("Not implemented yet")
}

#main_check_accuracy()
main <- function()
{
	test_univariate_zip()
	# TODO: Add some sort of test for the accuracy of the approximation?

	test_multivariate_zip_no_zeros()
	test_multivariate_zip_half_zeros()

	test_multivariate_zip_no_zeros_random_intercept()
	test_multivariate_zip_half_zeros_random_intercept()
	# TODO: Add a test for the random intercepts?
}
main()
