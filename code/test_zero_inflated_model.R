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

generate_multivariate_test_data <- function (mC, m, n, rho, vnu, sigma2_u)
{
	vy = rep(NA, sum(n))
	for (i in 1:m) {
		intercept = rnorm(1, 0, sigma2_u)
		for (j in 1:n[i]) {
			if (i == 1) {
				ind = j
			} else {
				ind = j + sum(n[1:(i-1)])
			}
			cat("ind", ind, "\n")

			if (runif(1) <= rho) {
				vy[ind] = rpois(1, exp(mC[ind,] %*% vnu) + intercept)
			} else {
				vy[ind] = 0
			}
		}
	}
	return(vy)
}

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
	result_var = zero_infl_var(univariate)

	expect_equal(result_var$a_lambda / result_var$b_lambda, expected_lambda, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

test_multivariate_zip_no_zeros <- function()
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	# FIXME: You have serious overflow issues
	m = 50
	n = rep(1, m)
	mC = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
	cat("mC", mC, "\n")
	expected_rho = 1
	expected_nu = c(1, 2)
	expected_sigma2_u = 0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mC, m, n, expected_rho, expected_nu, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mC, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate)

	expect_equal(as.vector(result_var$vnu), expected_nu, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

test_multivariate_zip_half_zeros <- function()
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	m = 100
	n = rep(1, m)
	mC = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
	expected_rho = .5
	expected_nu = c(1, 2)
	expected_sigma2_u = 0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mC, m, n, expected_rho, expected_nu, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mC, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate, trace=TRUE)

	expect_equal(as.vector(result_var$vnu), expected_nu, tolerance=1e-1)
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

	# TODO: Add a test for the random intercepts?
	test_multivariate_zip_no_zeros_random_intercept()
	test_multivariate_zip_half_zeros_random_intercept()
}
main()
