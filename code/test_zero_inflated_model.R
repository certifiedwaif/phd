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

generate_multivariate_test_data <- function (mX, mZ, m, n, rho, vnu, sigma2_u)
{
	mC = cbind(mX, mZ)
	vy = rep(NA, sum(n))
	for (i in 1:m) {
		intercept = rnorm(1, 0, sigma2_u)
		for (j in 1:n[i]) {
			if (i == 1) {
				ind = j
			} else {
				ind = j + sum(n[1:(i-1)])
			}
			#cat("ind", ind, "\n")

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
	mX = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
	cat("mX", mX, "\n")
	mZ = NULL
	expected_rho = 1
	expected_nu = c(1, 2)
	expected_sigma2_u = 0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mX, m, n, expected_rho, expected_nu, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, a_sigma, b_sigma)
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
	mX = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
	mZ = NULL
	expected_rho = .5
	expected_nu = c(1, 2)
	expected_sigma2_u = 0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mX, m, n, expected_rho, expected_nu, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate, trace=TRUE)

	expect_equal(as.vector(result_var$vnu), expected_nu, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

test_multivariate_zip_no_zeros_random_intercept <- function()
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	# FIXME: You have serious overflow issues
	m = 2
	n = c(25, 25)
	mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
	cat("mX", mX, "\n")
	# Indicator variables for groups
	v = c(rep(1, 25), rep(0, 25))
	mZ = matrix(cbind(v, 1-v), sum(n), 2)
	cat("mZ", mZ, "\n")
	expected_rho = 1
	expected_nu = c(1, 2, 1, 2)
	expected_sigma2_u = 10.0
	a_sigma = 1e5
	b_sigma = 1e5
	vy = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_nu, expected_sigma2_u)

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate)

	expect_equal(as.vector(result_var$vnu), expected_nu, tolerance=5e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=2e-1)
	print(str(result_var))
	expect_equal(result_var$b_sigma^2 / result_var$a_sigma, expected_sigma2_u, tolerance=1e-1)
}

test_multivariate_zip_half_zeros_random_intercept <- function()
{
	stop("Not implemented yet")
}

#main_check_accuracy()
main <- function()
{
	#test_univariate_zip()
	# TODO: Add some sort of test for the accuracy of the approximation?

	#test_multivariate_zip_no_zeros()
	#test_multivariate_zip_half_zeros()

	# TODO: Add a test for the random intercepts?
	test_multivariate_zip_no_zeros_random_intercept()
	#test_multivariate_zip_half_zeros_random_intercept()
}
main()
