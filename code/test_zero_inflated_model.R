# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC
source("zero_inflated_model.R")
require(testthat)

test_univariate_zip <- function()
{
	n = 10000
	expected_rho = .5
	expected_lambda = 1

	vx = generate_test_data(n, expected_rho, expected_lambda)

	a_lambda = 0.01
	b_lambda = 0.01

	univariate = create_univariate(vx, a_lambda, b_lambda)
	result_var = zero_infl_var(univariate)

	expect_equal(result_var$a_lambda / result_var$b_lambda, expected_lambda, tolerance=1e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=1e-1)
}

#main_check_accuracy()
main <- function()
{
	test_univariate_zip()
}
main()
