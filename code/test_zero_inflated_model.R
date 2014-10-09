# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC
source("zero_inflated_model.R")
source("rwmh.R")
require(testthat)
#require(Matrix)

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

generate_multivariate_test_data <- function (mX, mZ, m, n, rho, vbeta, sigma2_u, verbose=FALSE)
{
	if (is.null(mZ)) {
		mC = mX
	} else {
		mC = cbind(mX, mZ)
	}
	vy = rep(NA, sum(n))
	vu = rep(NA, length(n))
	if (verbose)
		cat("vy", vy, "\n")
	idx = 0
	for (i in 1:length(n)) {
		if (sigma2_u == 0)
			vu[i] = 0
		else
			vu[i] = rnorm(1, 0, sqrt(sigma2_u))

		for (j in 1:n[i]) {
			idx = idx + 1
			if (verbose)
				cat("idx", idx, "\n")

			if (runif(1) <= rho) {
				veta = mX[idx,] %*% vbeta + vu[i]
				if (verbose)
					cat("veta", veta, "\n")
				vy[idx] = rpois(1, exp(veta))
				if (verbose)
					cat("Generated vy[idx]", vy[idx], "\n")
			} else {
				vy[idx] = 0
			}
		}
	}
	if (verbose)
		cat("vy", vy, "\n")
	if(NA %in% vy)
		stop("NAs in vy")
	result = list(vy=vy, vu=vu)
	return(result)
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
  test_data = generate_multivariate_test_data(mX, NULL, m, n, expected_rho, expected_mu, expected_sigma2_u)
	vy = test_data$vy

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate, verbose=TRUE, method=approximation)

	print(result_var$vmu)
	expect_equal(as.vector(result_var$vmu), expected_mu, tolerance=2e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=2e-1)
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
	test_data = generate_multivariate_test_data(mX, NULL, m, n, expected_rho, expected_mu, expected_sigma2_u)
	vy = test_data$vy

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma)
	result_var = zero_infl_var(multivariate, verbose=TRUE, method=approximation)

	expect_equal(as.vector(result_var$vmu), expected_mu, tolerance=2e-1)
	expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=2e-1)
}

test_multivariate_zip_no_zeros_random_intercept <- function(approximation="gva")
{
	# Simulate data
	# Could we load test data from somewhere? I don't know that hardcoding the
	# test data into the source files is really the best idea.
	# FIXME: You have serious overflow issues
	m = 20
	ni = 10
	n = rep(ni,m)
	mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
	#print("mX=")
	#print(mX)
	#cat("dim(mX)", dim(mX), "\n")
	
	#v = c(rep(1, g), rep(0, g))
	# Indicator variables for groups
	
	#mZ = matrix(cbind(v, 1-v), sum(n), 2)
	#mZ <- matrix(0,sum(n),m)
	#count <- 0
	#for (i in 1:m) {
	#	mZ[count + (1:n[i]),i] <- 1
	#}
	#	count <- count + n[i]
	
	mZ <- kronecker(diag(1,m),rep(1,ni))
	
	#print("mZ=")
	#print(mZ)
	#cat("dim(mZ)", dim(mZ), "\n")
	
	expected_rho = 1
	expected_beta = c(2, 1)
	expected_sigma2_u = .5^2
	a_sigma = 1e-2
	b_sigma = 1e-2
	
	sigma2.beta <- 1.0E3
	
	tau = 1.0E2
	
	test_data = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=TRUE)
	vy = test_data$vy
	
	print(table(vy))
	
	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau)
	#result_var = zero_infl_var(multivariate, method="laplacian", verbose=TRUE)
	#expect_equal(as.vector(result_var$vmu[1:2]), expected_beta, tolerance=1e-1)
	#expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=2e-1)
	#print(str(result_var))
	#result_sigma2_u = (result_var$b_sigma / result_var$a_sigma)
	#expect_equal(result_sigma2_u, expected_sigma2_u, tolerance=3e-1)

	result_var = zero_infl_var(multivariate, method=approximation, verbose=TRUE)
	expect_equal(as.vector(result_var$vmu[1:2]), expected_beta, tolerance=1e-1)
	result_sigma2_u = (result_var$b_sigma / result_var$a_sigma)
	expect_equal(result_sigma2_u, expected_sigma2_u, tolerance=3e-1)
	#pdf("mult_no_zeroes_lower_bound.pdf")
	#plot(result_var$vlower_bound, type="l", xlab="Iterations", ylab="Lower bound")
  #dev.off()
}

test_multivariate_zip_half_zeros_random_intercept <- function(approximation="gva")
{
	m = 20
	ni = 10
	n = rep(ni,m)
	mX = matrix(as.vector(cbind(rep(1, sum(n)), runif(sum(n), -1, 1))), sum(n), 2)
	#print("mX=")
	#print(mX)
	#cat("dim(mX)", dim(mX), "\n")
	
	#v = c(rep(1, g), rep(0, g))
	# Indicator variables for groups
	
	#mZ = matrix(cbind(v, 1-v), sum(n), 2)
	#mZ <- matrix(0,sum(n),m)
	#count <- 0
	#for (i in 1:m) {
	#	mZ[count + (1:n[i]),i] <- 1
	#	count <- count + n[i]
	#}
	
	mZ <- kronecker(diag(1,m),rep(1,ni))
	
	#print("mZ=")
	#print(mZ)
	#cat("dim(mZ)", dim(mZ), "\n")
	
	expected_rho = 0.5
	expected_beta = c(2, 1)
	expected_sigma2_u = .5^2
	a_sigma = 1e-2
	b_sigma = 1e-2
	
	tau = 1.0E2
	
	sigma2.beta <- 1.0E3
	
	test_data = generate_multivariate_test_data(mX, mZ, m, n, expected_rho, expected_beta, expected_sigma2_u, verbose=TRUE)
	vy = test_data$vy
	
	#print(table(vy))
	#ans <- readline()

	# Test model fitting
	multivariate = create_multivariate(vy, mX, mZ, sigma2.beta, a_sigma, b_sigma, tau)
	#result_var = zero_infl_var(multivariate, method="laplacian", verbose=TRUE)
	#expect_equal(as.vector(result_var$vmu[1:2]), expected_beta, tolerance=1e-1)
	#expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=2e-1)
	#result_sigma2_u = (result_var$b_sigma / result_var$a_sigma)
	#expect_equal(result_sigma2_u, expected_sigma2_u, tolerance=3e-1)

	result_var = zero_infl_var(multivariate, method=approximation, verbose=TRUE)
	#expect_equal(as.vector(result_var$vmu[1:2]), expected_beta, tolerance=1e-1)
	#expect_equal(result_var$a_rho / (result_var$a_rho + result_var$b_rho), expected_rho, tolerance=2e-1)
  #result_sigma2_u = (result_var$b_sigma / result_var$a_sigma)
  #expect_equal(result_sigma2_u, expected_sigma2_u, tolerance=3e-1)
  pdf("mult_half_zeroes_lower_bound.pdf")
	plot(result_var$vlower_bound, type="l", xlab="Iterations", ylab="Lower bound")
  dev.off()
  return(result_var)
}

# Idea: Run each of the tests for convergence repeatedly.

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
}

#main()
