# variational_approximation_to_zero_inflated_model.R
source("zero_inflated_poisson_linear_model.R")

logit <- function(p) log(p/(1-p))
expit <- function(eta) 1/(1+exp(-eta))

zero_infl_mcmc.univariate <- function(iterations, vx, a, b)
{
	# Initialise
	n = length(vx)
	zero.set = which(vx == 0)
	nonzero.set = which(vx != 0)
	vlambda = rep(NA, iterations+1)
	vrho = rep(NA, iterations+1)
	vlambda[1] = mean(vx)
	vrho[1] = length(zero.set)/length(vx)
	veta = rep(NA, n)
	vr = rep(NA, n)
	# Build up a list of the zero observations. We'll only generate
	# r[i]'s for those. The non-zero observations will always have r[i] = 1.
	vr[nonzero.set] = 1
	# Iterate ----
	for (i in 1:iterations) {
		varg = exp(-vlambda[i] + logit(vrho[i]))
		veta = varg/(1 + varg)
		
		# The following code is correct, but 100 times slower as R is interpreted.
		# So we use the vectorised code instead.
		#for (j in 1:length(zero_observation_idx)) {
		#	r[zero_observation_idx[j]] = rbinom(1, 1, eta)
		#}
		
		vr[zero.set] = rbinom(length(zero.set), 1, veta) 
		vlambda[i+1] = rgamma(1, a + sum(vx), b + sum(vr))
		vrho[i+1] = rbeta(1, sum(vr) + 1, n - sum(vr) + 1)
	}
	return(list(vlambda=vlambda, vrho=vrho))
}

zero_infl_mcmc.multivariate <- function(multivariate)
{
	stop("Not implemented yet")
}

zero_infl_mcmc <- function(object)
{
	UseMethod("zero_infl_mcmc", object)
}

gamma_entropy <- function(alpha, beta)
{
	alpha - log(beta) + lgamma(alpha) - (alpha-1) * digamma(alpha)
}
	
beta_entropy <- function(alpha, beta)
{
	lbeta(alpha, beta) - (alpha-1)*digamma(alpha) - (beta-1)*digamma(beta) + (alpha+beta-2)*digamma(alpha+beta)
}	
	
calculate_lower_bound.univariate <- function(univariate)
{
	vx = univariate$vx
	vp = univariate$vp
	a_lambda = univariate$a_lambda
	b_lambda = univariate$b_lambda
	a_rho = univariate$a_rho
	b_rho = univariate$b_rho

	zero.set <- which(vx==0)
	
	E_lambda = a_lambda/b_lambda
	E_log_lambda = digamma(a_lambda) - log(b_lambda)	
	
	E_r = ifelse(vx == 0, vp, 1)
	E_xi_log_lambda_r = ifelse(vx == 0, 0, vx*E_log_lambda)
	
	E_log_rho = digamma(a_rho) - digamma(a_rho + b_rho)
	E_log_one_minus_rho = digamma(b_rho) - digamma(a_rho + b_rho)
	
	E_log_q_r = (vp[zero.set]*log(vp[zero.set]) + (1-vp[zero.set])*log(1-vp[zero.set]))
	E_log_q_lambda = -gamma_entropy(a_lambda, b_lambda)
	E_log_q_rho = -beta_entropy(a_rho, b_rho) # FIXME: Why is this entropy positive?
	
	result = a_lambda * log(b_lambda) + (a_lambda-1) * E_log_lambda - b_lambda * E_lambda - lgamma(a_lambda)
	result = result - E_lambda * sum(E_r)
	result = result + sum(E_xi_log_lambda_r) - sum(lgamma(vx+1))
	result = result + sum(E_r) * E_log_rho + sum(1 - E_r) * E_log_one_minus_rho
	result = result - sum(E_log_q_r) - E_log_q_lambda - E_log_q_rho
	
	return(result)
}

calculate_lower_bound.multivariate <- function(multivariate)
{
	# TODO: Re-use what you can from the univariate lower bound, rather than
  # duplicating code.
	# Take the lower bound returned by Poisson fit and combine it
	vy = multivariate$vy
	vp = multivariate$vp
	#a_lambda = multivariate$a_lambda
	#b_lambda = multivariate$b_lambda
	a_rho = multivariate$a_rho
	b_rho = multivariate$b_rho

	zero.set <- which(vy==0)
	
	#E_lambda = a_lambda/b_lambda
	#E_log_lambda = digamma(a_lambda) - log(b_lambda)	
	
	E_r = ifelse(vy == 0, vp, 1)
	#E_xi_log_lambda_r = ifelse(vy == 0, 0, vy*E_log_lambda)
	
	E_log_rho = digamma(a_rho) - digamma(a_rho + b_rho)
	E_log_one_minus_rho = digamma(b_rho) - digamma(a_rho + b_rho)
	
	E_log_q_r = (vp[zero.set]*log(vp[zero.set]) + (1-vp[zero.set])*log(1-vp[zero.set]))
	#E_log_q_lambda = -gamma_entropy(a_lambda, b_lambda)
	E_log_q_rho = -beta_entropy(a_rho, b_rho) # FIXME: Why is this entropy positive?
	
	result = 0
	#result = a_lambda * log(b_lambda) + (a_lambda-1) * E_log_lambda - b_lambda * E_lambda - lgamma(a_lambda)
	#result = result - E_lambda * sum(E_r)
	#result = result + sum(E_xi_log_lambda_r) - sum(lgamma(vx+1))
	result = result + sum(E_r) * E_log_rho + sum(1 - E_r) * E_log_one_minus_rho
	result = result - sum(E_log_q_r) #- E_log_q_lambda
	result = result - E_log_q_rho
	result = result + multivariate$f
	
	return(result)
}

calculate_lower_bound <- function(object)
{
	UseMethod("calculate_lower_bound", object)
}

expected_lambda.univariate <- function(univariate)
{
	univariate$a_lambda/univariate$b_lambda  
}

expected_lambda.multivariate <- function(multivariate)
{
	expected_lambda = exp(multivariate$vmu%*%multivariate$mC)
	return(expected_lambda)
}

expected_lambda <- function(object)
{
	UseMethod("expected_lambda", object)
}

create_univariate <- function(vx, a, b)
{
	# Initialise
	n = length(vx)
	vp = rep(1, n)
	a_lambda = a + sum(vx)
	b_lambda = b
	univariate = list(vx=vx, vp=vp, a_lambda=a_lambda, b_lambda=b_lambda)
	class(univariate) = "univariate"
	return(univariate)
}

zero_infl_var.univariate <- function(univariate, trace=FALSE, plot_lower_bound=FALSE)
{
	vx = univariate$vx
	a = univariate$a_lambda
	b = univariate$b_lambda
	n = length(univariate$vx)
	zero.set = which(univariate$vx == 0)
	nonzero.set = which(univariate$vx != 0)
	vlower_bound <- c()
	
	i = 0
	# Iterate ----
	while (i <= 2 || vlower_bound[i] > vlower_bound[i-1]) {
		i = i+1

		# Update parameter for q_lambda
		univariate$b_lambda = b + sum(univariate$vp)
		
		# Update parameters for q_rho
		univariate$a_rho = 1 + sum(univariate$vp)
		univariate$b_rho = n - sum(univariate$vp) + 1
		
		# Update parameters for q_vr
		univariate$vp[zero.set] = expit(-expected_lambda(univariate) + digamma(univariate$a_rho) - digamma(univariate$b_rho))

		#vlower_bound[i] <- calculate_lower_bound(vx, vp, a_lambda, b_lambda, a_rho, b_rho)
		vlower_bound[i] <- calculate_lower_bound(univariate)
		
		if (trace && i > 1)
			cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
					vlower_bound[i] - vlower_bound[i-1], " parameters ", "a_lambda", univariate$a_lambda,
					"b_lambda", univariate$b_lambda, "a_rho", univariate$a_rho, "b_rho", univariate$b_rho, "\n")
	}

	if (plot_lower_bound)
		plot(lower_bound_vector,type="l")

	params = list(a_lambda=univariate$a_lambda, b_lambda=univariate$b_lambda, a_rho=univariate$a_rho, b_rho=univariate$b_rho)
	return(params)
}

create_multivariate <- function(vy, mX, mZ, a_sigma, b_sigma)
{
	# Initialise
	n = length(vy)
	vp = rep(1, n)
	mC = cbind(mX, mZ)
	vnu = rep(0, ncol(mC))

	mLambda = diag(rep(1, ncol(mC)))
	mSigma.inv = diag(rep(b_sigma/a_sigma, ncol(mC)))
	prior = list(a_sigma=a_sigma, b_sigma=b_sigma)
	multivariate = list(vy=vy, mX=mX, mZ=mZ, mC=mC, vp=vp, vnu=vnu,
						a_sigma=a_sigma, b_sigma=b_sigma,
						mLambda=mLambda,
						prior=prior,
						mSigma.inv=mSigma.inv)
	class(multivariate) = "multivariate"
	return(multivariate)
}

zero_infl_var.multivariate <- function(m, trace=FALSE, plot_lower_bound=FALSE)
{
	# Initialise
	n = length(m$vy)
	zero.set = which(m$vy == 0)
	nonzero.set = which(m$vy != 0)
	vlower_bound <- c()
	
	i = 0
	# Iterate ----
	while (i <= 2 || vlower_bound[i] > vlower_bound[i-1]) {
		i = i+1

		# Update parameter for q_lambda
		# FIXME: We don't update, we maximise the log-likelihood using z_i = r_i y_i
		# and r_i
		# Maximise the Gaussian Variational Approximation using
		# Dr Ormerod's Poisson mixed model code
		fit1 = fit.Lap(m$vnu, m$vy, m$vp, m$mC, m$mSigma.inv, m$mLambda)
		fit2 = fit.GVA(fit1$vnu, fit1$mLambda, m$vy, m$vp, m$mC, m$mSigma.inv, "L-BFGS-B")
		m$vnu = fit2$vnu
		m$mLambda = fit2$mLambda
		m$f = fit2$res$value

		# Update parameters for q_sigma_u^2
		# a_sigma is fixed
		m$a_sigma = m$a_sigma + n/2
		vu = m$vnu[(ncol(m$mX)+1):ncol(m$mC)]
		# We know that mSigma = sigma_u^2 I. We should exploit this knowledge
		tr_mSigma = ncol(m$mZ) * m$prior$a_sigma/m$prior$b_sigma
		m$b_sigma = m$prior$b_sigma + sum(vu^2)/2 + (tr_mSigma)/2

		# Update parameters for q_rho
		m$a_rho = 1 + sum(m$vp)
		m$b_rho = n - sum(m$vp) + 1
		
		# Update parameters for q_vr
		#print(dim(m$mC))
		#print(dim(m$vnu))
		m$vp[zero.set] = expit(-exp(m$mC[zero.set,]%*%m$vnu) + digamma(m$a_rho) - digamma(m$b_rho))
		# FIXME: We get zeros in here sometimes, which plays havoc with the
		# lower bound calculation.
    
		#vlower_bound[i] <- calculate_lower_bound(vx, vp, a_lambda, b_lambda, a_rho, b_rho)
		vlower_bound[i] <- calculate_lower_bound(m)
		#print(m$vnu)
		#print(m$vp)
		#print(m$a_rho)
		#print(m$b_rho)
		#cat("End of iteration", i, "\n")
		
		if (trace && i > 1)
			cat("Iteration ", i, ": lower bound ", vlower_bound[i], " difference ",
					vlower_bound[i] - vlower_bound[i-1], " parameters ", "vnu", m$vnu,
					"a_rho", m$a_rho, "b_rho", m$b_rho, "\n")
	}

	if (plot_lower_bound)
		plot(lower_bound_vector,type="l")

	params = list(vnu=m$vnu, a_rho=m$a_rho, b_rho=m$b_rho)
	return(params)
}

zero_infl_var <- function(object, trace=FALSE, plot_lower_bound=FALSE)
{
	UseMethod("zero_infl_var", object)
}

# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
calculate_accuracy <- function(result_mcmc, result_var)
{
	density_mcmc_rho = density(result_mcmc$rho)
	integrand <- function(x)
	{
		fn = splinefun(density_mcmc_rho$x, density_mcmc_rho$y)
		return(abs(fn(x) - dbeta(x, result_var$a_rho, result_var$b_rho)))
	}
	integrate(integrand, min(density_mcmc_rho$x), max(density_mcmc_rho$x), subdivisions = length(density_mcmc_rho$x))
	
	density_mcmc_lambda = density(result_mcmc$lambda)
	integrand <- function(x)
	{
		fn = splinefun(density_mcmc_lambda$x, density_mcmc_lambda$y)
		return(abs(fn(x) - dgamma(x, result_var$a_lambda, result_var$b_lambda)))
	}
	result = integrate(integrand, min(density_mcmc_lambda$x), max(density_mcmc_lambda$x), subdivisions = length(density_mcmc_lambda$x))
	accuracy = 1 - .5 * result$value
	return(accuracy)
}

check_accuracy <- function(n, rho, lambda)
{
	x = generate_test_data(n, rho, lambda)
	
	a = 10
	b = 10
	
	# MCMC ----
	iterations = 1e6
	burnin = 1e3
	
	start = Sys.time()
	result_mcmc = zero_infl_mcmc(iterations+burnin, x, a, b)
	mcmc_runtime = Sys.time() - start
	# Throw away burn-in samples.
	# Brackets turn out to be incredibly important here!!!
	result_mcmc$rho = result_mcmc$rho[(burnin+1):(burnin+iterations+1)]
	result_mcmc$lambda = result_mcmc$lambda[(burnin+1):(burnin+iterations+1)]
	# 1000 iterations in 0.07461715 seconds.
	
	# Variational approximation ----
	start = Sys.time()
	result_var = zero_infl_var(x, a, b)
	var_runtime = Sys.time() - start
	# Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
	# or 200 iterations a second.
	var_lambda = result_var$a_lambda / result_var$b_lambda
	var_rho = result_var$a_rho / (result_var$a_rho + result_var$b_lambda)
	
	accuracy = calculate_accuracy(result_mcmc, result_var)
	return(list(n=n, rho=rho, lambda=lambda, accuracy=accuracy))
}

# Check accuracy of the approximation for a range of parameter values ----
main_check_accuracy <- function()
{
	n = 1000
	rho = .5
	lambda = 100
	
	for (rho in .1*1:9)
		for (lambda in seq(5, 10, 0.05))
			print(check_accuracy(n, rho, lambda))
}
