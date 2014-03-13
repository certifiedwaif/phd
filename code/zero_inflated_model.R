# variational_approximation_to_zero_inflated_model.R

generate_test_data = function (n, rho, lambda)
{
  vx = rep(NA, n)
  for (i in 1:n) {
    if (runif(1) >= rho) {
      vx[i] = rpois(1, lambda)
    } else {
      vx[i] = 0
    }
  }
  return(vx)
}

logit = function(p) log(p/(1-p))
expit = function(eta) 1/(1+exp(-eta))

zero_infl_mcmc = function(iterations, vx, a, b)
{
  # Initialise
  n = length(vx)
  vlambda = rep(NA, iterations+1)
  vrho = rep(NA, iterations+1)
  vlambda[1] = mean(vx)
  vrho[1] = sum(vx == 0)/length(vx)
  veta = rep(NA, n)
  vr = rep(NA, n)
  # Build up a list of the zero observations. We'll only generate
  # r[i]'s for those. The non-zero observations will always have r[i] = 1.
  vr[x != 0] = 1
  zero_observation_idx = which(vx == 0)
  # Iterate ----
  for (i in 1:iterations) {
    varg = exp(-vlambda[i] + logit(vrho[i]))
    veta = varg/(1 + varg)
    
    # The following code is correct, but 100 times slower as R is interpreted.
    # So we use the vectorised code instead.
    #for (j in 1:length(zero_observation_idx)) {
    #  r[zero_observation_idx[j]] = rbinom(1, 1, eta)
    #}
    
    vr[zero_observation_idx] = rbinom(length(zero_observation_idx), 1, veta) 
    vlambda[i+1] = rgamma(1, a + sum(vx), b + sum(vr))
    vrho[i+1] = rbeta(1, sum(vr) + 1, n - sum(vr) + 1)
  }
  return(list(vlambda=vlambda, vrho=vrho))
}

calculate_lower_bound = function(vx, vp, a_lambda, b_lambda, a_rho, b_rho)
{
  gamma_entropy = function(alpha, beta)
  {
    alpha - log(beta) + lgamma(alpha) - (alpha-1) * digamma(alpha)
  }
  
  beta_entropy = function(alpha, beta)
  {
    lbeta(alpha, beta) - (alpha-1)*digamma(alpha) - (beta-1)*digamma(beta) + (alpha+beta-2)*digamma(alpha+beta)
  }  
  
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

zero_infl_var = function(vx, a, b, verbose=FALSE, plot_lower_bound=FALSE)
{
	# Initialise
  n = length(vx)
  vp = rep(1, n)
  veta = rep(NA, n)
  a_lambda = a + sum(vx)
	zero.set = which(vx == 0)
	nonzero.set = which(vx != 0)

  vlower_bound <- c()
  
  lower_bound = -Inf
  last_lower_bound = -Inf
  i = 1
  # Iterate ----
  while (i == 1 || vlower_bound[iteration] > vlower_bound[iteration-1]) {
	# Update paramerer for q_lambda
    b_lambda = b + sum(vp)

	# Update parameters for q_rho
    a_rho = 1 + sum(vp)
    b_rho = n - sum(vp) + 1
    
	# Update parameters for q_vr
    vp[nonzero.set] = 1.0
    vp[zero.set] = expit(-a_lambda/b_lambda + digamma(a_rho) - digamma(b_rho))
    
    vlower_bound[i] <- calculate_lower_bound(vx, vp, a_lambda, b_lambda, a_rho, b_rho)
    
    params = list(a_lambda=a_lambda, b_lambda=b_lambda, a_rho=a_rho, b_rho=b_rho)
	if (verbose && i > 1)
		cat("Iteration ", iteration, ": lower bound ", vlower_bound, " difference ",
			vlower_bound[i] - vlower_bound[i-1], " parameters ", "a_lambda", a_lambda,
			"b_lambda", b_lambda, "a_rho", a_rho, "b_rho", b_rho, "\n")
    i = i + 1
  }
  if (plot_lower_bound)
  	plot(lower_bound_vector,type="l")
  return(params)
}

# Calculate accuracy ----
# Approximate the L1 norm between the variational approximation and
# the MCMC approximation
calculate_accuracy = function(result_mcmc, result_var)
{
  density_mcmc_rho = density(result_mcmc$rho)
  integrand = function(x)
  {
    fn = splinefun(density_mcmc_rho$x, density_mcmc_rho$y)
    return(abs(fn(x) - dbeta(x, result_var$a_rho, result_var$b_rho)))
  }
  integrate(integrand, min(density_mcmc_rho$x), max(density_mcmc_rho$x), subdivisions = length(density_mcmc_rho$x))
  
  density_mcmc_lambda = density(result_mcmc$lambda)
  integrand = function(x)
  {
    fn = splinefun(density_mcmc_lambda$x, density_mcmc_lambda$y)
    return(abs(fn(x) - dgamma(x, result_var$a_lambda, result_var$b_lambda)))
  }
  result = integrate(integrand, min(density_mcmc_lambda$x), max(density_mcmc_lambda$x), subdivisions = length(density_mcmc_lambda$x))
  accuracy = 1 - .5 * result$value
  return(accuracy)
}

check_accuracy = function(n, rho, lambda)
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
main_check_accuracy = function()
{
  n = 1000
  rho = .5
  lambda = 100
  
  for (rho in .1*1:9)
    for (lambda in seq(5, 10, 0.05))
      print(check_accuracy(n, rho, lambda))
}

# MCMC ----
#iterations = 1e6
#burnin = 1e3

#start = Sys.time()
#result_mcmc = zero_infl_mcmc(iterations+burnin, x, a, b)
#mcmc_runtime = Sys.time() - start
# Throw away burn-in samples.
# Brackets turn out to be incredibly important here!!!
#result_mcmc$rho = result_mcmc$rho[(burnin+1):(burnin+iterations+1)]
#result_mcmc$lambda = result_mcmc$lambda[(burnin+1):(burnin+iterations+1)]
# 1000 iterations in 0.07461715 seconds.

# Extension to multivariate case ----
# multivariate = function(y, X, Z)
# {
#   SigmaInv = solve(Sigma)
#   D_b = sum(t(y - B(1, beta, mu, Lambda)) %*% X)
#   D_Sigma = 0.5 * sum(t(vec(SigmaInv %*% (m[i] %*% t(mu[i]) + Lambda[i]) %*% SigmaInv)) %*% D)
#   D_mu[i] = t(y - B(1, mu[i], Lambda[i])) %*% Z[i] - t(mu[i]) %*% SigmaInv
#   D_Lambda[i] = .5 t(vec(Lambda[i]Inv - SigmaInv - t(Z[i]) %*% diag(B(2, beta, mu[i], Lambda[i])) %*% Z[i])) %*% D
# 
#   H_b_u[i] = -t(X[i]) %*% diag(B(2, beta, mu[i], Lambda[i])) %*% Z[i]
#   H_b_Lambda[i] = - .5 * t(X[i]) %*% diag(B(3, beta, mu[i], Lambda[i])) %*% O(Z[i]) %*% D
#   H_Sigma_mu[i] = t(D) %*% kron(SigmaInv %*% mu[i], SigmaInv)
#   H_mu_mu[i, i] = -t(Z[i]) %*% diag(B(2, beta, mu[i], Lambda[i])) %*% Z[i] - SigmaInv
#   H_mu_Lambda[i, i] = - .5 * t(Z[i] %*% diag(B(3, beta, mu[i], Lambda[i]))) %*% O(Z[i]) %*% D
#   LambdaInv[i] = solve(Lambda[i])
#   H_Lambda_Lambda[i, i] = -.25 * t(D) %*% (t(O(Z[i])) %*% diag(4, beta, mu[i], Lambda[i])) + 2(kron(LambdaInv[i], LambdaInv[i])) %*% D
#}
