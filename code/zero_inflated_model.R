# variational_approximation_to_zero_inflated_model.R

generate_test_data = function (n, rho, lambda)
{
  x = rep(NA, n)
  for (i in 1:n) {
    if (runif(1) >= rho) {
      x[i] = rpois(1, lambda)
    } else {
      x[i] = 0
    }
  }
  return(x)
}

logit = function(p) log(p/(1-p))
expit = function(eta) 1/(1+exp(-eta))

zero_infl_mcmc = function(iterations, x, a, b)
{
  # Initialise
  n = length(x)
  lambda = rep(NA, iterations+1)
  rho = rep(NA, iterations+1)
  lambda[1] = mean(x)
  rho[1] = sum(x == 0)/length(x)
  eta = rep(NA, n)
  r = rep(NA, n)
  # Build up a list of the zero observations. We'll only generate
  # r[i]'s for those. The non-zero observations will always have r[i] = 1.
  r[x != 0] = 1
  zero_observation_idx = which(x == 0)
  # Iterate ----
  for (i in 1:iterations) {
    arg = exp(-lambda[i] + logit(rho[i]))
    eta = arg/(1 + arg)
    
    # The following code is correct, but 100 times slower as R is interpreted.
    # So we use the vectorised code instead.
    #for (j in 1:length(zero_observation_idx)) {
    #  r[zero_observation_idx[j]] = rbinom(1, 1, eta)
    #}
    
    r[zero_observation_idx] = rbinom(length(zero_observation_idx), 1, eta) 
    lambda[i+1] = rgamma(1, a + sum(x), b + sum(r))
    rho[i+1] = rbeta(1, sum(r) + 1, n - sum(r) + 1)
  }
  return(list(lambda=lambda, rho=rho))
}

calculate_lower_bound = function(x, p, a_lambda, b_lambda, a_rho, b_rho)
{
  #browser()

  gamma_entropy = function(alpha, beta)
  {
    alpha - log(beta) + lgamma(alpha) - (alpha-1) * digamma(alpha)
  }
  
  beta_entropy = function(alpha, beta)
  {
    lbeta(alpha, beta) - (alpha-1)*digamma(alpha) - (beta-1)*digamma(beta) + (alpha+beta-2)*digamma(alpha+beta)
  }  
  
  zero.set <- which(x==0)
  
  E_lambda = a_lambda/b_lambda
  E_log_lambda = digamma(a_lambda) - log(b_lambda)  

  E_r = ifelse(x == 0, p, 1)
  E_xi_log_lambda_r = ifelse(x == 0, 0, x*E_log_lambda)
  

  
  E_log_rho = digamma(a_rho) - digamma(a_rho + b_rho)
  E_log_one_minus_rho = digamma(b_rho) - digamma(a_rho + b_rho)

  E_log_q_r = (p[zero.set]*log(p[zero.set]) + (1-p[zero.set])*log(1-p[zero.set]))
  E_log_q_lambda = -gamma_entropy(a_lambda, b_lambda)
  E_log_q_rho = -beta_entropy(a_rho, b_rho) # FIXME: Why is this entropy positive?
  
  result = a * log(b) + (a-1) * E_log_lambda - b * E_lambda - lgamma(a)
  result = result - E_lambda * sum(E_r)
  result = result + sum(E_xi_log_lambda_r) - sum(lgamma(x+1))
  result = result + sum(E_r) * E_log_rho + sum(1 - E_r) * E_log_one_minus_rho
  result = result - sum(E_log_q_r) - E_log_q_lambda - E_log_q_rho
  
  return(result)
}

zero_infl_var = function(x, a, b)
{
  #browser()
  n = length(x)
  #p = rep(sum(x == 0)/length(x), n)
  p = rep(1, n)
  eta = rep(NA, n)
  a_lambda = a + sum(x)
  
  lower_bound_vector<- c()
  
  lower_bound = -Inf
  last_lower_bound = -Inf
  iteration = 1
  # Iterate ----
  while (iteration == 1 || lower_bound > last_lower_bound) {
    last_lower_bound = lower_bound
    
    b_lambda = b + sum(p)
    a_rho = 1 + sum(p)
    b_rho = n - sum(p) + 1
    
    # TODO: You could vectorise the loop below.
    for (i in 1:n) {
      if (x[i] != 0)
        p[i] = 1.0
      else {
        eta[i] = - a_lambda/b_lambda + digamma(a_rho) - digamma(b_rho)
        p[i] = expit(eta[i])
      }
    }
    # TODO: Lower bound? ----
    lower_bound = calculate_lower_bound(x, p, a_lambda, b_lambda, a_rho, b_rho)
    
    lower_bound_vector[iteration] <- lower_bound
    
    if (iteration>2) {
    	plot(lower_bound_vector,type="l")
    }
    ans <- readline()
    
    
    #print(lower_bound)
    params = list(a_lambda=a_lambda, b_lambda=b_lambda, a_rho=a_rho, b_rho=b_rho)
    #print(params)
    cat("Iteration ", iteration, ": lower bound ", lower_bound, " difference ",
        lower_bound - last_lower_bound, " parameters ", "a_lambda", a_lambda,
        "b_lambda", b_lambda, "a_rho", a_rho, "b_rho", b_rho, "\n")
    iteration = iteration + 1
  }
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
  result_var = zero_infl_var(10, x, a, b)
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
    for (lambda in seq(5, 100, 5))
      print(check_accuracy(n, rho, lambda))
}
#main_check_accuracy()


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
multivariate = function(y, X, Z)
{
  SigmaInv = solve(Sigma)
  D_b = sum(t(y - B(1, beta, mu, Lambda)) %*% X)
  D_Sigma = 0.5 * sum(t(vec(SigmaInv %*% (m[i] %*% t(mu[i]) + Lambda[i]) %*% SigmaInv)) %*% D)
  D_mu[i] = t(y - B(1, mu[i], Lambda[i])) %*% Z[i] - t(mu[i]) %*% SigmaInv
  D_Lambda[i] = .5 t(vec(Lambda[i]Inv - SigmaInv - t(Z[i]) %*% diag(B(2, beta, mu[i], Lambda[i])) %*% Z[i])) %*% D

  H_b_u[i] = -t(X[i]) %*% diag(B(2, beta, mu[i], Lambda[i])) %*% Z[i]
  H_b_Lambda[i] = - .5 * t(X[i]) %*% diag(B(3, beta, mu[i], Lambda[i])) %*% O(Z[i]) %*% D
  H_Sigma_mu[i] = t(D) %*% kron(SigmaInv %*% mu[i], SigmaInv)
  H_mu[i]_mu[i] = -t(Z[i]) %*% diag(B(2, beta, mu[i], Lambda[i])) %*% Z[i] - SigmaInv
  H_mu[i]_Lambda[i] = - .5 * t(Z[i] %*% diag(B(3, beta, mu[i], Lambda[i]))) %*% O(Z[i]) %*% D
  Lambda[i]Inv = solve(Lambda[i])
  H_Lambda[i]_Lambda[i] = -.25 * t(D) %*% (t(O(Z[i])) %*% diag(4, beta, mu[i], Lambda[i])) + 2(kron(Lambda[i]Inv, Lambda[i]Inv)) %*% D
}