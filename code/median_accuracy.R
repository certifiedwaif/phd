# median_accuracy.R
source("test_zero_inflated_model.R")

# Graph of Var_q(theta) against Var(theta|y)
# How to get this?
# Run fits for a range of theta values?
var_fn = function(vbeta)
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
  approximation = "gva"
  result_var = zero_infl_var(multivariate, method=approximation, verbose=FALSE)
  return(result_var)
}
var_fn(c(1, 0))
var_fn(c(1, 1))
var_fn(c(1, 2))
var_fn(c(1, 3))
var_fn(c(1, 4))

# Graph of E_q(theta) against E(theta|y)

# For most important parameters: beta_1, beta_2