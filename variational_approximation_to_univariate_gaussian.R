# variational_approximation_to_univariate_gaussian.R
# Author: Mark Greenaway

# This is the very first model I've ever tried to fit using variational approximations

# Generate samples from true distribution x ~ N(25, 5^2)
true_mu = 2
true_sigma2 = 1^2
N = 100
EPSILON = 1e-20
x = rnorm(N, true_mu, true_sigma2)

mu_0 = 3
K_0 = 1/5
a_0 = 1
b_0 = 1
xbar = mean(x)

# These parameters are fixed
mu_N = (K_0 * mu_0 + N * xbar)/(K_0 + N)
a_N = a_0 + (N+1)/2

E_mu = function(mu_N)
{
  return(mu_N)
}

E_mu2 = function(mu_N, K_N)
{
  return(1/K_N + mu_N^2)
}

# Initialise
K_N = K_0
b_N = b_0
# Iteratively update until convergence
while (TRUE) {
  K_N_last = K_N
  b_N_last = b_N
  
  K_N = (K_0 + N) * a_N/b_N
  b_N = b_0 + K_0*(E_mu2(mu_N, K_N) + mu_0^2 - 2*E_mu(mu_N) * mu_0) + 1/2*(sum(x^2 + E_mu2(mu_N, K_N) - 2*E_mu(mu_N)*x))
  
  if (abs(K_N - K_N_last) < EPSILON)
    break
  
  if (abs(b_N - b_N_last) < EPSILON)
    break
  print(paste("K_N", K_N, "b_N", b_N))
  Sys.sleep(1)
}
# Compare the true distribution with the distribution that we fit
curve(dnorm(x, true_mu, true_sigma2), from=2-3, to=2+3)
