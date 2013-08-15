# variational_approximation_to_univariate_gaussian.R
# Author: Mark Greenaway

# This is the very first model I've ever tried to fit using variational approximations.
# This code is based entirely on the example from Section 21.5.1 on page 742 of
# Machine Learning: A Probabilistic Perspective

# Generate samples from true distribution x ~ N(2, 1^2) ----
true_mu = 2
true_sigma2 = 1^2
N = 100
x = rnorm(N, true_mu, true_sigma2)

mu_0 = 3
K_0 = 1/5
a_0 = 1
b_0 = 1
xbar = mean(x)

# These parameters are fixed ----
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

# Initialise ----
K_N = K_0
b_N = b_0
# Iteratively update until convergence ----
i = 1
likelihood = rep(NA, 1000)
while (TRUE) {
  K_N_last = K_N
  b_N_last = b_N
  
  K_N = (K_0 + N) * a_N/b_N
  b_N = b_0 + K_0*(E_mu2(mu_N, K_N) + mu_0^2 - 2*E_mu(mu_N) * mu_0) + 1/2*(sum(x^2 + E_mu2(mu_N, K_N) - 2*E_mu(mu_N)*x))
  #likelihood[i] = 1/2*log(1/K_N) + log(gamma(a_N)) - a_N*log(b_N)
  likelihood[i] = -1/2*log(K_N) + lgamma(a_N) - a_N*log(b_N)
  
  print(paste("Iteration", i, ":", "K_N", K_N, "b_N", b_N, "lik", likelihood[i]))
  i = i + 1
  if (isTRUE(all.equal(K_N, K_N_last)) && isTRUE(all.equal(b_N, b_N_last)))
    break
  
  if (i > 99) 
    break
}
# FIXME: The likelihood does not monotonically increase! The theory guarantees that it must,
# so something must be wrong - either in the maths or the code.

# Compare the true distribution with the distribution that we fit graphically ----
par(mfrow=c(2,2))
curve(dnorm(x, true_mu, true_sigma2), from=2-3, to=2+3)
lambda = rep(NA, N)
mu = rep(NA, N)
samples = rep(NA, N)
for(i in 1:N) {
  lambda[i] = rgamma(1, a_N, b_N)
  mu[i] = rnorm(1, mu_N, (K_N*lambda)^{-1})
  samples[i] = rnorm(1, mu[i], 1/lambda[i])
}
hist(lambda)
hist(samples, xlim=c(-1, 5))
hist(mu)

# Try John Ormerod's variant ----
# This code actually works correctly.
# Let X ~ N(mu, sigma2)
# mu ~ N(mu_mu, sigma_mu^2)
# sigma2 ~ IG(A,B)
# Initialise ----
mu_mu = 0
sigma2_mu = 10^8
A = 1/100
B = 1/100
n = 20
x = rnorm(n, 100, 15)
xbar = mean(x)
B_q_sigma2 = 1
# Iterate ----
sigma2_q_mu = (n*(A + n/2)/B_q_sigma2 + 1/sigma2_mu)^{-1}
mu_q_mu = (n*xbar*(A+n/2)/B_q_sigma2+mu_mu/sigma2_mu)*sigma2_q_mu
B_q_sigma2 = B + 1/2*(sum((x - mu_q_mu)^2)+n*sigma2_q_mu)
log_likelihood = 1/2-n/2*log(2*pi)+1/2*log(sigma2_q_mu/sigma2_mu) - ((mu_q_mu - mu_mu)^2 + sigma2_q_mu)/(2*sigma2_mu)
log_likelihood