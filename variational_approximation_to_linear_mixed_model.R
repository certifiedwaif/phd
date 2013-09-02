# variational_approximation_to_linear_mixed_model.R
# Author: Mark Greenaway
require(Matrix)

# Get Orthodont data

calcLogLik = function()
{
  logLik = 0
  logLik = logLik + 0.5*(p + sum(K))- 0.5*n*log(2*pi) - 0.5*p * log(sigma2_beta)
  logLik = logLik + 0.5 * log(det(Sigma_q_beta_u)) - (sum(mu_q_beta_u^2) + trace(Sigma_q_beta_u)) / (2*sigma2_beta)
  logLik = logLik + A_epsilon * log(B_epsilon) - (A_epsilon + 0.5 * n) * log(B_q_sigma2_epsilon)
  logLik = logLik + lgamma(A_epsilon + 0.5*n) - lgamma(A_epsilon)
  logLik = logLik + sum(A_u * log(B_u) - (A_u + 0.5*K) * log(B_q_sigma2_u) + lgamma(A_u + 0.5 * K) - lgamma(A_u))

  return(logLik)
}

trace = function(X)
{
  return(sum(diag(X)))
}

eye = function(p)
{
  return(diag(rep(1.0, p)))
}

# i varies from 1 to 27
# j varies from 1 to 4
# beta_i's all independently N(0, sigma2_beta), sigma2_u, sigma2_epsilon independently ~ IG(A, B)

X = matrix(c(1, 0,
             1, 1,
             1, 0,
             1, 1), 4, 2, byrow=TRUE)
Z = matrix(c(1, 0,
             1, 0,
             0, 1,
             0, 1), 4, 2, byrow=TRUE)
fit = lm(Y~X-1)
summary(fit)
Y = matrix(c(9, 8, 9, 8), 4, 1)
C = cbind(X, Z)
A_epsilon = 1/100
A_u = 1/100 # Quite sensitive, depends on data being standardised
B_epsilon = 1/100
B_u = 1/100
K = 2
# Should this really be five? What should it be?
# When I change it to 2, everything goes to hell!
# When I set it to 1e8, everything works like a charm. The prior was
# overwhelming the data, and funny things happened.
sigma2_beta = 1e8 
p = 2
n = 4
# Initialise ----
B_q_sigma2_epsilon = 1
B_q_sigma2_u = 1

# Debugging ideas: ----
# * Save lists of parameters
# * Visualise anything that you can
# * Go through your parameters and variables one by one and check that they're
# what they're supposed to be. That goes double for vague priors.

# Cycle ----
lastLogLik = -1e6
logLik = -1e6 + 1
paramList = list()
iteration = 1
#sink("linear_mixed_model_results.txt")
while (logLik > lastLogLik)
{
  lastLogLik = logLik
  Sigma_q_beta_u_inv = (A_epsilon + 0.5*n)/(B_q_sigma2_epsilon) * t(C) %*% C
  Sigma_q_beta_u_inv = Sigma_q_beta_u_inv + bdiag(eye(p) / sigma2_beta, (A_u + 0.5 * K)/(B_q_sigma2_u) * eye(K))
  print(Sigma_q_beta_u_inv)
  print(eigen(Sigma_q_beta_u_inv))
  Sigma_q_beta_u = solve(Sigma_q_beta_u_inv)
  print(Sigma_q_beta_u)
  print(eigen(Sigma_q_beta_u))
  mu_q_beta_u = (A_epsilon + 0.5 * n)/(B_q_sigma2_epsilon) * Sigma_q_beta_u %*% t(C) %*% Y
  #mu_q_beta_u
  B_q_sigma2_epsilon = B_epsilon + 0.5*(sum((Y - C %*% mu_q_beta_u)^2) + trace(t(C) %*% C %*% Sigma_q_beta_u))
  B_q_sigma2_u = B_u + 0.5 * (sum(mu_q_beta_u[3:4]^2) + trace(Sigma_q_beta_u[3:4, 3:4]))
  logLik = calcLogLik()
  
  params = list(Sigma_q_beta_u=Sigma_q_beta_u,
                mu_q_beta_u=mu_q_beta_u,
                B_q_sigma2_epsilon=B_q_sigma2_epsilon,
                B_q_sigma2_u=B_q_sigma2_u)
  # Build up list of parameters that we've fit.
  paramList = c(paramList, list(params))

  print(paste("Iteration:",  iteration))
  print(params)
  print(logLik)
  print(logLik - lastLogLik)
  iteration = iteration + 1
}
#sink()

# Compare against MCMC result?
