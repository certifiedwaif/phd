# variational_approximation_to_linear_mixed_model.R
# Author: Mark Greenaway
require(Matrix)

# Get Orthodont data

logLik = function()
{
  logLik = 0
  logLik = logLik + 0.5*(p + sum(K))- 0.5*n*log(2*pi) - 0.5*p * log(sigma2_beta)
  logLik = logLik + 0.5 * log(det(Sigma_q_beta_u)) - (sum(mu_q_beta^2) + trace(Sigma_q_beta)) / (2*sigma2_beta)
  logLik = logLik + A_epsilon * log(B_epsilon) - (A_epsilon + 0.5 * n) * log(B_q_sigma2_epsilon)
  logLik = logLik + lgamma(A_epsilon + 0.5*n) - lgamma(A_epsilon)
  logLik = logLik + sum(A_ul * log(B_ul) - (A_ul + 0.5*K_l) * log(B_q_sigma2_ul) + lgamma(A_ul + 0.2K_l) - lgamma(A_ul))

  return logLik
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


X = matrix(c(1, 1,
             1, 2,
             2, 3,
             4, 5), 4, 2)
Z = matrix(c(1, 0,
             1, 0,
             0, 1,
             0, 1), 4, 2)
Y = matrix(c(9, 8, 7, 6), 4, 1)
C = cbind(X, Z)
A_epsilon = 1
A_u = 1
B_epsilon = 1
B_u = 1
K = 2
sigma2_beta = 5
p = 2

# Initialise
B_q_sigma2_epsilon = 1
B_q_sigma2_u = 1

# Cycle
Sigma_q_beta_u_inv = (A_epsilon + 0.5*n)/(B_q_sigma2_epsilon) * t(C) %*% C
Sigma_q_beta_u_inv = Sigma_q_beta_u_inv + bdiag(1/sigma2_beta * eye(p), (A_u + 0.5 * K)/(B_q_sigma2_u) * eye(K))
Sigma_q_beta_u = solve(Sigma_q_beta_u_inv)
mu_q_beta_u = Sigma_q_beta_u %*% t(C) %*% Y
B_q_sigma2_epsilon = B_epsilon + 0.5*(sum((Y - C %*% mu_q_beta_u)^2) + trace(t(C) %*% C %*% Sigma_q_beta_u))
B_q_sigma2_ul = B_u + 0.5 * (sum(mu_q_u^2) + trace(Sigma_q_u))
print(logLik)