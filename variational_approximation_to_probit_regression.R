# variational_approximation_to_probit_regression.R

logLik = function(y, mu_q_a, mu_q_beta, X, Sigma_beta, p) {
  liks = y * log(pnorm(X %*% mu_q_beta)) + (1 - y)*log(1 - pnorm(X %*% mu_q_beta))
  lik= sum(liks) - .5*log(det(Sigma_beta %*% t(X) %*% X + diag(rep(1, p))))
  lik = lik -.5*t(mu_q_beta - mu_beta) %*% solve(Sigma_beta) %*% (mu_q_beta - mu_beta)
  return(sum(liks))
}

y = as.matrix(c(1, 0, 1, 1), byrow=TRUE)
X = matrix(c(1, 0, 
      1, 0,
      1, 1,
      1, 1), 4, 2, byrow=TRUE)
n = 4
p = 2
Sigma_beta = diag(rep(1, p))
mu_beta = as.matrix(c(0, 0), byrow=TRUE)
# Initialise ----
mu_q_a = rep(0, 4)
# Cycle ----
mu_q_beta = solve(t(X) %*% X + solve(Sigma_beta)) %*% (t(X) %*% mu_q_a + solve(Sigma_beta) %*% mu_beta)
mu_q_a = X %*% mu_q_beta + (dnorm(X %*% mu_q_beta))/(pnorm(X %*% mu_q_beta)^y * pnorm(X %*% mu_q_beta) - rep(1, n))^(rep(1, n) - y)
logLik(y, mu_q_a, mu_q_beta, X, Sigma_beta, p)
