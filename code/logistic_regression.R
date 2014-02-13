# variational_approximation_to_logistic_regression.R

# Test data
X = read.csv("binary.csv")
n = 400
p = 3
y = X[,1]
X = as.matrix(X[,2:4])
ones = c(rep(1, n))
# Priors
mu_beta = c(0, 0, 0)
Sigma_beta = matrix(c(1, 0, 0,
                      0, 1, 0,
                      0, 0, 1), 3, 3)

sigma = function(x)
{
  return(1 + exp(-x))^(-1)
}

lambda = function(xi)
{
  return(-tanh(xi/2)/(4*xi))
}

# Initialise
xi = rep(.1, n)
# Mean field update
stop = FALSE
ll = -Inf
iteration = 1
while (!stop) {
  Sigma_qbeta = solve(solve(Sigma_beta) - 2 * t(X) %*% diag(lambda(xi)) %*% X)
  mu_qbeta = Sigma_qbeta %*% (solve(Sigma_beta) %*% mu_beta + t(X) %*% (y - .5 * ones))
  xi = sqrt(diag(X %*% (Sigma_qbeta + mu_qbeta %*% t(mu_qbeta)) %*% t(X)))
  log_likelihood = function()
  {
    loglik = .5*log(det(Sigma_qbeta)/det(Sigma_beta)) + .5*t(mu_qbeta) %*% solve(Sigma_qbeta) %*% mu_qbeta - .5 * t(mu_beta) %*% solve(Sigma_beta) %*% mu_beta
    loglik = loglik + sum(xi/2 - log(1+exp(xi)) + (xi/4) * tanh(xi/2))
    return(loglik)
  }
  if (ll == log_likelihood()) {
    stop = TRUE
  }
  ll = log_likelihood()
  print(paste(iteration, " ", ll))
  iteration = iteration + 1
}
print(mu_qbeta)
