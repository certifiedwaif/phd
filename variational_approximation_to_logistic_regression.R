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
# Initialise
xi = rep(.1, n)
# Mean field update
A = function(xi)
{
  return(-tanh(xi/2)/(4*xi))
}
Sigma_qbeta = solve(solve(Sigma_beta) - 2 * t(X) %*% diag(A(xi)) %*% X)
mu_qbeta = Sigma_qbeta %*% (solve(Sigma_beta) %*% mu_beta + t(X) %*% (y - .5 * ones))
xi = sqrt(diag(X %*% (Sigma_qbeta + mu_qbeta %*% t(mu_qbeta)) %*% t(X)))
print(mu_qbeta)
