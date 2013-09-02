# variational_approximation_to_linear_mixed_model.R
# Author: Mark Greenaway

# Get Orthodont data
# i varies from 1 to 27
# j varies from 1 to 4
# beta_i's all independently N(0, sigma2_beta), sigma2_u, sigma2_epsilon independently ~ IG(A, B)


# My own test case ----
X = matrix(c(1, 0,
             1, 1,
             1, 0,
             1, 1), 4, 2, byrow=TRUE)
Z = matrix(c(1, 0,
             1, 0,
             0, 1,
             0, 1), 4, 2, byrow=TRUE)
Y = matrix(c(9, 8, 9, 8), 4, 1)
A_epsilon = 1/100
A_u = 1/100
B_epsilon = 1/100
B_u = 1/100
K = 2
sigma2_beta = 5
p = 2
n = 4

source("fit_linear_mixed_model.R")
fit_linear_mixed_model(Y, X, Z, A_epsilon, B_epsilon, A_u, B_u, K, sigma2_beta, p, n)

require(nlme)
data(Orthodont)
Y = Orthodont$distance
X = model.matrix(distance~age+factor(Sex)-1, Orthodont)
Z = model.matrix(distance~factor(Subject)-1, Orthodont)
n = 108
p = 3
dim(X)
dim(Z)
fit_linear_mixed_model(Y, X, Z, A_epsilon, B_epsilon, A_u, B_u, K, sigma2_beta, p, n)
# Compare against MCMC result?
# Or frequentist fit.
