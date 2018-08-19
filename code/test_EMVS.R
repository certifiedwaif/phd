# Linear regression with p>n variables
library(EMVS)

n = 100
p = 1000
X = matrix(rnorm(n * p), n, p)
beta = c(1, 2, 3, rep(0, p - 3))
Y = X[, 1] * beta[1] + X[, 2] * beta[2] + X[, 3] * beta[3] + rnorm(n)
Y = Y - mean(Y)
v0 = seq(0.1, 1, by = 0.1)
v1 = 1000
beta_init = rep(1, p)
a = b = 1
epsilon = 10^{
    -5
}

result = EMVS(Y, X, v0 = v0, v1 = v1, type = "betabinomial", beta_init = beta_init, 
    sigma_init = 1, epsilon = epsilon, a = a, b = b)

EMVSplot(result, "both", FALSE)

EMVSbest(result)
