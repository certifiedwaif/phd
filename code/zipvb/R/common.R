
logit <- function(p) log(p/(1-p))
expit <- function(eta) 1/(1+exp(-eta))
tr <- function(mX) sum(diag(mX))

gamma_entropy <- function(alpha, beta)
{
  alpha - log(beta) + lgamma(alpha) - (alpha - 1) * digamma(alpha)
}

beta_entropy <- function(alpha, beta)
{
  lbeta(alpha, beta) - (alpha - 1) * digamma(alpha) - (beta - 1) * digamma(beta) + (alpha + beta - 2) * digamma(alpha + beta)
} 

lgammap <- function(a, p)
{
  .25 * p * (p - 1) * log(pi) + sum(lgamma(a + .5 * (1 - 1:p)))
}
