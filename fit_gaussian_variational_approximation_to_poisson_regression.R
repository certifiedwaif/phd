# fit_gaussian_variational_approximation_to_poisson_regression.R
# Initialise mean field update parameters ----
# FIXME - this initialisation totally doesn't work
fit_poisson_regression = function(y=y, X=X, df=df, mubeta=mubeta, Sigmabeta=Sigmabeta, prior_sigma2=prior_sigma2)
{
  # Initialise muqbeta ----
  muqbeta = lm(log(y)~X[,2:p], data=df)$coef
  
  # Cycle mean field update equations until convergence ----
  last_loglik = -2e12
  loglik = -1e12
  while (last_loglik < loglik) {
    last_loglik = loglik
    wv = exp(X %*% muqbeta)
    W <- diag(c(wv), n, n)
    Sigmaqbeta <- solve(t(X) %*% W %*% X +  prior_sigma2^(-1) * diag(rep(1, p)))
    
    dmuqbeta = t(X) %*% y
    dmuqbeta = dmuqbeta - t(X) %*% exp(X %*% muqbeta + .5 * diag(X %*% Sigmaqbeta %*% t(X)))
    dmuqbeta = dmuqbeta - solve(Sigmabeta) %*% (muqbeta - mubeta)
    
    muqbeta <- muqbeta + Sigmaqbeta %*% dmuqbeta
    loglik = loglikelihood(y=y, X=X, muqbeta=muqbeta, Sigmaqbeta=Sigmaqbeta, mubeta=mubeta, Sigmabeta=Sigmabeta, prior_sigma2=prior_sigma2)
  }
}

# Log likelihood ----
loglikelihood = function(y=y, X=X, muqbeta=muqbeta, Sigmaqbeta=Sigmaqbeta, mubeta=mubeta, Sigmabeta=Sigmabeta, prior_sigma2=prior_sigma2)
{
  loglik = t(y) %*% X %*% muqbeta - t(rep(1, n)) %*% exp(X %*% muqbeta + .5 * diag(X %*% Sigmaqbeta %*% t(X)))
  loglik = loglik - .5 * t(mubeta - muqbeta) %*% solve(Sigmabeta) %*% (mubeta - muqbeta)
  loglik = loglik - .5 * p * log(2 * pi) - .5 * log(det(Sigmaqbeta)) - t(rep(1, n)) %*% lgamma(y + 1)
  
  return(loglik)
}