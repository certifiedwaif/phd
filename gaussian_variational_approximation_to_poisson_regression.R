# Example data ----
n = 272
p = 2
y = faithful$waiting
X = cbind(rep(1, n), faithful$eruptions)
  
# Priors ----
mubeta = rep(0, p)
Sigmabeta = diag(rep(1, p))
prior_sigma2 = 1e7

# Initialise mean field update parameters ----
# FIXME - this initialisation totally doesn't work
muqbeta = lm(log(waiting)~eruptions, data=faithful)$coef
#fit = glm(waiting~eruptions, family=poisson, data=faithful)
#muqbeta = coef(fit)

# Cycle mean field update equations until convergence ----
wv = exp(X %*% muqbeta)
W <- diag(c(wv), n, n)
Sigmaqbeta <- solve(t(X) %*% W %*% X +  prior_sigma2^(-1) * diag(rep(1, p)))

dmuqbeta = t(X) %*% y
# FIXME ----
dmuqbeta = dmuqbeta - t(X) %*% exp(X %*% muqbeta + .5 * diag(X %*% Sigmaqbeta %*% t(X)))
dmuqbeta = dmuqbeta - solve(Sigmabeta) %*% (muqbeta - mubeta)

muqbeta <- muqbeta + Sigmaqbeta %*% dmuqbeta

# Log likelihood ----
loglik = t(y) %*% X %*% muqbeta - t(rep(1, n)) %*% exp(X %*% muqbeta + .5 * diag(X %*% Sigmaqbeta %*% t(X)))
loglik = loglik - .5 * t(mubeta - muqbeta) %*% solve(Sigmabeta) %*% (mubeta - muqbeta)
loglik = loglik - .5 * p * log(2 * pi) - .5 * log(det(Sigmaqbeta)) - t(rep(1, n)) %*% lgamma(y + 1)
