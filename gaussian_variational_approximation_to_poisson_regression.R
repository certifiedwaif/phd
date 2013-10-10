source("fit_gaussian_variational_approximation_to_poisson_regression.R")

# Example data ----
df = faithful
n = 272
p = 2
y = faithful$waiting
X = cbind(rep(1, n), faithful$eruptions)
  
# Priors ----
mubeta = rep(0, p)
Sigmabeta = diag(rep(1, p))
prior_sigma2 = 1e7

fit_poisson_regression(y=y, X=X, df=df, mubeta=mubeta, Sigmabeta=Sigmabeta, prior_sigma2=prior_sigma2)
