
logit <- function(p) log(p/(1-p))
expit <- function(eta) 1/(1+exp(-eta))
tr <- function(mX) sum(diag(mX))
