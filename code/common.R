
logit <- function(p) log(p/(1-p))
expit <- function(eta) 1/(1+exp(-eta))
tr <- function(mX) sum(diag(mX))

my.diag <- function(mA) {
    k <- ncol(mA) 
    Dinds <- k*((1:k)-1) + (1:k)
    val <- mA[Dinds]
    return(val)
}

my.tr <- function(mA) {
	return( sum(my.diag(mA)) ) 
}
