library(gsl)

trapint <- function(xgrid, fgrid)
{
	ng <- length(xgrid)
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
	integ <- sum(xvec * fvec)/2
	return(integ)
}

n <- 50
p <- 30
alpha <- 1
a <- -3/4
b <- (n - p) / 2 - a - 2
c <- (n-1) / 2
R2 <- .5

sigma2 <- seq(0.1, 10, .01)
log_K <- c * log(n/2) + (b + 1) * log(1 - R2) - lgamma(c) - 0.5 * log(2*pi / n)
integrand <- exp(log_K-(n*alpha^2)/(2*sigma2) * -(c+1) * log(sigma2) - n/(2 * sigma2)) * hyperg_1F1(b + 1, c, (n * R2)/(2*sigma2))

trapint(sigma2, integrand)

exp(lgamma(c + 0.5) + (b + 1) * log(1 - R2) - lgamma(c) - 0.5 * log(pi) -n/2 * log(1 + alpha^2)) * hyperg_2F1(b + 1, c + 0.5, c, R2 / (1 + alpha^2))