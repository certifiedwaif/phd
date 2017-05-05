library(gsl)

trapint <- function(xgrid, fgrid)
{
	ng <- length(xgrid)
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
	integ <- sum(xvec * fvec)/2
	return(integ)
}

# Computer proof of \alpha|y
n <- seq(20, 100, 1)
p <- seq(10, 80, 1)
alpha <- seq(-10, 10, 1)
R2 <- seq(-0.9, 0.9, 0.1)
grid_df <- expand.grid(n = n, p = p, alpha = alpha, R2 = R2)
grid_df <- subset(grid_df, n >= p)

sigma2 <- seq(0.1, 1e5, 0.1)
for (i in 1:1000) {
	n <- grid_df[i, "n"]
	p <- grid_df[i, "p"]
	alpha <- grid_df[i, "alpha"]
	R2 <- grid_df[i, "R2"]

	a <- -3/4
	b <- (n - p) / 2 - a - 2
	c <- (n-1) / 2

	log_K <- c * log(n/2) + (b + 1) * log(1 - R2) - lgamma(c) - 0.5 * log(2*pi / n)
	integrand <- exp(-(n*alpha^2)/(2*sigma2) + -(c+3/2) * log(sigma2) -n/(2 * sigma2)) * hyperg_1F1(b + 1, c, (n * R2)/(2*sigma2))

	numeric <- exp(log_K) * trapint(sigma2, integrand)

	analytic <- exp(lgamma(c + 0.5) + (b + 1) * log(1 - R2) - lgamma(c) - 0.5 * log(pi) -n/2 * log(1 + alpha^2)) * hyperg_2F1(b + 1, c + 0.5, c, R2 / (1 + alpha^2))

	cat(numeric, analytic, abs((numeric - analytic) / numeric), "\n")
}

# Computer proof of E[g/(1 + g) | y]
n <- seq(20, 100, 1)
p <- seq(10, 80, 1)
R2 <- seq(-0.9, 0.9, 0.1)
grid_df <- expand.grid(n = n, p = p, R2 = R2)
grid_df <- subset(grid_df, n >= p)
g <- seq(0.1, 1e5, 0.1)

for (i in 1:1000) {
	n <- grid_df[i, "n"]
	p <- grid_df[i, "p"]
	R2 <- grid_df[i, "R2"]
	a <- -3/4
	b <- (n - p) / 2 - a - 2

	log_K <- (b + 1) * log(1-R2) - lbeta(p/2 + a + 1, b + 1)
	integrand <- exp((b + 1) * log(g) - log(1 + g) - n / 2 * log(1 + g * (1 - R2)))
	numeric <- exp(log_K) * trapint(g, integrand)

	analytic <- exp(log_K + lbeta(p / 2 + a + 1, b + 2)) * hyperg_2F1(n/2, b + 2, n / 2 + 1, R2)

	cat(numeric, analytic, abs((numeric - analytic) / numeric), "\n")
}

# Computer proof of E[g^2/(1 + g)^2 | y]
n <- seq(20, 100, 1)
p <- seq(10, 80, 1)
R2 <- seq(-0.9, 0.9, 0.1)
grid_df <- expand.grid(n = n, p = p, R2 = R2)
grid_df <- subset(grid_df, n >= p)
g <- seq(0.1, 1e5, 0.1)

for (i in 1:1000) {
	n <- grid_df[i, "n"]
	p <- grid_df[i, "p"]
	R2 <- grid_df[i, "R2"]
	a <- -3/4
	b <- (n - p) / 2 - a - 2

	log_K <- (b + 1) * log(1-R2) - lbeta(p/2 + a + 1, b + 1)
	integrand <- exp((b + 2) * log(g) -2 * log(1 + g) - n / 2 * log(1 + g * (1 - R2)))
	numeric <- exp(log_K) * trapint(g, integrand)

	analytic <- exp(log_K + lbeta(p / 2 + a + 1, b + 3)) * hyperg_2F1(n/2, b + 3, n / 2 + 2, R2)

	cat(numeric, analytic, abs((numeric - analytic) / numeric), "\n")
}
