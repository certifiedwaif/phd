library(correlation)
set.seed(as.numeric(Sys.time()))
n <- 100
p <- 20
sigma2 <- 1.
mX <- matrix(rnorm(n * p), n, p)
vy <- 5 * mX[, 1] - 10 * mX[, 3] + rnorm(n, 0., sigma2)
mX <- scale(mX)
vy <- (vy - mean(vy))
vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
K <- 20
result <- cva(vy, mX, K)
library(Matrix)
image(Matrix(result$bitstrings))
