# local_solutions.R
source("zero_inflated_model.R")

# Construct suitable mX, mZ and vy
vx <- rnorm(200)
mX <- cbind(1, vx)

# Construct mZ with three blocks
mZ <- cbind(rep(c(1, 0, 0), each=c(75, 75, 50)),
						rep(c(0, 1, 0), each=c(75, 75, 50)),
						rep(c(0, 0, 1), each=c(75, 75, 50)))
mZ <- mZ[1:200, 2:3]
vy <- rep(NA, 200)
true_intercept <- 3
true_beta <- 1.0
for (i in 1:200) {
	vy[i] <- rpois(1, lambda = true_intercept * mX[i, 1] + true_beta * mX[i, 2] + true_intercept * mZ[i, 1] + true_intercept * mZ[i, 2])
}

# Construct mult object
sigma2.beta <- 1e5
mult <- create_mult(vy, mX, mZ, sigma2.beta, m=3, blocksize=1, v=2)

# Start fits from a range of points on a 2 dimensional grid
for (init_intercept in seq(-3, 3, 1e-2)) {
	for (init_slope in seq(-3, 3, 1e-2)) {
		mult$vmu[1] <- init_intercept
		mult$vmu[2] <- init_slope
		fit <- zipvb(mult, method="laplace", glm_init=FALSE)
		print(fit$vmu)
	}
}
