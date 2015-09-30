#!/usr/bin/env Rscript
# local_solutions.R
library(optparse)
source("zero_inflated_model.R")

local_solutions <- function(approximation)
{
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

	counter <- 0
	errors <- 0
	vmus <- matrix(NA, 1e5, 4)

	# Start fits from a range of points on a 2 dimensional grid
	for (init_intercept in seq(-5, 5, 1e-1)) {
		for (init_slope in seq(-5, 5, 1e-1)) {
			cat("init_intercept", init_intercept, "\n")
			cat("init_slope", init_slope, "\n")
			cat("counter", counter, "\n")
			cat("errors", errors, "\n")
			mult$vmu[1] <- init_intercept
			mult$vmu[2] <- init_slope
			fit <- tryCatch(zipvb(mult, method=approximation, verbose=FALSE, glm_init=FALSE),
											error = function(e) NULL)
			if (is.null(fit)) {
				print("Caught error")
				errors <- errors + 1
			} else {
				print(str(fit))
				vmus[counter, ] <- fit$vmu
			}
			counter <- counter + 1
		}
	}

	# Bad point for GVA
	# > init_intercept
	# [1] -3
	# > init_slope
	# [1] -0.39
	# > 

	# mult$vmu[1] <- -3
	# mult$vmu[2] <- -0.17
	# fit <- zipvb(mult, method="gva", verbose=TRUE, glm_init=FALSE)
	
	cat("Counter", counter, "\n")
	cat("Errors", errors, "\n")
	print(vmus[1:(counter - 1), ])
}

main <- function()
{
  option_list <- list(make_option(c("-a", "--approximation"), default="gva"))
  opt <- parse_args(OptionParser(option_list=option_list))
  approximation <- opt$approximation
  local_solutions(approximation)
}
main()
