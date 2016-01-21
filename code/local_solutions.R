#!/usr/bin/env Rscript
# local_solutions.R
library(optparse)
library(zipvb)
source("zero_inflated_model.R")

local_solutions <- function(approximation)
{
	# Construct suitable mX, mZ and vy
	n <- 100
	vx <- rnorm(n)
	mX <- cbind(1, vx)

	# Construct mZ with three blocks
	mZ <- cbind(
		rep(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 1, 0, 0, 0, 0, 0, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 1, 0, 0, 0, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 1, 0, 0, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 0, 1, 0, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 0, 0, 1, 0, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 0, 0, 0, 1, 0, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 0, 0, 0, 0, 1, 0), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1), each=rep(10, 10)),
		rep(c(0, 0, 0, 0, 0, 0, 0, 0, 0, 1), each=rep(10, 10))
	)
	mZ <- mZ[1:n, 2:10]
	vy <- rep(NA, n)
	true_intercept <- 2
	true_beta <- 1.0
	for (i in 1:n) {
		lambda_param <- true_intercept * mX[i, 1] + true_beta * mX[i, 2]
		lambda_param <- lambda_param + sum(true_intercept * mZ[i, 1:9])
		vy[i] <- rpois(1, lambda = lambda_param)
	}

	# Construct mult object
	sigma2.beta <- 1e5
	mult <- create_mult(vy, mX, mZ, sigma2.beta, m=10, blocksize=1, v=2)

	counter <- 0
	error_count <- 0
	vmus <- matrix(NA, (10 * 1e2)^2, 2 + 9)

	# Start fits from a range of points on a 2 dimensional grid
	GRID_SIZE <- 1e-1
	error_locations <- c(-99, -99)
	for (init_intercept in seq(- 4.5, 5, GRID_SIZE)) {
		for (init_slope in seq(- 4.5, 5, GRID_SIZE)) {
			cat("init_intercept", init_intercept, "")
			cat("init_slope", init_slope, "")
			cat("counter", counter, "")
			cat("error_count", error_count, "")
			mult$vmu[1] <- init_intercept
			mult$vmu[2] <- init_slope
			fit <- tryCatch(zipvb(mult, method=approximation, verbose=FALSE, glm_init=FALSE),
											error = function(e) {
												print(e)
												return(NULL)
											})
			if (is.null(fit)) {
				print("Caught error")
				error_count <- error_count + 1
				error_locations <- rbind(error_locations,
																 c(init_intercept, init_slope))
			} else {
				# print(str(fit))
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
	cat("Error_count", error_count, "\n")
	return(list(vmus=vmus[1:(counter - 1), ],
							error_locations=error_locations[2:error_count, ]))
}

display_local_solutions <- function(error_locations)
{
	min_x <- min(x)
	max_x <- max(x)
	min_y <- min(y)
	max_y <- max(y)

	x <- seq(min_x, max_x, by=0.1)
	y <- seq(min_y, max_y, by=0.1)

	image_mat <- matrix(x, y, matrix(0, length(x), length(y)))
	for (i in 1:nrow(error_locations)) {
		err_x <- error_locations[i, 2]
		err_y <- error_locations[i, 3]
		image_mat[err_x, err_y] <- 1
	}

	image(x, y, image_mat)
}

main <- function()
{
  option_list <- list(make_option(c("-a", "--approximation"), default="gva"))
  opt <- parse_args(OptionParser(option_list=option_list))
  approximation <- opt$approximation
  result <- local_solutions(approximation)
  vmus <- result$vmus
  error_locations <- result$error_locations
  write.csv(vmus, file=sprintf("results/local_solutions_%s_vmus.csv", approximation))
  write.csv(error_locations, file=sprintf("results/local_solutions_%s_error_locations.csv", approximation))
}
main()
