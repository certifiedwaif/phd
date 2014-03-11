require(lme4)
source("PoissonLinearModel.R")
source("PoissonLinearModel_old.R")
require(testthat)

test_new <- function()
{
	vbeta = c(0.0, 0.0)
	vu = c(0.0, 0.0) 
	vy = c(1, 2, 3, 4)
	mat = matrix(c(1, 3, 1, 0,
					1, 4, 1, 0,
					1, 5, 0, 1,
					1, 7, 0, 1), 4, 4, byrow=TRUE)
	mX = mat[,1:2]
	mZ = mat[,3:4]
	
	mSigmaBeta.inv = matrix(0, 2, 2)
  diag(mSigmaBeta.inv) = rep(1.0, 2)
	mSigma = matrix(0, 2, 2)
	diag(mSigma) = rep(1.0, 2)
	rho = .2
	mSigma[2, 1] = rho
	mSigma[1, 2] = rho
	mSigma.inv = solve(mSigma)
	fit = fit.Lap(vbeta, vu, vy, mX, mZ, mSigmaBeta.inv, mSigma.inv, debug=TRUE)
	print(str(fit))
  # vg should equal (6, 35, 1, 5)
	expect_equal(fit$vmu, matrix(c(-0.05263632,
								  0.19422878,
								 -0.16430321,
								  0.10113962), 4, 1), tolerance=1e-6)
	
	test_data = data.frame(y=vy, x=mX[,2], Id=c(1, 1, 2, 2))
	print(test_data)
	lme4_fit = glmer(y~x+(1|Id), data=test_data, family=poisson)
	print(lme4_fit)
}

test_old <- function()
{
	vu = c(0.0, 0.0, 0.0, 0.0) 
	vy = c(1, 2, 3, 4)
	mZ = matrix(c(1, 3, 1, 0,
					1, 4, 1, 0,
					1, 5, 0, 1,
					1, 7, 0, 1), 4, 4, byrow=TRUE)
	
	mSigma = matrix(0, 4, 4)
	diag(mSigma) = rep(1.0, 4)
	rho = .2
	mSigma[4, 3] = rho
	mSigma[3, 4] = rho
	mSigma.inv = solve(mSigma)
	fit = fit.Lap.old(vu, vy, mZ, mSigma.inv)
	print(str(fit))
	expect_equal(fit$vmu, matrix(c(-0.05263632,
								  0.19422878,
								 -0.16430321,
								  0.10113962), 4, 1), tolerance=1e-6)
}

main <- function()
{
	test_old()
	test_new()
}

main()
