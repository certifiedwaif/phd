require(lme4)
source("PoissonLinearModel.R")
source("PoissonLinearModel_old.R")
require(testthat)

test_laplace_approximation <- function()
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
}

test_gaussian_approximation <- function()
{
  vbeta = c(0.0, 0.0)
  vu = c(0.0, 0.0) 
  vmu = c(vbeta, vu)
  vy = c(1, 2, 3, 4)
  mat = matrix(c(1, 3, 1, 0,
                 1, 4, 1, 0,
                 1, 5, 0, 1,
                 1, 7, 0, 1), 4, 4, byrow=TRUE)
  mX = mat[,1:2]
  mZ = mat[,3:4]
  
  mSigma = matrix(0, 4, 4)
  diag(mSigma) = rep(1.0, 4)
  rho = .2
  mSigma[4, 3] = rho
  mSigma[3, 4] = rho
  mSigma.inv = solve(mSigma)
  fit_lap = fit.Lap.old(vmu, vy, mat, mSigma.inv)
  fit_gaussian = fit.GVA(fit_lap$vmu, fit_lap$mLambda, vy, mat, mSigma.inv, "Nelder-Mead")  
  expect_equal(fit_gaussian$vmu, matrix(c(0.122079962,
                                          0.090590365,
                                          -0.008619054,
                                          0.420113211), 4, 1), tolerance=1e-6)
}

main <- function()
{
	test_laplace_approximation()
	test_gaussian_approximation()
}

main()
