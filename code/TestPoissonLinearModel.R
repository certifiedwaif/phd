require(lme4)
source("PoissonLinearModel.R")
require(testthat)

main = function()
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
	mSigma = matrix(0, 2, 2)
	diag(mSigma) = rep(1.0, 2)
	rho = .2
	mSigma[2, 1] = rho
	mSigma[1, 2] = rho
	mSigma.inv = solve(mSigma)
	fit = fit.Lap(vbeta, vu, vy, mX, mZ, mSigmaBeta.inv, mSigma.inv)
	print(str(fit))
	expect_equal(fit$vmu, matrix(c(-0.05263632,
								  0.19422878,
								 -0.16430321,
								  0.10113962), 4, 1), tolerance=1e-6)
	
	test_data = data.frame(y=vy, x=mX[,2], Id=c(1, 1, 2, 2))
	print(test_data)
	lme4_fit = glmer(y~x+(1|Id), data=test_data, family=poisson)
	print(lme4_fit)
}

main()
