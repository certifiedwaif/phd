source("PoissonLinearModel.Rs")

vmu = c(0.0, 0.0, 0.0, 0.0)
vy = c(1, 2, 3, 4)
mX = matrix(c(1, 3, 1, 0,
				1, 4, 1, 0,
				1, 5, 0, 1,
				1, 7, 0, 1), 4, 4, byrow=TRUE)
mSigma = matrix(0, 4, 4)
diag(mSigma) = rep(1.0, 4)
rho = .2
mSigma[4, 3] = rho
mSigma[3, 4] = rho
mSigma.inv = solve(mSigma)

fit.Lap(vmu, vy, mX, mSigma.inv)
