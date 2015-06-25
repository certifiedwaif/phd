require(Rcpp)
require(RcppEigen)

mR = matrix(c(2, 0, 0, 3), 2, 2)
sourceCpp(file="fastdiag.cpp")
print(fastinv(mR, 2, 0, 0, 0))
print("Done!")

mL = matrix(c(2, -.5, 0, 0, 0,
              -.5, 3, 0, 0, 0,
              0, 0, 1, 0, 0,
              0, 0, 0, 1, 0,
              0, 0, 0, 0, .8), 5, 5)
print(eigen(mL))
mR = t(chol(mL))
print(mR)
print(solve(mR))
print(fastinv(mR, 3, 1, 2, 0))
print("Done!")
print(mR%*%solve(mR))
print(t(mR)%*%fastinv(mR, 3, 1, 2, 0))

require(microbenchmark)
microbenchmark(solve(mR))
microbenchmark(fastinv(mR, 3, 1, 2, 0))
