# quadRegMFVB.R

set.seed(1)
nBatch <- 15
beta <- c(1.36, -1.1, 1)
sigma <- 0.04

x <- runif(nBatch)
X <- cbind(rep(1, nBatch), x, x^2)
y <- X %*% beta + sigma * rnorm(nBatch)

plot(x, y, bty="l", type="n", xlim=c(0, 1), ylim=c(1, 1.3))
points(x, y, col="dodgerblue", lwd=2)

XTX <- crossprod(X)
XTy <- crossprod(X, y)
yTy <- sum(y^2)

sigma.beta <- 1e5
A <- 1e5

E.q.recip.sigsq <- 1
E.q.recip.a <- 1
Cov.q.beta <- diag(3)
for (itnum in 1:1000) {
  Cov.q.beta <- solve(E.q.recip.sigsq * XTX + (1/sigma.beta^2) * diag(3))
  E.q.beta <- E.q.recip.sigsq * Cov.q.beta %*% XTy
  E.q.recip.a <- 1/(E.q.recip.sigsq + (1/A^2))
  E.q.recip.sigsq <- (nBatch + 1) / (2 * E.q.recip.a + yTy - 2 * sum(E.q.beta*XTy) + sum(diag(XTX %*% (Cov.q.beta + tcrossprod(E.q.beta)))))
}

ng <- 101
xg <- seq(0, 1, length=ng)
Xg <- cbind(rep(1, ng), xg, xg^2)
meang <- as.vector(Xg %*% E.q.beta)
seg <- sqrt(diag(Xg %*% Cov.q.beta %*% t(Xg)))
lowg <- meang - qnorm(0.975) * seg
uppg <- meang + qnorm(0.975) * seg
lines(xg, meang, col="blue", lwd=2)
lines(xg, lowg, col="blue", lwd=2, lty=2)
lines(xg, uppg, col="blue", lwd=2, lty=2)

nMax <- 500
for (n in (nBatch+1):nMax) {
  xnew <- runif(1)
  Xnew <- c(1, xnew, xnew^2)
  ynew <- sum(Xnew*beta) + sigma*rnorm(1)
  x <- c(x, xnew); y <- c(y, ynew)
  XTX <- XTX + tcrossprod(Xnew)
  XTy <- XTy + Xnew * ynew
  yTy <- yTy + sum(ynew^2)
  Cov.q.beta <- solve(E.q.recip.sigsq*XTX + (1 / sigma.beta^2) * diag(3))
  E.q.beta <- E.q.recip.sigsq * Cov.q.beta %*% XTy
  E.q.recip.a <- 1/(E.q.recip.sigsq + (1/A^2))
  E.q.recip.sigsq <- (nBatch + 1) / (2 * E.q.recip.a + yTy - 2 * sum(E.q.beta*XTy) + sum(diag(XTX %*% (Cov.q.beta + tcrossprod(E.q.beta)))))

  plot(x, y, bty="l", type="n", xlim=c(0, 1), ylim=c(1, 1.3))
  points(x, y, col="dodgerblue", lwd=2)
  text(.5, 1.28, paste("n=", n), col="green3", cex=3)
  ng <- 101
  xg <- seq(0, 1, length=ng)
  Xg <- cbind(rep(1, ng), xg, xg^2)
  meang <- as.vector(Xg %*% E.q.beta)
  seg <- sqrt(diag(Xg %*% Cov.q.beta %*% t(Xg)))
  lowg <- meang - qnorm(0.975) * seg
  uppg <- meang + qnorm(0.975) * seg
  lines(xg, meang, col="blue", lwd=2)
  lines(xg, lowg, col="blue", lwd=2, lty=2)
  lines(xg, uppg, col="blue", lwd=2, lty=2)
  cat("Hit enter to continue. \n")
  ans <- readline()
}