# tquadRegMFVB.R

set.seed(2)
nBatch <- 200
beta <- c(1.36, -1.1, 1)
sigma <- .04
nu <- 0.5
x <- runif(nBatch)
X <- cbind(rep(1, nBatch), x, x^2)
y <- X %*% beta + sigma * rt(nBatch, nu)
plot(x, y, bty="l", type="n")
points(x, y, col="dodgerblue", lwd=2)

sigma.beta <- 1e5
A <- 1e5
E.q.recip.sigsq <- 1
E.q.beta <- rep(0, 3)
Cov.q.beta <- diag(3)
for (itnum in 1:1000) {
  resid <- as.vector(y - X %*% E.q.beta)
  E.q.recip.b <- (nu + 1) / (nu + E.q.recip.sigsq * (resid^2 + diag(X %*% Cov.q.beta %*% t(X))))
  Cov.q.beta <- solve(E.q.recip.sigsq * crossprod(X, E.q.recip.b * X) + (1/sigma.beta^2) * diag(3))
  E.q.beta <- E.q.recip.sigsq * Cov.q.beta %*% crossprod(X, E.q.recip.b * y)
  E.q.recip.a <- 1/(E.q.recip.sigsq + (1/A^2))
  E.q.recip.sigsq <- (nBatch + 1) / (2 * E.q.recip.a + sum(E.q.recip.b * resid^2) + sum(diag(Cov.q.beta %*% crossprod(X, E.q.recip.b * X))))
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

# Resistant to outliers because of heavy tailed distribution

plot(0, 0, type="l", xlim=c(0, 1), ylim=range(c(lowg, uppg)), bty="l",
     xlab="x", ylab="y", main="zoomed view")
points(x, y, col="dodgerblue")
lines(xg, meang, col="blue", lwd=2)
lines(xg, lowg, col="blue", lwd=2, lty=2)
lines(xg, uppg, col="blue", lwd=2, lty=2)
lines(xg, as.vector(cbind(rep(1, ng), xg, xg^2) %*% beta), col="red")
legend("bottomright", legend=c("estimate", "truth"), lwd = rep(2, 1), col=c("blue", "red"))
