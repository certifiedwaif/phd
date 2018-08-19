trapint <- function(xgrid, fgrid) {
    ng <- length(xgrid)
    xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
    fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
    integ <- sum(xvec * fvec)/2
    return(integ)
}

vg <- 1:10
a <- 10
b <- 10
p <- 10
tau <- 0.5

mX_df <- read.csv("mX.csv", header = FALSE)
mX <- as.matrix(mX_df)
vy_df <- read.csv("vy.csv", header = FALSE)
vy <- as.vector(vy_df)
lm_fit <- lsfit(mX, vy, intercept = FALSE)
summ_fit <- summary(lm_fit)
vmu <- coef(lm_fit)
mXTX <- t(mX) %*% mX
mSigma <- diag(rep(1, ncol(mX)))

log.q.g <- (b - 0.5 * p) * log(vg) - (a + b + 2) * log(1 + vg)
log.q.g <- log.q.g + -0.5 * tau * (t(vmu) %*% mXTX %*% vmu + sum(diag(mXTX %*% mSigma)))/vg
log.q.g.til <- log.q.g - max(log.q.g)
tau.g <- trapint(vg, exp(log.q.g.til)/vg)/trapint(vg, exp(log.q.g.til))
