
greycode <- function(p) {
    A <- matrix(c(0, 1), 2, 1)
    if (p != 1) {
        for (i in 2:p) {
            P <- nrow(A)
            inds <- rev(1:P)
            A1 <- matrix(c(rep(0, P), rep(1, P)), 2 * P, 1)
            A2 <- rbind(A, matrix(A[inds, ], nrow(A), ncol(A)))
            A <- cbind(A1, A2)
        }
    }
    return(A)
}

ALL.bohning <- function(vy, mX, models, tau = 1e-05, MAXABSERR = 0.001) {
    MAXITER <- 1000
    
    n <- length(vy)
    p <- ncol(mX)
    
    mBeta <- matrix(0, nrow(models), p)
    
    for (j in 1:nrow(models)) {
        vgamma <- c(1, models[j, ])
        inds1 <- which(vgamma == 1)
        q <- length(inds1)
        mX1 <- matrix(mX[, inds1], n, q)
        vmu <- matrix(0, q, 1)
        mI <- diag(1, q)
        mSigma.inv <- 0.25 * t(mX1) %*% mX1 + tau * mI
        mA <- solve(mSigma.inv, t(mX1), tol = 1e-99)
        for (ITER in 1:MAXITER) {
            vmu.old <- vmu
            vxi <- mX1 %*% vmu
            vb <- 1/(1 + exp(-vxi))
            vc <- vy + 0.25 * vxi - vb
            vmu <- mA %*% vc
            err <- max(abs(vmu.old - vmu))
            # cat(ITER,err,'\n')
            if (err < MAXABSERR) {
                break
            }
        }
        mBeta[j, inds1] <- vmu
    }
    return(list(mBeta = mBeta))
}

library(compiler)

ALL.bohning2 <- cmpfun(ALL.bohning)


