library(tidyverse)
library(blma)
library(Matrix)
library(here)

# Generate the data
set.seed(1)

generate_data <- function(n = 50, K = 50) {
    p <- 12
    sigma2 <- 1
    mX <- matrix(0, n, p)
    mSigma_block <- matrix(c(1, 0.9, 0.9, 0.9, 1, 0.9, 0.9, 0.9, 1), 3, 3)
    mSigma <- as.matrix(bdiag(mSigma_block, mSigma_block, mSigma_block, mSigma_block))
    chol_mSigma <- chol(mSigma)
    for (i in 1:n) {
        mX[i, ] <- t(chol_mSigma) %*% rnorm(p)
    }
    mX <- scale(mX)
    vy <- 1.3 * mX[, 1] + 1.3 * mX[, 4] + 1.3 * mX[, 7] + 1.3 * mX[, 10] + rnorm(n, 
        0, sigma2)
    
    vy <- (vy - mean(vy))
    vy <- sqrt(n) * vy/sqrt(sum(vy^2))
    initial_gamma <- matrix(rbinom(K * p, 1, 0.5), K, p)
    return(list(n = n, p = p, vy = vy, mX = mX, K = K, initial_gamma = initial_gamma))
}

data <- generate_data(n = 50, K = 100)
n <- data$n
p <- data$p
vy <- data$vy
mX <- data$mX
K <- data$K
initial_gamma <- data$initial_gamma

# Fit models
blma_result_maruyama <- blma(vy, mX, "maruyama")
blma_result_bic <- blma(vy, mX, "BIC")
blma_result_ze <- blma(vy, mX, "ZE")
blma_result_g1 <- blma(vy, mX, "liang_g1")
blma_result_g2 <- blma(vy, mX, "liang_g2")
blma_result_g3 <- blma(vy, mX, "liang_g3")
blma_result_robust_bayarri1 <- blma(vy, mX, "robust_bayarri1")
blma_result_robust_bayarri2 <- blma(vy, mX, "robust_bayarri2")

source("~/Downloads/Marginal.Rs")

# vlogp_maruyama <-
vlogp_bic <- marginal.cake(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)
vlogp_ze <- marginal.ZE(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)
vlogp_g1 <- marginal.g.naive(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)
vlogp_g2 <- marginal.g.safe1(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)
vlogp_g3 <- marginal.g.safe2(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)
vlogp_robust_bayarri_1 <- marginal.robust.quad1(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)
vlogp_robust_bayarri_2 <- marginal.robust.naive4(blma_result_maruyama$vR2, n, blma_result_maruyama$vp_gamma)

sum(blma_result_bic$vlogp - vlogp_bic)
sum(blma_result_ze$vlogp - vlogp_ze)
sum(blma_result_g1$vlogp - vlogp_g1)
sum(blma_result_g2$vlogp - vlogp_g2)
sum(blma_result_g3$vlogp[2:4096] - vlogp_g3[2:4096])
sum(blma_result_robust_bayarri1$vlogp - vlogp_robust_bayarri_1)
sum(blma_result_robust_bayarri2$vlogp - vlogp_robust_bayarri_2)

# Timings
library(MASS)

mD <- UScrime
notlog <- c(2, ncol(UScrime))
mD[, -notlog] <- log(mD[, -notlog])

for (j in 1:ncol(mD)) {
    mD[, j] <- (mD[, j] - mean(mD[, j]))/sd(mD[, j])
}

varnames <- c("log(AGE)", "S", "log(ED)", "log(Ex0)", "log(Ex1)", "log(LF)", "log(M)", 
    "log(N)", "log(NW)", "log(U1)", "log(U2)", "log(W)", "log(X)", "log(prison)", 
    "log(time)")

y.t <- mD$y
X.f <- data.matrix(cbind(mD[1:15]))
colnames(X.f) <- varnames
UScrime_tbl <- blma::timings(y.t, X.f)
print(UScrime_tbl)
blma_result <- blma(y.t, X.f, "maruyama")
blma_result <- blma(y.t, X.f, "maruyama", "beta-binomial", c(1, 1))
blma_result <- blma(y.t, X.f, "maruyama", "bernoulli", rep(1/15, 15))
str(blma_result)

y.t <- mD$y
X.f <- data.matrix(cbind(mD[, 1:10]))
colnames(X.f) <- varnames
Z.f <- data.matrix(cbind(mD[, 11:15]))
blma_result <- blma_fixed(y.t, X.f, Z.f, "maruyama")

library(Ecdat)

dat = read.csv(here("../Kakadu.csv"))
vy <- as.vector(dat$income)
mX <- dat[, c(2:21, 23)]
mX <- model.matrix(~., data = mX)[, -1]
kakadu_tbl <- timings(vy, mX)
print(kakadu_tbl)

dat <- VietNamI
y.t <- as.vector(dat$lnhhexp)
X.f <- dat[, -2]
X.f <- model.matrix(~., data = X.f)[, -1]
VietNamI_tbl <- timings(y.t, X.f)
print(VietNamI_tbl)

# 
sum.na <- function(x) {
    sum(is.na(x))
}
apply(X, 1, sumNA)
apply(X, 2, sumNA)

inds = which(apply(X, 2, sum.na) == 0)
mX = X[, inds]
i <- 1
vy = Y[, i]

full_fit <- lm(vy ~ mX)
summ <- summary(full_fit)

true_model <- summ$coefficients[summ$coefficients[, 4] < 0.05, ]
p <- ncol(mX)
K <- 20
initial_gamma <- matrix(rbinom(K * p, 1, 0.5), K, p)
cva_result <- cva(initial_gamma, vy, mX, K)
