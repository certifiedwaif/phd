library(tidyverse)
library(correlation)
library(Matrix)

# Generate the data
set.seed(1)

generate_data <- function(n=50, K=50)
{
  p <- 12
  sigma2 <- 1.
  mX <- matrix(0, n, p)
  mSigma_block <- matrix(c(1, .9, .9,
                           .9, 1, .9,
                           .9, .9, 1), 3, 3)
  mSigma <- as.matrix(bdiag(mSigma_block, mSigma_block, mSigma_block,
                            mSigma_block))
  chol_mSigma <- chol(mSigma)
  for (i in 1:n) {
    mX[i, ] <- t(chol_mSigma) %*% rnorm(p)
  }
  mX <- scale(mX)
  vy <- 1.3 * mX[, 1] + 1.3 * mX[, 4] + 1.3 * mX[, 7] + 1.3 * mX[, 10] + rnorm(n, 0., sigma2)

  vy <- (vy - mean(vy))
  vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
  initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
  return(list(n=n, p=p, vy=vy, mX=mX, K=K, initial_gamma=initial_gamma))
}

data <- generate_data(n=50, K=100)
n <- data$n
p <- data$p
vy <- data$vy
mX <- data$mX
K <- data$K
initial_gamma <- data$initial_gamma

# Fit models
corr_result_maruyama <- correlation::all_correlations_mX(vy, mX, "maruyama")
corr_result_bic <- correlation::all_correlations_mX(vy, mX, "BIC")
corr_result_ze <- correlation::all_correlations_mX(vy, mX, "ZE")
corr_result_g1 <- correlation::all_correlations_mX(vy, mX, "liang_g1")
corr_result_g2 <- correlation::all_correlations_mX(vy, mX, "liang_g2")
corr_result_g3 <- correlation::all_correlations_mX(vy, mX, "liang_g3")
corr_result_robust_bayarri1 <- correlation::all_correlations_mX(vy, mX, "robust_bayarri1")
corr_result_robust_bayarri2 <- correlation::all_correlations_mX(vy, mX, "robust_bayarri2")

source("~/Downloads/Marginal.Rs")

#vlogp_maruyama <-
vlogp_bic <- marginal.cake(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)
vlogp_ze <- marginal.ZE(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)
vlogp_g1 <- marginal.g.naive(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)
vlogp_g2 <- marginal.g.safe1(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)
vlogp_g3 <- marginal.g.safe2(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)
vlogp_robust_bayarri_1 <- marginal.robust.quad1(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)
vlogp_robust_bayarri_2 <- marginal.robust.naive4(corr_result_maruyama$vR2, n, corr_result_maruyama$vp_gamma)

sum(corr_result_bic$vlogp - vlogp_bic)
sum(corr_result_ze$vlogp - vlogp_ze)
sum(corr_result_g1$vlogp - vlogp_g1)
sum(corr_result_g2$vlogp - vlogp_g2)
sum(corr_result_g3$vlogp[2:4096] - vlogp_g3[2:4096])
sum(corr_result_robust_bayarri1$vlogp - vlogp_robust_bayarri_1)
sum(corr_result_robust_bayarri2$vlogp - vlogp_robust_bayarri_2)

# Timings
library(MASS)

mD <- UScrime
notlog <- c(2,ncol(UScrime))
mD[,-notlog] <- log(mD[,-notlog])

for (j in 1:ncol(mD)) {
  mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
}

varnames <- c(
  "log(AGE)",
  "S",
  "log(ED)",
  "log(Ex0)",
  "log(Ex1)",
  "log(LF)",
  "log(M)",
  "log(N)",
  "log(NW)",
  "log(U1)",
  "log(U2)",
  "log(W)",
  "log(X)",
  "log(prison)",
  "log(time)")

y.t <- mD$y
X.f <- data.matrix(cbind(mD[1:15]))
colnames(X.f) <- varnames 
UScrime_tbl <- correlation::timings(y.t, X.f)
print(UScrime_tbl)

library(Ecdat)

dat = read.csv("Kakadu.csv")
vy <- as.vector(dat$income)
mX <- dat[,c(2:21,23)]  
mX <- model.matrix(~.,data=mX)[,-1]
kakadu_tbl <- timings(vy, mX)
print(kakadu_tbl)

dat <- VietNamI
y.t <- as.vector(dat$lnhhexp)
X.f <- dat[,-2]  
X.f <- model.matrix(~.,data=X.f)[,-1]
VietNamI_tbl <- timings(vy, mX)
print(VietNamI_tbl)
