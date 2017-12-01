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
timings <- function(vy, mX)
{
  
  
  a1 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "maruyama")
  b1 <- proc.time()[3]
  t1 <- b1 - a1
  cat("maruyama", t1, "\n")

  a2 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "BIC")
  b2 <- proc.time()[3]
  t2 <- b2 - a2
  cat("BIC", t2, "\n")
  
  a3 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "ZE")
  b3 <- proc.time()[3]
  t3 <- b3 - a3
  cat("ZE", t3, "\n")
  
  a4 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "liang_g1")
  b4 <- proc.time()[3]
  t4 <- b4 - a4
  cat("liang_g1", t4, "\n")
  
  a5 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "liang_g2")
  b5 <- proc.time()[3]
  t5 <- b5 - a5
  cat("liang_g2", t5, "\n")
  
  a6 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "liang_g3")
  b6 <- proc.time()[3]
  t6 <- b6 - a6
  cat("liang_g3", t6, "\n")
  
  a7 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "robust_bayarri1")
  b7 <- proc.time()[3]
  t7 <- b7 - a7
  cat("robust_bayarri1", t7, "\n")
  
  a8 <- proc.time()[3]
  correlation::all_correlations_mX(vy, mX, "robust_bayarri2")
  b8 <- proc.time()[3]
  t8 <- b8 - a8
  cat("robust_bayarri2", t8, "\n")
}


timings2 <- function(y, X)
{
  # Package, prior
  packages <- c("BLMA",
                "BLMA",
                "BAS",
                "BAS",
                "BVS",
                "BMS",
                "BLMA",
                "BLMA",
                "BAS",
                "BLMA",
                "BLMA",
                "BLMA",
                "BVS",
                "BLMA")

  priors <- c("BIC",
              "ZE",
              "hyper_g",
              "hyper_g_laplace",
              "hyper_g",
              "g",
              "liang_g1",
              "liang_g2",
              "hyper_g_n",
              "liang_g3",
              "quad",
              "approx",
              "Robust",
              "Robust")

  tbl <- tibble(package=packages, prior=priors)
  tbl$time <- map2_dbl(tbl$package, tbl$prior, function(package, prior.val) {
    start_time <- proc.time()[3]
    if (package == "BAS") {
      library(BAS)
      bas.lm(y~X, prior=prior.val, model=uniform())
    }
    if (package == "BVS") {
      library(BayesVarSelect)
      Bvs(formula="y~.", data=data.frame(y=y, X=X), prior.betas=prior.val, prior.models="Constant",
          time.test=FALSE, n.keep=50000)
    }
    if (package == "BMS") {
      library(BMS)
      bms(cbind(y, X), nmodel=50000, mcmc="enumerate", g="hyper=3", mprior="uniform")
    }
    if (package == "BLMA") {
      library(correlation)
      all_correlations_mX(vy, mX, prior.val)
    }
    end_time <- proc.time()[3]
    end_time - start_time
  })
  tbl
}


library(MASS)
dat <- UScrime
vy <- UScrime$y
mX <- as.matrix(UScrime[,-16])
tbl <- timings2(vy, mX)

library(Ecdat)

dat = read.csv("Kakadu.csv")
vy <- as.vector(dat$income)
mX <- dat[,c(2:21,23)]  
mX <- model.matrix(~.,data=mX)[,-1]
timings(vy, mX)

dat <- VietNamI
y.t <- as.vector(dat$lnhhexp)
X.f <- dat[,-2]  
X.f <- model.matrix(~.,data=X.f)[,-1]
timings(vy, mX)
