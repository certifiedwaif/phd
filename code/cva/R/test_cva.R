library(correlation)
library(Matrix)
library(EMVS)

# Test 1
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


binary_to_model <- function(binary_vec)
{
  acc <- 0
  mul <- 1
  for (i in 1:length(binary_vec)) {
    acc <- acc + mul * binary_vec[i]
    mul <- mul * 2
  }
  return(acc)
}


log_p <- function(n, p, vR2, vp_gamma)
{
  R2 <- R2[2:length(R2)]
  p_gamma <- p_gamma[2:length(p_gamma)]
  a <- 1
  b <- p
  return(-n / 2 * log(1 - R2) - p_gamma / 2 * log(n) + lbeta(a + p_gamma, b + p - p_gamma))
}


posterior_percentage <- rep(NA, 100)
for (i in 1:1) {
  cat("i", i, "\n")
  data <- generate_data(n=50, K=100)
  n <- data$n
  p <- data$p
  vy <- data$vy
  mX <- data$mX
  K <- data$K
  initial_gamma <- data$initial_gamma

  cva_result <- cva(initial_gamma, vy, mX, K, 1.0)
  corr_result <- correlation::all_correlations_mX(vy, mX, "maruyama", 1, TRUE)
  R2 <- corr_result$vR2
  p_gamma <- corr_result$vp_gamma
  cva_models <- apply(cva_result$models, 1, binary_to_model)

  vlog_p <- log_p(n, p, R2, p_gamma)
  vlog_p.til <- vlog_p - max(vlog_p)
  vp <- exp(vlog_p.til)/sum(exp(vlog_p.til))
  plot(vlog_p.til, pch=21, xlab="Model Index", ylab="Log Posterior Model Probability", col="black", bg="black")
  points(cva_models, vlog_p.til[cva_models], pch=21, col="red", bg="red")
  fact <- rep("Non-CVA", length(vlog_p.til))
  fact[cva_models] <- "CVA"
  library(tidyverse)
  dat <- tibble(model=1:length(vlog_p.til), vlog_p=vlog_p.til, cva=fact)
  pdf("cva_low_dimensional.pdf")
  ggplot(dat, aes(x=model, y=vlog_p, color=cva)) + geom_point() + xlab("Model Index") + ylab("Log Posterior Model Probability") + labs(color = "Model type")
  dev.off()
  posterior_percentage[i] <- sum(vp[cva_models])
}
mean(posterior_percentage)

# Compare with EMVS
v0=seq(0.001,0.2,by=0.005)
v1=1.0E3
beta_init=rep(1,p)
a=1
b=1
epsilon=1.0E-5

EMVS_result=EMVS(vy,mX,v0=v0,v1=v1,type="betabinomial",beta_init=beta_init,sigma_init=1,epsilon=epsilon,a=a,b=b)

EMVSplot(EMVS_result,"both",FALSE)
EMVSbest(EMVS_result)
EMVSsummary(EMVS_result)

x <- diff(vlog_p)
n <- length(x)

# Why do we care about the zero crossings?
# Modes
zero_crossings <- c()
for (i in 1:(n-1)) {
  if (x[i] > 0 & x[i + 1] < 0) {
    zero_crossings <- c(zero_crossings, i)
  }
}


plot_traj_probs <- function(cva_result)
{
  traj_prob <- log(cva_result$trajectory_probs)
  plot(1:ncol(traj_prob), traj_prob[1, ], ylim=c(min(traj_prob), max(traj_prob)), type="l",
       xlab="Iteration", ylab="Posterior Probability")
  for (k in 1:nrow(traj_prob)) {
    lines(1:ncol(traj_prob), traj_prob[k, ])
  }
}
plot_traj_probs(cva_result)


posterior_percentages <- function(K, lambda=1.)
{
  n_sims <- 1e4
  cva_results <- list()
  global_mode <- rep(NA, n_sims)
  posterior_percentage <- rep(NA, n_sims)
  for (i in 1:n_sims) {
    data <- generate_data()
    n <- data$n
    p <- data$p
    vy <- data$vy
    mX <- data$mX
    K <- data$K
    initial_gamma <- data$initial_gamma

    cva_results[[i]] <- cva(initial_gamma, vy, mX, K, lambda=lambda)
    corr_result <- correlation::all_correlations_mX(vy, mX)
    R2 <- corr_result$vR2
    p_gamma <- corr_result$vp_gamma
    vp <- exp(log_p(n, p, R2, p_gamma))
    cva_models <- apply(cva_results[[i]]$models, 1, binary_to_model)
    global_mode[i] <- ifelse(any(cva_models == which.max(vp)), 1, 0)
    posterior_percentage[i] <- sum(vp[cva_models])/sum(vp)
  }

  return(list(median_posterior_percentage=median(posterior_percentage),
              mean_global_mode=mean(global_mode)))
}

for (K in c(20, 50, 100, 200)) {
  for (lambda in 0:3) {
    cat(K, lambda, "\n")
    cat(str(posterior_percentages(K, lambda)), "\n")
  }
}

covariate_selection_probs <- apply(cva_result$models, 2, sum)/K

# Distance between covariation inclusion probabilities from CVA and EMVS
# The probabilities seem to be on different scales.
apply(EMVS_result$prob_inclusion, 1, function(x) { sqrt(sum((covariate_selection_probs - x)^2)) })

v0=exp(seq(log(0.05),log(5),,20))
v1=1000

beta_init=rep(0,p)

a=1
b=p
epsilon=10E-3
sigma_init=1

res.emvs <- EMVS(vy, mX, v0=v0,v1=v1,type="betabinomial",beta_init=beta_init,sigma_init=1,epsilon=epsilon,a=a,b=b)

# Test 2
K <- 2
n <- 1e3
p <- 2
sigma2 <- 1.
mX <- matrix(rnorm(n * p), n, p)
vy <- 0 * mX[, 1] + 1 * mX[, 2] + rnorm(n, 0., sigma2)
mX <- scale(mX)
vy <- (vy - mean(vy))
vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
fit <- lm(vy~mX[, 1] + mX[, 2])
#initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
initial_gamma <- matrix(c(1, 0,
                          0, 1), K, p)
cva_result <- cva(initial_gamma, vy, mX, K, lambda=10.)
image(Matrix(cva_result$models))

# Replicate tests from Veronika Rockova's paper
vbeta_0 <- c(1.3, 0, 0, 1.3, 0, 0, 1.3, 0, 0, 1.3, 0, 0)
n <- 100
p <- length(vbeta_0)
mX <- matrix(rnorm(n * p), n, p)
vy <- mX %*% vbeta_0 + rnorm(n)
mX <- scale(mX)
vy <- (vy - mean(vy))
vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
K <- 20
initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
result <- cva(initial_gamma, vy, mX, K, lambda=10.)
v0=seq(0.1,1,by=0.1)
v1=1000
beta_init=rep(1,p)
a=b=1
epsilon=10^{-5}

# Test 3
library(EMVS)
result2 <- EMVS(vy,mX,v0=v0,v1=v1,type="betabinomial",beta_init=beta_init,sigma_init=1,epsilon=epsilon,a=a,b=b)
EMVSplot(result2,"both",FALSE)
EMVSbest(result2)
library(Matrix)
image(Matrix(result$models))

vbeta_0 <- c(1.3, 0, 0, 1.3, 0, 0, 1.3, 0, 0)
n <- 100
p <- length(vbeta_0)
mX <- matrix(rnorm(n * p), n, p)
vy <- mX %*% vbeta_0 + rnorm(n)
mX <- scale(mX)
vy <- (vy - mean(vy))
vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
K <- 20
initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
result <- cva(initial_gamma, vy, mX, K)
v0=seq(0.1,1,by=0.1)
v1=1000
beta_init=rep(1,p)
a=b=1
epsilon=10^{-5}
library(EMVS)
result2 <- EMVS(vy,mX,v0=v0,v1=v1,type="betabinomial",beta_init=beta_init,sigma_init=1,epsilon=epsilon,a=a,b=b)
EMVSplot(result2,"both",FALSE)
EMVSbest(result2)
library(Matrix)
image(Matrix(result$models))

# Test the new interface
corr_result1 <- correlation::all_correlations_mX(vy, mX, "maruyama")
corr_result2 <- correlation::all_correlations_mX(vy, mX, "BIC")
corr_result3 <- correlation::all_correlations_mX(vy, mX, "ZE")
corr_result4 <- correlation::all_correlations_mX(vy, mX, "liang_g1")
corr_result5 <- correlation::all_correlations_mX(vy, mX, "liang_g2")
corr_result6 <- correlation::all_correlations_mX(vy, mX, "liang_g3")
corr_result7 <- correlation::all_correlations_mX(vy, mX, "robust_bayarri1")
corr_result8 <- correlation::all_correlations_mX(vy, mX, "robust_bayarri2")
summary(corr_result1$vlogp)
summary(corr_result2$vlogp)
summary(corr_result3$vlogp)
summary(corr_result4$vlogp)
summary(corr_result5$vlogp)
summary(corr_result6$vlogp)
summary(corr_result7$vlogp)
summary(corr_result8$vlogp)
