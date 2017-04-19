library(correlation)
set.seed(as.numeric(Sys.time()))
n <- 50
p <- 12
sigma2 <- 1.
mX <- matrix(rnorm(n * p), n, p)
vy <- 5 * mX[, 1] - 10 * mX[, 3] + 0 * mX[, 12] + rnorm(n, 0., sigma2)
mX <- scale(mX)
vy <- (vy - mean(vy))
vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
K <- 20
initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
now <- Sys.time()
cva_result <- cva(initial_gamma, vy, mX, K)
cat(Sys.time() - now, "\n")
library(Matrix)
image(Matrix(cva_result$models))
now <- Sys.time()
corr_result <- correlation::all_correlations_mX(vy, mX)
cat(Sys.time() - now, "\n")
summary(corr_result$vR2)
summary(lm(vy~mX[, 1] + mX[, 3]-1))

R2 <- corr_result$vR2
p_gamma <- corr_result$vp_gamma
R2 <- R2[2:length(R2)]
p_gamma <- p_gamma[2:length(p_gamma)]
a <- 1
b <- p
log_p <- -n / 2 * log(1 - R2) - p_gamma / 2 * log(n) + lbeta(a + p_gamma, b + p - p_gamma)
plot(exp(log_p), pch=24)

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
cva_models <- apply(cva_result$models, 1, binary_to_model)
points(cva_models, exp(log_p)[cva_models], pch=24, col="red")

library(correlation)
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
library(Matrix)
image(Matrix(cva_result$models))

# Replicate tests from Veronika Rockova's paper
library(correlation)
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

traj_prob <- log(cva_result$trajectory_probs)
plot(1:ncol(traj_prob), traj_prob[1, ], ylim=c(min(traj_prob), max(traj_prob)), type="l")
for (k in 1:nrow(traj_prob)) {
  lines(1:ncol(traj_prob), traj_prob[k, ])
}
