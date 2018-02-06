#!/usr/bin/env Rscript
# application.R

# library(sqldf)
# source("accuracy.R")

# # Load the data set.
# long12_analyse <- read.csv(file="data/long12_analyse.csv", header=TRUE)
# long12_analyse <- long12_analyse[order(long12_analyse$Id, long12_analyse$time), ]
# # Drop individuals that we don't observe at all three time points.
# Ids <- sqldf("select Id, count(*) from long12_analyse group by Id having count(*) = 3")
# long12_subset <- subset(long12_analyse, Id %in% Ids$Id)
# # Construct vy, mX and mZ matrices.
# vy <- long12_subset$mintot
# mX <- subset(long12_subset, select=-c(mintot, Id))
# mX <- cbind(1, subset(long12_subset, select=time))
# # Scale time, otherwise everything overlows and fails horribly
# mX$time <- scale(mX$time)
# m <- dim(Ids)[1]
# ni <- 3
# mZ <- matrix(0, ni * m, m * 2)
# i <- 1:(ni * m)
# j <- rep(1:m, each=ni)
# for (idx in 1:length(i)) {
# 	mZ[i[idx], (j[idx] - 1) * 2 + 1] <- 1
# 	mZ[i[idx], (j[idx] - 1) * 2 + 2] <- mX$time[idx]
# }
# # threshold <- 300
# # vy[vy > threshold] <- threshold
# # vy <- vy / 10C
# # Construct mult object
# # mult <- create_mult(vy, mX, mZ[, 3:ncol(mZ)], 1e-5, m=125, blocksize=2)
# # There were two zeroes in a row, so I interpolated.
# vy[194] <- 30
# mult <- create_mult(vy, as.matrix(mX), mZ[, 3:ncol(mZ)], 1e5,
# 										m=125, blocksize=2)
# # Initialise vmu to something not too far off, but without giving the whole game away.
# mult$vmu[1:2] <- c(log(1 + mean(vy)), 1)
# #mult$mLambda <- solve(crossprod(mult$mC))
# # mult$vmu <- c(coef(fit))
# # VB fit.
# # TODO: This overflows! Fix. Look at the likelihood, and see where it's
# # coming from.

# # Laplace works because there's no 0.5 * sigma2 argument to the exponential
# now <- Sys.time()
# fit <- zipvb(mult, method="gva", verbose=FALSE)
# print(Sys.time() - now)

# now <- Sys.time()
# fit2 <- zipvb(mult, method="gva2", verbose=FALSE)
# print(Sys.time() - now)
# # MCMC fit.
# # Compare.

# I might not be able to use that data set. Instead, I could use the roach data set
# from ARM.
# IPM_BASELINE_R2.csv        RoachCounts.csv
# IPM_BASELINE_R2_032006.csv roachdata.csv

# Roaches ------

library(zipvb)
setwd("~/Dropbox/phd/code/")
source("accuracy.R")
library(sqldf)

roachdata <- read.csv(file="../ARM_Data/roaches/roachdata.csv")
RoachCounts <- read.csv(file="../ARM_Data/roaches/RoachCounts.csv")
IPM_BASELINE_R2 <- read.csv(file="../ARM_Data/roaches/IPM_BASELINE_R2.csv")
IPM_BASELINE_R2_032006 <- read.csv(file="../ARM_Data/roaches/IPM_BASELINE_R2.csv")
roaches <- cbind(RoachCounts, IPM_BASELINE_R2_032006)
# Must exclude entries where second observation is NA
# Get first five from each building.
roaches <- sqldf("select * from roaches where hasround2 = 1")
# roaches_new <- data.frame()
# for (building in 1:9) {
# 	roaches_new <- rbind(roaches_new, sqldf(sprintf("select * from roaches where building = %d limit 5", building)))
# }
# roaches <- roaches_new
# Convert from wide to long
roaches_long <- matrix(NA, nrow(roaches) * 2, 9)
colnames(roaches_long) <- c("time", "treatment", "senior", "building", "stories", "hasround2", "roachsum", "trapmiss", "trapdays")
roaches_long <- as.data.frame(roaches_long)
for (row_idx in 1:nrow(roaches)) {
	roaches_long[row_idx * 2 - 1, "time"] <- 1
	roaches_long[row_idx * 2 - 1, "roachsum"] <- roaches[row_idx, "roachsum1"]
	roaches_long[row_idx * 2 - 1, "trapmiss"] <- roaches[row_idx, "trapmiss1"]
	roaches_long[row_idx * 2 - 1, "trapdays"] <- roaches[row_idx, "trapdays1"]
	roaches_long[row_idx * 2 - 1, 2:6] <- roaches[row_idx, 7:11]

	roaches_long[row_idx * 2, "time"] <- 2
	roaches_long[row_idx * 2, "roachsum"] <- roaches[row_idx, "roachsum2"]
	roaches_long[row_idx * 2, "trapmiss"] <- roaches[row_idx, "trapmiss2"]
	roaches_long[row_idx * 2, "trapdays"] <- roaches[row_idx, "trapdays2"]
	roaches_long[row_idx * 2, 2:6] <- roaches[row_idx, 7:11]
}

#fit <- glm(roachsum~time+treatment+senior, data=roaches_long, family=poisson())
#summary(fit)
# Need to fix NAs
model_mat <- model.matrix(roachsum~time+time:treatment+senior+factor(building), data=roaches_long)

# Construct vy, mX, and mZ
#mX <- model.matrix(~time+time:treatment+senior, data=roaches_long)
mX <- model.matrix(~time+time:treatment, data=roaches_long)
mZ <- model.matrix(~factor(building), data=roaches_long)
mZ <- mZ[, 2:ncol(mZ)]

vy <- round(with(roaches_long, roachsum / trapdays))

mult <- create_mult(vy, mX, mZ, 1e5, m=13, blocksize=1, v=2)
kappa(mult$mC)
# Condition number for this matrix was very high. It's okay now.
# Why does this work with GVA NR, but go crazy with GVA2?
# Need more rounds of GVA NR to get mLambda in the ballpark.
# 2 is too few. 5 is enough.

now <- Sys.time()
fit1 <- zipvb(mult, method="laplace", verbose=FALSE, glm_init=FALSE)
cat("Laplace", Sys.time() - now, "\n")

now <- Sys.time()
fit2 <- zipvb(mult, method="gva", verbose=FALSE, glm_init=FALSE)
cat("GVA", Sys.time() - now, "\n")

now <- Sys.time()
fit3 <- zipvb(mult, method="gva2", verbose=FALSE, glm_init=FALSE)
cat("GVA inv. param", Sys.time() - now, "\n")

now <- Sys.time()
fit4 <- zipvb(mult, method="gva_nr", verbose=FALSE, glm_init=FALSE)
cat("GVA fixed point", Sys.time() - now, "\n")

# Check accuracy
# save <- TRUE
save <- FALSE
if (save) {
  mcmc_result <- mcmc(mult, p=3, iterations=1e6, warmup=1e5, mc.cores = 1)
  mcmc_samples <- mcmc_result$mcmc_samples
  fit <- mcmc_result$fit
  print(fit)
  save(mult, mcmc_result, mcmc_samples, fit, file="data/accuracy_application_2017_06_15.RData")
  #save(mult, mcmc_samples, fit, file="data/accuracy_application_2015_11_20.RData")
  # save(mult, mcmc_samples, fit, allKnots, file="/tmp/accuracy_spline_2015_05_19.RData")
} else {
  load(file="data/accuracy_application_2017_06_15.RData")
  # load(file="data/accuracy_application_2015_11_20.RData")
  # load(file="/tmp/accuracy_spline_2015_05_19.RData")
}

var_accuracy <- calculate_accuracies("application", mult, mcmc_samples, fit1, "Laplace", plot_flag=TRUE)
var_accuracy2 <- calculate_accuracies("application", mult, mcmc_samples, fit2, "GVA", plot_flag=TRUE)
var_accuracy3 <- calculate_accuracies("application", mult, mcmc_samples, fit3, "GVA inv par.", plot_flag=TRUE)
var_accuracy4 <- calculate_accuracies("application", mult, mcmc_samples, fit4, "GVA fixed point", plot_flag=TRUE)

# Police stops ----
library(zipvb)
setwd("~/Dropbox/phd/code/")
source("accuracy.R")

frisk <- read.table("frisk_with_noise.dat",skip=6,header=TRUE)
# Crime should be a factor
C <- model.matrix(stops~factor(eth)+factor(crime)+factor(precinct), data=frisk)
vy <- frisk$stops
mX <- C[, 1:6]
mZ <- C[, 7:ncol(C)]
mult <- create_mult(vy, mX, mZ, 1e5, m=75, blocksize=1, v=2)
mult$prior$a_rho <- 3

now <- Sys.time()
fit1 <- zipvb(mult, method="laplace", verbose=FALSE)
cat("Laplace", Sys.time() - now, "\n")

now <- Sys.time()
fit2 <- zipvb(mult, method="gva", verbose=FALSE)
cat("GVA", Sys.time() - now, "\n")

now <- Sys.time()
fit3 <- zipvb(mult, method="gva2", verbose=FALSE)
cat("GVA2", Sys.time() - now, "\n")

now <- Sys.time()
fit4 <- zipvb(mult, method="gva_nr", verbose=FALSE)
cat("GVA NR", Sys.time() - now, "\n")

save <- FALSE
if (save) {
  mcmc_result <- mcmc(mult, p=6, iterations=1e6, warmup=1e5, mc.cores = 1)
  mcmc_samples <- mcmc_result$mcmc_samples
  fit <- mcmc_result$fit
  print(fit)
  save(mult, mcmc_result, mcmc_samples, fit, file="data/accuracy_application2_2017_03_06.RData")
  # save(mult, mcmc_samples, fit, allKnots, file="/tmp/accuracy_spline_2015_05_19.RData")
} else {
  load(file="data/accuracy_application2_2017_03_06.RData")
  # load(file="/tmp/accuracy_spline_2015_05_19.RData")
}

var_accuracy <- calculate_accuracies("application2", mult, mcmc_samples, fit1, "Laplace", plot_flag=TRUE)
var_accuracy2 <- calculate_accuracies("application2", mult, mcmc_samples, fit2, "GVA", plot_flag=TRUE)
var_accuracy3 <- calculate_accuracies("application2", mult, mcmc_samples, fit3, "GVA inv par.", plot_flag=TRUE)
var_accuracy4 <- calculate_accuracies("application2", mult, mcmc_samples, fit4, "GVA fixed point", plot_flag=TRUE)

# > var_accuracy$vbeta_accuracy
# [1] 74.37920 96.78242 99.17660 97.05468
# > mean(var_accuracy$vu_accuracy)
# [1] 78.72966
# > var_accuracy$sigma2_vu_accuracy
# [1] 44.44708
# > var_accuracy$rho_accuracy
# [1] 98.30838

# Owls ----
library(zipvb)
setwd("~/Dropbox/phd/code/")
source("accuracy.R")
load("Owls.Rdata")
# 27 nests
matrix_owls <- model.matrix(SiblingNegotiation~factor(FoodTreatment)+ArrivalTime+factor(Nest), family=poisson, data=Owls)
vy <- Owls$SiblingNegotiation
mZ <- matrix_owls[, c(1, 4:29)]
mX <- matrix_owls[, 2:3]
mult <- create_mult(vy, mX, mZ, 1e5, m=28, blocksize=1, v=2)

now <- Sys.time()
fit1 <- zipvb(mult, method="laplace", verbose=FALSE, glm_init=TRUE)
cat("Laplace", Sys.time() - now, "\n")

now <- Sys.time()
fit2 <- zipvb(mult, method="gva", verbose=FALSE, glm_init=TRUE)
cat("GVA", Sys.time() - now, "\n")

options(safe_exp=TRUE)
options(threshold=2)
now <- Sys.time()
fit3 <- zipvb(mult, method="gva2", verbose=FALSE, glm_init=TRUE)
cat("GVA inv param", Sys.time() - now, "\n")

now <- Sys.time()
fit4 <- zipvb(mult, method="gva_nr", verbose=FALSE, glm_init=TRUE)
cat("GVA fixed point", Sys.time() - now, "\n")

save <- FALSE
if (save) {
  mcmc_result <- mcmc(mult, p=2, iterations=1e5, warmup=1e4, mc.cores = 1)
  mcmc_samples <- mcmc_result$mcmc_samples
  fit <- mcmc_result$fit
  print(fit)
  save(mult, mcmc_result, mcmc_samples, fit, file="data/accuracy_application_owls_2017_06_13.RData")
} else {
  load("data/accuracy_application_owls_2017_06_13.RData")
}
apply(mcmc_samples$vbeta, 2, mean)
apply(mcmc_samples$vu, 2, mean)
fit2$vmu
apply(mcmc_samples$vbeta, 2, sd)
apply(mcmc_samples$vu, 2, sd)
sqrt(diag(fit2$mLambda))

fit1 <- zipvb(mult, method="laplace", verbose=FALSE, glm_init=TRUE)
cat("Laplace", Sys.time() - now, "\n")

now <- Sys.time()
var_accuracy1 <- calculate_accuracies("application_owls", mult, mcmc_samples, fit1, "laplace", plot_flag=TRUE)
cat("Laplace", Sys.time() - now, "\n")

now <- Sys.time()
var_accuracy2 <- calculate_accuracies("application_owls", mult, mcmc_samples, fit2, "GVA", plot_flag=TRUE)
cat("GVA", Sys.time() - now, "\n")

now <- Sys.time()
var_accuracy3 <- calculate_accuracies("application_owls", mult, mcmc_samples, fit3, "GVA inv param", plot_flag=TRUE)
cat("GVA inv param", Sys.time() - now, "\n")

now <- Sys.time()
var_accuracy4 <- calculate_accuracies("application_owls", mult, mcmc_samples, fit4, "GVA fixed point", plot_flag=TRUE)
cat("GVA fixed point", Sys.time() - now, "\n")

# Biochemists ----
library(zipvb)
setwd("~/Dropbox/phd/code/")
source("accuracy.R")
library(pscl)
#bioChemists <- bioChemists[bioChemists$art != 0, ]
bioChemists <- pscl::bioChemists[1:915, ]
#matrix_bioChemists <- model.matrix(art~fem+mar+kid5+phd+ment, data=bioChemists)
matrix_bioChemists <- model.matrix(art~fem+mar+kid5+phd+ment, data=bioChemists)
vy <- bioChemists$art
mZ <- NULL
mX <- matrix_bioChemists[, 1:5]
mult <- create_mult(vy, mX, mZ, 1e5, mPsi=NULL, m=0, blocksize=1, v=2)
options(safe_exp=FALSE)
options(threshold=2)

now <- Sys.time()
fit1 <- zipvb(mult, method="laplace", verbose=TRUE, glm_init=TRUE)
cat("Laplace", Sys.time() - now, "\n")

now <- Sys.time()
fit2 <- zipvb(mult, method="gva", verbose=TRUE, glm_init=TRUE)
cat("GVA", Sys.time() - now, "\n")

now <- Sys.time()
fit3 <- zipvb(mult, method="gva2", verbose=TRUE, glm_init=TRUE)
cat("GVA inv param", Sys.time() - now, "\n")

now <- Sys.time()
fit4 <- zipvb(mult, method="gva_nr", verbose=TRUE, glm_init=TRUE)
cat("GVA fixed point", Sys.time() - now, "\n")

save <- FALSE
if (save) {
  mcmc_result <- mcmc(mult, p=5, iterations=1e6, warmup=1e5, mc.cores = 1,
                      stan_file="multivariate_zip_fixed.stan")
  mcmc_samples <- mcmc_result$mcmc_samples
  fit <- mcmc_result$fit
  save(mult, mcmc_result, mcmc_samples, fit, file="data/accuracy_application_biochemists_2017_06_13.RData")
} else {
  load("data/accuracy_application_biochemists_2017_06_13.RData")
}

apply(mcmc_samples$vbeta, 2, mean)
apply(mcmc_samples$vu, 2, mean)
fit1$vmu
apply(mcmc_samples$vbeta, 2, sd)
apply(mcmc_samples$vu, 2, sd)
sqrt(diag(fit1$mLambda))

var_accuracy1 <- calculate_accuracies("application_biochemists", mult, mcmc_samples, fit1, "laplace", plot_flag=TRUE)
var_accuracy2 <- calculate_accuracies("application_biochemists", mult, mcmc_samples, fit2, "GVA", plot_flag=TRUE)
var_accuracy3 <- calculate_accuracies("application_biochemists", mult, mcmc_samples, fit3, "GVA inv par.", plot_flag=TRUE)
var_accuracy4 <- calculate_accuracies("application_biochemists", mult, mcmc_samples, fit4, "GVA fixed point", plot_flag=TRUE)

# Epilepsy
library(robustbase)
data(epilepsy)

n <- 59
m <- 5

vy <- rep(0, n * m)
mX <- matrix(0, n * m, 2)
mZ_m <- model.matrix(~factor(epilepsy$ID))
mZ <- matrix(0, n * m, n)

for (i in 1:59) {
  for (j in 1:5) {
    idx <- (i - 1) * m + j
    if (j == 1) vy[idx] <- epilepsy[i, "Base"]
    if (j == 2) vy[idx] <- epilepsy[i, "Y1"]
    if (j == 3) vy[idx] <- epilepsy[i, "Y2"]
    if (j == 4) vy[idx] <- epilepsy[i, "Y3"]
    if (j == 5) vy[idx] <- epilepsy[i, "Y4"]
    mX[idx, 1] <- 1
    if (epilepsy[i, "Trt"] == "progabide") mX[idx, 2] <- 1
    mZ[idx, ] <- mZ_m[i, ]
  }
}
mult <- create_mult(vy, mX, mZ, 1e5, m=60, blocksize=1, v=2)
start <- Sys.time()
fit_epilepsy <- zipvb(mult, method="gva2", verbose=FALSE, glm_init=FALSE)
Sys.time() - start