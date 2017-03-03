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
library(zipvb)
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

fit <- glm(roachsum~time+treatment+senior, data=roaches_long, family=poisson())
summary(fit)
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
# Condition number for this matrix is very high.
# Why does this work with GVA NR, but go crazy with GVA2?
# Need more rounds of GVA NR to get mLambda in the ballpark.
# 2 is too few. 5 is enough.

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
cat("GVA_NR", Sys.time() - now, "\n")

# Check accuracy
# save <- TRUE
save <- FALSE
if (save) {
  mcmc_result <- mcmc(mult, p=3, iterations=1e5, warmup=1e4, mc.cores = 1)
  mcmc_samples <- mcmc_result$mcmc_samples
  fit <- mcmc_result$fit
  print(fit)
  save(mult, mcmc_samples, fit, file="data/accuracy_application_2015_11_20.RData")
  # save(mult, mcmc_samples, fit, allKnots, file="/tmp/accuracy_spline_2015_05_19.RData")
} else {
  load(file="data/accuracy_application_2015_11_20.RData")
  # load(file="/tmp/accuracy_spline_2015_05_19.RData")
}

var_accuracy <- calculate_accuracies("application", mult, mcmc_samples, fit1, "Laplace", plot_flag=TRUE)
# var_accuracy2 <- calculate_accuracies("application", mult, mcmc_samples, fit2, "GVA", plot_flag=TRUE)
var_accuracy3 <- calculate_accuracies("application", mult, mcmc_samples, fit3, "GVA", plot_flag=TRUE)
# var_accuracy4 <- calculate_accuracies("application", mult, mcmc_samples, fit4, "GVA", plot_flag=TRUE)

# Police stops example
	
library(arm) # for display() function
frisk <- read.table("frisk_with_noise.dat",skip=6,header=TRUE)
model.matrix(, data=frisk)
display(fit.1)
C <- model.matrix(stops~factor(precinct)+factor(eth)+crime, data=frisk)
vy <- frisk$stops
mZ <- C[, 2:75]
mX <- C[, c(1, 76:78)]
mult <- create_mult(vy, mX, mZ, 1e5, m=75, blocksize=1, v=2)
zipvb(mult, method="gva2", verbose=TRUE)

save <- FALSE
if (save) {
  mcmc_result <- mcmc(mult, p=4, iterations=1e5, warmup=1e4, mc.cores = 1)
  mcmc_samples <- mcmc_result$mcmc_samples
  fit <- mcmc_result$fit
  print(fit)
  save(mult, mcmc_samples, fit, file="data/accuracy_application_2015_11_20.RData")
  # save(mult, mcmc_samples, fit, allKnots, file="/tmp/accuracy_spline_2015_05_19.RData")
} else {
  load(file="data/accuracy_application_2015_11_20.RData")
  # load(file="/tmp/accuracy_spline_2015_05_19.RData")
}

var_accuracy <- calculate_accuracies("application2", mult, mcmc_samples, fit, "GVA2", plot_flag=TRUE)
# var_accuracy2 <- calculate_accuracies("application", mult, mcmc_samples, fit2, "GVA", plot_flag=TRUE)
# var_accuracy3 <- calculate_accuracies("application", mult, mcmc_samples, fit3, "GVA", plot_flag=TRUE)
# var_accuracy4 <- calculate_accuracies("application", mult, mcmc_samples, fit4, "GVA", plot_flag=TRUE)

# > var_accuracy$vbeta_accuracy
# [1] 74.37920 96.78242 99.17660 97.05468
# > mean(var_accuracy$vu_accuracy)
# [1] 78.72966
# > var_accuracy$sigma2_vu_accuracy
# [1] 44.44708
# > var_accuracy$rho_accuracy
# [1] 98.30838
