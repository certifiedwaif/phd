# application.R

library(sqldf)
source("zero_inflated_model.R")

# Load the data set.
long12_analyse <- read.csv(file="data/long12_analyse.csv", header=TRUE)
long12_analyse <- long12_analyse[order(long12_analyse$Id, long12_analyse$time), ]
# Drop individuals that we don't observe at all three time points.
Ids <- sqldf("select Id, count(*) from long12_analyse group by Id having count(*) = 3")
long12_subset <- subset(long12_analyse, Id %in% Ids$Id)
# Construct vy, mX and mZ matrices.
vy <- long12_subset$mintot
mX <- subset(long12_subset, select=-c(mintot, Id))
mX <- cbind(1, subset(long12_subset, select=time))
# Scale time, otherwise everything overlows and fails horribly
mX$time <- scale(mX$time)
m <- dim(Ids)[1]
ni <- 3
mZ <- matrix(0, ni * m, m * 2)
i <- 1:(ni * m)
j <- rep(1:m, each=ni)
for (idx in 1:length(i)) {
	mZ[i[idx], (j[idx] - 1) * 2 + 1] <- 1
	mZ[i[idx], (j[idx] - 1) * 2 + 2] <- mX$time[idx]
}
# threshold <- 300
# vy[vy > threshold] <- threshold
# vy <- vy / 10C
# Construct mult object
# mult <- create_mult(vy, mX, mZ[, 3:ncol(mZ)], 1e-5, m=125, blocksize=2)
# There were two zeroes in a row, so I interpolated.
vy[194] <- 30
mult <- create_mult(vy, as.matrix(mX), mZ[, 3:ncol(mZ)], 1e5,
										m=125, blocksize=2)
# Initialise vmu to something not too far off, but without giving the whole game away.
mult$vmu[1:2] <- c(log(1 + mean(vy)), 1)
#mult$mLambda <- solve(crossprod(mult$mC))
# mult$vmu <- c(coef(fit))
# VB fit.
# TODO: This overflows! Fix. Look at the likelihood, and see where it's
# coming from.

# Laplace works because there's no 0.5 * sigma2 argument to the exponential
now <- Sys.time()
fit <- zipvb(mult, method="gva", verbose=FALSE)
print(Sys.time() - now)

now <- Sys.time()
fit2 <- zipvb(mult, method="gva2", verbose=FALSE)
print(Sys.time() - now)
# MCMC fit.
# Compare.

# I might not be able to use that data set. Instead, I could use the roach data set
# from ARM.
roaches <- read.csv(file="../ARM_Data/roaches/roachdata.csv")
fit <- glm (y ~ roach1 + treatment + senior, family=poisson, offset=log(exposure2), data=roaches)
summary(fit)
# Construct vy, mX, and mZ
# Maybe there's a better way to do these things
colnames(roaches)[2] <- "roach2"
colnames(roaches)
# Doing things this way will not work for our case. This data set is wide, and you need it to
# be long.
roaches_long <- matrix(NA, nrow(roaches) * 2, ncol(roaches) - 2)
idx <- 1
roaches_mat <- as.matrix(roaches)
for (i in 1:nrow(roaches_mat)) {
	roaches_long[idx, 1] <- 0
	roaches_long[idx + 1, 1] <- roaches_mat[i, "exposure2"]
	roaches_long[idx, 2] <- roaches_mat[i, "roach1"]
	roaches_long[idx, 3:4] <- roaches_mat[i, 4:5]
	roaches_long[idx + 1, 2] <- roaches_mat[i, "roach2"]
	roaches_long[idx + 1, 3:4] <- roaches_mat[i, 4:5]
	idx <- idx + 2
}
colnames(roaches_long) <- c("time", "roaches", "treatment", "senior")
roaches_long_df <- as.data.frame(roaches_long)
