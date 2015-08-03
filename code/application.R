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
m <- dim(Ids)[1]
ni <- 3
mZ <- matrix(0, ni * m, m * 2)
i <- 1:(ni * m)
j <- rep(1:m, each=ni)
for (idx in 1:length(i)) {
	mZ[i[idx], (j[idx] - 1) * 2 + 1] <- 1
	mZ[i[idx], (j[idx] - 1) * 2 + 2] <- mX$time[idx]
}
mX2 <- model.matrix(~., data=mX)
# Construct mult object
mult <- create_mult(vy, mX2, mZ, 1e-5, m=125, blocksize=2)
# VB fit.
fit <- zipvb(mult)
# MCMC fit.
# Compare.
