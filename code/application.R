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
mX <- subset(long12_subset, select=-mintot)
n <- dim(Ids)[1]
m <- 3
mZ <- matrix(0, m * n, n * 2)
i <- 1:(m * n)
j <- rep(1:n, each=m)
for (idx in 1:length(i)) {
	mZ[i[idx], (j[idx] - 1) * 2 + 1] <- 1
	mZ[i[idx], (j[idx] - 1) * 2 + 2] <- mX$time[idx]
}
# VB fit.
# MCMC fit.
# Compare.
