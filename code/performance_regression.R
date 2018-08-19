rm(list = ls())


doBVS <- FALSE
doMFVB <- TRUE
doLASSO <- TRUE
doMCP <- TRUE
doSCAD <- TRUE
doVARBVS <- TRUE
doEMVS <- TRUE
doBMS <- TRUE

# rm(list = ls())

dataset <- "communities.Rdata"
print(dataset)



load("comData.Rdata")

sum.na <- function(x) {
    sum(is.na(x))
}
inds <- which(apply(X, 2, sum.na) == 0)
X2 <- X[, inds]
X3 <- X2[, !colnames(X2) %in% c("ownHousQrange", "rentQrange")]

y <- Y[, 18]
inds <- which(is.na(y))

vy <- y[-inds]
mX <- X3[-inds, ]
mX.til <- cbind(1, mX)

n <- length(vy)
p <- ncol(mX)


mult <- sqrt(n/(n - 1))
mX <- mX
for (j in 1:p) {
    mX[, j] = mult * (mX[, j] - mean(mX[, j]))/sd(mX[, j])
}
vy <- mult * (vy - mean(vy))/sd(vy)

################################################################################ 


# ans <- readline()

res.lm <- lm(vy ~ mX)
inds <- which(summary(res.lm)$coef[, 4] < 0.1) - 1

library(correlation)

a0 <- proc.time()[3]
res.fast <- all_correlations_mX(vy, mX[, inds], intercept_col = 0, bIntercept = FALSE, 
    bCentre = FALSE, cores = 4)
b0 <- proc.time()[3]
t0 <- b0 - a0
cat(t0, "\n")
cat((2^length(inds))/t0, "\n")
