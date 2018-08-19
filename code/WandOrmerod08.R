########## R script: WandOrmerod08 ##########

# Code corresonding to the appendix of the paper: Wand, M.P. & Ormerod, J.T.
# (2008).  `On semiparametric regression with O'Sullivan penalised splines'.

# Last changed: 29 APR 2010

# Direct scatterplot smoothing with user choice of smoothing parameter
# --------------------------------------------------------------------

# Obtain scatterplot data corresponding to environmental data from the R package
# `lattice'. Set up plotting grid, knots and smoothing parameter:

library(lattice)
attach(environmental)
x <- radiation
y <- ozone^(1/3)
a <- 0
b <- 350
xg <- seq(a, b, length = 101)
numIntKnots <- 20
lambda <- 1000

# Set up the design matrix and related quantities:

library(splines)
intKnots <- quantile(unique(x), seq(0, 1, length = (numIntKnots + 2))[-c(1, (numIntKnots + 
    2))])
names(intKnots) <- NULL
B <- bs(x, knots = intKnots, degree = 3, Boundary.knots = c(a, b), intercept = TRUE)
BTB <- crossprod(B, B)
BTy <- crossprod(B, y)

# Create the Omega matrix:

formOmega <- function(a, b, intKnots) {
    allKnots <- c(rep(a, 4), intKnots, rep(b, 4))
    K <- length(intKnots)
    L <- 3 * (K + 8)
    xtilde <- (rep(allKnots, each = 3)[-c(1, (L - 1), L)] + rep(allKnots, each = 3)[-c(1, 
        2, L)])/2
    wts <- rep(diff(allKnots), each = 3) * rep(c(1, 4, 1)/6, K + 7)
    Bdd <- spline.des(allKnots, xtilde, derivs = rep(2, length(xtilde)), outer.ok = TRUE)$design
    Omega <- t(Bdd * wts) %*% Bdd
    return(Omega)
}

Omega <- formOmega(a, b, intKnots)

# Obtain the coefficients:

nuHat <- solve(BTB + lambda * Omega, BTy)

# For large K the following alternative Cholesky-based approach can be
# considerably faster (O(K), because B'B+ lambda*Omega is banded diagonal):

cholFac <- chol(BTB + lambda * Omega)
nuHat <- backsolve(cholFac, forwardsolve(t(cholFac), BTy))

# Display the fit:

Bg <- bs(xg, knots = intKnots, degree = 3, Boundary.knots = c(a, b), intercept = TRUE)
fhatg <- Bg %*% nuHat
par(mfrow = c(1, 2))
plot(x, y, xlim = range(xg), bty = "l", type = "n", xlab = "radiation", ylab = "cuberoot of ozone", 
    main = "(a) direct fit; user 
     choice of smooth. par.")
lines(xg, fhatg, lwd = 2)
points(x, y, lwd = 2)

# Mixed model scatterplot smoothing with REML choice of smoothing parameter
# -------------------------------------------------------------------------

# Obtain the spectral decomposition of $\\bOmega$:

eigOmega <- eigen(Omega)

# Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$:

indsZ <- 1:(numIntKnots + 2)
UZ <- eigOmega$vectors[, indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))

# Perform stability check:

indsX <- (numIntKnots + 3):(numIntKnots + 4)
UX <- eigOmega$vectors[, indsX]
L <- cbind(UX, LZ)
stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))
if (sum(stabCheck^2) > 1.0001 * (numIntKnots + 2)) print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")

# Form the X and Z matrices:

X <- cbind(rep(1, length(x)), x)
Z <- B %*% LZ

# Fit using lme() with REML choice of smoothing parameter:

library(nlme)
group <- rep(1, length(x))
gpData <- groupedData(y ~ x | group, data = data.frame(x, y))
fit <- lme(y ~ -1 + X, random = pdIdent(~-1 + Z), data = gpData)

# Extract coefficients and plot scatterplot smooth over a grid:

betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)
Zg <- Bg %*% LZ
fhatgREML <- betaHat[1] + betaHat[2] * xg + Zg %*% uHat
plot(x, y, xlim = range(xg), bty = "l", type = "n", xlab = "radiation", ylab = "cuberoot of ozone", 
    main = "(b) mixed model fit; 
      REML choice of smooth. par.")
lines(xg, fhatgREML, lwd = 2)
points(x, y, lwd = 2)

cat("Hit Enter to continue\n")
ans <- readline()

# Fitting an additive mixed model -------------------------------

# The spinal bone mineral density data of Bachrach et al (1999) are not publicly
# available. Therefore we will illustrate fitting of additive mixed models using
# simulated data.  For simplicity we will use two ethnicity categories rather
# than four.

# Generated data:

set.seed(394600)
m <- 230
nVals <- sample(1:4, m, replace = TRUE)
betaVal <- 0.1
sigU <- 0.25
sigEps <- 0.05
f <- function(x) return(1 + pnorm((2 * x - 36)/5)/2)
U <- rnorm(m, 0, sigU)
age <- NULL
ethnicity <- NULL
Uvals <- NULL
idNum <- NULL
for (i in 1:m) {
    idNum <- c(idNum, rep(i, nVals[i]))
    stt <- runif(1, 8, 28 - (nVals[i] - 1))
    age <- c(age, seq(stt, by = 1, length = nVals[i]))
    xCurr <- sample(c(0, 1), 1)
    ethnicity <- c(ethnicity, rep(xCurr, nVals[i]))
    Uvals <- c(Uvals, rep(U[i], nVals[i]))
}

epsVals <- rnorm(sum(nVals), 0, sigEps)
SBMD <- f(age) + betaVal * ethnicity + Uvals + epsVals

# Set up basic variables for the spline component:

a <- 8
b <- 28
numIntKnots <- 15
intKnots <- quantile(unique(age), seq(0, 1, length = (numIntKnots + 2))[-c(1, (numIntKnots + 
    2))])

# Obtain the spline component of the Z matrix:

B <- bs(age, knots = intKnots, degree = 3, Boundary.knots = c(a, b), intercept = TRUE)
Omega <- formOmega(a, b, intKnots)
eigOmega <- eigen(Omega)
indsZ <- 1:(numIntKnots + 2)
UZ <- eigOmega$vectors[, indsZ]
LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
ZSpline <- B %*% LZ

# Obtain the X matrix.

X <- cbind(rep(1, length(SBMD)), age, ethnicity)

# Set up variables required for fitting via lme().  Note that the random
# intercept is taken care of via the tree identification numbers variable
# `idNum', and that explicit formation of the random effect contribution to the Z
# matrix is not required.

groupVec <- factor(rep(1, length(SBMD)))
ZBlock <- list(list(groupVec = pdIdent(~ZSpline - 1)), list(idNum = pdIdent(~1)))
ZBlock <- unlist(ZBlock, recursive = FALSE)
dataFr <- groupedData(SBMD ~ ethnicity | groupVec, data = data.frame(SBMD, X, ZSpline, 
    idNum))
fit <- lme(SBMD ~ -1 + X, data = dataFr, random = ZBlock)
betaHat <- fit$coef$fixed
uHat <- unlist(fit$coef$random)
uSplineHat <- uHat[1:ncol(ZSpline)]

# Plot the data and fitted curve estimates together:

ng <- 101
ageg <- seq(a, b, length = ng)
Bg <- bs(ageg, knots = intKnots, degree = 3, Boundary.knots = c(a, b), intercept = TRUE)
ZgSpline <- Bg %*% LZ

plotMatrix0 <- cbind(rep(1, ng), ageg, rep(0, ng), ZgSpline)
fhatgREML <- plotMatrix0 %*% c(betaHat, uSplineHat)

xLabs <- paste("ethnicity =", as.character(ethnicity))
pobj <- xyplot(SBMD ~ age | xLabs, groups = idNum, xlab = "age (years)", ylab = "spinal bone mineral density", 
    subscripts = TRUE, panel = function(x, y, subscripts, groups) {
        panel.grid()
        panel.superpose(x, y, subscripts, groups, type = "b", col = "grey60", pch = 16)
        panelInd <- any(ethnicity[subscripts] == 1)
        panel.xyplot(ageg, fhatgREML + panelInd * betaHat[3], lwd = 3, type = "l", 
            col = "black")
    })
print(pobj)

# Print approximate 95% confidence intervals for key parameters:

print(intervals(fit))

########## End of WandOrmerod08 ##########
