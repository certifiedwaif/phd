# splines.R
library(lattice)
library(splines)

# Code corresonding to the appendix of the paper:
# Wand,  M.P. & Ormerod,  J.T. (2008).
# `On semiparametric regression with O'Sullivan
#  penalised splines'.

# Last changed: 29 APR 2010

# Direct scatterplot smoothing with user choice of smoothing parameter
# --------------------------------------------------------------------

# Obtain scatterplot data corresponding to environmental
# data from the R package `lattice'. Set up plotting 
# grid,  knots and smoothing parameter:

# Outline:
# Need to generate spline data

fit_spline <- function(x,  y)
{
  # x could be anything you like
  # y could be anything you like
  a <- min(x)
  b <- max(x)
  #x <- radiation
  #y <- ozone^(1/3)
  #a <- 0 
  #b <- 350
  xg <- seq(a, b, length=101)
  numIntKnots <- 20
  lambda <-  0.001
  
  # Set up the design matrix and related quantities:
  
  intKnots <- quantile(unique(x), seq(0, 1, length=
                                       (numIntKnots + 2))[-c(1, (numIntKnots + 2))])
  names(intKnots) <- NULL
  B <- bs(x, knots=intKnots, degree=3, Boundary.knots=c(a, b), intercept=TRUE)
  BTB <- crossprod(B, B)
  BTy <- crossprod(B, y)   
  
  # Create the Omega matrix:
  
  formOmega <- function(a, b, intKnots)
  {
    allKnots <- c(rep(a, 4), intKnots, rep(b, 4)) 
    K <- length(intKnots)
    L <- 3 * (K + 8)
    xtilde <- (rep(allKnots, each=3)[-c(1, (L - 1), L)] + 
                 rep(allKnots, each=3)[-c(1, 2, L)]) / 2
    wts <- rep(diff(allKnots), each=3) * rep(c(1, 4, 1)/6, K + 7)
    Bdd <- spline.des(allKnots, xtilde, derivs=rep(2, length(xtilde)),
                      outer.ok=TRUE)$design  
    Omega <- t(Bdd * wts) %*% Bdd     
    return(Omega)
  }
  
  Omega <- formOmega(a, b, intKnots)
  
  # Obtain the coefficients:
  
  nuHat <- solve(BTB + lambda*Omega, BTy)
  
  # For large K the following alternative Cholesky-based
  # approach can be considerably faster (O(K),  because 
  # B'B+ lambda*Omega is banded diagonal):
  
  cholFac <- chol(BTB + lambda * Omega)
  nuHat <- backsolve(cholFac, forwardsolve(t(cholFac), BTy)) 
  
  # Display the fit:
  
  Bg <- bs(xg, knots=intKnots, degree=3,
           Boundary.knots=c(a, b), intercept=TRUE)
  fhatg <- Bg %*% nuHat
  
  par(mfrow=c(1, 2))
  plot(x, y, xlim=range(xg), bty="l", type="n", xlab="radiation",
       ylab="cuberoot of ozone", main="(a) direct fit; user 
       choice of smooth. par.")
  lines(xg, fhatg, lwd=2)   
  points(x, y, lwd=2)
  
  # Mixed model scatterplot smoothing with REML choice of smoothing parameter
  # -------------------------------------------------------------------------
  
  # Obtain the spectral decomposition of $\\bOmega$:
  
  eigOmega <- eigen(Omega)
  
  # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$:
  
  indsZ <- 1:(numIntKnots+2)
  UZ <- eigOmega$vectors[, indsZ]
  LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
  
  # Perform stability check:   
  
  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX <- eigOmega$vectors[, indsX]   
  L <- cbind( UX,  LZ )
  stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))          
  if (sum(stabCheck^2) > 1.0001 * (numIntKnots + 2))
    print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")
  
  # Form the X and Z matrices:
  
  X <- cbind(rep(1, length(x)), x)
  Z <- B %*% LZ
  
  # Fit using lme() with REML choice of smoothing parameter:
  
  library(nlme)
  group <- rep(1, length(x))
  gpData <- groupedData(y~x|group, data=data.frame(x, y))
  fit <- lme(y~-1+X, random=pdIdent(~-1+Z), data=gpData)
  
  # Extract coefficients and plot scatterplot smooth over a grid:
  
  betaHat <- fit$coef$fixed
  uHat <- unlist(fit$coef$random)
  Zg <- Bg %*% LZ
  fhatgREML <- betaHat[1] + betaHat[2] * xg + Zg %*% uHat
  plot(x, y, xlim=range(xg), bty="l", type="n", xlab="radiation",
       ylab="cuberoot of ozone", main="(b) mixed model fit; 
       REML choice of smooth. par.")
  lines(xg, fhatgREML, lwd=2)   
  points(x, y, lwd=2)
  
  return(list(X=X, Z=Z))
}

########## R function: ZOSull ##########

# For creation of O'Sullivan-type Z matrices.

# Last changed: 03 SEP 2007

ZOSull <- function(x, range.x, intKnots, drv=0, stability_check=FALSE)
{
  if (drv > 2) stop("splines not smooth enough for more than 2 derivatives")

  # Set defaults for `range.x' and `intKnots'

  if (missing(range.x))
    range.x <- c(1.05 * min(x) - 0.05 * max(x), 1.05 * max(x) - 0.05 * min(x))

  if (missing(intKnots))
  {
    numIntKnots <- min(length(unique(x)), 35)
    intKnots <- quantile(unique(x), seq(0, 1, length=(numIntKnots + 2))[-c(1, (numIntKnots + 2))])
  }
  numIntKnots <- length(intKnots) 

  # Obtain the penalty matrix.

  allKnots <- c(rep(range.x[1],  4),  intKnots,  rep(range.x[2],  4)) 
  K <- length(intKnots)
  L <- 3 * (K + 8)
  xtilde <- (rep(allKnots, each=3)[-c(1, (L - 1), L)] +  
            rep(allKnots, each=3)[-c(1, 2, L)]) / 2

  wts <- rep(diff(allKnots), each=3) * rep(c(1, 4, 1) / 6, K + 7)
  Bdd <- spline.des(allKnots, xtilde, derivs=rep(2, length(xtilde)),
                   outer.ok=TRUE)$design  
  Omega     <- t(Bdd * wts) %*% Bdd     

  # Use the spectral decomposition of Omega to obtain Z.

  eigOmega <- eigen(Omega)
  indsZ <- 1:(numIntKnots + 2)
  UZ <- eigOmega$vectors[, indsZ]
  LZ <- t(t(UZ) / sqrt(eigOmega$values[indsZ]))

  if (stability_check) {
     # Perform stability check.   

     indsX <- (numIntKnots + 3):(numIntKnots + 4)
     UX <- eigOmega$vectors[, indsX]   
     L <- cbind(UX, LZ)
     stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))          
     if (sum(stabCheck ^ 2) > 1.0001 * (numIntKnots + 2))
         stop("WARNING: NUMERICAL INSTABILITY ARISING\\
                FROM SPECTRAL DECOMPOSITION")
  }

  # Obtain B and post-multiply by LZ matrix to get Z.

  B <- spline.des(allKnots, x, derivs=rep(drv, length(x)), outer.ok=TRUE)$design  

  Z <- B %*% LZ

  # Add the `range.x' and 'intKnots' as attributes
  # of the return object.

  attr(Z, "range.x") <- range.x
  attr(Z, "intKnots") <- intKnots

  # Return Z matrix with 2 attributes.

  return(Z)
}

formOmega <- function(a,b,intKnots)
{
  allKnots <- c(rep(a,2),intKnots,rep(b,4)) 
  u <- unique(allKnots)
  R <- length(u)
  xtilde <- c()
  wts <- c()
  for (i in 1:(R-1)) {
    A <- u[i]
    B <- u[i+1]
    d <- B - A
    wt <- d * c(1, 4, 1) / 6
    xtilde <- c(xtilde, A, (A + B) / 2, B)
    wts <- c(wts, wt)
  }
   
  Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),outer.ok=TRUE)$design  
  Omega <- t(Bdd * wts) %*% Bdd     
  
  return(list(Omega=Omega, allKnots=allKnots))
}

ZOSull2 <- function(x, range.x, intKnots, allKnots=NULL, Omega=NULL, drv=0, stability_check=FALSE)
{
  # Use the spectral decomposition of Omega to obtain Z.
  numIntKnots <- length(intKnots)
  eigOmega <- eigen(Omega)
  indsZ <- 1:(numIntKnots + 2)
  UZ <- eigOmega$vectors[, indsZ]
  LZ <- t(t(UZ) / sqrt(eigOmega$values[indsZ]))

  if (stability_check) {
     # Perform stability check.   

     indsX <- (numIntKnots + 3):(numIntKnots + 4)
     UX <- eigOmega$vectors[, indsX]   
     L <- cbind(UX, LZ)
     stabCheck <- t(crossprod(L, t(crossprod(L, Omega))))          
     if (sum(stabCheck ^ 2) > 1.0001 * (numIntKnots + 2))
         stop("WARNING: NUMERICAL INSTABILITY ARISING\\
                FROM SPECTRAL DECOMPOSITION")
  }

  # Obtain B and post-multiply by LZ matrix to get Z.

  B <- spline.des(allKnots, x, derivs=rep(drv, length(x)), outer.ok=TRUE)$design  

  #Z <- B %*% LZ
  Z <- B

  # Add the `range.x' and 'intKnots' as attributes
  # of the return object.

  attr(Z, "range.x") <- range.x
  attr(Z, "intKnots") <- intKnots

  # Return Z matrix with 2 attributes.

  return(Z)
}
