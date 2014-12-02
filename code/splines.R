# splines.R

# Code corresonding to the appendix of the paper:
# Wand, M.P. & Ormerod, J.T. (2008).
# `On semiparametric regression with O'Sullivan
#  penalised splines'.

# Last changed: 29 APR 2010

# Direct scatterplot smoothing with user choice of smoothing parameter
# --------------------------------------------------------------------

# Obtain scatterplot data corresponding to environmental
# data from the R package `lattice'. Set up plotting 
# grid, knots and smoothing parameter:

# Outline:
# Need to generate spline data

fit_spline = function(x, y)
{
  library(lattice) ; attach(environmental)
  # x could be anything you like
  # y could be anything you like
  a = min(x)
  b = max(x)
  #x <- radiation ; y <- ozone^(1/3)
  #a <- 0 ; b <- 350 ;
  xg <- seq(a,b,length=101)
  numIntKnots <- 20 ; lambda <-  1000
  
  # Set up the design matrix and related quantities:
  
  library(splines)
  intKnots <- quantile(unique(x),seq(0,1,length=
                                       (numIntKnots+2))[-c(1,(numIntKnots+2))])
  names(intKnots) <- NULL
  B <- bs(x,knots=intKnots,degree=3,Boundary.knots=c(a,b),intercept=TRUE)
  BTB <- crossprod(B,B) ; BTy <- crossprod(B,y)   
  
  # Create the Omega matrix:
  
  formOmega <- function(a,b,intKnots)
  {
    allKnots <- c(rep(a,4),intKnots,rep(b,4)) 
    K <- length(intKnots) ; L <- 3*(K+8)
    xtilde <- (rep(allKnots,each=3)[-c(1,(L-1),L)]+ 
                 rep(allKnots,each=3)[-c(1,2,L)])/2
    wts <- rep(diff(allKnots),each=3)*rep(c(1,4,1)/6,K+7)
    Bdd <- spline.des(allKnots,xtilde,derivs=rep(2,length(xtilde)),
                      outer.ok=TRUE)$design  
    Omega     <- t(Bdd*wts)%*%Bdd     
    return(Omega)
  }
  
  Omega <- formOmega(a,b,intKnots)
  
  # Obtain the coefficients:
  
  nuHat <- solve(BTB+lambda*Omega,BTy)
  
  # For large K the following alternative Cholesky-based
  # approach can be considerably faster (O(K), because 
  # B'B+ lambda*Omega is banded diagonal):
  
  cholFac <- chol(BTB+lambda*Omega)
  nuHat <- backsolve(cholFac,forwardsolve(t(cholFac),BTy)) 
  
  # Display the fit:
  
  Bg <- bs(xg,knots=intKnots,degree=3,
           Boundary.knots=c(a,b),intercept=TRUE)
  fhatg <- Bg%*%nuHat
  par(mfrow=c(1,2))
  plot(x,y,xlim=range(xg),bty="l",type="n",xlab="radiation",
       ylab="cuberoot of ozone",main="(a) direct fit; user 
       choice of smooth. par.")
  lines(xg,fhatg,lwd=2)   
  points(x,y,lwd=2)
  
  # Mixed model scatterplot smoothing with REML choice of smoothing parameter
  # -------------------------------------------------------------------------
  
  # Obtain the spectral decomposition of $\\bOmega$:
  
  eigOmega <- eigen(Omega)
  
  # Obtain the matrix for linear transformation of $\\bB$ to $\\bZ$:
  
  indsZ <- 1:(numIntKnots+2)
  UZ <- eigOmega$vectors[,indsZ]
  LZ <- t(t(UZ)/sqrt(eigOmega$values[indsZ]))
  
  # Perform stability check:   
  
  indsX <- (numIntKnots+3):(numIntKnots+4)
  UX <- eigOmega$vectors[,indsX]   
  L <- cbind( UX, LZ )
  stabCheck <- t(crossprod(L,t(crossprod(L,Omega))))          
  if (sum(stabCheck^2) > 1.0001*(numIntKnots+2))
    print("WARNING: NUMERICAL INSTABILITY ARISING FROM SPECTRAL DECOMPOSITION")
  
  # Form the X and Z matrices:
  
  X <- cbind(rep(1,length(x)),x)
  Z <- B%*%LZ
  
  return(list(X=X, Z=Z))
  # Fit using lme() with REML choice of smoothing parameter:
  
  library(nlme)
  group <- rep(1,length(x))
  gpData <- groupedData(y~x|group,data=data.frame(x,y))
  fit <- lme(y~-1+X,random=pdIdent(~-1+Z),data=gpData)
  
  # Extract coefficients and plot scatterplot smooth over a grid:
  
  betaHat <- fit$coef$fixed
  uHat <- unlist(fit$coef$random)
  Zg <- Bg%*%LZ
  fhatgREML <- betaHat[1] + betaHat[2]*xg + Zg%*%uHat
  plot(x,y,xlim=range(xg),bty="l",type="n",xlab="radiation",
       ylab="cuberoot of ozone",main="(b) mixed model fit; 
       REML choice of smooth. par.")
  lines(xg,fhatgREML,lwd=2)   
  points(x,y,lwd=2)
}

m = 50
n = rep(1, m)
mX = matrix(as.vector(cbind(rep(1, m), runif(m, -1, 1))), m, 2)
mZ = NULL
expected_rho = 1
expected_mu = c(2, 1)
expected_sigma2_u = 0
sigma2.beta = 1e5
a_sigma = 1e5
b_sigma = 1e5
test_data = generate_multivariate_test_data(mX, NULL, m, n, expected_rho, expected_mu, expected_sigma2_u)
vy = test_data$vy
vy = vy ^ 3

result = fit_spline(mX[,2], vy)
