###############################################################################

rawdata <- matrix(c(1, 53,
                    1, 57,
                    1, 58,
                    1, 63,
                    0, 66,
                    0, 67,
                    0, 67,
                    0, 67,
                    0, 68,
                    0, 69,
                    0, 70,
                    0, 70,
                    1, 70,
                    1, 70,
                    0, 72,
                    0, 73,
                    0, 75,
                    1, 75,
                    0, 76,
                    0, 76,
                    0, 78,
                    0, 79,
                    0, 81),23,2,byrow=T)

w <- rawdata[,1]
s <- rawdata[,2]

vy <- w
n  <- length(vy)
mX <- cbind(rep(1,n),s-mean(s))
d  <- ncol(mX)
mSigma     <- diag(1.0E4,d)
mSigma.inv <- diag(1.0E-4,d)

res.glm <- glm(vy~-1+mX,family=binomial(link=probit))
summary(res.glm)

###############################################################################

trapint2D <- function(xgrid, ygrid, fgrid) 
{
  ng1   <- length(xgrid)
  ng2   <- length(ygrid)
  area <- 0   
  for (i in 2:ng1) {
    for (j in 2:ng2) {
      area <- area + (xgrid[i] - xgrid[i-1])*(ygrid[j] - ygrid[j-1])*(fgrid[i,j] + fgrid[i-1,j] + fgrid[i,j-1] + fgrid[i-1,j-1])/4
    }
  }
  return(area)
}

djointLMM <- function(u,v,vy,mX,mSigma.inv) {
  log.f <- matrix(0,length(u),length(v))
  for (i in 1:length(u)) {
    for (j in 1:length(v)) {
      vtheta <- matrix(c(u[i],v[j]))
      veta <- mX%*%vtheta
      log.f[i,j] <- sum(vy*pnorm(veta, log.p=TRUE) + (1-vy)*pnorm(veta, lower.tail=FALSE, log.p=TRUE)) - 0.5*t(vtheta)%*%mSigma.inv%*%vtheta
    }
  }
  f <- exp(log.f - max(log.f))
  area <- trapint2D(u,v,f) 
  return(f/area)
}

u <- seq(-3,0.5,,100)
v <- seq(-0.7,0.1,,100)

f <- djointLMM(u,v,vy,mX,mSigma.inv)
levels <- pretty(range(f),12)

contour(u,v,f,levels=levels,lwd=2,main="Exact Posterior Density",cex.main=3)

ans <- readline()

###############################################################################

library(mvtnorm)

fast.f <- function(mTheta,vy,mX,mSigma) 
{
  N <- nrow(mTheta)
  n <- length(vy)
  mY <- matrix(vy,n,N)
  vEta <- mX%*%t(mTheta)
  
  #vy*pnorm(veta, log.p=TRUE) + (1-vy)*pnorm(veta, lower.tail=FALSE, log.p=TRUE)
  log.vp <- matrix(1,1,n)%*%(vy*pnorm(vEta, log.p=TRUE) + (1-vy)*pnorm(vEta, lower.tail=FALSE, log.p=TRUE)) + dmvnorm(mTheta,sigma=mSigma,log=TRUE)
  return(log.vp)
} 

###############################################################################

ImportanceSampling <- function(N,vmu,mLambda,nu,vy,mX,mSigma,EPS) 
{
  MAXITER <- 10000000
  log.vq <- c()
  log.vp <- c()
  for (ITER in 1:MAXITER) 
  {
    mZ <- rmvt(N, delta=as.vector(vmu), sigma=mLambda, df = nu, type = "shifted")
    log.vq.new <- dmvt(mZ, delta=as.vector(vmu), sigma=mLambda, df = nu, log=TRUE, type = "shifted")
    log.vp.new <- as.vector(fast.f(mZ,vy,mX,mSigma)) 
    log.vq <- c(log.vq,log.vq.new)
    log.vp <- c(log.vp,log.vp.new)				
    vw <- exp(log.vp - log.vq)
    N <- length(vw)
    I.hat <- mean(vw)
    se <- sqrt(var(vw)/length(vw))
    print(c(ITER,N,I.hat,se))
    if (se<EPS) {
      break;
    }
  }
  return(list(I.hat=I.hat,se=se,vw=vw))
}

fast.f2 <- function(mult, mTheta,vy,vr,mX,mSigma) 
{
  N <- nrow(mTheta)
  n <- length(vy)
  mY <- matrix(vy,n,N)
  vEta <- mX%*%t(mTheta)
  
  #vy*pnorm(veta, log.p=TRUE) + (1-vy)*pnorm(veta, lower.tail=FALSE, log.p=TRUE)
  #log.vp <- matrix(1,1,n)%*%(vy*pnorm(vEta, log.p=TRUE) + (1-vy)*pnorm(vEta, lower.tail=FALSE, log.p=TRUE)) + dmvnorm(mTheta,sigma=mSigma,log=TRUE)

  # Log-likelihood pertaining to vbeta and vu
  veta = mC%*%vtheta  
  log.vp <- t(vy*vr)%*%veta - t(vr)%*%exp(veta) - sum(lgamma(vy+1))
  
  # Prior
  mSigma.vbeta = solve(mult$mSigma.beta.inv)
  mSigma.vu = solve(mult$mSigma.u.inv)
  log.vp <- log.vp + dmvnorm(vtheta, blockDiag(mSigma.vbeta, mSigma.vu), log=TRUE)
  
  return(log.vp)
} 

###############################################################################

NormalisedImportanceSampling <- function(N,vmu,mLambda,nu,vy,mX,mSigma) 
{
  log.vq <- c()
  log.vp <- c()
  
  mZ <- rmvt(N, delta=as.vector(vmu), sigma=mLambda, df = nu, type = "shifted")
  log.vq <- dmvt(mZ, delta=as.vector(vmu), sigma=mLambda, df = nu, log=TRUE, type = "shifted")
  log.vp <- as.vector(fast.f(mZ,vy,mX,mSigma)) 
  
  vw <- exp(log.vp - log.vq)
  vw <- vw/sum(vw)
  
  EXPECTED.VAL <- apply(t(mZ)*vw,1,sum)
  
  return(list(EXPECTED.VAL=EXPECTED.VAL))
}

###############################################################################

LaplaceApproxPosterior <- function(vy,mX,mSigma.inv) 
{
  n <- length(vy)
  d <- ncol(mX)
  vtheta <- matrix(0,d,1)
  for (ITER in 1:1000) 
  {
    #vmu <- 1/(1+exp(-mX%*%vtheta))
    vmu <- pnorm(mX%*%vtheta)
    vg <- t(mX)%*%(vy - vmu) - mSigma.inv%*%vtheta
    mH <- - t(mX*as.vector(vmu*(1 - vmu)))%*%mX - mSigma.inv
    vtheta <- vtheta - solve(mH)%*%vg
    err <- max(abs(vg))
    if (err<1.0E-6) {
      break;
    }
    print(c(ITER,err))
  }
  return(list(vtheta=vtheta,mLambda=solve(-mH)))
}


###############################################################################

res.lap <- LaplaceApproxPosterior(vy,mX,mSigma.inv) 

vmu     <- res.lap$vtheta
mLambda <- res.lap$mLambda
nu      <- 8
EPS     <- 1.0E-12
N       <- 1000

res.is <- ImportanceSampling(N,vmu,mLambda,nu,vy,mX,mSigma,EPS) 

plot(density(res.is$vw))


res.nis <- NormalisedImportanceSampling(100*length(res.is$vw),vmu,mLambda,nu,vy,mX,mSigma) 
