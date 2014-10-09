# John Ormerod's Poisson regression code, slightly modified to be used by my
# zero-inflated model code.
###############################################################################

source("CalculateB.R")

###############################################################################

f.lap <- function(vmu,vy,vr,mC,mSigma.inv,mLambda) 
{       
  d <- length(vmu)
  veta <- mC%*%vmu
  mSigma <- solve(mSigma.inv)
  mDiag <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
	f <- sum(vy*vr*veta - vr*exp(veta+0.5*mDiag)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu - 0.5*sum(diag(mLambda%*%mSigma))
  return(f)
}

###############################################################################

norm <- function(v) sqrt(sum(v^2))

###############################################################################

vg.lap <- function(vmu,vy,vr,mC,mSigma.inv,mLambda) 
{       
    vg <- t(mC)%*%(vy - vr*exp(mC%*%vmu)) - mSigma.inv%*%vmu
    
    return(vg)
}

###############################################################################

mH.lap <- function(vmu,vy,vr,mC,mSigma.inv,mLambda) 
{
    #mDiag <- sapply(1:ncol(mC), function(i) sum(mC[i,] * mLambda[, i] * mC[, i]))
    vw <- exp(mC%*%vmu); dim(vw) <- NULL
    mH <- -t(mC*vw)%*%mC - mSigma.inv
    return(mH)
}

###############################################################################

fit.Lap <- function(vmu,vy,vr,mC,mSigma.inv,mLambda) 
{
    MAXITER <- 100
    
    for (ITER in 1:MAXITER) {
        f  <- f.lap(vmu,vy,vr,mC,mSigma.inv,mLambda)
        vg <- vg.lap(vmu,vy,vr,mC,mSigma.inv,mLambda)
        mH <- mH.lap(vmu,vy,vr,mC,mSigma.inv,mLambda)
        mLambda <- solve(-mH,tol=1.0E-99)
        vmu <- vmu + mLambda%*%vg
        #print(c(ITER,f,max(abs(vg))))
        #cat("iterations", ITER, "f", f, "max(abs(vg))", max(abs(vg)))
        if (max(abs(vg))<1.0E-8) {
            break;
        }
    } 

    f <- f.lap(vmu,vy,vr,mC,mSigma.inv,mLambda) + 0.5*log(det(mLambda%*%mSigma.inv))
    return(list(vmu=vmu,mLambda=mLambda,f=f))
}

###############################################################################

f.G <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
{
  d <- length(vmu)
  
  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  vB0 <- B0.fun("POISSON",vmu.til,vsigma2.til,gh) 
  
  f <- sum(vr*(vy*vmu.til - vB0)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu 
  f <- f - 0.5*d*log(2*pi) + 0.5*log(det(mSigma.inv)) 
  return(f)
}

###############################################################################

f.GVA <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
   d <- ncol(mC)  
   vmu <- vtheta[1:d]
   mR[Rinds] <- vtheta[(1+d):length(vtheta)]
   mR[Dinds] <- exp(mR[Dinds]) 
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
    }   
  mLambda <- mR%*%t(mR)   
   
   f <- sum(log(diag(mR))) + f.G(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
   f <- f + 0.5*d*log(2*pi) + 0.5*d
   
   if (!is.finite(f)) {
       f <- -1.0E16
   }
   return(f)
}

###############################################################################

vg.G <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
{
    vg <- t(mC)%*%(vr*(vy - vB1)) - mSigma.inv%*%vmu     
    return(vg)
}

###############################################################################

mH.G <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2) 
{
    vw <-  vB2; dim(vw) <- NULL
    mH <- -t(mC*vw)%*%mC - mSigma.inv
    return(mH)    
}

###############################################################################

vg.GVA.approx <- function(vmu,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  n <- length(vy); P <- length(vmu); d <- ncol(mC)
	eps <- 1.0E-6
	f <- f.GVA(vmu,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	vg.approx <- matrix(0,P,1)
	for (i in 1:P) {
	   vmup <- vmu 
	   vmup[i] <- vmu[i] + eps
	   fp <- f.GVA(vmup,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	   
	   vmum <- vmu 
	   vmum[i] <- vmu[i] - eps
	   fm <- f.GVA(vmum,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	   
	   vg.approx[i] <- (fp - fm)/(2*eps)
	}
	return(vg.approx)
}

###############################################################################

vg.GVA <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  d <- ncol(mC)
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]

  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
  }    
  mLambda <- mR%*%t(mR)   

  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  vg <- 0*vmu
  vg[1:d] <- vg.G(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 

  mLambda.inv <- solve(mR%*%t(mR))
  mH <- mH.G(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2)
  dmLambda <- (mLambda.inv + mH)%*%mR
  
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
 
  return(vg)
}

###############################################################################

fit.GVA <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-12)
{
	gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))
		
  d <- length(vmu)
  Dinds <- d*((1:d)-1)+(1:d)
          
  mR <- t(chol(mLambda + diag(1.0E-8,d)))
  mR[Dinds] <- log(mR[Dinds])
  Rinds <- which(lower.tri(mR,diag=TRUE))
	vmu <- c(vmu,mR[Rinds])
  P <- length(vmu)
  lower_constraint <- rep(-Inf, length(vmu))
  
  if (method=="L-BFGS-B") {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
  } else {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
  }
  res <- optim(par=vmu, fn=f.GVA, gr=vg.GVA,
      method=method,lower=lower_constraint, upper=Inf, control=controls,
      vy=vy,vr=vr,mC=mC,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds)        
      
  vtheta <- res$par 
  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):P]
  mR[Dinds] <- exp(mR[Dinds])  
  mLambda <- mR%*%t(mR)
     
  return(list(res=res,vmu=vmu,mLambda=mLambda))
}
###############################################################################

f.G_new <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
{
  d <- length(vmu)
  
  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  vB0 <- B0.fun("POISSON",vmu.til,vsigma2.til,gh) 
  
  f <- sum(vr*(vy*vmu.til - vB0)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu 
  f <- f - 0.5*d*log(2*pi) + 0.5*log(det(mSigma.inv)) 
  return(f)
}

###############################################################################

f.GVA_new <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  d <- ncol(mC)  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]
  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
  }   
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  
  f <- -sum(log(diag(mR))) + f.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
  f <- f + 0.5*d*log(2*pi) + 0.5*d
  
  if (!is.finite(f)) {
    f <- -1.0E16
  }
  return(f)
}

###############################################################################

vg.G_new <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
{
  vg <- t(mC)%*%(vr*(vy - vB1)) - mSigma.inv%*%vmu     
  #cat("vg.G", vg, "\n")
  return(vg)
}

###############################################################################

mH.G_new <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2) 
{
  vw <-  vB2; dim(vw) <- NULL
  mH <- -t(mC*vw)%*%mC - mSigma.inv
  return(mH)    
}

###############################################################################

vg.GVA.approx_new <- function(vmu,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  n <- length(vy); P <- length(vmu); d <- ncol(mC)
  eps <- 1.0E-6
  f <- f.GVA_new(vmu,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
  vg.approx <- matrix(0,P,1)
  for (i in 1:P) {
    vmup <- vmu 
    vmup[i] <- vmu[i] + eps
    fp <- f.GVA_new(vmup,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
    
    vmum <- vmu 
    vmum[i] <- vmu[i] - eps
    fm <- f.GVA_new(vmum,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
    
    vg.approx[i] <- (fp - fm)/(2*eps)
  }
  return(vg.approx)
}

###############################################################################

vg.GVA_new <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  d <- ncol(mC)
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]
  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
  }    
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  
  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  vg <- 0*vmu
  vg[1:d] <- vg.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
  
  mLambda.inv <- mR%*%t(mR)
  mH <- mH.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2)
  dmLambda <- -solve(mR, tol=1.0E-99)%*%(mLambda.inv + mH)%*%mLambda
  
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
  
  return(vg)
}

###############################################################################

fit.GVA_new <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-12)
{
  #library(statmod)
  #N <- 15
  #gh  <- gauss.quad(N,kind="hermite")
  gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))
  
  d <- length(vmu)
  Dinds <- d*((1:d)-1)+(1:d)
  
  mR <- t(chol(solve(mLambda, tol=1.0E-99) + diag(1.0E-8,d)))
  cat("mR", mR, "\n")
  mR[Dinds] <- log(mR[Dinds])
  Rinds <- which(lower.tri(mR,diag=TRUE))
  vmu <- c(vmu,mR[Rinds])
  P <- length(vmu)
  lower_constraint <- rep(-Inf, length(vmu))
  
  if (method=="L-BFGS-B") {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
  } else {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
  }
  res <- optim(par=vmu, fn=f.GVA_new, gr=vg.GVA_new,
               method=method,lower=lower_constraint, upper=Inf, control=controls,
               vy=vy,vr=vr,mC=mC,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds)        
  
  vtheta <- res$par 
  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):P]
  mR[Dinds] <- exp(mR[Dinds])  
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  
  return(list(res=res,vmu=vmu,mLambda=mLambda))
}

###############################################################################

vg.G_nr <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
{
  vg <- t(mC)%*%(vr*(vy - vB1)) - mSigma.inv%*%vmu     
  return(vg)
}

###############################################################################

mH.G_nr <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2) 
{
  mH <- -t(mC*as.vector(vr*vB2))%*%mC - mSigma.inv 
  return(mH)
}

###############################################################################

fit.GVA_nr <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-12)
{
 
  MAXITER <- 1000
  TOL <- reltol
  
  for (ITER in 1:MAXITER) 
  {
      vmu.old <- vmu
      
      # calculate B1
      # calculate B2
      vmu.til <- mC%*%vmu
      vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
      res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
      vB1 <- res.B12$vB1
      vB2 <- res.B12$vB2      
    
      vg <- vg.G_nr(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
      mH <- mH.G_nr(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2) 
      
      mLambda <- solve(-mH,tol=1.0E-99)
      vmu <- vmu + mLambda%*%vg
        
      err <- max(abs(vmu - vmu.old)) 
      cat(ITER,err,"\n")
      if (err<TOL) {
        break;
      }
  } 
  
  return(list(vmu=vmu,mLambda=mLambda))
}
