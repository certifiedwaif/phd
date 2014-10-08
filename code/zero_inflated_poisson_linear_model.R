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
	# diag(mC mLambda mC^T)_i i = mC[i, ] . mLambda[, i] . mC[, i] 
    #f <- sum(vy*veta - vr*exp(veta+0.5*diag(mC%*%mLambda%*%t(mC)))) - 0.5*t(vmu)%*%mSigma.inv%*%vmu - 0.5*sum(diag(mLambda%*%mSigma))
	#mDiag <-  diag(mC%*%mLambda%*%t(mC))
  	mDiag <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
	f <- sum(vy*vr*veta - vr*exp(veta+0.5*mDiag)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu - 0.5*sum(diag(mLambda%*%mSigma))
    return(f)
}

###############################################################################

norm <- function(v) sqrt(sum(v^2))

###############################################################################

vg.lap <- function(vmu,vy,vr,mC,mSigma.inv,mLambda) 
{       
	#cat("vg.lap vmu", vmu, "\n")
	#cat("vg.lap vy", head(vy), "\n")
	#cat("vg.lap vr", head(vr), "\n")
	#cat("vg.lap mC", head(mC), "\n")
	#mDiag <- sapply(1:ncol(mC), function(i) sum(mC[i,] * mLambda[, i] * mC[, i]))
	#mDiag2 <-  diag(mC%*%mLambda%*%t(mC))
    vg <- t(mC)%*%(vy - vr*exp(mC%*%vmu)) - mSigma.inv%*%vmu
    
    #print(mC%*%vmu)
    
    
	#cat("vg.lap vg", vg, "norm", norm(vg), "\n")
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

    #vmu <- vmu
    #controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=1.0e-15)
    #res <- optim(par=vmu, fn=f.lap, gr=vg.lap,
    #    method="BFGS",lower=-Inf, upper=Inf, control=controls,
    #    vy=vy,mC=mC,mSigma.inv) 
        
    #vmu <- res$par
    #mLambda <- solve(-mH.lap(vmu,vy,mC,mSigma.inv),tol=1.0E-99)
    
    f <- f.lap(vmu,vy,vr,mC,mSigma.inv,mLambda) + 0.5*log(det(mLambda%*%mSigma.inv))
    #print(f)
    return(list(vmu=vmu,mLambda=mLambda,f=f))
}

###############################################################################

f.G <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
{
	#print(vmu)
	#print(mSigma.inv)
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
	# Old parameterisation
  mLambda <- mR%*%t(mR)   
	# New parameterisation
   #mLambda <- solve(mR%*%t(mR))
  #mLambda <- solve(t(mR)%*%mR)
   
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
	#cat("vg.G", vg, "\n")
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
	   
	   #vmupp <- vmu 
	   #vmupp[i] <- vmu[i] + 2*eps
	   #fpp <- f.GVA(vmupp,vy,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	   
	   #vmumm <- vmu 
	   #vmumm[i] <- vmu[i] - 2*eps
	   #fmm <- f.GVA(vmumm,vy,mC,mSigma.inv,gh,mR,Rinds,Dinds)  
	   
	   #vg.approx[i] <- (-fpp + 8*fp - 8*fm + fmm)/(12*eps)
	}
	return(vg.approx)
}

###############################################################################

vg.GVA <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
	#cat("vtheta", vtheta, "\n")
	#cat("ncol(mC)", ncol(mC), "\n")
  d <- ncol(mC)
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]

  #diag(mR) <- abs(diag(mR))
  #cat("mR", mR, "\n")
  #diag(mR) <- abs(diag(mR)) # FIXME: This may be the wrong way to solve this. Get outside the
                            # allowable solution space on every second iteration, diagonal entries
                            # of mR become negative.
  #cat("mR after diag. abs", mR, "\n")

  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
  }    
  # Old parameterisation
  mLambda <- mR%*%t(mR)   
  # New parameterisation
  #mLambda <- solve(mR%*%t(mR),tol=1.0E-99)
  #mLambda <- solve(t(mR)%*%mR)

  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  vg <- 0*vmu
  vg[1:d] <- vg.G(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 

  # Old parameterisation
  mLambda.inv <- solve(mR%*%t(mR))
  # New parameterisation
  #mLambda.inv <- t(mR)%*%mR
  mH <- mH.G(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2)
  # Old parameterisation
  dmLambda <- (mLambda.inv + mH)%*%mR
  # New parameterisation
  #dmLambda <- -solve(mR)%*%(mLambda.inv + mH)%*%mLambda
  #dmLambda <- -solve(t(mR))%*%(mLambda.inv + mH)%*%mLambda
  
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  # Worth a try ...
  #dmLambda[Dinds] <- log(dmLambda[Dinds])
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
 
  cat("vtheta.GVA vtheta", round(vtheta[(1+d):length(vtheta)], 2), "norm", norm(vtheta[(1+d):length(vtheta)]), "\n")
  cat("vg.GVA vg", round(vg[(1+d):length(vtheta)], 2), "norm", norm(vg[(1+d):length(vtheta)]), "\n")

  return(vg)
}

###############################################################################

fit.GVA <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-12)
{
	#library(statmod)
	#N <- 15
	#gh  <- gauss.quad(N,kind="hermite")
	gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))
		
  d <- length(vmu)
	#cat("vmu", vmu, "\n")
	#cat("d", d, "\n")
  Dinds <- d*((1:d)-1)+(1:d)
          
  mR <- t(chol(mLambda + diag(1.0E-8,d)))
  cat("mR", mR, "\n")
  mR[Dinds] <- log(mR[Dinds])
  Rinds <- which(lower.tri(mR,diag=TRUE))
	#Rinds <- which(upper.tri(mR,diag=TRUE))
  #cat("Rinds", Rinds, "\n")
	#cat("mR", mR, "\n")
	vmu <- c(vmu,mR[Rinds])
  #cat("vmu", vmu, "\n")
  P <- length(vmu)
  lower_constraint <- rep(-Inf, length(vmu))
  #lower_constraint[d+Dinds] <- 0
  #cat("lower_constraint", lower_constraint, "\n")
  
  if (method=="L-BFGS-B") {
    controls <- list(maxit=1000,trace=1,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
  } else {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
  }
  res <- optim(par=vmu, fn=f.GVA, gr=vg.GVA,
      method=method,lower=lower_constraint, upper=Inf, control=controls,
      vy=vy,vr=vr,mC=mC,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds)        
      
  vtheta <- res$par 
  
  #f <- f.GVA(vmu,vy,mC,CTC,CTy,p,K,mR,Rinds,Dinds)
  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):P]
  mR[Dinds] <- exp(mR[Dinds])  
  # Old parameterisation
  mLambda <- mR%*%t(mR)
  # New parameterisation
  #mLambda <- solve(t(mR)%*%mR)
     
  return(list(res=res,vmu=vmu,mLambda=mLambda))
}
###############################################################################

f.G_new <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
{
  #print(vmu)
  #print(mSigma.inv)
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
  # Old parameterisation
  #mLambda <- mR%*%t(mR)   
  # New parameterisation
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  #mLambda <- solve(t(mR)%*%mR)
  
  f <- -sum(log(diag(mR))) + f.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
  #f <- 0.5*sum(log(det(mLambda))) + f.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,gh) 
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
    
    #vmupp <- vmu 
    #vmupp[i] <- vmu[i] + 2*eps
    #fpp <- f.GVA(vmupp,vy,mC,mSigma.inv,gh,mR,Rinds,Dinds)
    
    #vmumm <- vmu 
    #vmumm[i] <- vmu[i] - 2*eps
    #fmm <- f.GVA(vmumm,vy,mC,mSigma.inv,gh,mR,Rinds,Dinds)  
    
    #vg.approx[i] <- (-fpp + 8*fp - 8*fm + fmm)/(12*eps)
  }
  return(vg.approx)
}

###############################################################################

vg.GVA_new <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  #cat("vtheta", vtheta, "\n")
  #cat("ncol(mC)", ncol(mC), "\n")
  d <- ncol(mC)
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]
  #diag(mR) <- abs(diag(mR))
  #cat("mR", mR, "\n")
  #diag(mR) <- abs(diag(mR)) # FIXME: This may be the wrong way to solve this. Get outside the
  # allowable solution space on every second iteration, diagonal entries
  # of mR become negative.
  #cat("mR after diag. abs", mR, "\n")
  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
  }    
  # Old parameterisation
  #mLambda <- mR%*%t(mR)   
  # New parameterisation
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  
  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  vg <- 0*vmu
  vg[1:d] <- vg.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
  
  # Old parameterisation
  #mLambda.inv <- solve(mR%*%t(mR))
  # New parameterisation
  mLambda.inv <- mR%*%t(mR)
  #mLambda.inv <- t(mR)%*%mR
  mH <- mH.G_new(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2)
  # Old parameterisation
  #dmLambda <- (mLambda.inv + mH)%*%mR
  # New parameterisation
  #dmLambda <- -solve(mR)%*%(mLambda.inv + mH)%*%mLambda
  # FIXME: I think this derivative must be wrong. Check again.
  dmLambda <- -solve(mR, tol=1.0E-99)%*%(mLambda.inv + mH)%*%mLambda
  #dmLambda <- -mLambda%*%(mLambda.inv + mH)%*%mLambda%*%mR
  #browser()
  #dmLambda <- -solve(t(mR))%*%(mLambda.inv + mH)%*%mLambda
  
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  # Worth a try ...
  #dmLambda[Dinds] <- log(dmLambda[Dinds])
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
  
  #cat("vtheta.GVA vtheta", vtheta, "norm", norm(vtheta), "\n")
  #cat("vg.GVA vg[mR]", vg, "norm", norm(vg), "\n")
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
  #cat("vmu", vmu, "\n")
  #cat("d", d, "\n")
  Dinds <- d*((1:d)-1)+(1:d)
  
  mR <- t(chol(solve(mLambda, tol=1.0E-99) + diag(1.0E-8,d)))
  cat("mR", mR, "\n")
  mR[Dinds] <- log(mR[Dinds])
  Rinds <- which(lower.tri(mR,diag=TRUE))
  #Rinds <- which(upper.tri(mR,diag=TRUE))
  cat("Rinds", Rinds, "\n")
  cat("mR", mR, "\n")
  vmu <- c(vmu,mR[Rinds])
  cat("vmu", vmu, "\n")
  P <- length(vmu)
  lower_constraint <- rep(-Inf, length(vmu))
  #lower_constraint[d+Dinds] <- 0
  #cat("lower_constraint", lower_constraint, "\n")
  
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
  
  #f <- f.GVA(vmu,vy,mC,CTC,CTy,p,K,mR,Rinds,Dinds)
  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):P]
  mR[Dinds] <- exp(mR[Dinds])  
  # Old parameterisation
  #mLambda <- mR%*%t(mR)
  # New parameterisation
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  #mLambda <- solve(t(mR)%*%mR)
  
  return(list(res=res,vmu=vmu,mLambda=mLambda))
}

# NR scratch
# Re-use one of the existing f's
# Implement new vg
# Is there any point to using the Cholesky factor rather than just storing the
# upper or lower triangular matrix?
###############################################################################

vg.G_nr <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
{
  vg <- mLambda %*% t(mC)%*%(vr*(vy - vB1)) - mSigma.inv%*%vmu     
  #cat("vg.G", vg, "\n")
  return(vg)
}

###############################################################################

vg.GVA_nr <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  #cat("vtheta", vtheta, "\n")
  #cat("ncol(mC)", ncol(mC), "\n")
  d <- ncol(mC)
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]
  
  #diag(mR) <- abs(diag(mR))
  #cat("mR", mR, "\n")
  #diag(mR) <- abs(diag(mR)) # FIXME: This may be the wrong way to solve this. Get outside the
  # allowable solution space on every second iteration, diagonal entries
  # of mR become negative.
  #cat("mR after diag. abs", mR, "\n")
  
  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
  }    
  # Old parameterisation
  mLambda <- mR%*%t(mR)   
  # New parameterisation
  #mLambda <- solve(mR%*%t(mR),tol=1.0E-99)
  #mLambda <- solve(t(mR)%*%mR)
  
  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  vg <- 0*vmu
  vg[1:d] <- vg.G(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
  
  # Old parameterisation
  mLambda.inv <- solve(mR%*%t(mR))
  # New parameterisation
  #mLambda.inv <- t(mR)%*%mR
  # Old parameterisation
  #dmLambda <- (mLambda.inv + mH)%*%mR
  cat("dim(mC)", dim(mC), "\n")
  cat("dim(t(mC))", dim(t(mC)), "\n")
  cat("vB2", vB2, "\n")
  cat("class(vB2)", class(vB2), "\n")
  cat("dim(vB2)", dim(vB2))
  cat("dim(diag(vB2))", dim(diag(as.vector(vB2), nrow = length(vB2))), "\n")
  cat("length(vB2)", length(vB2), "\n")
  dmLambda <- solve(t(mC)%*%diag(as.vector(vB2), nrow = length(vB2))%*%mC + mSigma.inv) - mLambda
  
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  # Worth a try ...
  #dmLambda[Dinds] <- log(dmLambda[Dinds])
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
  
  cat("vtheta.GVA vtheta", round(vtheta[(1+d):length(vtheta)], 2), "norm", norm(vtheta[(1+d):length(vtheta)]), "\n")
  cat("vg.GVA vg", round(vg[(1+d):length(vtheta)], 2), "norm", norm(vg[(1+d):length(vtheta)]), "\n")
  
  return(vg)
}

###############################################################################

fit.GVA_nr <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-12)
{
  #library(statmod)
  #N <- 15
  #gh  <- gauss.quad(N,kind="hermite")
  gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))
  
  d <- length(vmu)
  #cat("vmu", vmu, "\n")
  #cat("d", d, "\n")
  Dinds <- d*((1:d)-1)+(1:d)
  
  mR <- t(chol(mLambda + diag(1.0E-8,d)))
  cat("mR", mR, "\n")
  mR[Dinds] <- log(mR[Dinds])
  Rinds <- which(lower.tri(mR,diag=TRUE))
  #Rinds <- which(upper.tri(mR,diag=TRUE))
  cat("Rinds", Rinds, "\n")
  cat("mR", mR, "\n")
  vmu <- c(vmu,mR[Rinds])
  cat("vmu", vmu, "\n")
  P <- length(vmu)
  lower_constraint <- rep(-Inf, length(vmu))
  #lower_constraint[d+Dinds] <- 0
  cat("lower_constraint", lower_constraint, "\n")
  
  if (method=="L-BFGS-B") {
    controls <- list(maxit=1000,trace=1,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
  } else {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
  }
  res <- optim(par=vmu, fn=f.GVA, gr=vg.GVA_nr,
               method=method,lower=lower_constraint, upper=Inf, control=controls,
               vy=vy,vr=vr,mC=mC,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds)        
  
  vtheta <- res$par 
  
  #f <- f.GVA(vmu,vy,mC,CTC,CTy,p,K,mR,Rinds,Dinds)
  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):P]
  mR[Dinds] <- exp(mR[Dinds])  
  # Old parameterisation
  mLambda <- mR%*%t(mR)
  # New parameterisation
  #mLambda <- solve(t(mR)%*%mR)
  
  return(list(res=res,vmu=vmu,mLambda=mLambda))
}
