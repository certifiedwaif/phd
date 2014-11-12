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
    vg <- t(mC)%*%(vr*vy - vr*exp(mC%*%vmu)) - mSigma.inv%*%vmu
    
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
  
  f <- sum(vr*(vy*vmu.til - vB0)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu - 0.5*tr(mSigma.inv%*%mLambda)
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
    #mH <- -t(mC*vw)%*%mC - mSigma.inv
    mH <- -t(mC*vr*vw)%*%(mC*vr) - mSigma.inv
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
  
  f <- sum(vr*(vy*vmu.til - vB0)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu - 0.5*tr(mSigma.inv%*%mLambda)
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
  mLambda.inv = mR%*%t(mR)
  mLambda <- solve(mLambda.inv, tol=1.0E-99)
  
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
  #mH <- -t(mC*vw)%*%mC - mSigma.inv
  mH <- -t(mC*vr*vw)%*%(mC*vr) - mSigma.inv
  #mH <- -t(mC*(vr*vw))%*%(mC)    - mSigma.inv
  return(mH)    
}

###############################################################################

vg.GVA.approx_new <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  n <- length(vy); P <- length(vtheta); d <- ncol(mC)
  eps <- 1.0E-6
  f <- f.GVA_new(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
  vg.approx <- matrix(0,P,1)
  for (i in 1:P) {
    vthetap <- vtheta
    vthetap[i] <- vtheta[i] + eps
    fp <- f.GVA_new(vthetap,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
    
    vthetam <- vtheta
    vthetam[i] <- vtheta[i] - eps
    fm <- f.GVA_new(vthetam,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)

    
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
  # mR is lower triangular. Can you rewrite this using forward solves and
  # backsolves?
  dmLambda <- -solve(mR, tol=1.0E-99)%*%(mLambda.inv + mH)%*%mLambda  
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]

  #res <- vg.GVA.approx_new(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
 
  #print("john test")

  #print("mLambda.inv")
  #print(mLambda.inv)
  
  #print("mH")
  #print(mH)
  
  #print("mLambda")
  #print(mLambda)

  #print("mSigma.inv")
  #print(mSigma.inv)

  #print("diff")
  #print(mLambda - solve(mLambda.inv,tol=1.0E-99))


  #print("mR")
  #print(mR)


  #print(res)
  #bit <- res[(1+d):length(vtheta)]
  
  #print(cbind(bit,dmLambda[Rinds],bit-dmLambda[Rinds]))    

  #ans <- readline()



  # Check with numeric derivative
  #h = 1e-8
  #dmLambda2 = matrix(0, nrow(mLambda), ncol(mLambda))
  #for (i in 1:length(Rinds)) {
  #  mR_high = mR
  #  mR_high[Rinds[i]] = mR_high[Rinds[i]] + h
  #  mR_low = mR
  #  mR_low[Rinds[i]] = mR_low[Rinds[i]] - h
    #mLambda_high = solve(mR_high%*%t(mR_high), tol=1e-99)
    #mLambda_low = solve(mR_low%*%t(mR_low), tol=1e-99)
    #dmLambda2[Rinds[i]] = (-0.5*log(det(mLambda_high)) + f.G_new(vmu,mLambda_high,vy,vr,mC,mSigma.inv,gh)  + 0.5*log(det(mLambda_low)) - f.G_new(vmu,mLambda_low,vy,vr,mC,mSigma.inv,gh))/(2*h)
    #vtheta_high = c(vmu, mR_high[Rinds])
    #vtheta_low = c(vmu, mR_low[Rinds])
    #dmLambda2[Rinds[i]] = (f.GVA_new(vtheta_high,vy,vr,mC,mSigma.inv,gh,mR_high,Rinds,Dinds) - f.GVA_new(vtheta_low,vy,vr,mC,mSigma.inv,gh,mR_low,Rinds,Dinds))/(2*h)
  #  mLambda_high = solve(mR_high%*%t(mR_high), tol=1e-99)
  #  mLambda_low = solve(mR_low%*%t(mR_low), tol=1e-99)
  #  dmLambda2[Rinds[i]] = (f.G_new(vmu,mLambda_high,vy,vr,mC,mSigma.inv,gh) - f.G_new(vmu,mLambda_low,vy,vr,mC,mSigma.inv,gh))/(2*h)
  #}  
  
  #cat("dmLambda", dmLambda, "\n")
  #cat("dmLambda2", dmLambda2, "\n")
  #cat("dmLambda diff", dmLambda - dmLambda2, "\n")
  #for (i in 1:ncol(dmLambda))
  #  for (j in 1:nrow(dmLambda))
  #    cat(j, i, dmLambda[i, j], dmLambda2[i, j], dmLambda[i, j] - dmLambda2[i, j], "\n")
  #cat("\n")
  
  printMatrix = function(mat) {
    for (i in 1:nrow(mat)) {
      for (j in 1:ncol(mat)) {
        cat(mat[i, j], " ")
      }
      cat("\n")
    }
  }
  #printMatrix(dmLambda)
  #printMatrix(dmLambda2)
  #printMatrix(dmLambda-dmLambda2)
  
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
  #vg[(1+d):length(vtheta)] <- dmLambda2[Rinds]    
  
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
  cat("mSigma.inv", mSigma.inv, "\n")
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
