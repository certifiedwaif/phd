# John Ormerod's Poisson regression code, slightly modified to be used by my
# zero-inflated model code.
###############################################################################

source("CalculateB.R")

###############################################################################

f.lap <- function(vnu,vy,vr,mC,mSigma.inv,mLambda) 
{       
    d <- length(vnu)
    veta <- mC%*%vnu
	# FIXME: This is awful, we shouldn't need to invert the whole matrix. There's
	# a lot of room for improvement here.
	mSigma <- solve(mSigma.inv)
	# diag(mC mLambda mC^T)_i i = mC[i, ] . mLambda[, i] . mC[, i] 
    #f <- sum(vy*veta - vr*exp(veta+0.5*diag(mC%*%mLambda%*%t(mC)))) - 0.5*t(vnu)%*%mSigma.inv%*%vnu - 0.5*sum(diag(mLambda%*%mSigma))
	# Are these really equivalent?
	mDiag <- sapply(1:ncol(mC), function(i) sum(mC[i,] * mLambda[, i] * mC[, i]))
	#mDiag2 <-  diag(mC%*%mLambda%*%t(mC))
	#print(mDiag)
	#print(mDiag2)
	f <- sum(vy*veta - vr*exp(veta+0.5*mDiag)) - 0.5*t(vnu)%*%mSigma.inv%*%vnu - 0.5*sum(diag(mLambda%*%mSigma))
    return(f)
}

###############################################################################

norm <- function(v) sqrt(sum(v^2))

###############################################################################

vg.lap <- function(vnu,vy,vr,mC,mSigma.inv,mLambda) 
{       
	#cat("vg.lap vnu", vnu, "\n")
	#cat("vg.lap vy", head(vy), "\n")
	#cat("vg.lap vr", head(vr), "\n")
	#cat("vg.lap mC", head(mC), "\n")
	#mDiag <- sapply(1:ncol(mC), function(i) sum(mC[i,] * mLambda[, i] * mC[, i]))
	#mDiag2 <-  diag(mC%*%mLambda%*%t(mC))
    vg <- t(mC)%*%(vy - vr*exp(mC%*%vnu)) - mSigma.inv%*%vnu
	#cat("vg.lap vg", vg, "norm", norm(vg), "\n")
    return(vg)
}

###############################################################################

mH.lap <- function(vnu,vy,vr,mC,mSigma.inv,mLambda) 
{
    #mDiag <- sapply(1:ncol(mC), function(i) sum(mC[i,] * mLambda[, i] * mC[, i]))
    vw <- exp(mC%*%vnu); dim(vw) <- NULL
    mH <- -t(mC*vw)%*%mC - mSigma.inv
    return(mH)
}

###############################################################################

fit.Lap <- function(vnu,vy,vr,mC,mSigma.inv,mLambda) 
{
    MAXITER <- 100
    
    for (ITER in 1:MAXITER) {
        f  <- f.lap(vnu,vy,vr,mC,mSigma.inv,mLambda)
        vg <- vg.lap(vnu,vy,vr,mC,mSigma.inv,mLambda)
        mH <- mH.lap(vnu,vy,vr,mC,mSigma.inv,mLambda)
        mLambda <- solve(-mH,tol=1.0E-99)
        vnu <- vnu + mLambda%*%vg
        print(c(ITER,f,max(abs(vg))))
        if (max(abs(vg))<1.0E-8) {
            break;
        }
    } 

    #vnu <- vnu
    #controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=1.0e-15)
    #res <- optim(par=vnu, fn=f.lap, gr=vg.lap,
    #    method="BFGS",lower=-Inf, upper=Inf, control=controls,
    #    vy=vy,mC=mC,mSigma.inv) 
        
    #vnu <- res$par
    #mLambda <- solve(-mH.lap(vnu,vy,mC,mSigma.inv),tol=1.0E-99)
    
	# FIXME: Should the mLambda on the next line be something else?
    f <- f.lap(vnu,vy,vr,mC,mSigma.inv,mLambda) + 0.5*log(det(mLambda%*%mSigma.inv))
    #print(f)
    return(list(vnu=vnu,mLambda=mLambda,f=f))
}

###############################################################################

f.G <- function(vnu,mLambda,vy,vr,mC,mSigma.inv,gh) 
{
	#print(vnu)
	#print(mSigma.inv)
    d <- length(vnu)
    
    vnu.til     <- mC%*%vnu
    vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
    vB0 <- B0.fun("POISSON",vnu.til,vsigma2.til,gh) 
    
    f <- sum(vr*(vy*vnu.til - vB0)) - 0.5*t(vnu)%*%mSigma.inv%*%vnu 
    f <- f - 0.5*d*log(2*pi) + 0.5*log(det(mSigma.inv)) 
    return(f)
}

###############################################################################

f.GVA <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
   d <- ncol(mC)  
   vnu <- vtheta[1:d]
   mR[Rinds] <- vtheta[(1+d):length(vtheta)]
   mR[Dinds] <- exp(mR[Dinds]) 
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
    }   
   mLambda <- mR%*%t(mR)   
   
   f <- sum(log(diag(mR))) + f.G(vnu,mLambda,vy,vr,mC,mSigma.inv,gh) 
   f <- f + 0.5*d*log(2*pi) + 0.5*d
   
   if (!is.finite(f)) {
       f <- -1.0E16
   }
   return(f)
}

###############################################################################

vg.G <- function(vnu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
{
    vg <- t(mC)%*%(vr*(vy - vB1)) - mSigma.inv%*%vnu     
	#cat("vg.G", vg, "\n")
    return(vg)
}

###############################################################################

mH.G <- function(vnu,mLambda,vy,vr,mC,mSigma.inv,vB2) 
{
    vw <-  vB2; dim(vw) <- NULL
    mH <- -t(mC*vw)%*%mC - mSigma.inv
    return(mH)    
}

###############################################################################

vg.GVA.approx <- function(vnu,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
    n <- length(vy); P <- length(vnu); d <- ncol(mC)
	eps <- 1.0E-6
	f <- f.GVA(vnu,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	vg.approx <- matrix(0,P,1)
	for (i in 1:P) {
	   vnup <- vnu 
	   vnup[i] <- vnu[i] + eps
	   fp <- f.GVA(vnup,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	   
	   vnum <- vnu 
	   vnum[i] <- vnu[i] - eps
	   fm <- f.GVA(vnum,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	   
	   vg.approx[i] <- (fp - fm)/(2*eps)
	   
	   #vnupp <- vnu 
	   #vnupp[i] <- vnu[i] + 2*eps
	   #fpp <- f.GVA(vnupp,vy,mC,mSigma.inv,gh,mR,Rinds,Dinds)
	   
	   #vnumm <- vnu 
	   #vnumm[i] <- vnu[i] - 2*eps
	   #fmm <- f.GVA(vnumm,vy,mC,mSigma.inv,gh,mR,Rinds,Dinds)  
	   
	   #vg.approx[i] <- (-fpp + 8*fp - 8*fm + fmm)/(12*eps)
	}
	return(vg.approx)
}

###############################################################################

# FIXME: I think there's an error in this routine. The derivatives become
# very large indeed.
vg.GVA <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
	#cat("vtheta", vtheta, "\n")
	#cat("ncol(mC)", ncol(mC), "\n")
    d <- ncol(mC)
    vnu <- vtheta[1:d]
    mR[Rinds] <- vtheta[(1+d):length(vtheta)]
    mR[Dinds] <- exp(mR[Dinds]) 
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
    }    
    mLambda <- mR%*%t(mR)   
  
    vnu.til     <- mC%*%vnu
    vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
    res.B12 <- B12.fun("POISSON",vnu.til,vsigma2.til,gh)
    vB1 <- res.B12$vB1
    vB2 <- res.B12$vB2
    
    vg <- 0*vnu
    vg[1:d] <- vg.G(vnu,mLambda,vy,vr,mC,mSigma.inv,vB1) 

    mLambda.inv <- solve(mLambda,tol=1.0E-99)
    mH <- mH.G(vnu,mLambda,vy,vr,mC,mSigma.inv,vB2)
    dmLambda <- (mLambda.inv + mH)%*%mR

    dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
    # FIXME: Something broken in this function
    vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
   
	cat("vtheta.GVA vtheta", vtheta, "norm", norm(vtheta), "\n")
	cat("vg.GVA vg", vg, "norm", norm(vg), "\n")
    return(vg)
}

###############################################################################

fit.GVA <- function(vnu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-8)
{
	#library(statmod)
	#N <- 15
	#gh  <- gauss.quad(N,kind="hermite")
	gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))
		
    d <- length(vnu)
	#cat("vnu", vnu, "\n")
	#cat("d", d, "\n")
    Dinds <- d*((1:d)-1)+(1:d)
          
    mR <- t(chol(mLambda + diag(1.0E-8,d)))
    mR[Dinds] <- log(mR[Dinds])
    Rinds <- which(lower.tri(mR,diag=TRUE))
    vnu <- c(vnu,mR[Rinds])
    P <- length(vnu)
    
    if (method=="L-BFGS-B") {
        controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
    } else if (method=="Nelder-Mead") {
         controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
    } else {
        controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
    }
    res <- optim(par=vnu, fn=f.GVA, gr=vg.GVA,
        method=method,lower=-Inf, upper=Inf, control=controls,
        vy=vy,vr=vr,mC=mC,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds)        
        
    vtheta <- res$par 
    
    #f <- f.GVA(vnu,vy,mC,CTC,CTy,p,K,mR,Rinds,Dinds)
    
    vnu <- vtheta[1:d]
    mR[Rinds] <- vtheta[(1+d):P]
    mR[Dinds] <- exp(mR[Dinds])  
    mLambda <- mR%*%t(mR)
       
    return(list(res=res,vnu=vnu,mLambda=mLambda))
}

###############################################################################
