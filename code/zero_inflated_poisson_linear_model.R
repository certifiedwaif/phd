# John Ormerod's Poisson regression code, slightly modified to be used by my
# zero-inflated model code.
###############################################################################

source("CalculateB.R")

###############################################################################

f.lap <- function(vtheta,vy,vr,mX,mSigma.inv,mLambda) 
{       
    d <- length(vtheta)
    veta <- mX%*%vtheta
    f <- sum(vy*veta - vr*exp(veta+.5*diag(mX%*%mLambda%*%t(mX)))) - 0.5*t(vtheta)%*%mSigma.inv%*%vtheta - 0.5*tr(mLambda%*%Dhat)
    return(f)
}

###############################################################################

vg.lap <- function(vtheta,vy,vr,mX,mSigma.inv,mLambda) 
{       
    vg <- t(mX)%*%(vy - vr*exp(mX%*%vtheta+.5*diag(mX%*%mLambda*t(mX)))) - mSigma.inv%*%vtheta
    return(vg)
}

###############################################################################

mH.lap <- function(vtheta,vy,vr,mX,mSigma.inv,mLambda) 
{
    vw <- exp(mX%*%vtheta+.5*diag(mX%*%mLambda*t(mX))); dim(vw) <- NULL
    mH <- -t(mX*vw)%*%mX - mSigma.inv
    return(mH)
}

###############################################################################

fit.Lap <- function(vmu,vy,vr,mX,mSigma.inv,mLambda) 
{
    MAXITER <- 100
    
    for (ITER in 1:MAXITER) {
        f  <- f.lap(vmu,vy,vr,mX,mSigma.inv,mLambda)
        vg <- vg.lap(vmu,vy,vr,mX,mSigma.inv,mLambda)
        mH <- mH.lap(vmu,vy,vr,mX,mSigma.inv,mLambda)
        mLambda <- solve(-mH,tol=1.0E-99)
        vmu <- vmu + mLambda%*%vg
        #print(c(ITER,f,max(abs(vg))))
        if (max(abs(vg))<1.0E-8) {
            break;
        }
    } 

    #vtheta <- vmu
    #controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=1.0e-15)
    #res <- optim(par=vtheta, fn=f.lap, gr=vg.lap,
    #    method="BFGS",lower=-Inf, upper=Inf, control=controls,
    #    vy=vy,mX=mX,mSigma.inv) 
        
    #vtheta <- res$par
    #mLambda <- solve(-mH.lap(vtheta,vy,mX,mSigma.inv),tol=1.0E-99)
    
    f <- f.lap(vmu,vy,vr,mX,mSigma.inv,mLambda) + 0.5*log(det(mLambda%*%mSigma.inv))
    print(f)
    return(list(vmu=vmu,mLambda=mLambda,f=f))
}

###############################################################################

f.G <- function(vmu,mLambda,vy,mX,mSigma.inv,gh) 
{
    d <- length(vmu)
    
    vmu.til     <- mX%*%vmu
    vsigma2.til <- diag(mX%*%mLambda%*%t(mX))
    vB0 <- B0.fun("POISSON",vmu.til,vsigma2.til,gh) 
    
    f <- sum(vy*vmu.til - vB0) - 0.5*t(vmu)%*%mSigma.inv%*%vmu 
    f <- f - 0.5*d*log(2*pi) + 0.5*log(det(mSigma.inv)) 
    return(f)
}

###############################################################################

f.GVA <- function(vtheta,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
{
   d <- ncol(mX)  
   vmu <- vtheta[1:d]
   mR[Rinds] <- vtheta[(1+d):length(vtheta)]
   mR[Dinds] <- exp(mR[Dinds]) 
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
    }   
   mLambda <- mR%*%t(mR)   
   
   f <- sum(log(diag(mR))) + f.G(vmu,mLambda,vy,mX,mSigma.inv,gh) 
   f <- f + 0.5*d*log(2*pi) + 0.5*d
   
   if (!is.finite(f)) {
       f <- -1.0E16
   }
   return(f)
}

###############################################################################

vg.G <- function(vmu,mLambda,vy,mX,mSigma.inv,vB1) 
{
    vg <- t(mX)%*%(vy - vB1) - mSigma.inv%*%vmu     
    return(vg)
}

###############################################################################

mH.G <- function(vmu,mLambda,vy,mX,mSigma.inv,vB2) 
{
    vw <-  vB2; dim(vw) <- NULL
    mH <- -t(mX*vw)%*%mX - mSigma.inv
    return(mH)    
}

###############################################################################

vg.GVA.approx <- function(vtheta,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
{
    n <- length(vy); P <- length(vtheta); d <- ncol(mX)
		eps <- 1.0E-6
		f <- f.GVA(vtheta,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
		vg.approx <- matrix(0,P,1)
		for (i in 1:P) {
		   vthetap <- vtheta 
		   vthetap[i] <- vtheta[i] + eps
		   fp <- f.GVA(vthetap,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
		   
		   vthetam <- vtheta 
		   vthetam[i] <- vtheta[i] - eps
		   fm <- f.GVA(vthetam,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
		   
		   vg.approx[i] <- (fp - fm)/(2*eps)
		   
		   #vthetapp <- vtheta 
		   #vthetapp[i] <- vtheta[i] + 2*eps
		   #fpp <- f.GVA(vthetapp,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
		   
		   #vthetamm <- vtheta 
		   #vthetamm[i] <- vtheta[i] - 2*eps
		   #fmm <- f.GVA(vthetamm,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)  
		   
		   #vg.approx[i] <- (-fpp + 8*fp - 8*fm + fmm)/(12*eps)
		}
		return(vg.approx)
}

###############################################################################

vg.GVA <- function(vtheta,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
{
    d <- length(vmu)
    vmu <- vtheta[1:d]
    mR[Rinds] <- vtheta[(1+d):length(vtheta)]
    mR[Dinds] <- exp(mR[Dinds]) 
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
    }    
    mLambda <- mR%*%t(mR)   
  
    vmu.til     <- mX%*%vmu
    vsigma2.til <- diag(mX%*%mLambda%*%t(mX))
    res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
    vB1 <- res.B12$vB1
    vB2 <- res.B12$vB2
    
    vg <- 0*vtheta
    vg[1:d] <- vg.G(vmu,mLambda,vy,mX,mSigma.inv,vB1) 

    mLambda.inv <- solve(mLambda,tol=1.0E-99)
    mH <- mH.G(vmu,mLambda,vy,mX,mSigma.inv,vB2)
    dmLambda <- (mLambda.inv + mH)%*%mR

    #count <- d+1
    #for (j in 1:d) {
    #    for (i in j:d) {
    #       if (i==j) {
    #           vg[count] <- dmLambda[i,i]*mR[i,i]
    #       } else {
    #           vg[count] <- dmLambda[i,j]
    #       }
    #       count <- count + 1
    #    }      
    #}
    
    dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
    vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
   
    #vg.approx <- vg.GVA.approx(vtheta,vy,mX,mSigma.inv,gh,mR,Rinds,Dinds)
    #print(cbind(vg,vg.approx,abs(vg - vg.approx),abs(vg - vg.approx)/abs(vg),vg/vg.approx))
    #ans <- readline()            
              
    return(vg)
}

###############################################################################

fit.GVA <- function(vmu,mLambda,vy,mX,mSigma.inv,method,reltol=1.0e-8)
{
		#library(statmod)
		#N <- 15
		#gh  <- gauss.quad(N,kind="hermite")
		gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))    
		
    d <- length(vmu)
    Dinds <- d*((1:d)-1)+(1:d) 		
          
    mR <- t(chol(mLambda + diag(1.0E-8,d)))
    mR[Dinds] <- log(mR[Dinds]) 
    Rinds <- which(lower.tri(mR,diag=TRUE))
    vtheta <- c(vmu,mR[Rinds])
    P <- length(vtheta)
    
    if (method=="L-BFGS-B") {
        controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
    } else if (method=="Nelder-Mead") {
         controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
    } else {
        controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
    }
    res <- optim(par=vtheta, fn=f.GVA, gr=vg.GVA,
        method=method,lower=-Inf, upper=Inf, control=controls,
        vy=vy,mX=mX,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds)        
        
    vtheta <- res$par 
    
    #f <- f.GVA(vtheta,vy,mC,CTC,CTy,p,K,mR,Rinds,Dinds)
    
    vmu <- vtheta[1:d]
    mR[Rinds] <- vtheta[(1+d):P]
    mR[Dinds] <- exp(mR[Dinds])  
    mLambda <- mR%*%t(mR)
       
    return(list(res=res,vmu=vmu,mLambda=mLambda))
}

###############################################################################
