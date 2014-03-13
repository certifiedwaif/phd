###############################################################################

source("CalculateB.R")
require(Matrix)

###############################################################################

f.lap <- function(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv) 
{       
    vtheta <- c(vbeta, vu)
    d <- length(vtheta)
    veta <- cbind(mX, mZ)%*%vtheta
    f <- sum(vy*veta - exp(veta)) - 0.5*(t(vbeta)%*%mSigmaBeta.inv%*%vbeta + t(vu)%*%mSigma.inv%*%vu)
    return(f)
}

###############################################################################

vg.lap <- function(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv) 
{       
  vtheta <- c(vbeta, vu)
  vg <- t(cbind(mX, mZ))%*%(vy - exp(cbind(mX, mZ)%*%vtheta)) - rbind(mSigmaBeta.inv%*%vbeta, mSigma.inv%*%vu)
  return(vg)
}

###############################################################################

mH.lap <- function(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv) 
{
  mW <- diag(as.vector(exp(mX%*%vbeta + mZ%*%vu)))
  mH <- rbind(cbind(-t(mX)%*%mW%*%mX - mSigmaBeta.inv, -t(mX)%*%mW%*%mZ),
              cbind(-t(mZ)%*%mW%*%mX, -t(mZ)%*%mW%*%mZ - mSigma.inv))
  return(mH)
}

###############################################################################

fit.Lap <- function(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv, debug=FALSE) 
{
    MAXITER <- 100
    vmu = c(vbeta, vu)
    
    for (ITER in 1:MAXITER) {
        f  <- f.lap(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv) 
        vg <- vg.lap(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv) 
        mH <- mH.lap(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv)
        mLambda <- solve(-mH,tol=1.0E-99)
        vmu <- vmu + mLambda%*%vg
        vbeta <- vmu[1:2]
        vu <- vmu[3:4]
        if (debug)
            print(c(ITER,f,max(abs(vg))))
        if (max(abs(vg))<1.0E-8) {
            break;
        }
    } 

    #vu <- vmu
    #controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=1.0e-15)
    #res <- optim(par=vu, fn=f.lap, gr=vg.lap,
    #    method="BFGS",lower=-Inf, upper=Inf, control=controls,
    #    vy=vy,mZ=mZ,mSigma.inv) 
        
    #vu <- res$par
    #mLambda <- solve(-mH.lap(vu,vy,mZ,mSigma.inv),tol=1.0E-99)

	# FIXME: The term at the end involving determinants is probably wrong.
    f <- f.lap(vbeta,vu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv) + 0.5*log(det(mLambda%*%bdiag(mSigmaBeta.inv, mSigma.inv)))
    if (debug)
		print(f)
    return(list(vmu=vmu,mLambda=mLambda,f=f))
}

###############################################################################

f.G <- function(vmu,mLambda,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh) 
# Calculate the main part of the lower bound.
{
    d1 <- ncol(mX)
	d2 <- ncol(mZ)
	vbeta <- vmu[1:d1]
	vu <- vmu[(d1+1):(d1+d2)]
    
    vmu.til     <- mX%*%vbeta+mZ%*%vu
    vsigma2.til <- diag(cbind(mX,mZ)%*%mLambda%*%t(cbind(mX,mZ)))
    vB0 <- B0.fun("POISSON",vmu.til,vsigma2.til,gh) 
    
    f <- sum(vy*vmu.til - vB0) - 0.5*(t(vbeta)%*%mSigmaBeta.inv%*%vbeta + t(vu)%*%mSigma.inv%*%vu)
    f <- f - 0.5*(d1+d2)*log(2*pi) + 0.5*(log(det(mSigmaBeta.inv)) + log(det(mSigma.inv)) )
    return(f)
}

###############################################################################

f.GVA <- function(vmu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh,mR,Rinds,Dinds)
# Calculate lower bound.
{
	# Unpack vbeta, vu and mR from vu
	d1 <- ncol(mX)
   d2 <- ncol(mZ)  
	vbeta = vmu[1:d1]
   vu <- vmu[(d1+1):d2]
   mR[Rinds] <- vmu[(1+d1+d2):length(vmu)]
   mR[Dinds] <- exp(mR[Dinds]) 
	# Threshold entries in the diagonal of mR at 10,000
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
    }   
   mLambda <- mR%*%t(mR)   
   
	# Calculate f
   f <- sum(log(diag(mR))) + f.G(vmu,mLambda,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh) 
   f <- f + 0.5*(d1+d2)*log(2*pi) + 0.5*(d1+d2)
   
   if (!is.finite(f)) {
       f <- -1.0E16
   }
   return(f)
}

###############################################################################

vg.G <- function(vmu,mLambda,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,vB1) 
{
    vg <- t(cbind(mX,mZ))%*%(vy - vB1) - mSigmaBeta.inv%*%vbeta - mSigma.inv%*%vu     
    return(vg)
}

###############################################################################

mH.G <- function(vmu,mLambda,vy,mZ,mSigma.inv,vB2) 
{
    vw <-  vB2; dim(vw) <- NULL
    mH <- -t(cbind(mX,mZ)*vw)%*%cbind(mX, mZ) - rbind(mSigmaBeta.inv, mSigma.inv)
    return(mH)    
}

###############################################################################

vg.GVA.approx <- function(vmu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh,mR,Rinds,Dinds)
{
    n <- length(vy); P <- length(vu); d <- ncol(mZ)
		# Question: Why do we do this epsilon averaging trick?
		eps <- 1.0E-6
		f <- f.GVA(vmu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh,mR,Rinds,Dinds)
		vg.approx <- matrix(0,P,1)
		for (i in 1:P) {
		   vup <- vu 
		   vup[i] <- vu[i] + eps
		   fp <- f.GVA(vup,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh,mR,Rinds,Dinds)
		   
		   vum <- vu 
		   vum[i] <- vu[i] - eps
		   fm <- f.GVA(vum,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh,mR,Rinds,Dinds)
		   
		   vg.approx[i] <- (fp - fm)/(2*eps)
		   
		   #vthetapp <- vtheta 
		   #vthetapp[i] <- vtheta[i] + 2*eps
		   #fpp <- f.GVA(vthetapp,vy,mZ,mSigma.inv,gh,mR,Rinds,Dinds)
		   
		   #vthetamm <- vtheta 
		   #vthetamm[i] <- vtheta[i] - 2*eps
		   #fmm <- f.GVA(vthetamm,vy,mZ,mSigma.inv,gh,mR,Rinds,Dinds)  
		   
		   #vg.approx[i] <- (-fpp + 8*fp - 8*fm + fmm)/(12*eps)
		}
		return(vg.approx)
}

###############################################################################

vg.GVA <- function(vmu,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,gh,mR,Rinds,Dinds)
{
	# Unpack vbeta, vu and mR from vmu
	d1 <- ncol(mX)
    d2 <- ncol(mZ)
	vbeta <- vmu[1:d1]
    vu <- vmu[(d1+1):(d1+d2)]
    mR[Rinds] <- vmu[(1+d1+d2):length(vmu)]

    mR[Dinds] <- exp(mR[Dinds]) 
	# Threshold entries of mR at 1,000
    for (i in 1:length(Dinds)) {
        mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
    }    
    mLambda <- mR%*%t(mR)   
  
    vmu.til     <- mZ%*%vmu
    vsigma2.til <- diag(mZ%*%mLambda%*%t(mZ))
    res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
    vB1 <- res.B12$vB1
    vB2 <- res.B12$vB2
    
	# Set the derivative of vmu in the first d entries of vg
    vg <- 0*vu
    vg[1:(d1+d2)] <- vg.G(vmu,mLambda,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,vB1) 

    mLambda.inv <- solve(mLambda,tol=1.0E-99)
    mH <- mH.G(vmu,mLambda,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,vB2)
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
    
	# Set the derivative of mR in the last entries of vg
    dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
    vg[(1+d1+d2):length(vu)] <- dmLambda[Rinds]    
   
    #vg.approx <- vg.GVA.approx(vu,vy,mZ,mSigma.inv,gh,mR,Rinds,Dinds)
    #print(cbind(vg,vg.approx,abs(vg - vg.approx),abs(vg - vg.approx)/abs(vg),vg/vg.approx))
    #ans <- readline()            
              
    return(vg)
}

###############################################################################

fit.GVA <- function(vbeta,vu,mLambda,vy,mX,mZ,mSigmaBeta.inv,mSigma.inv,method,reltol=1.0e-8)
{
	#library(statmod)
	#N <- 15
	#gh  <- gauss.quad(N,kind="hermite")
	gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))    
	vmu <- c(vbeta, vu)
    d <- length(vmu)
    # Calculate the diagonal indices
    Dinds <- d*((1:d)-1)+(1:d) 		
          
    mR <- t(chol(mLambda + diag(1.0E-8,d)))
    # Question: Why do we only take the log of the diagonal entries?
    mR[Dinds] <- log(mR[Dinds]) 
    Rinds <- which(lower.tri(mR,diag=TRUE))
    vu <- c(vmu,mR[Rinds])
    P <- length(vu)
    
    if (method=="L-BFGS-B") {
        controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
    } else if (method=="Nelder-Mead") {
         controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
    } else {
        controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
    }
    res <- optim(par=vu, fn=f.GVA, gr=vg.GVA,
        method=method,lower=-Inf, upper=Inf, control=controls,
        vy=vy,mX=mX,mZ=mZ,mSigmaBeta.inv=mSigmaBeta.inv,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,
		Rinds=Rinds,Dinds=Dinds)        
        
    vu <- res$par 
    
    #f <- f.GVA(vu,vy,mC,CTC,CTy,p,K,mR,Rinds,Dinds)
    
    vmu <- vu[1:d]
    mR[Rinds] <- vu[(1+d):P]
    mR[Dinds] <- exp(mR[Dinds])  
    mLambda <- mR%*%t(mR)
       
    return(list(res=res,vmu=vmu,mLambda=mLambda))
}

###############################################################################
