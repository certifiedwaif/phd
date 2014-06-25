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

ImportanceSample <- function(mult)
{
	with(mult, {
    # FIXME: WTF should nu be?
    nu <- ncol(mX) + ncol(mZ)
    N <- nrow(mX)
		mZ <- rmvt(N, delta=as.vector(vmu), sigma=mLambda, df = nu, type = "shifted")
		log.vq <- dmvt(mZ, delta=as.vector(vmu), sigma=mLambda, df = nu, log=TRUE, type = "shifted")
		log.vp <- as.vector(fast.f.zip(mult)) 
		#log.vq <- c(log.vq,log.vq.new)
		#log.vp <- c(log.vp,log.vp.new)				
		vw <- exp(log.vp.new - log.vq.new)
		N <- length(vw)
		I.hat <- mean(vw)
		se <- sqrt(var(vw)/length(vw))
		#print(c(ITER,N,I.hat,se))
		#if (se<EPS) {
		#  break;
		#}
    list(I.hat=I.hat, se=se)
	})
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

fast.f.zip <- function(mult, vtheta,vr) 
{
  cat("fast.f.zip", vtheta, "\n")
  
  #mult$vr <- vr
  
  log.vp <- with(mult,{
    
    #N <- nrow(mTheta)
    n <- length(vy)
    #mY <- matrix(vy,n,N)
    
    #vy*pnorm(veta, log.p=TRUE) + (1-vy)*pnorm(veta, lower.tail=FALSE, log.p=TRUE)
    #log.vp <- matrix(1,1,n)%*%(vy*pnorm(vEta, log.p=TRUE) + (1-vy)*pnorm(vEta, lower.tail=FALSE, log.p=TRUE)) + dmvnorm(mTheta,sigma=mSigma,log=TRUE)
  
    # Log-likelihood pertaining to vbeta and vu
    veta = mC%*%as.vector(vtheta)
    log.vp <- t(vy*vr)%*%veta - t(vr)%*%exp(veta) - sum(lgamma(vy+1))
    cat("fast.f.zip: log.vp", log.vp, "\n")
    # Prior
    vtheta = as.vector(vtheta)
    #print(length(vmu))
    #print(dim(mSigma.beta.inv))
    #print(dim(mSigma.u.inv))
    mSigma.vbeta = solve(mSigma.beta.inv)
    mSigma.vu = solve(mSigma.u.inv)
    #print(dim(blockDiag(mSigma.vbeta, mSigma.vu)))
    log.vp <- log.vp + dmvnorm(vtheta, mean=rep(0, length(vmu)), sigma=blockDiag(mSigma.vbeta, mSigma.vu), log=TRUE)
    cat("fast.f.zip: log.vp 2 ", log.vp, "\n")
    
    log.vp
  })
  cat("fast.f.zip: Returning", log.vp, "\n")
  return(log.vp)
} 

###############################################################################


my.diag <- function(mA) {
    k <- ncol(mA) 
    Dinds <- k*((1:k)-1) + (1:k)
    val <- mA[Dinds]
    return(val)
}

my.tr <- function(mA) {
	return( sum(my.diag(mA)) ) 
}

mcmc <- function(mult, iterations=1e3)
{
	# Initialise with Laplacian approximation
  #print(str(mult))
	lap_approx = with(mult, {
    mSigma.inv = blockDiag(mSigma.beta.inv, mSigma.u.inv)
    fit.Lap(vmu, vy, vp, mC, mSigma.inv, mLambda)
  })
  print(str(lap_approx))
  mult$vmu = lap_approx$vmu
  mult$mLambda = lap_approx$mLambda
  result = with(mult,
	{
    # Initialise and set constants
	  ITERATIONS = iterations
	  n = length(vy)
	  m = ncol(mZ)
	  zero.set = which(vy == 0)
	  nnz <- length(zero.set)
    vr = matrix(1, nrow = length(vy), ncol = ITERATIONS)
    vr[!zero.set,] <- 1
    vr[zero.set,]  <- rbinom(nnz,1,0.5) 
    vr[,1] = rep(1, length(vy))
    u_idx = (ncol(mX)+1):ncol(mC)
    vnu = matrix(NA, nrow = length(vmu), ncol = ITERATIONS)
    vnu[,1] = rep(0, length(lap_approx$vmu))
    veta = matrix(NA, nrow = length(vy), ncol = ITERATIONS)
    rho = rep(NA, ITERATIONS)
    sigma2_u = rep(NA, ITERATIONS)
    
    d <- length(mult$vmu)
    mR <- chol(((2.38^2)/d)*mult$mLambda)
    
    # Iterate
    # TODO: To reduce memory consumption, implement thinning
    for (i in 2:ITERATIONS) {
      vnu[,i] <- RandomWalkMetropolisHastings(mult, vnu[,i-1],mR,vr[,i])
  		rho[i] <- rbeta(1, a_rho + sum(vr[,i]), b_rho + n - sum(vr[,i]))
      # FIXME: This is only needed on the zero set vy == 0
  		veta[zero.set,i] <- -exp(mC[zero.set,]%*%as.vector(vnu[,i])) + logit(rho[i])
  		
  		print(nnz)
  		print(length(zero.set))
  		
  		val <- c()
  		for (j in 1:nnz) {
  			val[j] <- rbinom(1, 1, expit(veta[zero.set[j],i]))
  		}
  		
  		print(val)
  		print(vr[zero.set,i] )
  		
  		print(n)
  		
  		vr[zero.set,i] <- val
  		sigma2_u[i] <- 1/rgamma(1, a_sigma + 0.5*m, b_sigma + 0.5*sum(vnu[u_idx, i]^2))
  		#sigma2_u[i] <- 1/rgamma(1, a_sigma + 0.5*m, b_sigma + 0.5*(n-1)*var(vnu[u_idx, i]) + 0.5*tr(mLambda[u_idx, u_idx]))
  		
  		cat("\n")
		cat("i=",i,"\n")
		cat("\n")
  		
    }
    result = list(vnu=vnu, rho=rho, vr=vr, sigma2_u=sigma2_u, vy=vy)

    result
	})
  print(result)  
}

###############################################################################

RandomWalkMetropolisHastings <- function(mult, vtheta,mR,vr)
{
  with(mult, {
    mSigma.beta = solve(mSigma.beta.inv)
    mSigma.u = solve(mSigma.u.inv)
    #mult.new = mult
    
    d <- length(vtheta)
    #vnu_new = vtheta + rmvnorm(1, mean=0*vtheta, sigma=((2.38^2)/d)*mLambda)
    vnu_new = vtheta + t(mR)%*%matrix(rnorm(d))
    cat("RandomWalkMetropolisHastings: vnu_new", vnu_new, "\n")
    ratio = min(1, exp(fast.f.zip(mult, vnu_new,vr)-fast.f.zip(mult, vtheta,vr)))
    cat("RandomWalkMetropolisHastings: fast.f.zip(mult, vtheta,vr)", fast.f.zip(mult, vtheta,vr), "\n")
    cat("RandomWalkMetropolisHastings: fast.f.zip(mult, vnu_new,vr)", fast.f.zip(mult, vnu_new,vr), "\n")
    cat("RandomWalkMetropolisHastings: vtheta", vtheta, "\n")
    cat("RandomWalkMetropolisHastings: vnu_new", vnu_new, "\n")
    cat("RandomWalkMetropolisHastings: ratio", ratio, "\n")
    if (runif(1) < ratio) {
      vnu_new
    } else {
      vtheta
    }
  })
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

rwmh_main <- function()
{
  u <- seq(-3,0.5,,100)
  v <- seq(-0.7,0.1,,100)
  
  f <- djointLMM(u,v,vy,mX,mSigma.inv)
  levels <- pretty(range(f),12)
  
  contour(u,v,f,levels=levels,lwd=2,main="Exact Posterior Density",cex.main=3)
  
  ans <- readline()
  
  res.lap <- LaplaceApproxPosterior(vy,mX,mSigma.inv) 
  
  vmu     <- res.lap$vtheta
  mLambda <- res.lap$mLambda
  nu      <- 8
  EPS     <- 1.0E-12
  N       <- 1000
  
  res.is <- ImportanceSampling(N,vmu,mLambda,nu,vy,mX,mSigma,EPS) 
  
  plot(density(res.is$vw))
  
  res.nis <- NormalisedImportanceSampling(100*length(res.is$vw),vmu,mLambda,nu,vy,mX,mSigma) 
}
