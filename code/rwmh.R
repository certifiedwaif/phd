
###############################################################################

library(mvtnorm)

###############################################################################

fast.f.zip <- function(mult, vtheta, vr) 
{
  #cat("fast.f.zip", vtheta, "\n")
  
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
    #cat("fast.f.zip: log.vp", log.vp, "\n")
    # Prior
    vtheta = as.vector(vtheta)
    #print(length(vmu))
    #print(dim(mSigma.beta.inv))
    #print(dim(mSigma.u.inv))
    mSigma.vbeta = solve(mSigma.beta.inv)
    #mSigma.vu = solve(mSigma.u.inv)
    #print(dim(blockDiag(mSigma.vbeta, mSigma.vu)))
    # Break into fixed effects part and random effects part
    #log.vp <- log.vp + dmvnorm(vtheta, mean=rep(0, length(vmu)), sigma=blockDiag(mSigma.vbeta, mSigma.vu), log=TRUE)
    fixed_idx = 1:ncol(mX)
    u_idx = (ncol(mX) + 1):ncol(mC)
    log.vp <- log.vp + dmvnorm(vtheta[fixed_idx], mean=rep(0, ncol(mX)), sigma=mSigma.vbeta, log=TRUE)
    log.vp <- log.vp + dmvnorm(vtheta[u_idx], mean=rep(0, ncol(mZ)), sigma=mSigma.vu, log=TRUE)
    #cat("fast.f.zip: log.vp 2 ", log.vp, "\n")
    
    log.vp
  })
  #cat("fast.f.zip: Returning", log.vp, "\n")
  return(log.vp)
} 

###############################################################################

mcmc <- function(mult, iterations=1e3, burnin=round(iterations/10), thinning=10)
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
    vr = matrix(NA, nrow = length(vy), ncol = ITERATIONS)
	  vr[,1] = rep(1, length(vy))
    vr[zero.set,1]  <- rbinom(nnz,1,0.5) 
    u_idx = (ncol(mX)+1):ncol(mC)
    vnu = matrix(NA, nrow = length(vmu), ncol = ITERATIONS)
    vnu[,1] = rep(0, length(lap_approx$vmu))
	accept = rep(NA, ITERATIONS)
    veta = matrix(NA, nrow = length(vy), ncol = ITERATIONS)
    rho = rep(NA, ITERATIONS)
    sigma2_u = rep(NA, ITERATIONS)
    
    d <- length(mult$vmu)
    mR <- chol(mult$mLambda)
    
    # Iterate
    # TODO: To reduce memory consumption, implement thinning
    for (i in 2:ITERATIONS) {
		ret <- RandomWalkMetropolisHastings(mult, vnu[,i-1],mR,vr[,i-1])
      vnu[,i] <- ret[[1]]
		accept[i] <- ret[[2]]
  		rho[i] <- rbeta(1, prior$a_rho + sum(vr[,i-1]), prior$b_rho + n - sum(vr[,i-1]))
      # FIXME: This is only needed on the zero set vy == 0
  		veta[zero.set,i] <- -exp(mC[zero.set,]%*%as.vector(vnu[,i])) + logit(rho[i])
  		
  		#print(nnz)
  		#print(length(zero.set))
  		
  		val <- c()
  		for (j in 1:nnz) {
  			val[j] <- rbinom(1, 1, expit(veta[zero.set[j],i]))
  		}
  		
  		#print(val)
  		#print(vr[zero.set,i] )
  		
  		#print(n)
  		vr[,i] <- 1
  		vr[zero.set,i] <- val
  		sigma2_u[i] <- 1/rgamma(1, a_sigma + 0.5*m, b_sigma + 0.5*sum(vnu[u_idx, i]^2))
      mult$mSigma.vu = diag(sigma2_u[i], ncol(mZ))
  		#sigma2_u[i] <- 1/rgamma(1, a_sigma + 0.5*m, b_sigma + 0.5*(n-1)*var(vnu[u_idx, i]) + 0.5*tr(mLambda[u_idx, u_idx]))
  		# FIXME: Where do you update mSigma.inv?
  		#cat("\n")
		if ((i %% 1000) == 0)
      cat("i=",i,"\n")
		#cat("\n")
  		
    }
    idx = burnin:iterations
    idx = idx[idx %% thinning == 0]
    result = list(vnu=vnu[,idx], accept=accept[idx], rho=rho[idx], vr=vr[idx],
                  sigma2_u=sigma2_u[idx], vy=vy)

    result
	})
  print(result)  
}

###############################################################################

RandomWalkMetropolisHastings <- function(mult, vtheta, mR, vr)
{
  with(mult, {
    # TODO: Can we avoid doing these matrix inversions on every iteration?
    # Answer: No, you can't, at least for mSigma.u. sigma_u^2 is changing every
    # iteration.
    mSigma.beta = solve(mSigma.beta.inv)
    mSigma.u = solve(mSigma.u.inv)
    #mult.new = mult
    
    d <- length(vtheta)
    #vnu_new = vtheta + rmvnorm(1, mean=0*vtheta, sigma=((2.38^2)/d)*mLambda)
    vnu_new = vtheta + (mR)%*%matrix(sqrt((2.38^2)/d)*rnorm(d))
    #cat("RandomWalkMetropolisHastings: vnu_new", vnu_new, "\n")
    ratio = min(1, exp(fast.f.zip(mult, vnu_new,vr)-fast.f.zip(mult, vtheta,vr)))
    #cat("RandomWalkMetropolisHastings: fast.f.zip(mult, vtheta,vr)", fast.f.zip(mult, vtheta,vr), "\n")
    #cat("RandomWalkMetropolisHastings: fast.f.zip(mult, vnu_new,vr)", fast.f.zip(mult, vnu_new,vr), "\n")
    #cat("RandomWalkMetropolisHastings: vtheta", vtheta, "\n")
    #cat("RandomWalkMetropolisHastings: vnu_new", vnu_new, "\n")
    #cat("RandomWalkMetropolisHastings: ratio", ratio, "\n")
    if (runif(1) < ratio) {
      return(list(vnu_new, 1))
    } else {
      return(list(vtheta, 0))
    }
  })
}

###############################################################################


