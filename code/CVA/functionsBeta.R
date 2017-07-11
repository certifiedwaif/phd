
beta.calculations <- function(vy,mX,inds)
{
	n = length(vy)
	p = ncol(mX)
	q = length(inds)
	
	res = normalize(vy,mX) 
	vy = res$vy
	mX = res$mX	
	XTy = t(mX)%*%vy
	XTX = t(mX)%*%mX
	
	XTX.inv = solve(XTX[inds,inds])
	vbeta.hat = XTX.inv%*%XTy[inds]
	R2 = as.vector(t(XTy[inds])%*%solve(XTX[inds,inds])%*%XTy[inds]/n)
	sigma2.hat = 1 - R2
	
	return(list(XTX.inv=XTX.inv,vbeta.hat=vbeta.hat,R2=R2,sigma2.hat=sigma2.hat))
}

posterior.betaGivenY.MonteCarlo <- function(lbeta,vy,mX,inds,N)
{
	n = length(vy)
	p = ncol(mX)
	q = length(inds)
	
	post = posterior.MC.nonBeta(n,p,R2,N) 
	vu = post$vu
	
	calc = beta.calculations(vy,mX,inds)
	vbeta.hat = calc$vbeta.hat 
	XTX.inv   = calc$XTX.inv
	R2        = calc$R2
	
	mBeta.marg.samples = matrix(0,q,N)
	for (k in 1:N) {
		vmu = vu[k]*vbeta.hat 
		mSigma = n*vu[k]*(1 - vu[k]*R2)*XTX.inv/(n-1)
		for (j in 1:q) {
			mBeta.marg.samples[j,k] = rnorm(1,vmu[j],sqrt(mSigma[j,j]))
		} 
	}
	
	lx <- vector("list", q)
	ly <- vector("list", q)
	
	for (j in 1:q) {
		res = density(mBeta.marg.samples[j,])
		if (is.null(lbeta)) {
			lx[[j]] = res$x
			ly[[j]] = res$y
		} else {
			lx[[j]] = lbeta[[j]]
			y.hat = spline(x=res$x,y=res$y,xout=lbeta[[j]])$y
			y.hat = y.hat/trapint(x,y.hat)	
			ly[[j]] = y.hat
		}
	}
	
	return(list(mBeta.marg.samples=mBeta.marg.samples,lx=lx,ly=ly))
}
	
	
posterior.betaGivenY.RaoBlackwell <- function(lbeta,vy,mX,inds,N)
{
	n = length(vy)
	p = ncol(mX)
	q = length(inds)
	
	post = posterior.MC.nonBeta(n,p,R2,N) 
	vu = post$vu
	
	calc = beta.calculations(vy,mX,inds)
	vbeta.hat = calc$vbeta.hat 
	XTX.inv   = calc$XTX.inv
	R2        = calc$R2
	
	lx <- vector("list", q)
	ly <- vector("list", q)
	
	for (j in 1:q) {
		lx[[j]] = lbeta[[j]]
		ly[[j]] = 0*lbeta[[j]]
		for (k in 1:N) {
			vmu = vu[k]*vbeta.hat 
			mSigma = n*vu[k]*(1 - vu[k]*R2)*XTX.inv/(n-1)
			ly[[j]] = ly[[j]] + (1/N)*dnorm(lx[[j]],vmu[j],sqrt(mSigma[j,j]))
		} 
	}
	
	return(list(lx=lx,ly=ly))
}
	
	
posterior.betaGivenY.approx <- function(lbeta,vy,mX,inds)
{
	n = length(vy)
	p = ncol(mX)
	q = length(inds)
	
	calc = beta.calculations(vy,mX,inds)
	vbeta.hat = calc$vbeta.hat 
	XTX.inv   = calc$XTX.inv
	R2        = calc$R2
	
	a = -0.75
	b = (n - q - 5)/2 - a
	c = (n - 1)/2
	d = q/2 + a
	
	M1 = (b + 1)*hyperg_2F1(d+1, 1, c+1, R2, give=FALSE, strict=TRUE)/c
	
	mu = M1*vbeta.hat
	Sigma = M1*(1-M1*R2)*XTX.inv 
	
	lx <- vector("list", q)
	ly <- vector("list", q)
	for (j in 1:q) {
		lx[[j]] = lbeta[[j]]
		ly[[j]] = dnorm(lx[[j]], mu[j],sqrt(Sigma[j,j]))
	}
	
	return(list(lx=lx,ly=ly))
}
		
posterior.betaGivenY.MLE <- function(lbeta,vy,mX,inds,N)
{
	n = length(vy)
	p = ncol(mX)
	q = length(inds)
	
	calc = beta.calculations(vy,mX,inds)
	vbeta.hat = calc$vbeta.hat 
	XTX.inv   = calc$XTX.inv
	R2        = calc$R2
	
	mu = vbeta.hat 
	Sigma = (1-R2)*XTX.inv
	
	lx <- vector("list", q)
	ly <- vector("list", q)
	
	for (j in 1:q) {
		lx[[j]] = lbeta[[j]]
		ly[[j]] = dnorm(lx[[j]], mu[j],sqrt(Sigma[j,j]))
	}
	
	return(list(lx=lx,ly=ly))
}

 