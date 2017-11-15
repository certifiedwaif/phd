
library(correlation)
library(gsl)
library(mvtnorm)
library(latex2exp)

source("functionsMix.Rs")
source("functionsPlot.Rs")
source("functionsBeta.Rs")

####################################################################################################

trapint <- function(xgrid, fgrid) 
{
	ng <- length(xgrid)
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
	integ <- sum(xvec * fvec)/2
	return(integ)
}

####################################################################################################

caclulate.constants <- function(n,p) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a

	return(list(a=a,b=b,c=c,d=d))
}

####################################################################################################

posterior.MC.nonBeta <- function(n,p,R2,N) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d
 
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)
	vu = vg/(1 + vg)
	vsigma2 = 1/rgamma(N,c,0.5*n*(1 - R2*vu))
	valpha = rnorm(N,0,sqrt(vsigma2/n))
	
    return(list(vg=vg,vu=vu,vsigma2=vsigma2,valpha=valpha))		
}

####################################################################################################

posterior.gGivenY.MonteCarlo <- function(g,n,p,R2,N) 
{
	post = posterior.MC.nonBeta(n,p,R2,N) 
	vg = post$vg
	
	den = density(vg)

	y = den$y[den$x>0]
	x = den$x[den$x>0]
	y = y/trapint(x,y)	
	
	if (is.null(g)) {
		g     = x
		y.hat = y
	} else {
		#y[y==0] = min(y[y>0]) 
		#log.y = log(y)
		#log.y.hat = spline(x=x,y=log.y,xout=g)$y
		#y.hat = exp(log.y.hat)
		y.hat = spline(x=x,y=y,xout=g,method="natural")$y
		y.hat[y.hat<0] = 0
		y.hat = y.hat/trapint(g,y.hat)	
		#plot(x,y,type="l")
		#lines(g,y.hat,col="red")
		#ans <- readline()
	}
	
	return(list(samples=vg,x=g,y=y.hat))
}

posterior.gGivenY.exact <- function(g,n,p,R2) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d
	
	if (is.null(g)) {
		res = posterior.gGivenY.MonteCarlo(g,n,p,R2,N=1000)
		g = res$x
	}
	
	log.f = 0
	log.f = log.f + (b + 1)*log(1 - R2) 
	log.f = log.f + b*log(g) 
	log.f = log.f - (d + b + 2)*log(1 + g*(1 - R2)) 
	log.f = log.f - lbeta(d+1,b+1)
	
	return(list(x=g,y=exp(log.f)))
}

####################################################################################################

logit = function(x) {
	return( log(x/(1-x)) )
}

expit = function(x) {
	return( 1/(1+exp(-x)) )
}

posterior.uGivenY.MonteCarlo <- function(u,n,p,R2,N) 
{
	post = posterior.MC.nonBeta(n,p,R2,N) 
	vu = post$vu
	den = density(vu)
	
	cond = (den$x>=0)&(den$x<1)
	y = den$y[cond]
	x = den$x[cond]
	y = y/trapint(x,y)	
	
	if (is.null(u)) {
		u     = x
		y.hat = y
	} else {
		log.y = log(1+y)
		log.y.hat = spline(x=x,y=log.y,xout=u,method="natural")$y
		
		#print(u)
		
		#par(mfrow=c(1,2))
		#plot(x,log.y,type="l",col="black",xlim=range(u))
		#lines(u,log.y.hat,col="red")
		
		y.hat = exp(log.y.hat) - 1
		y.hat[y.hat<0] = 0
		y.hat = y.hat/trapint(u,y.hat)	
		
		#plot(x,y,type="l",col="black",xlim=range(u))
		#lines(u,y.hat,col="red")
		
		#print(y.hat)
		#ans <- readline()
	}
	
	return(list(samples=vu,x=u,y=y.hat))
}


posterior.uGivenY.exact <- function(u,n,p,R2) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d
	
	#if (is.null(u)) {
	#	res = posterior.uGivenY.MonteCarlo(u,n,p,R2,N=1000)
	#	u = res$x
	#}
	
	log.f = 0
	log.f = log.f +  (b + 1)*log(1 - R2) 
	log.f = log.f - lbeta(d+1,b+1) 
	log.f = log.f + b*log(u) 
	log.f = log.f + d*log(1 - u) 
	log.f = log.f - (d + b + 2)*log(1 - u*R2)
	
	return(list(x=u,y=exp(log.f)))
}

####################################################################################################

posterior.sigma2GivenY.MonteCarlo <- function(sigma2,n,p,R2,N) 
{
	post = posterior.MC.nonBeta(n,p,R2,N) 
	vsigma2 = post$vsigma2
	den = density(vsigma2)
	
	y = den$y[den$x>0]
	x = den$x[den$x>0]
	y = y/trapint(x,y)	
	
	if (is.null(sigma2)) {
		sigma2 = x
		y.hat  = y
	} else {
		#y[y==0] = min(y[y>0]) 
		#log.y = log(y)
		#log.y.hat = spline(x=x,y=log.y,xout=sigma2)$y
		#y.hat = exp(log.y.hat)
		y.hat = spline(x=x,y=y,xout=sigma2,method="natural")$y
		y.hat[y.hat<0] = 0
		y.hat = y.hat/trapint(sigma2,y.hat)	
		#plot(x,y,type="l")
		#lines(sigma2,y.hat,col="red")
		#ans <- readline()
	}
	
	return(list(samples=vsigma2,x=sigma2,y=y.hat))
}


posterior.sigma2GivenY.exact <- function(sigma2,n,p,R2) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d

	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(sigma2,n,p,R2,N=1000)
		sigma2 = res$x
	}

	log.f = c*log(n/2)
	log.f = log.f + (b + 1)*log(1 - R2)
	log.f = log.f - lgamma(c)
	log.f = log.f - (c + 1)*log(sigma2)
	log.f = log.f - 0.5*n*(1-R2)/sigma2
	log.f = log.f + log( hyperg_1F1(d + 1, c, -0.5*n*R2/sigma2, give=FALSE, strict=TRUE) )

	return(list(x=sigma2,y=exp(log.f)))
}

dinvgamma <- function(x,a,b) {
	log.f = a*log(b) - lgamma(a) - (a + 1)*log(x) - b/x
	return(exp(log.f))
}

posterior.sigma2GivenY.RaoBlackwell <- function(sigma2,n,p,R2,N) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d
	
	post = posterior.MC.nonBeta(n,p,R2,N)
	vu = post$vu
	
	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(sigma2,n,p,R2,N=1000)
		sigma2 = res$x
	}

	vs = rep(c,N)
	vt = 0.5*n*(1 - R2*vu)

	y = 0*sigma2
	for (i in 1:N) {
		y = y + (1/N)*dinvgamma(sigma2,vs[i],vt[i])
	}
	
	return(list(vs=vs,vt=vt,x=sigma2,y=y))
}

posterior.sigma2GivenY.approx <- function(sigma2,n,p,R2) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d

	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(sigma2,n,p,R2,N=1000)
		sigma2 = res$x
	}
	
	M1 = (b + 1)*hyperg_2F1(d+1, 1, c+1, R2, give=FALSE, strict=TRUE)/c
	
	s = c
	t = 0.5*n*(1 - M1*R2)
	y = dinvgamma(sigma2,s,t)
	
	return(list(x=sigma2,y=y,s=s,t=t))
}

####################################################################################################

posterior.alphaGivenY.MonteCarlo <- function(alpha,n,p,R2,N) 
{
	post = posterior.MC.nonBeta(n,p,R2,N) 
	valpha = post$valpha
	den = density(valpha)
	
	if (is.null(alpha)) {
		alpha = den$x
		y.hat = den$y
	} else {
		y.hat = spline(x=den$x,y=den$y,xout=alpha)$y
		y.hat = y.hat/trapint(x,y.hat)	
	}
	
	return(list(samples=valpha,x=alpha,y=y.hat))
}

log.f1.fun <- function(cval,b,R2,d,n,alpha) 
{
	log.f1 = 0
	log.f1 = log.f1 + lgamma(cval+0.5)   
	log.f1 = log.f1 + (b+1)*log(1 - R2)
	log.f1 = log.f1 - lgamma(cval)
	log.f1 = log.f1 - 0.5*log(pi)
	log.f1 = log.f1 - 0.5*n*log(1 + alpha^2)
	log.f1 = log.f1 + log( hyperg_2F1(b+1, cval+0.5, cval, R2/(1 + alpha^2), give=FALSE, strict=TRUE))
	
	return(log.f1)
}

log.f2.fun <- function(cval,b,R2,d,n,alpha) 
{
	log.f2 = 0
	log.f2 = log.f2 + lgamma(cval + 0.5)   
	log.f2 = log.f2 - lgamma(cval)
	log.f2 = log.f2 - 0.5*log(pi)
	log.f2 = log.f2 + (b+1)*log(1 - R2)
	log.f2 = log.f2 - (b+3/2)*log(alpha^2 + 1 - R2)
	log.f2 = log.f2 - (d+1)*log(1+alpha^2)
	log.f2 = log.f2 + log( hyperg_2F1(d+1, -0.5, cval, R2/(1 + alpha^2), give=FALSE, strict=TRUE))
	
	return(log.f2)
}

posterior.alphaGivenY.exact <- function(alpha,n,p,R2) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d
		
	if (is.null(alpha)) {
		res = posterior.alphaGivenY.MonteCarlo(alpha,n,p,R2,N=1000)
		alpha = res$x
	}
	
	if (n<=100) 
	{
		log.f = log.f1.fun(c,b,R2,d,n,alpha) 
		if (any(is.nan(log.f))) {
			log.f = log.f2.fun(c,b,R2,d,n,alpha) 
		}
	} 
	else 
	{ 
		log.f = log.f2.fun(c,b,R2,d,n,alpha) 
		if (any(is.nan(log.f))) {
			log.f = log.f1.fun(c,b,R2,d,n,alpha) 
		}
	}
	
	#print(log.f1)
	#print(log.f2)
	#ans <- readline()
	#log.f = ifelse(is.nan(log.f1),log.f2,log.f1)

	return(list(x=alpha,y=exp(log.f)))
}


posterior.alphaGivenY.RaoBlackwell <- function(alpha,n,p,R2,N) 
{
	post = posterior.MC.nonBeta(n,p,R2,N) 
	vsigma2 = post$vsigma2
	
	if (is.null(alpha)) {
		res = posterior.sigma2GivenY.MonteCarlo(alpha,n,p,R2,N=1000)
		alpha = res$x
	}
			
	vmu    = rep(0,N)
	vsigma = sqrt(vsigma2/n)
	
	y = 0*alpha
	for (i in 1:N) {
		y = y + (1/N)*dnorm(alpha, vmu[i], vsigma[i]) 
	}

	return(list(mu=vmu,sigma=vsigma,x=alpha,y=y))
}


posterior.alphaGivenY.approx <- function(alpha,n,p,R2,N) 
{
	cons = caclulate.constants(n,p)
	a = cons$a; b=cons$b; c=cons$c; d=cons$d
	
	if (is.null(alpha)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		alpha = res$x
	}
	
	M1 = (b + 1)*hyperg_2F1(d+1, 1, c+1, R2, give=FALSE, strict=TRUE)/c
	
	y = dnorm(alpha,0,sqrt((1 - M1*R2)/(n-1)))
	
	return(list(x=alpha,y=y))
}

####################################################################################################

normalize <- function(y,X) 
{	
	n = length(y)
	p = ncol(X)
	
	mu.y     = mean(y)
	sigma2.y = (n-1)*var(y)/n
	vy       = (y - mu.y)/sqrt(sigma2.y)

	# Normalise covariates
	mX <- matrix(0,n,p)
	mu.x <- c()
	sigma2.x <- c()
	for (j in 1:p) 
	{
		mu.x[j]      <- mean(X[,j])
		sigma2.x[j]  <- (n-1)*var(X[,j])/n
		mX[,j] <- (X[,j] - mu.x[j])/sqrt(sigma2.x[j])
	} 
	
	return(list(vy=vy,mX=mX,mu.y=mu.y,sigma2.y=sigma2.y,mu.x=mu.x,sigma2.x=sigma2.x))
}

calculate.model.posteriors <- function(vy,mX,CHECK) 
{	
	res = normalize(vy,mX) 
	vy = res$vy
	mX = res$mX

	n = length(vy)
	p = ncol(mX)
	
	vR2 = correlation::all_correlations_mX(vy, mX, "maruyama")$vR2
	M = length(vR2)
	
	mGamma = graycode(p,0)
	vq = mGamma%*%matrix(1,p,1)
	

	if (CHECK) {
		XTy = t(mX)%*%vy
		XTX = t(mX)%*%mX
		
		# Check Graycode order
		verr = c()
		for (i in 1:min(c(M,1000))) {
			inds = which(mGamma[i,]==1)
			if (length(inds)>2) {
				R2 = t(XTy[inds])%*%solve(XTX[inds,inds])%*%XTy[inds]/n
				verr <- c(verr,abs(R2-vR2[i]))
			}
		}
		print("The maximum error is:")
		print(max(verr))
	}
	
	vBIC = n*log(1 - vR2) + vq*log(n) 
	
	a  = -0.75
	vb = 0.5*(n - vq - 5) - a
	c  = 0.5*(n - 1)
	vd = 0.5*vq + a
	
	log.vp = -(vb+1)*log(1 - vR2) + lbeta(vd+1,vb+1) - lbeta(a+1,vb+1)
	vZE <- -2*log.vp
	
	model.post <- exp(log.vp - max(log.vp))/sum(exp(log.vp - max(log.vp)))
	variable.post <- t(mGamma)%*%model.post
	
	return(list(vR2=vR2,mGamma=mGamma,vq=vq,vBIC=vBIC,vZE=vZE,log.vp=log.vp,model.post=model.post,variable.post=variable.post))
}

