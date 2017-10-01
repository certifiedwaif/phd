
library(gsl)

posterior.gGivenY.MonteCarlo <- function(n,p,R2,N) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)

	den = density(vg)
	
	return(list(samples=vg,x=den$x,y=den$y))
}

posterior.gGivenY.exact <- function(g,n,p,R2) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	if (is.null(g)) {
		res = posterior.gGivenY.MonteCarlo(n,p,R2,N=1000)
		g = res$x
		g = g[g>0] 
	}
	
	log.f = 0
	log.f = log.f + (b + 1)*log(1 - R2) 
	log.f = log.f + b*log(g) 
	log.f = log.f - (d + b + 2)*log(1 + g*(1 - R2)) 
	log.f = log.f - lbeta(d+1,b+1)
	
	return(list(x=g,y=exp(log.f)))
}


posterior.uGivenY.MonteCarlo <- function(n,p,R2,N) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)
	vu = vg/(1 + vg)
	
	den = density(vu)
	
	return(list(samples=vu,x=den$x,y=den$y))
}


posterior.uGivenY.exact <- function(u,n,p,R2) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	if (is.null(u)) {
		res = posterior.uGivenY.MonteCarlo(n,p,R2,N=1000)
		u = res$x
		u = u[(u>0)&(u<1)] 
	}
	
	log.f = 0
	log.f = log.f +  (b + 1)*log(1 - R2) 
	log.f = log.f - lbeta(d+1,b+1) 
	log.f = log.f + b*log(u) 
	log.f = log.f + d*log(1 - u) 
	log.f = log.f - (d + b + 2)*log(1 - u*R2)
	
	return(list(x=u,y=exp(log.f)))
}



posterior.sigma2GivenY.MonteCarlo <- function(n,p,R2,N) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)
	vu = vg/(1 + vg)
	vsigma2.post = 1/rgamma(N,c,0.5*n*(1 - R2*vu))
	
	den = density(vsigma2.post)
	
	return(list(samples=vsigma2.post,x=den$x,y=den$y))
}


posterior.sigma2GivenY.exact <- function(sigma2,n,p,R2) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a

	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		sigma2 = res$x
		sigma2 = sigma2[sigma2>0] 
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
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		sigma2 = res$x
		sigma2 = sigma2[sigma2>0] 
	}
	
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)
	vu = vg/(1 + vg)

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
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		sigma2 = res$x
		sigma2 = sigma2[sigma2>0] 
	}
	
	M1 = (b + 1)*hyperg_2F1(d+1, 1, c+1, R2, give=FALSE, strict=TRUE)/c
	
	s = c
	t = 0.5*n*(1 - M1*R2)
	y = dinvgamma(sigma2,s,t)
	
	return(list(x=sigma2,y=y,s=c,t=t))
}



posterior.alphaGivenY.MonteCarlo <- function(n,p,R2,N) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)
	vu = vg/(1 + vg)
	vsigma2.post = 1/rgamma(N,c,0.5*n*(1 - R2*vu))
	valpha.post = rnorm(N,0,sqrt(vsigma2.post/n))
	
	den = density(valpha.post)
	
	return(list(samples=valpha.post,x=den$x,y=den$y))
}


posterior.alphaGivenY.exact <- function(alpha,n,p,R2) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
		
	if (is.null(alpha)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		alpha = res$x
	}
	
	
	log.f1 = 0
	log.f1 = log.f1 + lgamma(c+0.5)   
	log.f1 = log.f1 + (b+1)*log(1 - R2)
	log.f1 = log.f1 - lgamma(c)
	log.f1 = log.f1 - 0.5*log(pi)
	log.f1 = log.f1 - 0.5*n*log(1 + alpha^2)
	log.f1 = log.f1 + log( hyperg_2F1(b+1, c+0.5, c, R2/(1 + alpha^2), give=FALSE, strict=TRUE))
 
	
	log.f2 = 0
	log.f2 = log.f2 + lgamma(c + 0.5)   
	log.f2 = log.f2 - lgamma(c)
	log.f2 = log.f2 - 0.5*log(pi)
	log.f2 = log.f2 + (b+1)*log(1 - R2)
	log.f2 = log.f2 - (b+3/2)*log(alpha^2 + 1 - R2)
	log.f2 = log.f2 - (d+1)*log(1+alpha^2)
	log.f2 = log.f2 + log( hyperg_2F1(d+1, -0.5, c, R2/(1 + alpha^2), give=FALSE, strict=TRUE))
	
	log.f = ifelse(is.nan(log.f1),log.f2,log.f1)

	return(list(x=alpha,y=exp(log.f)))
}


posterior.alphaGivenY.RaoBlackwell <- function(alpha,n,p,R2,N) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	if (is.null(alpha)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		alpha = res$x
	}
	
		
	vB = rbeta(N,b+1,d+1)
	vg.tilde = vB/(1 - vB)
	vg = vg.tilde/(1 - R2)
	vu = vg/(1 + vg)
	vsigma2.post = 1/rgamma(N,c,0.5*n*(1 - R2*vu))
	
	vmu    = rep(0,N)
	vsigma = sqrt(vsigma2.post/n)
	
	y = 0*alpha
	for (i in 1:N) {
		y = y + (1/N)*dnorm(alpha, vmu[i], vsigma[i]) 
	}

	return(list(mu=vmu,sigma=vsigma,x=alpha,y=y))
}


posterior.alphaGivenY.approx <- function(alpha,n,p,R2,N) 
{
	a = -0.75
	b = 0.5*(n - p - 5) - a
	c = 0.5*(n - 1)
	d = 0.5*p + a
	
	if (is.null(alpha)) {
		res = posterior.sigma2GivenY.MonteCarlo(n,p,R2,N=1000)
		alpha = res$x
	}
	
	M1 = (b + 1)*hyperg_2F1(d+1, 1, c+1, R2, give=FALSE, strict=TRUE)/c
	
	y = dnorm(alpha,0,sqrt((1 - M1*R2)/(n-1)))
	
	return(list(x=alpha,y=y))
}


trapint <- function(xgrid, fgrid) 
{
	ng <- length(xgrid)
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
	integ <- sum(xvec * fvec)/2
	return(integ)
}



