
mix.posterior.gGivenY <- function(g,n,vq,vR2,vw,normalize=TRUE,methods=c("exact","MC"),N=1000) 
{
	if (normalize) {
		vw = vw/sum(vw)
	}	
	M = length(vw)
	ord = order(vw,decreasing=TRUE)
	
	vq  = vq[ord]
	vR2 = vR2[ord]
	vw  = vw[ord] 
		
	if (is.null(g)) {
		res = posterior.gGivenY.MonteCarlo(g,n,vq[1],vR2[1],N=10000)
		g = res$x
	}
	
	lx <- vector("list", length(methods))
	ly <- vector("list", length(methods))
	
	for (m in 1:length(methods)) 
	{
		lx[[m]] = g
		ly[[m]] = 0*g
		for (i in 1:M) {
			if ("exact"%in%methods) {
				res = posterior.gGivenY.exact(g,n,vq[i],vR2[i])
			} 
			if ("MC"%in%methods) {
				res = posterior.gGivenY.MonteCarlo(g,n,vq[i],vR2[i],N)
			}	
			ly[[m]] = ly[[m]] + vw[i]*res$y
		}
	}
	
	return(list(lx=lx,ly=ly))
}

mix.posterior.uGivenY <- function(u,n,vq,vR2,vw,normalize=TRUE,methods=c("exact","MC"),N=1000) 
{
	if (normalize) {
		vw = vw/sum(vw)
	}	
	M = length(vw)
	ord = order(vw,decreasing=TRUE)
	
	vq  = vq[ord]
	vR2 = vR2[ord]
	vw  = vw[ord] 
	
	if (is.null(u)) {
		res = posterior.uGivenY.MonteCarlo(u,n,vq[1],vR2[1],N=10000)
		u = res$x
	}
	
	lx <- vector("list", length(methods))
	ly <- vector("list", length(methods))
	
	for (m in 1:length(methods)) 
	{
		lx[[m]] = u
		ly[[m]] = 0*u
		for (i in 1:M) {
			if ("exact"%in%methods) {
				res = posterior.uGivenY.exact(u,n,vq[i],vR2[i])
			} 
			if ("MC"%in%methods) {
				res = posterior.uGivenY.MonteCarlo(u,n,vq[i],vR2[i],N)
			}	
			ly[[m]] = ly[[m]] + vw[i]*res$y
		}
	}
	
	return(list(lx=lx,ly=ly))
}

mix.posterior.sigma2GivenY <- function(sigma2,n,vq,vR2,vw,normalize=TRUE,methods=c("exact","MC","RB","delta"),N=1000) 
{
	if (normalize) {
		vw = vw/sum(vw)
	}	
	M = length(vw)
	ord = order(vw,decreasing=TRUE)
	
	vq  = vq[ord]
	vR2 = vR2[ord]
	vw  = vw[ord] 
	
	if (is.null(sigma2)) {
		res = posterior.sigma2GivenY.MonteCarlo(sigma2,n,vq[1],vR2[1],N=10000)
		sigma2 = res$x
	}
	
	lx <- vector("list", length(methods))
	ly <- vector("list", length(methods))
	
	for (m in 1:length(methods)) 
	{
		lx[[m]] = sigma2
		ly[[m]] = 0*sigma2
		for (i in 1:M) {
			if ("exact"%in%methods) {
				res = posterior.sigma2GivenY.exact(sigma2,n,vq[i],vR2[i])
			} 
			if ("MC"%in%methods) {
				res = posterior.sigma2GivenY.MonteCarlo(sigma2,n,vq[i],vR2[i],N)
			}	
			if ("RB"%in%methods) {
				res = posterior.sigma2GivenY.RaoBlackwell(sigma2,n,vq[i],vR2[i],N)
			}
			if ("delta"%in%methods) {
				res = posterior.sigma2GivenY.approx(sigma2,n,vq[i],vR2[i])
			}		
			ly[[m]] = ly[[m]] + vw[i]*res$y
		}
	}
	
	return(list(lx=lx,ly=ly))
}


mix.posterior.alphaGivenY <- function(alpha,n,vq,vR2,vw,normalize=TRUE,methods=c("exact","MC","RB","delta"),N=1000) 
{
	if (normalize) {
		vw = vw/sum(vw)
	}	
	M = length(vw)
	ord = order(vw,decreasing=TRUE)
	
	vq  = vq[ord]
	vR2 = vR2[ord]
	vw  = vw[ord] 
	
	if (is.null(alpha)) {
		res = posterior.alphaGivenY.MonteCarlo(alpha,n,vq[1],vR2[1],N=10000)
		sigma2 = res$x
	}
	
	lx <- vector("list", length(methods))
	ly <- vector("list", length(methods))
	
	for (m in 1:length(methods)) 
	{
		lx[[m]] = alpha
		ly[[m]] = 0*alpha
		for (i in 1:M) {
			if ("exact"%in%methods) {
				res = posterior.alphaGivenY.exact(alpha,n,vq[i],vR2[i])
			} 
			if ("MC"%in%methods) {
				res = posterior.alphaGivenY.MonteCarlo(alpha,n,vq[i],vR2[i],N)
			}	
			if ("RB"%in%methods) {
				res = posterior.alphaGivenY.RaoBlackwell(alpha,n,vq[i],vR2[i],N)
			}
			if ("delta"%in%methods) {
				res = posterior.alphaGivenY.approx(alpha,n,vq[i],vR2[i])
			}		
			ly[[m]] = ly[[m]] + vw[i]*res$y
		}
	}
	
	return(list(lx=lx,ly=ly))
}
