
plot.gGivenY <- function(g,n,p,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
{	
	den1 = posterior.gGivenY.MonteCarlo(g,n,p,R2,N=1000)
	
	if (is.null(g)) {
		g = den1$x 
	}
	
	den2 = posterior.gGivenY.exact(g,n,p,R2) 
	
	xlim = range(g)
	ylim = range(c(den1$y,den2$y))
	
	plot(NA,type="n", xlim=xlim,ylim=ylim,ylab="density",xlab="g",cex.lab=1.5,
		cex.main=1.5,main=paste("g|y   (n=",n,", p=",p,", R2=",round(R2,3),")",sep=""))
	lines(den1,col="red",lwd=2,lty=2)
	lines(den2,col="black",lwd=2,lty=1)

	if (LEGEND) {
		legend("topright", legend = c("Exact","Samples"), col = c("black","red"), lwd=c(2,2), lty=c(1,2), cex=1.5, bty="n")
	}
	
	if (CHECK) {
		print(trapint(den1$x,den1$y))
		print(trapint(den2$x,den2$y))
	}
}


plot.uGivenY <- function(u,n,p,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
{
	den1 = posterior.uGivenY.MonteCarlo(u,n,p,R2,N=1000)
	
	if (is.null(u)) {
		u = den1$x 
	}
	
	den2 = posterior.uGivenY.exact(u,n,p,R2) 
	
	xlim = range(u)
	ylim = range(c(den1$y,den2$y))
	
	plot(NA,type="n", xlim=xlim,ylim=ylim,ylab="density",xlab="u",cex.lab=1.5,
		cex.main=1.5,main=paste("u|y   (n=",n,", p=",p,", R2=",round(R2,3),")",sep=""))
	lines(den1,col="red",lwd=2,lty=2)
	lines(den2,col="black",lwd=2,lty=1)
	
	if (LEGEND) {
		legend("topleft", legend = c("Exact","Samples"), col = c("black","red"), lwd=c(2,2), lty=c(1,2), cex=1.5, bty="n")
	}
	
	if (CHECK) {
		print(trapint(den1$x,den1$y))
		print(trapint(den2$x,den2$y))
	}
}


plot.sigma2GivenY <- function(sigma2,n,p,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
{	
	den1 = posterior.sigma2GivenY.MonteCarlo(sigma2,n,p,R2,N=1000)
	
	if (is.null(sigma2)) {
		sigma2 = den1$x
	}
	
	den2 = posterior.sigma2GivenY.exact(sigma2,n,p,R2) 
	den3 = posterior.sigma2GivenY.RaoBlackwell(sigma2,n,p,R2,N=100)
	den4 = posterior.sigma2GivenY.approx(sigma2,n,p,R2)
	
	xlim = range(sigma2)
	ylim = range(c(den1$y,den3$y,den4$y))
	
	plot(NA,type="n", xlim=xlim,ylim=ylim,ylab="density",xlab=expression(sigma^2),cex.lab=1.5,
		cex.main=1.5,main=TeX(paste("$\\sigma^2|y$   (n=",n,", p=",p,", R2=",round(R2,3),")",sep="")))
	lines(den1,col="red",lwd=1,lty=1)
	lines(den2,col="black",lwd=3,lty=3)
	lines(den3,col="blue",lwd=3,lty=3)
	lines(den4,col="green",lwd=3,lty=3)
	
	if (LEGEND) {
		legend("topright", 
			legend = c("Exact","Samples","RB","Delta"), 
			col = c("black","red","blue","green"), 
			lwd=c(3,1,3,3), 
			lty=c(3,1,3,3), 
			cex=1.5,
			bty="n")
	}
	
	if (CHECK) {
		print(trapint(den1$x,den1$y))
		print(trapint(den2$x,den2$y))
		print(trapint(den3$x,den3$y))
		print(trapint(den4$x,den4$y))
	}
}

plot.alphaGivenY <- function(alpha,n,p,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
{	
	den1 = posterior.alphaGivenY.MonteCarlo(alpha,n,p,R2,N=1000)
	
	if (is.null(alpha)) {
		alpha = den1$x # den1$samples
	}
	
	den2 = posterior.alphaGivenY.exact(alpha,n,p,R2)
	den3 = posterior.alphaGivenY.RaoBlackwell(alpha,n,p,R2,N=100)
	den4 = posterior.alphaGivenY.approx(alpha,n,p,R2)
	
	xlim = range(alpha)
	ylim = range(c(den1$y,den3$y,den4$y))
	
	plot(NA,type="n", xlim=xlim,ylim=ylim,ylab="density",xlab=expression(alpha),cex.lab=1.5,
		cex.main=1.5,main=TeX(paste("$\\alpha |y$   (n=",n,", p=",p,", R2=",round(R2,3),")",sep="")))
	lines(den1,col="red",lwd=1,lty=1)
	lines(den2,col="black",lwd=3,lty=3)
	lines(den3,col="blue",lwd=3,lty=3)
	lines(den4,col="green",lwd=3,lty=3)
	
	if (LEGEND) {
		legend("topright", 
			legend = c("Exact","Samples","RB","Delta"), 
			col = c("black","red","blue","green"), 
			lwd=c(3,1,3,3), 
			lty=c(3,1,3,3), 
			cex=1.5,
			bty="n")
	}
	
	if (CHECK) {
		print(trapint(den1$x,den1$y))
		print(trapint(den2$x,den2$y))
		print(trapint(den3$x,den3$y))
		print(trapint(den4$x,den4$y))
	}
}

