
rm(list = ls())

source("functions.Rs")

library(ISLR)
Hitters <- na.omit(Hitters)

# Get y vector and X matrix
y.t <- Hitters$Salary
X.f <- model.matrix(Salary~.,Hitters)[,-1]

n = length(y.t)
p = ncol(X.f)

res = calculate.model.posteriors(y.t,X.f,CHECK=TRUE) 
ord = order(res$log.vp,decreasing=TRUE)

TOPMODELS = 10

best = which.min(res$vZE)
tab = cbind( res$mGamma, round(100*res$model.post,2), round(res$vR2,3), res$vq, round(res$vZE - res$vZE[best],2), round(res$vBIC - res$vBIC[best],2) )[ord[1:TOPMODELS],]

print(tab)
print(res$variable.post)

##################################

# Plot posteriors for best model

q  = res$vq[best]
R2 = res$vR2[best]

##################################

pdf("NotBetaHitters.pdf",width=12,height=8)
par(mfrow=c(2,2))
plot.gGivenY(g=NULL,n,q,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
plot.uGivenY(u=NULL,n,q,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
plot.sigma2GivenY(sigma2=NULL,n,q,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
plot.alphaGivenY(alpha=NULL,n,q,R2,N=1000,LEGEND=TRUE,CHECK=TRUE)
dev.off()

##################################

attach(normalize(y.t,X.f),warn.conflicts=FALSE)
 
inds = which(res$mGamma[best,]==1)
res2 = posterior.betaGivenY.MonteCarlo(lbeta=NULL,vy,mX,inds,N=1000)
res3 = posterior.betaGivenY.RaoBlackwell(res2$lx,vy,mX,inds,N=100)
res4 = posterior.betaGivenY.approx(res2$lx,vy,mX,inds)
res5 = posterior.betaGivenY.MLE(res2$lx,vy,mX,inds)

pdf("BetaHitters.pdf",width=12,height=8)

par(mfrow=c(2,4))
for (j in 1:length(inds)) {

	J = inds[j]
	xlim = range(res2$lx[[j]])
	ylim = range(c(res2$ly[[j]],res3$ly[[j]],res4$ly[[j]],res5$ly[[j]]))
	
	plot(NA,type="n",xlim=xlim,ylim=ylim,main=colnames(X.f)[J],xlab="",ylab="Posterior density",cex.lab=1.5,cex.main=1.5,lwd=2)

	lines(res2$lx[[j]],res2$ly[[j]],col="red",lwd=2)
	lines(res3$lx[[j]],res3$ly[[j]],col="blue",lwd=2)
	lines(res4$lx[[j]],res4$ly[[j]],col="green",lwd=2)
	lines(res5$lx[[j]],res5$ly[[j]],col="purple",lwd=2)
	
	if (j==1) {
		legend("topleft", 
			legend = c("Samples","RB","Delta","MLE"), 
			col = c("red","blue","green","purple"), 
			lwd=c(2,2,2,2), 
			lty=c(1,1,1,1), 
			cex=1.5,
			bty="n")
	}
	
	
}

dev.off()

TOPMODELS = 1000

vq.a  = res$vq[ord[1:TOPMODELS]]
vR2.a = res$vR2[ord[1:TOPMODELS]]
vw.a  = res$model.post[ord[1:TOPMODELS]]

source("functions.Rs")

g = seq(0.001,500,,500)
u = seq(0.8,0.9999,,300)
sigma2 = seq(0.001,0.8,,1000)
alpha = seq(-0.3,0.3,,1000)

res.gGivenY = mix.posterior.gGivenY(g,n,vq=vq.a,vR2=vR2.a,vw=vw.a,normalize=TRUE,methods=c("exact"),N=100) 
res.uGivenY = mix.posterior.uGivenY(u,n,vq=vq.a,vR2=vR2.a,vw=vw.a,normalize=TRUE,methods=c("exact"),N=100) 
res.sigma2GivenY = mix.posterior.sigma2GivenY(sigma2,n,vq=vq.a,vR2=vR2.a,vw=vw.a,normalize=TRUE,methods=c("exact","RB","delta"),N=100)
res.alphaGivenY = mix.posterior.alphaGivenY(alpha,n,vq=vq.a,vR2=vR2.a,vw=vw.a,normalize=TRUE,methods=c("exact","RB","delta"),N=100) 

pdf("NotBetaMixHitters.pdf",width=12,height=8)
par(mfrow=c(2,2))

plot(g,res.gGivenY$ly[[1]],type="l",xlim=range(g))
#lines(g,res.gGivenY$ly[[2]],col="red")

plot(u,res.uGivenY$ly[[1]],type="l",xlim=range(u))
#lines(u,res.uGivenY$ly[[2]],col="red")

plot(sigma2,res.sigma2GivenY$ly[[1]],type="l",xlim=range(sigma2))
#lines(sigma2,res.sigma2GivenY$ly[[2]],col="red")
lines(sigma2,res.sigma2GivenY$ly[[2]],col="blue")
lines(sigma2,res.sigma2GivenY$ly[[3]],col="green")

plot(alpha,res.alphaGivenY$ly[[1]],type="l",xlim=range(alpha))
#lines(alpha,res.alphaGivenY$ly[[2]],col="red")
lines(alpha,res.alphaGivenY$ly[[2]],col="blue")
lines(alpha,res.alphaGivenY$ly[[3]],col="green")

dev.off()

XTy = t(mX)%*%vy
XTX = t(mX)%*%mX

 
lbeta =  vector("list", p)

G = 301
for (j in 1:p) {
	lbeta[[j]] = seq(-2,2,,G)
}


TOPMODELS = 1000

vw.b = res$model.post[ord[1:TOPMODELS]]
vw.b = vw.b/sum(vw.b)

mBeta.post   = matrix(0,p,G)
mBeta.MLE    = matrix(0,p,G)
mBeta.approx = matrix(0,p,G)

for (top in 1:TOPMODELS) 
{	
	inds = which(res$mGamma[ord[top],]==1)
	q = length(inds)
	
	if (q!=0) {
		XTX.inv = solve(XTX[inds,inds])
		vbeta.hat = XTX.inv%*%XTy[inds]
		R2 = t(XTy[inds])%*%solve(XTX[inds,inds])%*%XTy[inds]/n
		sigma2.hat = 1 - R2
		
		for (J in 1:length(inds)) 
		{
			j = inds[J]
			
			N = 100
			post = posterior.MC.nonBeta(n,q,R2,N) 
			vu = post$vu
			for (i in 1:N) {
				mBeta.post[j,] = mBeta.post[j,] + vw.b[top]*(1/N)*dnorm(lbeta[[j]], vu[i]*vbeta.hat[J], sqrt(n*vu[i]*(1 - vu[i]*R2)*XTX.inv[J,J]/(n-1)))
			}
			
			cons = caclulate.constants(n,p)
			a = cons$a; b=cons$b; c=cons$c; d=cons$d
			M1 = (b + 1)*hyperg_2F1(d+1, 1, c+1, R2, give=FALSE, strict=TRUE)/c
			mBeta.approx[j,] = mBeta.approx[j,] + vw.b[top]*dnorm(lbeta[[j]], M1*vbeta.hat[J],sqrt(M1*(1-M1*R2)*XTX.inv[J,J]))
			
			mBeta.MLE[j,]    = mBeta.MLE[j,] + vw.b[top]*dnorm(lbeta[[j]],vbeta.hat[J],sqrt(sigma2.hat*XTX.inv[J,J]))
		}
	
		par(mfrow=c(4,5))
		
		ylim = range(c(mBeta.post,mBeta.MLE,mBeta.approx))
		for (j in 1:p) {
			plot(lbeta[[j]],mBeta.post[j,],type="l",ylim=ylim,main=colnames(X.f)[j],
				xlab="",ylab="Posterior density",cex.lab=1.5,cex.main=1.5,lwd=2)
			lines(lbeta[[j]],mBeta.MLE[j,],col="red",lwd=2)
			lines(lbeta[[j]],mBeta.approx[j,],col="green",lwd=2)
		}
	}
}	


pdf("BetaMixHitters.pdf",width=12,height=8)
par(mfrow=c(4,5))

ylim = range(c(mBeta.post,mBeta.MLE,mBeta.approx))
for (j in 1:p) {	
	plot(lbeta[[j]],mBeta.post[j,],type="l",ylim=ylim,main=colnames(X.f)[j],
		xlab="",ylab="Posterior density",cex.lab=1.5,cex.main=1.5,lwd=2)
	lines(lbeta[[j]],mBeta.MLE[j,],col="red",lwd=2)
	lines(lbeta[[j]],mBeta.approx[j,],col="green",lwd=2)
	lines(c(0,0),c(0,ylim[2]*(1-res$variable.post[j])),col="black",lwd=3)
}

dev.off()

