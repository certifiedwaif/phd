###############################################################################

library(glmmAK)
data(toenail)

vy <- c()
mX <- c()
mZ <- c()
vn <- c()
ids <- unique(toenail$idnr) 					
m <- length(ids)		
for (i in 1:m) {		
	  ind <- (toenail$idnr==ids[i])
	  ni  <- sum(ind)
	  vyi <- toenail$infect[ind]
	  mXi <- cbind(rep(1,ni),toenail$trt[ind], toenail$time[ind], toenail$time[ind]*toenail$trt[ind])
	  mZi <- matrix(0,ni,m)
	  mZi[,i] <- 1
	  vy <- c(vy,vyi)
	  mX <- rbind(mX,mXi)
	  mZ <- rbind(mZ,mZi)
	  vn[i] <- ni
}

p <- ncol(mX)
vy <- matrix(vy)
mC <- cbind(mX,mZ)

###############################################################################

# Set initial values

sigma2.p <- 0.001
sigma2.u <- 1
sigma2.beta <- 1.0E8
A <- 0.01
B <- 0.01

vbeta <- rep(0,p)
vu    <- rnorm(m,0,4)

N <- 50000

mBeta    <- matrix(0,N,p)
mU       <- matrix(0,N,m)
mSigma.u <- matrix(0,N,1)


for (ITER in 1:N) 
{
    vbeta.prop <- rnorm(p,vbeta,sqrt(sigma2.p))
    vu.prop    <- rnorm(m,vu,sqrt(sigma2.p))
    
    vnu      <- matrix(c(vbeta,vu))
    vnu.prop <- matrix(c(vbeta.prop,vu.prop))
    
    log.RT <- sum(vy*mC%*%vnu.prop - log(1 + exp(mC%*%vnu.prop))) - 0.5*sum(vbeta.prop^2)/sigma2.beta - 0.5*sum(vu.prop^2)/sigma2.u 
    log.RB <- sum(vy*mC%*%vnu      - log(1 + exp(mC%*%vnu)))      - 0.5*sum(vbeta^2)/sigma2.beta      - 0.5*sum(vu^2)/sigma2.u 
    R <- exp(log.RT - log.RB)
    
    U <- runif(1)
    if (U<R) {
        vbeta <- vbeta.prop
        vu    <- vu.prop
    }
    
    sigma2.u <- 1/rgamma(1,A + m/2, B + 0.5*sum(vu^2))
    
		mBeta[ITER,] <- vbeta
		mU[ITER,] <- vu
		mSigma.u[ITER,] <- sigma2.u
}



###############################################################################

# Plot results

pdf("GibbsMHExample.pdf",width=12)

mai.dft <- c(1.02,0.82,0.82,0.42)
op    <- par(mai=mai.dft)
pcnt1  <- 0.40
pcnt2  <- 0.40
pcnt3  <- 0.40
pcnt4  <- 0.10
mai.new <- c(pcnt1,pcnt2,pcnt3,pcnt4)*mai.dft				
par(mfrow=c(2,5),mai=mai.new)

plot(mSigma.u,type="l",lwd=2,main=expression(sigma[u]^2),cex.main=2,xlab="",ylab="")
plot(mBeta[,1],type="l",lwd=2,main=expression(beta[0]),cex.main=2,xlab="",ylab="")
plot(mBeta[,2],type="l",lwd=2,main=expression(beta[1]),cex.main=2,xlab="",ylab="")
plot(mBeta[,3],type="l",lwd=2,main=expression(beta[2]),cex.main=2,xlab="",ylab="")
plot(mBeta[,4],type="l",lwd=2,main=expression(beta[3]),cex.main=2,xlab="",ylab="")

dens <- density(mSigma.u)
plot(dens,lwd=2)

for (i in 1:4) {
		dens <- density(mBeta[,i])
		plot(dens,lwd=2)
}

dev.off()

###############################################################################

# Plot results

pdf("GibbsMHExample2.pdf",width=12)

mai.dft <- c(1.02,0.82,0.82,0.42)
op    <- par(mai=mai.dft)
pcnt1  <- 0.40
pcnt2  <- 0.40
pcnt3  <- 0.40
pcnt4  <- 0.10
mai.new <- c(pcnt1,pcnt2,pcnt3,pcnt4)*mai.dft				
par(mfrow=c(2,5),mai=mai.new)

inds <- 25001:50000

plot(mSigma.u[inds],type="l",lwd=2,main=expression(sigma[u]^2),cex.main=2,xlab="",ylab="")
plot(mBeta[inds,1],type="l",lwd=2,main=expression(beta[0]),cex.main=2,xlab="",ylab="")
plot(mBeta[inds,2],type="l",lwd=2,main=expression(beta[1]),cex.main=2,xlab="",ylab="")
plot(mBeta[inds,3],type="l",lwd=2,main=expression(beta[2]),cex.main=2,xlab="",ylab="")
plot(mBeta[inds,4],type="l",lwd=2,main=expression(beta[3]),cex.main=2,xlab="",ylab="")

dens <- density(mSigma.u[inds])
plot(dens,lwd=2)

for (i in 1:4) {
		dens <- density(mBeta[inds,i])
		plot(dens,lwd=2)
}

dev.off()

###############################################################################

