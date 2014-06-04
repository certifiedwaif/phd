###############################################################################

# Simulate Example

set.seed(1)

k <- 20   # Number of cities
n <- 100   # Number of trials per city

# Simulate true values for vp
alpha <- 1
beta  <- 1
a <- 0.2
b <- 0.8
vp.true <- a + (b-a)*rbeta(k,alpha,beta)

vy <- rbinom(k,n,vp.true)

vp.hat <- vy/n
vz <- log(vp.hat/(1 - vp.hat))
vsigma2 <- 1/(n*vp.hat*(1 - vp.hat))

###############################################################################

# Starting values for mu and vpsi
mu   <- 0
vpsi <- rep(0,k)

# Number of Gibbs samples
N <- 1000

# Matrices to store chain
mMu  <- matrix(0,N,1)
mPsi <- matrix(0,N,k)
mP   <- matrix(0,N,k)

###############################################################################

# Main loop for Gibbs sampling

for (ITER in 1:N) {
   # Sample from mu
   mu   <- rnorm(1,mean(vpsi),1/sqrt(k))
   
   # Sample from vpis
   vd2  <- 1/(1 + 1/vsigma2)
   ve   <- vd2*(vz/vsigma2 + mu)
   vpsi <- rnorm(k,ve,sqrt(vd2))
   
   # Transformed vpsi to vp scale
   vp   <- 1/(1 + exp(-vpsi))
   
   # Store results
   mMu[ITER,]  <- mu
   mPsi[ITER,] <- vpsi
   mP[ITER,] <- vp
}

###############################################################################

# Plot results

pdf("GibbsExample.pdf",width=12)

mai.dft <- c(1.02,0.82,0.82,0.42)
op    <- par(mai=mai.dft)
pcnt1  <- 0.40
pcnt2  <- 0.40
pcnt3  <- 0.40
pcnt4  <- 0.10
mai.new <- c(pcnt1,pcnt2,pcnt3,pcnt4)*mai.dft				
par(mfrow=c(2,4),mai=mai.new)

plot(mMu,type="l",lwd=2,main=expression(mu),cex.main=2,xlab="",ylab="")

plot(mP[,1],type="l",lwd=2,main=expression(p[1]),cex.main=2,xlab="",ylab="")
lines(1:N,rep(vp.true[1],N),col="red",lwd=3)
plot(mP[,2],type="l",lwd=2,main=expression(p[2]),cex.main=2,xlab="",ylab="")
lines(1:N,rep(vp.true[2],N),col="red",lwd=3)
plot(mP[,3],type="l",lwd=2,main=expression(p[3]),cex.main=2,xlab="",ylab="")
lines(1:N,rep(vp.true[3],N),col="red",lwd=3)

dens <- density(mMu)
plot(dens,lwd=2)

for (i in 1:3) {
		dens <- density(mP[,i])
		plot(dens,lwd=2)
		lines(rep(vp.true[i],2),range(dens$y),col="red",lwd=3)
}

dev.off()

###############################################################################

