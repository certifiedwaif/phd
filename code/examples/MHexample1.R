

MH.Cauchy <- function(X0,b,N) 
{
    vX <- c()
    X <- X0
    for (i in 1:N) {
        Y <- rnorm(1,X,b)
        r <- min(c(1, (1 + X^2)/(1 + Y^2) ))
        if (runif(1)<r) { X <- Y; }  
        vX <- c(vX,X)
    }
    return(vX)
}


N <- 1000

b1 <- 0.1
b2 <- 1
b3 <- 10

X0 <- 0

vX1 <- MH.Cauchy(X0,b1,N) 
vX2 <- MH.Cauchy(X0,b2,N) 
vX3 <- MH.Cauchy(X0,b3,N) 

###############################################################################

pdf("MHexample1.pdf",width=12)

mai.dft <- c(1.02,0.82,0.82,0.42)
op    <- par(mai=mai.dft)
pcnt1  <- 0.40
pcnt2  <- 0.40
pcnt3  <- 0.40
pcnt4  <- 0.10
mai.new <- c(pcnt1,pcnt2,pcnt3,pcnt4)*mai.dft				
par(mfrow=c(2,3),mai=mai.new)

plot(vX1,type="l",lwd=2,main="b = 0.1",cex.main=2,xlab="",ylab="")
plot(vX2,type="l",lwd=2,main="b = 1"  ,cex.main=2,xlab="",ylab="")
plot(vX3,type="l",lwd=2,main="b = 10" ,cex.main=2,xlab="",ylab="")

range.x <- range(c(vX1,vX2,vX3))
x <- seq(range.x[1],range.x[2],,1000)
f <- dt(x,1)

plot(density(vX1),type="l",lwd=2,main="b = 0.1",cex.main=2,xlab="",ylab="")
lines(x,f,lwd=2,col="red")
plot(density(vX2),type="l",lwd=2,main="b = 1"  ,cex.main=2,xlab="",ylab="")
lines(x,f,lwd=2,col="red")
plot(density(vX3),type="l",lwd=2,main="b = 10" ,cex.main=2,xlab="",ylab="")
lines(x,f,lwd=2,col="red")

dev.off()

###############################################################################