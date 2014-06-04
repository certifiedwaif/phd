
N <- 1000

X.curr <- 1

X <- c()

alpha <- 2.43

for (i in 1:N) {
   Y <- rgamma(1,floor(alpha),floor(alpha)/alpha)
   ratio <-  ( (Y/X.curr)*exp((X.curr - Y)/alpha) )^(alpha - floor(alpha))
   r <- min(c(1,ratio))
   if (runif(1)<r) {
       X.curr <- Y
       X <- c(X,Y)
   }
}

x <- seq(min(X),max(X),,1000)
f <- dgamma(x,alpha,1)

pdf("GammaMH.pdf",width=12)

mai.dft <- c(1.02,0.82,0.82,0.42)
op    <- par(mai=mai.dft)
pcnt1  <- 0.40
pcnt2  <- 0.40
pcnt3  <- 0.40
pcnt4  <- 0.10
mai.new <- c(pcnt1,pcnt2,pcnt3,pcnt4)*mai.dft				
par(mfrow=c(1,2),mai=mai.new)

plot(x,f,type="l",lwd=3)
lines(density(X),lwd=3,col="red")

plot(X,type="l")

dev.off()

