
set.seed(1)

###############################################################################

x <- seq(-5,5,,1000)
f <- dt(x,3)
plot(x,f,type="l",ylim=c(0,0.4))
h <- abs(x)
lines(x,f*h,col="blue")

g1 <- dt(x,1)
lines(x,g1,col="green3")

g2 <- dnorm(x)
lines(x,g2,col="red")


###############################################################################

MAXITER <- 1000000
N <- 1000000

print("Classical Monte Carlo")
h <- c()
for (ITER in 1:MAXITER) 
{
		x <- rt(N,3)
		h <- c(h,abs(x))
		I1.hat <- mean(h)
	  se <- sqrt(var(h)/length(h))
		if ((se)<(5.0E-4)) {
				break;
		}
		print(c(length(h),I1.hat,se))
}

###############################################################################

print("Importance Sampling with t3 proposal")

N <- 100000
h <- c()
w <- c()
for (ITER in 1:MAXITER) 
{
		x <- rt(N,1)
		f <- dt(x,3)
		g <- dt(x,1)
		w <- c(w,f/g)
		h <- c(h,abs(x))
		I2.hat <- mean(w*h)
	  se <- sqrt(var(w*h)/length(h))
		if ((se)<(5.0E-4)) {
				break;
		}
		print(c(length(h),I2.hat,se))
}

###############################################################################

print("Importance Sampling with standard normal proposal")

N <- 1000000
h <- c()
w <- c()
for (ITER in 1:100) 
{
		x <- rnorm(N)
		f <- dt(x,3)
		g <- dnorm(x)
		w <- c(w,f/g)
		h <- c(h,abs(x))
		I3.hat <- mean(w*h)
	  se <- sqrt(var(w*h)/length(h))
		if ((se)<(5.0E-4)) {
				break;
		}
		print(c(length(h),I3.hat,se))
}

