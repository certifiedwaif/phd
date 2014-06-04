
set.seed(1)

# Suppose that we want to calculate P(Z>3) given by:
1 - pnorm(3)

MAXITER <- 1000000
N <- 10000000

print("Classical Monte Carlo")
h <- c()
for (ITER in 1:MAXITER) 
{
		x <- rnorm(N)
		h <- c(h, as.integer(x > 3))
		I1.hat <- mean(h)
	  se <- sqrt(var(h)/length(h))
		if ((se)<(5.0E-6)) {
				break;
		}
		print(c(length(h),I1.hat,se))
}


x <- seq(3,10,,1000)
f <- dnorm(x)
plot(x,f,type="l")
g <- dnorm(x,4,1)
lines(x,g,col="red")


print("Importance Sampling")

N <- 100000
h <- c()
w <- c()
for (ITER in 1:MAXITER) 
{
		x <- rnorm(N,4,1)
		f <- dnorm(x)
		g <- dnorm(x,4,1)
		w <- c(w,f/g)
		h <- c(h, as.integer(x > 3))
		I2.hat <- mean(w*h)
	  se <- sqrt(var(w*h)/length(h))
		if ((se)<(1.0E-8)) {
				break;
		}
		print(c(length(h),I2.hat,se))
}

