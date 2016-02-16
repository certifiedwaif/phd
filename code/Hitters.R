
rm(list = ls())
set.seed(1)

##########################################################################3

library(ISLR)
Hitters <- na.omit(Hitters)


# Get y vector and X matrix
y.t <- Hitters$Salary
X.f <- model.matrix(Salary~.,Hitters)[,-1]

# Get dimensions
n <- nrow(X.f)
p <- ncol(X.f)

# Normalise
y.n <- (y.t - mean(y.t))/sd(y.t)
X.n <- matrix(0,n,p)
for (j in 1:p) {
	X.n[,j] <- (X.f[,j] - mean(X.f[,j]))/sd(X.f[,j])
} 

##########################################################################

# Simulate data (based on real data so we know truth)

# Indices of true model
inds.true <- c(1, 2, 6, 10, 11, 13, 15, 16)

q <- length(inds.true)
res.lm <- lm(y.n~-1+X.n[,inds.true])
beta.true <- res.lm$coef
y.pred <- X.n[,inds.true]%*%beta.true
sigma2.pred <- sum((y.n - y.pred)^2)/(n-q)
y.sim <- y.pred + rnorm(n,0,sqrt(sigma2.pred))
y.sim <- y.sim - mean(y.sim)
y.n <- (y.sim - mean(y.sim))/sd(y.sim)

vgamma <- rep(0,p)
vgamma[inds.true] <- 1

##############################################

write.table(y.n, sep=",", file="vy.csv", row.names = FALSE, col.names = FALSE)
write.table(X.n, sep=",", file="mX.csv", row.names = FALSE, col.names = FALSE)

source("ZE.R")

library(lineprof)

if (FALSE) {
	print("Exact using lm")
	a3 <- proc.time()[3]
	res.exact1 <- ZE.exact.slow(y.n,X.n,FALSE)
	b3 <- proc.time()[3]     
	print(b3-a3)
	write.table(res.exact1$vR2, sep=",", file="Hitters_exact1.csv", row.names = FALSE, col.names = FALSE)
}

if (TRUE) {
	print("Exact not using lm")
	a3 <- proc.time()[3]
	#library(lineprof)
	#l1 <- lineprof( 
	res.exact2 <- ZE.exact(y.n,X.n,FALSE,TRUE)
	#, interval = 0.01)
	b3 <- proc.time()[3]     
	print(b3-a3)
	
	vlog.ZE <- res.exact2$vlog.ZE
	vlog.ZE <- vlog.ZE - max(vlog.ZE)
	vprobs.temp <- exp(vlog.ZE)/sum(exp(vlog.ZE))
	vp.post <- vprobs.temp%*%res.exact2$A
	
	ord <- order(vprobs.temp,decreasing=TRUE)
	print("The top 10 models are:")
	for (i in 1:10) {
		cat(vprobs.temp[ord[i]],which(res.exact2$A[ord[i],]==1),"\n")
	}
	#shine(l1)
	write.table(res.exact2$vR2, sep=",", file="Hitters_exact2.csv", row.names = FALSE, col.names = FALSE)
}
 

 

a3 <- proc.time()[3]
#l2 <- lineprof( 
	res.exact3 <- ZE.exact.fast(y.n,X.n,FALSE)
#, interval = 0.01)
b3 <- proc.time()[3]     
print(b3-a3)

cat("Fitted ", length(res.exact3$vlog.ZE)/(b3-a3), " models per second\n")

vlog.ZE <- res.exact3$vlog.ZE
vlog.ZE <- vlog.ZE - max(vlog.ZE)
vprobs.temp <- exp(vlog.ZE)/sum(exp(vlog.ZE))
vp.post <- vprobs.temp%*%res.exact3$mA

ord <- order(vprobs.temp,decreasing=TRUE)
print("The top 10 models are:")
for (i in 1:10) {
	cat(vprobs.temp[ord[i]],which(res.exact3$mA[ord[i],]==1),"\n")
}

print(max( res.exact3$vlog.ZE ))

 
#shine(l2)
#max(abs(res.exact3$vlog.ZE - res.exact2$vlog.ZE))

 