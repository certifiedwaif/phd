################################################################################

set.seed(11)

rm(list = ls())

dataset <- "communities.Rdata"
print(dataset)
 
load("comData.Rdata")

################################################################################

# Data preparation

sum.na <- function(x) {  sum(is.na(x)); }
inds <- which(apply(X,2,sum.na)==0)
X2 <- X[,inds]
X3 <- X2[,!colnames(X2)%in%c("ownHousQrange","rentQrange")]

y <- Y[,18]
inds <- which(is.na(y))

vy <- y[-inds]
mX <- X3[-inds,]
mX.til <- cbind(1,mX) 

n <- length(vy)
p <- ncol(mX)

mult <- sqrt(n/(n-1))
mX <- mX
for (j in 1:p) {
	mX[,j] = mult*(mX[,j] - mean(mX[,j]))/sd(mX[,j])	
}
vy <- mult*(vy - mean(vy))/sd(vy)

################################################################################

library(correlation)

res.lm <- lm(vy~mX)
inds <- which( summary( res.lm )$coef[,4] < 0.05) - 1

a0 <- proc.time()[3]
res.fast <- all_correlations_mX(vy, mX[,inds], intercept_col = 0, bIntercept = FALSE, bCentre = FALSE, cores = 4)
b0 <- proc.time()[3]    
t0 <- b0-a0
cat(t0,"\n")
#cat((2^length(inds))/t0,"\n")

################################################################################

digitsBase <- function (x, base = 2, ndigits = 1 + floor(1e-09 + log(max(x), base))) 
{
	if (any(x < 0)) 
		stop("'x' must be non-negative integers")
	if (any(x != trunc(x))) 
		stop("'x' must be integer-valued")
	r <- matrix(0, nrow = ndigits, ncol = length(x))
	if (ndigits >= 1) { 
		for (i in ndigits:1) { 
			r[ndigits-i+1, ] <- x%%base
			if (i > 1) {
				x <- x%/%base
			}
		}
	}
	return(r)
}

valToProb <- function(x) {
	x <- x - max(x)
	return( exp(x)/sum(exp(x)))
}

#################################################################################

# Iteratively go through all combinations of for a subset of variables of size K.
# Hill climb over these sets using BIC

# Size of the subset to search over
K <- 18

MAXITER <- 500

my.all_correlations_mX <- function(vy, mX) 
{
	n <- length(vy)
	p <- ncol(mX)
	mM <- graycode(p,0)
	XTX <- t(mX)%*%mX
	XTy <- t(mX)%*%vy
	vr2 <- rep(0,nrow(mM))
	for (i in 1:nrow(mM)) 
	{
		inds <- which(mM[i,]==1)
		q <- length(inds)
		if (q>0) {
			vb <- matrix(XTy[inds],q,1)
			mA <- matrix(XTX[inds,inds],q,q)
			vr2[i] <- t(vb)%*%solve(mA,vb)/n	
		} 
	}
	return(vr2)
}



my.all_correlations_mX_mZ <- function(vy, mX,mZ) 
{
	n  <- length(vy)
	p0 <- ncol(mX)
	p  <- ncol(mZ)
	mM <- graycode(p,p0)
	mC <- cbind(mX,mZ)
	CTC <- t(mC)%*%mC
	CTy <- t(mC)%*%vy
	vr2 <- rep(0,nrow(mM))
	for (i in 1:nrow(mM)) 
	{
		inds <- which(mM[i,]==1)
		q <- length(inds)
		if (q>0) {
			vb <- matrix(CTy[inds],q,1)
			mA <- matrix(CTC[inds,inds],q,q)
			vr2[i] <- t(vb)%*%solve(mA,vb)/n	
		}
	}
	return(vr2)
}

bestInds <- c()
vbest <- c()

for (ITER in 1:MAXITER) {

	inds <- sample(p)[1:K]
	fixedInds <- setdiff( bestInds, inds )
	nFixed <- length(fixedInds)
	
	cat("inds=",inds,"\n")
	
	if (nFixed==0) {
		res <- all_correlations_mX(vy, mX[,inds], intercept_col = 0, bIntercept = FALSE, bCentre = FALSE, cores = 1)
		res[1] <- 0

		res2 <- my.all_correlations_mX(vy, mX[,inds]) 
		
		res <- res2
		
		err <-  max(abs(res2 - res))
		#cat(ITER,err,"\n")
		
		#ans <- readline()
	} else { 
		res <- all_correlations_mX_mZ(vy, mX[,fixedInds], mX[,inds], intercept_col = 0, bIntercept = FALSE, bCentre = FALSE, cores = 1)
		#res[1] <- 0
		
		res2 <- my.all_correlations_mX_mZ(vy, mX[,fixedInds], mX[,inds]) 
	
		res <- res2
	
		err <-  max(abs(res2 - res))
		#cat(ITER,err,"\n")
		
		#ans <- readline()
	}
	
	mM <- graycode(K,0)
	colnames(mM) <- colnames(mX)[inds]
	vq <- apply(mM,1,sum) + nFixed
	
	# Compute the BIC for all models
	log.vp <-  -0.5*n*log(1 - res) - 0.5*vq*log(n)
	
	# Reset the current best model
	bestInds <- c(fixedInds, inds[ which( mM[ which.max(log.vp),]==1) ] )
	
	#cat("fixedInds=",fixedInds,"\n")
	#cat("bestInds=",bestInds,"\n")
	
	#ans <- readline()
	
	vbest[ITER] <- max(log.vp)
	if (ITER>1) {
		plot(vbest,type="b")
	}
	
	cat(ITER,max(vbest),"\n")
}


################################################################################
 