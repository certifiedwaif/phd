

#rm(list = ls())

##########################################################################

# Choose a simulation setting

SETTING <- 1

if (SETTING==1) 
{
	library(ISLR)
	Hitters <- na.omit(Hitters)
	
	# Get y vector and X matrix
	y.t <- Hitters$Salary
	X.f <- model.matrix(Salary~.,Hitters)[,-1] 
	varnames <- colnames(X.f)
} 

if (SETTING==2) 
{
	# Read the bodyfat data
	dat  = read.table(file="bodyfat.txt",header=TRUE)
		
	# delete a number of obs with errors or otherwise extreme
	s.i  = c(39,42,48,96,76,182,31,86) 
	dat2 = dat[-s.i,-1]
	dat2$Weight  = round(0.45359237*dat2$Weight,2) # convert lbs into kg
	
	# Get y vector and X matrix
	y.t <- matrix(dat2$Bodyfat)
	X.f <- as.matrix(dat2[,-1]) # note: includes intercept

	varnames <- colnames(X.f)
}

if (SETTING==3) 
{
	library(ISLR)
	Wage <- na.omit(Wage) 
		
	# Get y vector and X matrix
	y.t <- Wage$wage
	X.f <- model.matrix(wage~.,Wage)[,-1]
	
	# Remove some columns (I think because of a lack of information from memory)
	X.f <- X.f[,-which(colnames(X.f)%in%c("sex2. Female",
	"region2. Middle Atlantic",
	"region3. East North Central",
	"region4. West North Central",     
	"region5. South Atlantic", 
	"region6. East South Central", 
	"region7. West South Central",           
	"region8. Mountain",
	"region9. Pacific"))]
 
	varnames <- colnames(X.f)
}

if (SETTING==4) 
{
	library(ISLR)
	College <- na.omit(College)
	
	# Get y vector and X matrix
	y.t <- College$Grad.Rate
	X.f <- model.matrix(Grad.Rate~.,College)[,-1]
 
	varnames <- colnames(X.f)
}

if (SETTING==5) 
{
	# Famous example, used in most papers.

	library(MASS)
	
	mD <- UScrime
	notlog <- c(2,ncol(UScrime))
	mD[,-notlog] <- log(mD[,-notlog])
	
	for (j in 1:ncol(mD)) {
		mD[,j] <- (mD[,j] - mean(mD[,j]))/sd(mD[,j])
	}
	
	varnames <- c(
	"log(AGE)",
	"S",
	"log(ED)",
	"log(Ex0)",
	"log(Ex1)",
	"log(LF)",
	"log(M)",
	"log(N)",
	"log(NW)",
	"log(U1)",
	"log(U2)",
	"log(W)",
	"log(X)",
	"log(prison)",
	"log(time)")
	
	y.t <- mD$y
	X.f <- data.matrix(cbind(mD[1:15]))
	colnames(X.f) <- varnames 
}

n <- nrow(X.f)
p <- ncol(X.f)

# Normalise
y.n <- (y.t - mean(y.t))/sd(y.t)
X.n <- matrix(0,n,p)
for (j in 1:p) {
	X.n[,j] <- (X.f[,j] - mean(X.f[,j]))/sd(X.f[,j])
} 

vy <- y.n
mX <- X.n

################################################################################

# Trapezoidal integration over a potentially irregular grid

trapint <- function(xgrid, fgrid) 
{ 
	ng <- length(xgrid) 
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)] 
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng] 
	integ <- sum(xvec * fvec)/2 
	return(integ) 
}   

################################################################################

# Solve a*x*x + b*x + c = 0
# Probably a numerically bad way of doing this see Num Rec in C

solveQuad <- function(a,b,c) 
{
	disc <- b*b - 4*a*c
	#cat("disc=",disc,"\n")
	if (disc<0) {
		val <- c(NA)
	} else {
		val1 <-(-b + sqrt(disc))/(2*a)
		val2 <-(-b - sqrt(disc))/(2*a) 
		val <- c(val1,val2)
	}
	return( val )
}

################################################################################

# The log the q-density for g

log.qg <- function(x,A,B,C) {
	return( A*log(x) + B*log(1 + x) - C/x )	
}


################################################################################

# The log the q-density for h (where h is the inverse of g)

log.qh <- function(x,U,B,C) {
	# cat(x, U, B, C, "\n")
	return( U*log(x) + B*log(1 + x) - C*x )	
}

################################################################################

# Calculate the mode of the q-density for g

mode.qg <- function(A,B,C) {
	res <- solveQuad(A+B,A+C,C) 
	if (length(res)>1) {
		res <- res[res>0]
	}
	return(res)
}

################################################################################

# Calculate the mode of the q-density for h (where h is the inverse of g)

mode.qh <- function(U,B,C) {
	res <- solveQuad(C,C-B-U,-U) 
	if (length(res)>1) {
		res <- res[res>0]
	}
	# cat("res ", res, "\n")
	return(res)
}

################################################################################

# Calculate the Hessian of the q-density for g evaluated at x

hessian.qg <- function(x,A,B,C) 
{
	x2 <- x*x
	x3 <- x2*x
	xp12 <- (1+x)*(1+x)
	val <-  -A/x2 - B/xp12 - 2*C/x3
	return(val)
}

################################################################################

# Calculate the Hessian of the q-density for h evaluated at x (where h is the inverse of g)

hessian.qh <- function(x,U,B,C) 
{	
	x2 <- x*x
	xp12 <- (1+x)*(1+x)
	return( -U/x2 - B/xp12 )	
}

################################################################################

# Calculate the Laplace approximation of the normalizing constant for the
# q-density for g

Z.g.Laplace <- function(A,B,C)
{
	res <- mode.qg(A,B,C)
	g.hat <- res[1]
	qg.hat <- exp(log.qg(g.hat,A,B,C))
	sigma2.inv <- -hessian.qg(g.hat,A,B,C) 
	sigma2 <- 1/sigma2.inv
	return( sqrt(2*pi*sigma2)*qg.hat )
}

################################################################################

# Calculate the Laplace approximation of the normalizing constant for the
# q-density for h (where h is the inverse of g)

Z.h.Laplace <- function(U,B,C)
{
	res <- mode.qh(U,B,C)
	h.hat <- res[1]
	qh.hat <- exp(log.qh(h.hat,U,B,C))
	sigma2.inv <- -hessian.qh(h.hat,U,B,C) 
	sigma2 <- 1/sigma2.inv
	return( sqrt(2*pi*sigma2)*qh.hat )
}

################################################################################

# Approximate the expected value of g with respect to q(g) using the mode

E.g.Laplace <- function(A,B,C) {
	res <- mode.qg(A,B,C)
	g.hat <- res[1]
	return(g.hat)
}

################################################################################

# Approximate the expected value of h with respect to q(h) using the mode

E.h.Laplace <- function(U,B,C) {
	res <- mode.qh(U,B,C)
	h.hat <- res[1]
	return(h.hat)
}

################################################################################

# Approximate the mth moment of g using the Fully-Exponential-Laplace method

moment.g.FullyExponentialLaplace <- function(m,A,B,C) {
	val1 <- Z.g.Laplace(A+m,B,C)
	val2 <- Z.g.Laplace(A,B,C)
	return(val1/val2)
}

################################################################################

# Approximate the mth moment of h using the Fully-Exponential-Laplace method

moment.h.FullyExponentialLaplace <- function(m,U,B,C) {
	val1 <- Z.h.Laplace(U+m,B,C)
	val2 <- Z.h.Laplace(U,B,C)
	return(val1/val2)
}

################################################################################

# Approximate the normalizing constant for the q-density for g using trapezoidal
# integration

Z.g.trapint <- function(A,B,C,N=10000,PLOT=FALSE)
{
	TOL <- 1.0E-14
	TOR <- 1.0E-5
	
	# Find the mode of g
	res <- mode.qg(A,B,C)
	g.hat <- res[1]
	
	# Evaluate q(g) at the mode
	log.qg.max <- log.qg(g.hat,A,B,C)
	qg.hat <- exp(log.qg.max)
	
	# Use a normal approximation at the mode
	sigma2.inv <- -hessian.qg(g.hat,A,B,C) 
	sigma2 <- 1/sigma2.inv
	delta   <- 0.5*sqrt(sigma2)
	
	# Loop using steps of size delta to the right until the difference in the 
	# log-densities between the maximum and the current location is smaller 
	# than TOR
	g.curr  <- g.hat + delta
	log.qg.curr <- log.qg(g.curr,A,B,C)
	#count <- 0	
	#cat(count,g.hat,g.curr,log.qg.curr,log.qg.max,log.qg.curr - log.qg.max,log(TOL),"\n")
	while ( (log.qg.curr - log.qg.max) > log(TOR) ) {
		#cat(count,g.hat,g.curr,log.qg.curr,log.qg.max,log.qg.curr - log.qg.max,log(TOL),"\n")
		#count <- count + 1
		g.curr <- g.curr + delta
		log.qg.curr <- log.qg(g.curr,A,B,C)
	}
	R <- g.curr
	
	# Loop using steps of size delta to the left until the difference in the 
	# log-densities between the maximum and the current location is smaller 
	# than TOL
	delta <- g.hat/4
	g.curr  <- g.hat/2
	log.qg.curr <- log.qg(g.curr,A,B,C)
	while ( (log.qg.curr - log.qg.max) > log(TOL) ) {
		while ((g.curr - delta)<=0) {
			# Reduce the step size if necessary
			delta <- delta/5
		}
		g.curr <- g.curr - delta
		log.qg.curr <- log.qg(g.curr,A,B,C)
	}
	L <- g.curr
	
	#print(L)
	
	# Calculate a grid between L and R with N points
	vg <- seq(L,R,,N)
	log.qg.vg <- log.qg(vg,A,B,C)
	
	# Use trapezoidal integration
	intVal <- trapint(vg,exp(log.qg.vg))
	
	# Plot the points for diagnostic purposes
	if (PLOT) {
		plot(vg,exp(log.qg.vg),type="l")
		points(g.hat,exp(log.qg.max))
	}
	
	return(list(intVal=intVal,vg=vg,log.qg.vg=log.qg.vg))
}

################################################################################

Z.h.trapint <- function(U,B,C,N=10000,PLOT=FALSE)
{
	TOL <- 1.0E-5
	L <- 1.0E-3
	res <- mode.qh(U,B,C)
	h.hat <- res[1]
	log.qh.max <- log.qh(h.hat,U,B,C)
	log.qh.L <- log.qh(L,U,B,C)
	delta   <- h.hat
	h.curr  <- h.hat
	log.qh.curr <- log.qh(h.curr,U,B,C)
	while ( (log.qh.curr - log.qh.max)>log(TOL) ) {
		h.curr <- h.curr + delta
		log.qh.curr <- log.qh(h.curr,U,B,C)
	}
	vh <- seq(L,h.curr,,N)
	log.qh.vh <- log.qh(vh,U,B,C)
	intVal <- trapint(vh,exp(log.qh.vh))
	
	if (PLOT) {
		plot(vh,exp(log.qh.vh),type="l")
		points(h.hat,exp(log.qh.max))
	}

	return(intVal)
}

################################################################################

# Calculate the mth moment of q(g) using trapezoidal integration

moment.g.trapint <- function(m,A,B,C,N=1000000)
{
	val1 <- Z.g.trapint(A+m,B,C,N)$intVal
	val2 <- Z.g.trapint(A,B,C,N)$intVal
	return(val1/val2)
}

################################################################################

# Calculate the mth moment of q(h) using trapezoidal integration

moment.h.trapint <- function(m,U,B,C,N=1000000)
{
	val1 <- Z.h.trapint(U+m,B,C,N)
	val2 <- Z.h.trapint(U,B,C,N)
	return(val1/val2)
}

################################################################################

library(fAsianOptions)

# Calculate the normalizing constant using the exact result

Z.g.exact <- function(A,B,C) 
{
	nu <- A + 1
	mu <- B + 1
	beta <- C
	
	# Check the conditions under which the result holds
	if ( ((1-mu)>nu)&(nu>0) ) {
		#print("fine")
	} else {
		#print("not fine")
	}
	
	val1 <- 0.5*(nu-1)*log(beta) + lgamma(1 - mu - nu) + 0.5*beta
	val2 <- Re( whittakerW(beta, 0.5*(nu-1) + mu, -0.5*nu) )
	
	return( list(val=exp(val1)*val2,val1=val1,val2=val2) )
}

################################################################################

# Calculate the mth moment of q(g) using the exact result

moment.g.exact <- function(m,A,B,C)
{
	res1 <- Z.g.exact(A+m,B,C)
	res2 <- Z.g.exact(A,B,C)
	return(list(val=res1$val/res2$val,res1=res1,res2=res2))
}

################################################################################

# Use the plug in approximation for tau

tau.plugin  <- function(a,n,p,R2)
{
	b <- (n-p)/2 - a - 2
	A <- b - p/2
	B <- -(n-p)/2
	
	tau <- (1 - R2)*(1 + (0.5*p + a + 1)/b) 
	return(tau)
}

################################################################################

# Apply the update for tau_g using trapezoidal integration for ITER iterations

tau.g.trapint <- function(a,n,p,R2,ITER,N=1000)
{
	tau <- tau.plugin(a,n,p,R2)
	
	b <- (n-p)/2 - a - 2
	A <- b - p/2
	B <- -(n-p)/2
	
	for (i in 1:ITER) {
		C <- 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau)
		tau <- moment.g.trapint(m=-1,A,B,C,N)
	}
	
	return(tau)		
}

################################################################################

# Use the Laplace approximation for tau_g

tau.g.Laplace <- function(a,n,p,R2)
{
	A <- 2*a + p
	res <- A*(1 - R2)/(n - A) 
	return(res)
}

################################################################################

# Apply the update for tau_g using FEL for ITER iterations

tau.g.FullyExponentialLaplace <- function(a,n,p,R2,ITER)
{
	tau <- tau.g.Laplace(a,n,p,R2)
	
	b <- (n-p)/2 - a - 2
	A <- b - p/2
	B <- -(n-p)/2
	
	U <-  -(A+B+2)
	
	for (i in 1:ITER) {
		C <- 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau)
		# cat(i,tau,C,"\n")
		tau <- moment.h.FullyExponentialLaplace(m=1,U,B,C)
		# cat(i,tau,C,"\n")
		if (is.na(tau)) {
			tau <- tau.plugin(a,n,p,R2)
			break;
		}
	}
	
	return(tau)		
}

################################################################################

# Apply the update for tau_g using Laplace's method for ITER iterations

tau.g.IterativeLaplace  <- function(a,n,p,R2,ITER) {

	tau <- tau.plugin(a,n,p,R2)
	
	b <- (n-p)/2 - a - 2
	A <- b - p/2
	B <- -(n-p)/2
	
	U <-  -(A+B+2)
	
	for (i in 1:ITER) {
	C <- 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau)
	tau <- E.h.Laplace(U,B,C)
	}
	
	return(tau)	
}	


####################################################################

source("ZE.R")

# Perform the fully Bayesian analysis
res <- ZE.exact.Rcpp(vy,mX,LARGEP=FALSE) 

if (SETTING == 1) {
	write.table(vy, file = "Hitters_vy.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(mX, file = "Hitters_mX.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vR2, file = "Hitters_vR2.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vlog.ZE, file = "Hitters_vlog_ZE.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$mA, file = "Hitters_mGraycode.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 2) {
	write.table(vy, file = "bodyfat_vy.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(mX, file = "bodyfat_mX.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vR2, file = "bodyfat_vR2.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vlog.ZE, file = "bodyfat_vlog_ZE.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$mA, file = "bodyfat_mGraycode.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 3) {
	write.table(vy, file = "Wage_vy.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(mX, file = "Wage_mX.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vR2, file = "Wage_vR2.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vlog.ZE, file = "Wage_vlog_ZE.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$mA, file = "Wage_mGraycode.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 4) {
	write.table(vy, file = "GradRate_vy.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(mX, file = "GradRate_mX.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vR2, file = "GradRate_vR2.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vlog.ZE, file = "GradRate_vlog_ZE.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$mA, file = "GradRate_mGraycode.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 5) {
	write.table(vy, file = "USCrime_vy.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(mX, file = "USCrime_mX.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vR2, file = "USCrime_vR2.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$vlog.ZE, file = "USCrime_vlog_ZE.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(res$mA, file = "USCrime_mGraycode.csv", col.names=FALSE, row.names = FALSE, sep=",")
}


cat("vy sum", sum(vy), "\n")
cat("vR2 sum", sum(res$vR2), "\n")
for (col_idx in 1:ncol(mX)) {
	col <- mX[, col_idx]
	cat(sum(mX[, col_idx]), " ")
}
cat("\n")
logpy <- res$vlog.ZE
logpy.til <- logpy - max(logpy)

# Calculate the marginal variable inclusion probability
vp <- exp(logpy.til)/sum(exp(logpy.til))
 
a <- -0.75

velbo <- c(-Inf)
for (i in 1:length(res$vR2)) 
{
	# Note res$mA contains the gray-code
	p <- sum(res$mA[i,])
	b <- (n - p)/2 - a - 2 
	
	# Calculate
	A <-  n/2 - p - a - 2
	B <- -(n-p)/2
	U <- -(A+B+2)
	
	# Calculate tau_g
	tau.g <- tau.g.FullyExponentialLaplace(a,n,p,res$vR2[i],ITER=20)
	
	# Calculate the constant C (needed to calculate the normalizing constant for q(g)	
	C <- 0.5*n*res$vR2[i]/((1 + tau.g)*(1 - res$vR2[i] + tau.g)) + 0.5*p/(1 + tau.g)
	
	# Calculate the
	Z <- Z.g.trapint(A,B,C,N=1000,PLOT=FALSE)$intVal    # Z.h.Laplace(U,B,C)
	
	# Calculate the lower bound on the log-likelihood	
	velbo[i] <- 0.5*p - 0.5*n*log(2*pi) - lbeta(a+1,b+1)  - 0.5*n*log(1 + tau.g - res$vR2[i]) 
	velbo[i] <- velbo[i] - 0.5*(n+p)*log(0.5*(n+p))+ lgamma(0.5*(n+p)) + C*tau.g + log(Z)  + 0.5*(n-p)*log(1 + tau.g)
	
	# cat(i,velbo[i],tau.g,C,Z,"\n")
		
	# If there is an error stop here and have a look
	if (is.na(Z)) {
		print("error! press escape and have a look")
		ans <- readline()
	}

}

velbo[is.na(velbo)] <- -Inf

logqy.til <- velbo - max(velbo)
vq <- exp(logqy.til)/sum(exp(logqy.til))

# Calculate the variable inclusion probabilities
vw1 <- round( t(vp)%*%res$mA, 3)
vw2 <- round( t(vq)%*%res$mA, 3)
cat(vw1)
cat(vw2)

mW <- rbind(vw1,vw2)
colnames(mW) <- varnames

# Put the results in a table
tab <- cbind( round(100*vp,3),round(100*vq,3),res$mA )
colnames(tab) <- c("prob p","prob q",varnames)
ord <- order(vp,decreasing=TRUE)

# Print the table for the top 20 models
print( tab[ord[1:20],] ) 

if (SETTING == 1) {
	write.table(vw1, file = "Hitters_vw1.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(vw2, file = "Hitters_vw2.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 2) {
	write.table(vw1, file = "bodyfat_vw1.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(vw2, file = "bodyfat_vw2.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 3) {
	write.table(vw1, file = "Wage_vw1.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(vw2, file = "Wage_vw2.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 4) {
	write.table(vw1, file = "GradRate_vw1.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(vw2, file = "GradRate_vw2.csv", col.names=FALSE, row.names = FALSE, sep=",")
}

if (SETTING == 5) {
	write.table(vw1, file = "USCrime_vw1.csv", col.names=FALSE, row.names = FALSE, sep=",")
	write.table(vw2, file = "USCrime_vw2.csv", col.names=FALSE, row.names = FALSE, sep=",")
}
