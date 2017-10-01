

rm(list = ls())

##########################################################################3

SETTING <- 4

if (SETTING==1) {
	library(ISLR)
	Hitters <- na.omit(Hitters)
	
	# Get y vector and X matrix
	y.t <- Hitters$Salary
	X.f <- model.matrix(Salary~.,Hitters)[,-1]
 
	varnames <- colnames(X.f)
} 

if (SETTING==2) {
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

if (SETTING==3) {
	library(ISLR)
	Wage <- na.omit(Wage) 
	
	
	# Get y vector and X matrix
	y.t <- Wage$wage
	X.f <- model.matrix(wage~.,Wage)[,-1]
	
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

if (SETTING==4) {
	library(ISLR)
	College <- na.omit(College)
	
	# Get y vector and X matrix
	y.t <- College$Grad.Rate
	X.f <- model.matrix(Grad.Rate~.,College)[,-1]
 
	varnames <- colnames(X.f)
}

if (SETTING==5) {

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

mR <- cor(mX) - diag(1,p)

print(range(eigen(mR)$values))

vm <- matrix(0,p,1)

for (i in 1:100) 
{	
vm <- vm + t(mX)%*%(vy - mX%*%vm)/n	
}

print("done")

print(cbind(lm(vy~mX)$coef[-1],vm))

################################################################################

library(correlation)

a0 <- proc.time()[3]
res.fast <- all_correlations_mX(vy, mX, intercept_col = 0, bIntercept = FALSE, bCentre = FALSE)$vR2
b0 <- proc.time()[3]    
t0 <- b0-a0
cat(t0,"\n")
cat((2^p)/t0,"\n")

mA <- graycode(p,0)

vlog.BIC <- -0.5*n*log(1 - res.fast) - 0.5*apply(mA,1,sum)*log(n) 
vlog.BIC.til <- vlog.BIC - max(vlog.BIC)
vp <- exp(vlog.BIC.til)/sum(exp(vlog.BIC.til))
pip.BIC <- rev( t(mA)%*%vp)


ans <- readline()

################################################################################

source("mfvb.r")

a1 <- proc.time()[3]
vlambda <- seq(0,12,,25)
nTrials <- 1
nFolds  <- 5

library(cvTools)
cvSets <- cvFolds(n, K = nFolds, R = nTrials)

vmse <- matrix(0,nTrials,length(vlambda))
for (trial in 1:nTrials) 
{
	for (fold in 1:nFolds) 
	{
		inds.test <- cvSets$subsets[cvSets$which==fold,trial]
		inds.train <- (1:n)[-inds.test ]
		n.test   <- length(inds.test)
		n.train <- n - n.test
		
		vy.test  <- vy[inds.test]
		vy.train <- vy[inds.train]
		mX.test  <- mX[inds.test,]
		mX.train <- mX[inds.train,]
		
		
		for (i in 1:length(vlambda)) 
		{
			lambda <- vlambda[i]
			res <- mfvb(vy.train,mX.train,tau=1,lambda)
			vmse[trial,i] <- mean(  (vy.test - mX.test%*%(res$vw*res$vm))^2 )	
		}
	}

}

for (trial in 1:nTrials) 
{
	if (trial==1) {
		plot(vlambda,vmse[trial,],type="l",ylim=range(vmse))
	} else {
		lines(vlambda,vmse[trial,],col=trial)
	}
}

mse <-  apply( vmse,2,median)
lines(vlambda,mse,col=1,lwd=3)

	
res <- mfvb(vy,mX,tau=0.01,lambda=vlambda[which.min(mse)])
b1 <- proc.time()[3]    
t1 <- b1-a1
cat(t1,"\n")

#############

if (FALSE) {
	
	source("ZE.R")
	
	a2 <- proc.time()[3]
	res.ZE <- ZE.exact.Rcpp(vy,mX,LARGEP=FALSE)
	
	vlog.ZE.til <- res.ZE$vlog.ZE - max(res.ZE$vlog.ZE)
	vp <- exp(vlog.ZE.til)/sum(exp(vlog.ZE.til))
	pip.ZE <- rev( t(res.ZE$mA)%*%vp)
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")	
	
	XTX <- t(mX)%*%mX/n
	vlogdet <- c()
	for (i in 1:nrow(res.ZE$mA)) {
		vgamma <- res.ZE$mA[i,]
		inds <- which(vgamma==1)
		if (length(inds)==0) {
			vlogdet[i] <- 0
		}
		if (length(inds)==1) {
			vlogdet[i] <- log(XTX[inds,inds])
		}
		if (length(inds)>1) {
			vlogdet[i] <- log(det(XTX[inds,inds]))
		}
	}
	
    a3 <- proc.time()[3]
    res.BIC <- ZE.exact.Rcpp(vy,mX,LARGEP=FALSE)
	vlog.BIC <- -0.5*n*log(1 - res.ZE$vR2) - 0.5*apply(res.ZE$mA,1,sum)*log(n) # - 0.5*vlogdet
	vlog.BIC.til <- vlog.BIC - max(vlog.BIC)
	vp.BIC <- exp(vlog.BIC.til)/sum(exp(vlog.BIC.til))
	pip.BIC <- rev( t(res.ZE$mA)%*%vp.BIC)
	b3 <- proc.time()[3]   
	t3 <- b3-a3
	cat(t3,"\n")	
}


################################################################################

doBMS <- TRUE
if (doBMS) {
	print("BMS")
	library(BMS)
	
	start.value <- which(res$vw>0.5)
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), burn = 100, iter = 100000, nmodel = 10000, mcmc = "bd",
	g="hyper=3", mprior = "random", mprior.size = NA, user.int = TRUE,
	start.value = start.value, g.stats = TRUE, logfile = TRUE, logstep = 100, 
	force.full.ols = FALSE, fixed.reg=numeric(0))	
	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	postmean   <- coef.bms[,2]
	prob0      <- coef.bms[,1]
	ind        <- coef.bms[,5]
	vgamma.hat      <- rep(0,p-1)
	vgamma.hat[ind] <- round(prob0)
}

do.mombf <- FALSE
if (do.mombf) {
	library(mombf)
	
	a5 <- proc.time()[3]
	res.mombf <- modelSelection(vy, mX, center=TRUE, scale=TRUE, niter=10000, thinning=1,
	burnin=100, family='normal', priorCoef=momprior(tau=0.348),
	priorDelta=modelbbprior(alpha.p=1,beta.p=1),
	priorVar=igprior(alpha=.01,lambda=.01),
	priorSkew=momprior(tau=0.348), deltaini=rep(FALSE,ncol(mX)),
	initSearch='none', method='auto', optimMethod='LMA', B=10^5, verbose=TRUE)
	b5 <- proc.time()[3] 
	t5 <- b3-a5   
	cat(t5,"\n")
}

doBVS <- TRUE
if (doBVS) {
	print("BayesVarSel")
	library( BayesVarSel )
	
	source("Bvs.r")
	
	
	a6 <- proc.time()[3]
	res.bvs <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="FLS",
										prior.models="Constant",time.test= FALSE, priorprobs=NULL, n.keep=200)
	b6 <- proc.time()[3]  
	t6 <- b6-a6      
	print(t6) 
	
	
	a7 <- proc.time()[3]
	res.bvs2 <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="Robust",
											prior.models="Constant",time.test= FALSE, priorprobs=NULL, n.keep=200)
	b7 <- proc.time()[3]   
	t7 <- b7-a7       
	print(t7) 
	
	
	#a3 <- proc.time()[3]
	#res.bvs3 <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="Liang",prior.models="Constant",time.test= FALSE, #priorprobs=NULL)
	#b3 <- proc.time()[3]     
	#print(b3-a3) 
}

 
doVARBVS <- TRUE
if (doVARBVS) 
{
	print("varbvs")
	library(varbvs)
	
	a8 <- proc.time()[3]
	
	# Fit the variable selection model.
	fit <- varbvs(mX,Z=NULL,vy,logodds = seq(-3,-1,0.1))
	print(summary(fit))
	
	# Compute posterior inclusion probabilities (pip) and posterior mean
	# regression coefficients (beta) averaged over the hyperparameter settings.
	w <- normalizelogweights(fit$logw);
	pip <- fit$alpha %*% c(w)
	beta.hat <- fit$mu %*% c(w)
	
	# Compute the posterior mean estimate of hyperparameter sa.
	sa <- sum(fit$sa * w)
	
	# Compare estimated outcomes against observed outcomes.
	y.fit <- predict(fit,mX,Z=NULL)
	print(cor(vy,y.fit))
	
	b8 <- proc.time()[3] 
	t8 <- b8 - a8
	cat("varbvs: ",t8,"\n") 
}

doEMVS <- TRUE
if (doEMVS) {
	print("EMVS")
	library(EMVS)
	
	a7 <- proc.time()[3]
	
	beta_init=as.numeric(res$vw>0.5)
	
	v0=exp(seq(log(0.0001),log(10),,200))
	v1=1000
	a=b=1
	epsilon=10E-3
	sigma_init=1   
	
	res.emvs <- EMVS(vy, mX, v0=v0,v1=v1,type="betabinomial",beta_init=beta_init,sigma_init=1,epsilon=epsilon,a=a,b=b)
	
	posts = res.emvs$log_g_function
	best  = which.max(posts)
	beta.emvs = res.emvs$betas[best, ]
	gamma.emvs <- res.emvs$prob_inclusion[best, ]
	vy.emvs <- mX%*%beta.emvs
	
	b7 <- proc.time()[3] 
	t7 <- b7 - a7
	cat("EMVS: ",t7,"\n") 
}


doCVA <- TRUE
if (doCVA) {
	print("CVA")
	library(correlation)

	binary_to_model <- function(binary_vec)
	{
	  acc <- 0
	  mul <- 1
	  for (i in 1:length(binary_vec)) {
	    acc <- acc + mul * binary_vec[i]
	    mul <- mul * 2
	  }
	  return(acc)
	}

	log_p <- function(n, p, vR2, vp_gamma)
	{
	  R2 <- R2[2:length(R2)]
	  p_gamma <- p_gamma[2:length(p_gamma)]
	  a <- 1
	  b <- p
	  return(-n / 2 * log(1 - R2) - p_gamma / 2 * log(n) + lbeta(a + p_gamma, b + p - p_gamma))
	}

  corr_result <- correlation::all_correlations_mX(vy, mX, 1, TRUE)
	set.seed(1)
  K <- 400
  posterior_percentage <- rep(NA, 10)
  t8 <- rep(NA, 10)
  for (i in 1:10) {
	  initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)

		a8 <- proc.time()[3]

	  cva_result <- cva(initial_gamma, vy, mX, K, 1.0)

		# posts = res.emvs$log_g_function
		# best  = which.max(posts)
		# beta.emvs = res.emvs$betas[best, ]
		# gamma.emvs <- res.emvs$prob_inclusion[best, ]
		# vy.emvs <- mX%*%beta.emvs
		
		b8 <- proc.time()[3] 
		t8[i] <- b8 - a8
		cat("CVA: ",t8[i],"\n")

	  R2 <- corr_result$vR2
	  p_gamma <- corr_result$vp_gamma
	  cva_models <- apply(cva_result$models, 1, binary_to_model)

	  vlog_p <- log_p(n, p, R2, p_gamma)
	  vlog_p.til <- vlog_p - max(vlog_p)
	  vp <- exp(vlog_p.til)/sum(exp(vlog_p.til))

	  #plot(vlog_p.til, pch=21, xlab="Model Index", ylab="Log Posterior Model Probability", col="black", bg="black")
	  #points(cva_models, vlog_p.til[cva_models], pch=21, col="red", bg="red")
	  posterior_percentage[i] <- sum(vp[cva_models])
	}
	print(mean(t8))
	print(mean(posterior_percentage))

	CVA_inclprob <- apply(vp[cva_models] * cva_result$models, 2, sum)/sum(vp[cva_models])
	CVA_inclprob <- apply(cva_result$models, 2, sum)/K
}


tab <- round(100*cbind(pip.BIC,res.bvs$inclprob,res.bvs2$inclprob,prob0,pip,gamma.emvs, CVA_inclprob),1)
colnames(tab) <- c("BIC","BVS-FLS","BVS-Robust","BMS","varbvs","EMVS","CVA")

#tab <- round(100*cbind(pip.ZE,pip.BIC,prob0,pip,res$vw,gamma.emvs),1)
#colnames(tab) <- c("ZE","BIC","BMS","varbvs","MFVB","EMVS")

print(tab )



