

rm(list = ls())

##########################################################################3

SETTING <- 5

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

source("~/Dropbox/phd/code/functions.Rs")

res = normalize(vy,mX) 
vy = res$vy
mX = res$mX

################################################################################

library(BAS)

doBAS = TRUE
if (doBAS) {

	print("BAS")
	
	a2 <- proc.time()[3]
	crime.bic =  bas.lm(vy~mX, data=UScrime,prior="BIC",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	a2 <- proc.time()[3]
	crime.aic =  bas.lm(vy~mX, data=UScrime,prior="AIC",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	
	a2 <- proc.time()[3]
	crime.g =  bas.lm(vy~mX, data=UScrime,prior="g-prior",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	
	a2 <- proc.time()[3]
	crime.hg =  bas.lm(vy~mX, data=UScrime,prior="hyper-g",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	a2 <- proc.time()[3]
	crime.hgl =  bas.lm(vy~mX, data=UScrime,prior="hyper-g-laplace",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	
	a2 <- proc.time()[3]
	crime.hgn =  bas.lm(vy~mX, data=UScrime,prior="hyper-g-n",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	a2 <- proc.time()[3]
	crime.ZSN =  bas.lm(vy~mX, data=UScrime,prior="ZS-null",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	a2 <- proc.time()[3]
	crime.ZSF =  bas.lm(vy~mX, data=UScrime,prior="ZS-full",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	a2 <- proc.time()[3]
	crime.EBL =  bas.lm(vy~mX, data=UScrime,prior="EB-local",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	a2 <- proc.time()[3]
	crime.EBG =  bas.lm(vy~mX, data=UScrime,prior="EB-global",modelprior=uniform(), initprobs="uniform") 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	tab = cbind(crime.bic$probne0,crime.aic$probne0,crime.g$probne0,crime.hg$probne0,crime.hgl$probne0,crime.hgn$probne0,crime.ZSN$probne0,crime.ZSF$probne0,crime.EBL$probne0,crime.EBG$probne0)
	colnames(tab) = c("BIC","AIC","G","HG","HGL","HGN","ZSN","ZNF","EBL","EBG")
	
	print(round(100*tab,1))
	
}

################################################################################

doBVS <- TRUE
if (doBVS) {
	print("BayesVarSel")
	library( BayesVarSel )
	
	source("Bvs.r")
	
	a6 <- proc.time()[3]
	bvs.rob <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="Robust",prior.models="Constant",time.test= FALSE, priorprobs=NULL,n.keep=10000)
	b6 <- proc.time()[3]  
	t6 <- b6-a6      
	print(t6) 
	
	a7 <- proc.time()[3]
	bvs.liang <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="Liangetal",prior.models="Constant", priorprobs=NULL,n.keep=10000)
	b7 <- proc.time()[3]   
	t7 <- b7-a7       
	print(t7) 
	
	a7 <- proc.time()[3]
	bvs.gz <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="gZellner",prior.models="Constant", priorprobs=NULL,n.keep=10000)
	b7 <- proc.time()[3]   
	t7 <- b7-a7       
	print(t7) 
	
	a7 <- proc.time()[3]
	bvs.zs <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="ZellnerSiow",prior.models="Constant", priorprobs=NULL,n.keep=10000)
	b7 <- proc.time()[3]   
	t7 <- b7-a7       
	print(t7) 
	
	a7 <- proc.time()[3]
	bvs.fls <- my.Bvs(formula="vy~.",fixed.cov=c("Intercept"),data=data.frame(vy=vy,mX=mX),prior.betas="FLS",prior.models="Constant", priorprobs=NULL,n.keep=10000)
	b7 <- proc.time()[3]   
	t7 <- b7-a7       
	print(t7) 
	

	tab2 = cbind(bvs.rob$inclprob,bvs.liang$inclprob,bvs.gz$inclprob,bvs.zs$inclprob,bvs.fls$inclprob)
	colnames(tab2) = c("Robust","Liang","GZ","ZS","FLS")
	
	print(round(100*tab2,1))
}

################################################################################

doBMS <- TRUE
if (doBMS) {
	print("BMS")
	library(BMS)
	
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="UIP", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.UIP   <- coef.bms[,1]
	
	#######################
	
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="BRIC", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.BRIC   <- coef.bms[,1]
	
	
	#######################
	
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="RIC", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.RIC   <- coef.bms[,1]
	
	#######################
	
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="HQ", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.HQ   <- coef.bms[,1]
	
	
	#######################
	
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="EBL", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.EBL   <- coef.bms[,1]
	
	#######################
	
	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="hyper=3", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.hyper <- coef.bms[,1]
	
	#######################

	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="hyper=UIP", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.HUIP <- coef.bms[,1]	
	
	#######################

	a4 <- proc.time()[3]
	res.bms <- bms(cbind(vy,mX), nmodel = 10000, mcmc = "enumerate", g="hyper=BRIC", g.stats = TRUE, mprior="uniform")	
	b4 <- proc.time()[3]   
	t4 <- b4-a4 
	cat(t4,"\n")
	
	fit        <- predict(res.bms, exact=TRUE)
	coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	prob.HBRIC <- coef.bms[,1]	
	
	
	tab3 = cbind(prob.UIP,prob.BRIC,prob.RIC,prob.HQ,prob.EBL,prob.hyper,prob.HUIP,prob.HBRIC)
	colnames(tab3) = c("UIP","BRIC","RIC","HQ","EBL","hyper","HUIP","HBRIC")
	
	print(round(100*tab3,1))

}


################################################################################

library(correlation)

if (TRUE) {

	print("AllLinearModels")
		
	calcVarProbs <- function(score,mGamma) {
		model.prob = exp(score - max(score))/sum(exp(score - max(score)))
		var.prob = t(mGamma)%*%model.prob
		return(var.prob)
	}
	
	###################
	
	a2 <- proc.time()[3]

	vR2 = correlation::all_correlations_mX(vy, mX)
	M = length(vR2)
	mGamma = graycode(p,0)
	vq = mGamma%*%matrix(1,p,1)
	
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")	
	
	###################
	
	a2 <- proc.time()[3]
	vBIC = n*log(1 - vR2) + vq*log(n) 
	var.prob1 = calcVarProbs(-0.5*vBIC,mGamma) 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")	
	
	###################
	
	a2 <- proc.time()[3]
	a  = -0.75
	vb = 0.5*(n - vq - 5) - a
	c  = 0.5*(n - 1)
	vd = 0.5*vq + a
	
	log.vp = -(vb+1)*log(1 - vR2) + lbeta(vd+1,vb+1) - lbeta(a+1,vb+1)
	vZE <- -2*log.vp
	
	var.prob2 = calcVarProbs(log.vp,mGamma) 
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	###################	
	
	log.hyperg_2F1 = function(b,c,x) {
		val = 0 
		val = val + log(c-1) 
		val = val + (1-c)*log(x) 
		val = val + (c-b-1)*log(1-x) 
		val = val + lbeta(c-1,b-c+1) 
		val = val + pbeta(x,shape1=(c-1),shape2=(b-c+1),log=TRUE)
		return(val)
	}
	
	log.hyperg_2F1.naive = function(b,c,x) {
		val = log( hyperg_2F1( b, 1, c, x, give=FALSE, strict=TRUE) )
		return(val)
	}	
	
	a2 <- proc.time()[3]
	
	a = 3
	log.vp.g = log(a - 2) - log(vq + a - 2) + log( hyperg_2F1(0.5*(n-1), 1, 0.5*(vq+a), vR2, give=FALSE, strict=TRUE) )
	var.prob3 = calcVarProbs(log.vp.g,mGamma) 
	
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
 
	###################	
	
	a2 <- proc.time()[3]
	
	a = 3
	log.vp.g2 = log(a - 2) - log(vq + a - 2) + log.hyperg_2F1( 0.5*(n-1), 0.5*(vq+a), vR2 )
	var.prob4 = calcVarProbs(log.vp.g2,mGamma) 
	
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	###################	
	
	#library(appell)
	#M = nrow(res$mGamma)
	#log.vp.gprior3 = rep(0,M)
	#for (i in 2:M) {
	#	log.vp.gprior3[i] = log(a - 2) - log(res$vq[i] + a - 2) + log( Re(appellf1(1, 0.5*a, 0.5*(n-1), 0.5*(res$vq[i]+a), (1 - 1/n), res$vR2[i] ,hyp2f1 = "forrey")$val)  )
	#}
	#log.vp.gprior3[1] =0
	#model.prob3 = exp(log.vp.gprior3 - max(log.vp.gprior3))/sum(exp(log.vp.gprior3 - max(log.vp.gprior3)))
	#var.prob3 = t(mGamma)%*%model.prob3

	###########################

	a2 <- proc.time()[3]
	
	log.vp.gprior5 = rep(0,M)
	log.vp.gprior5[-1] = log.vp.gprior5[-1] - 0.5*vq[-1]*log(n+1)
	log.vp.gprior5[-1] = log.vp.gprior5[-1] + 0.5*vq[-1]*log(vq[-1]+1)
	log.vp.gprior5[-1] = log.vp.gprior5[-1] - 0.5*(n - 1)*log(vR2[-1])
	log.vp.gprior5[-1] = log.vp.gprior5[-1] - log(vq[-1]+1)
	log.vp.gprior5[-1] = log.vp.gprior5[-1] + log( hyperg_2F1( 0.5*(vq[-1]+1), 0.5*(n-1), 0.5*(vq[-1]+3), (1-1/vR2[-1])*(vq[-1]+1)/(n+1), give=FALSE, strict=TRUE) )
	
	model.prob5 = exp(log.vp.gprior5 - max(log.vp.gprior5))/sum(exp(log.vp.gprior5 - max(log.vp.gprior5)))
	var.prob5 = t(mGamma)%*%model.prob5
	
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	###########################	
	
	a2 <- proc.time()[3]
	
	vL = (1 + n)/(1 + vq) - 1
	vsigma2 = 1 - vR2
	vz = vR2/(1 + vL*vsigma2)
	
	log.vp.gprior6 = rep(0,M)
	log.vp.gprior6[-1] = log.vp.gprior6[-1] + 0.5*(n - vq[-1] - 1)*log( n + 1 )
	log.vp.gprior6[-1] = log.vp.gprior6[-1] - 0.5*(n - vq[-1] - 1)*log( vq[-1] + 1)
	log.vp.gprior6[-1] = log.vp.gprior6[-1] - 0.5*(n - 1)*log(1 + vL[-1]*vsigma2[-1])
	log.vp.gprior6[-1] = log.vp.gprior6[-1] - log (vq[-1] + 1)
	log.vp.gprior6[-1] = log.vp.gprior6[-1] + log.hyperg_2F1( 0.5*(n-1), 0.5*(vq[-1]+3), vz[-1] )
	
	model.prob6 = exp(log.vp.gprior6 - max(log.vp.gprior6))/sum(exp(log.vp.gprior6 - max(log.vp.gprior6)))
	var.prob6 = t(mGamma)%*%model.prob6
	
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")
	
	###########################	
	
	a2 <- proc.time()[3]
	
	log.vp.gprior7 = rep(0,M)
	log.vp.gprior7[-1] = log.vp.gprior7[-1] + 0.5*(n - vq[-1] - 1)*log( n + 1 )
	log.vp.gprior7[-1] = log.vp.gprior7[-1] - 0.5*(n - vq[-1] - 1)*log( vq[-1] + 1)
	log.vp.gprior7[-1] = log.vp.gprior7[-1] - 0.5*(n - 1)*log(1 + vL[-1]*vsigma2[-1])
	log.vp.gprior7[-1] = log.vp.gprior7[-1] - log (vq[-1] + 1)
	log.vp.gprior7[-1] = log.vp.gprior7[-1] + log( hyperg_2F1( 0.5*(n-1), 1, 0.5*(vq[-1] + 3), vz[-1], give=FALSE, strict=TRUE) )
 
	model.prob7 = exp(log.vp.gprior7 - max(log.vp.gprior7))/sum(exp(log.vp.gprior7 - max(log.vp.gprior7)))
	var.prob7 = t(mGamma)%*%model.prob7
	
	b2 <- proc.time()[3]    
	t2 <- b2-a2
	cat(t2,"\n")	
	
	###########################	

	tab4 <- cbind(var.prob1,var.prob2,var.prob3,var.prob4,var.prob5,var.prob6,var.prob7)
	colnames(tab4) = c("BIC","ZE","g.naive","g.safe","Robust.naive","Robust","Robust.safe")
	print(round(100*tab4,1))
	
	
	
}





#dyn.load("ZS.approx.null.np.dll")

#LogBF_ZS_null_vect <- function(r2curr, n, dim) {
#	nmodels =length(r2curr)
#	res <- .C("R_LogBF_ZS_null_vect",r2curr=r2curr, n=as.integer(n), dim=as.integer(dim), nmodels=as.integer(nmodels), logmarg=as.double(rep(0,nmodels)))
#	return(res$logmarg)
#}

#vals = LogBF_ZS_null_vect(res$vR2, n, res$vq)
#model.prob8 = exp(vals - max(vals))/sum(exp(vals - max(vals)))
#var.prob8 = t(res$mGamma)%*%model.prob8



