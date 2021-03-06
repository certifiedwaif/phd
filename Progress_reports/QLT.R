################################################################################
 
doBAS <- FALSE
doBMS <- TRUE
doEMVS <- TRUE
doVARBVS <- TRUE

library(ncvreg)
if (doBAS) {
	library(BAS)
}
if (doBMS) {	
	library(BMS)
}
library(Matrix)
library(EMVS)
library(varbvs)
 
################################################################################

dataset <- "QTL"
print(dataset)

response <- read.table("phe_simulat.csv",header=FALSE,sep=",")
covariates <- read.table("gen_simulat.csv",header=FALSE,sep=",")

n <- nrow(response)
P <- ncol(covariates)

vbeta <- c()
mX <- matrix(0,n,7381)
count <- 1
for (i in 1:P) {
    for (j in i:P) {
        if (i==j) {
            mX[,count] <- covariates[,i] 
        } else {
            mX[,count] <- covariates[,i]*covariates[,j]
        }
        
        vbeta[count] <- 0
        
        if ((i==1)&(j==1))     { vbeta[count] <- 4.47; }
        if ((i==21)&(j==21))   { vbeta[count] <- 3.16; }
        if ((i==31)&(j==31))   { vbeta[count] <- 2.24; }
        if ((i==51)&(j==51))   { vbeta[count] <- 1.58; }
        if ((i==71)&(j==71))   { vbeta[count] <- 1.58; }
        if ((i==91)&(j==91))   { vbeta[count] <- 1.10; }
        if ((i==101)&(j==101)) { vbeta[count] <- 1.10; }
        if ((i==111)&(j==111)) { vbeta[count] <- 0.77; }
        if ((i==121)&(j==121)) { vbeta[count] <- 0.77; }
        if ((i==1)&(j==11))    { vbeta[count] <- 1.00; }
        if ((i==2)&(j==119))   { vbeta[count] <- 3.87; }
        if ((i==10)&(j==91))   { vbeta[count] <- 1.30; }
        if ((i==15)&(j==75))   { vbeta[count] <- 1.73; }
        if ((i==20)&(j==46))   { vbeta[count] <- 1.00; }
        if ((i==21)&(j==22))   { vbeta[count] <- 1.00; }
        if ((i==26)&(j==91))   { vbeta[count] <- 1.00; }
        if ((i==41)&(j==61))   { vbeta[count] <- 0.71; }
        if ((i==56)&(j==91))   { vbeta[count] <- 3.16; }
        if ((i==65)&(j==85))   { vbeta[count] <- 2.24; }
        if ((i==86)&(j==96))   { vbeta[count] <- 0.89; }
        if ((i==101)&(j==105)) { vbeta[count] <- 1.00; }
        if ((i==111)&(j==121)) { vbeta[count] <- 2.24; }
       
        count <- count + 1
    }
}

mX.intercept <- cbind(1,mX)
p <- ncol(mX.intercept)

################################################################################

ebic.ncvreg <- function(vy,mX,penalty)
{
	res <- ncvreg(mX,vy,penalty=penalty)
	ebic <- c()
	for (i in 1:ncol(res$beta)) {
		dof  <- sum(res$beta[,i]!=0)
		rss  <- res$loss[i]
		ebic[i] <- log(rss/n) + (dof/n)*(log(n) + 2*log(p))
	}
	best <- which.min(ebic)
	return(list(beta=res$beta[,best],rss=res$loss[best]))
}
################################################################################

transfer <- function(mX,vbeta){ #mX without intcpt
  vmean <- apply(mX,2,mean)
  vsd   <- apply(mX,2,sd)
  new.vbeta   <- vbeta / vsd
  new.intcpt  <- -(vmean / vsd) %*% vbeta
  return (c(new.intcpt,new.vbeta)) 
}

################################################################################

EMVSbest.my <- function (result) 
{
    posts = result$log_g_function
    posts[is.finite(posts) == F] = NaN
    which <- which.max(posts)
    logpost <- posts[which]
    print("Best Model Found")
    list(log_g_function = logpost, indices = (1:dim(result$betas)[2])[result$prob_inclusion[which, 
        ] > 0.5], beta.hat = result$betas[which, ], gamma.hat <- result$prob_inclusion[which, ])
}
################################################################################

CalcSelectionScores <- function(vgamma,vgamma.hat) 
{
	pos <- which(vgamma==1)
	neg <- which(vgamma==0)

	TP <- sum(vgamma.hat[pos])
	TN <- sum(1 - vgamma.hat[neg])
	FP <- sum(vgamma.hat[neg])
	FN <- sum(1 - vgamma.hat[pos])
	
	sensitivity <- TP/length(pos)
	specificity <- TN/length(neg)
	precision <- TP/sum(vgamma.hat)
	recall <- TP/(TP + FN)
	accuracy <- (TP + TN)/(TP + TN + FP + FN)
		
	F1 <- 2*precision*recall/(precision+recall)		
	BER <- 0.5*(FP/(FP+TN) + FN/(FN+TP))
	
	return(list(TP=TP,TN=TN,FP=FP,FN=FN,
		sensitivity=sensitivity,
		specificity=specificity,
		precision=precision,
		recall=recall,
		accuracy=accuracy,
		F1=F1,
		BER=BER))
}

 
################################################################################

TRIALS <- 5
sigma2.true <- 20
percent <- 0.7

mX.til <- cbind(1,mX)

t.lasso <- rep(0,TRIALS)
t.mcp   <- rep(0,TRIALS)
t.scad  <- rep(0,TRIALS)
t.emvs  <- rep(0,TRIALS)
t.bas   <- rep(0,TRIALS)
t.vb.screening <- rep(0,TRIALS)
t.vb    <- rep(0,TRIALS)
t.varbvs <- rep(0,TRIALS)
t.bms <- rep(0,TRIALS)
t.cva <- rep(0,TRIALS)

SCORES.lasso <- NULL
SCORES.mcp   <- NULL
SCORES.scad  <- NULL
SCORES.emvs  <- NULL
SCORES.bas   <- NULL
SCORES.bms   <- NULL
SCORES.vb.screening <- NULL
SCORES.vb    <- NULL
SCORES.varbvs    <- NULL
SCORES.bms    <- NULL
SCORES.screen <- NULL
SCORES.HOLP   <- NULL
SCORES.both <- NULL
SCORES.cva <- NULL

MSE.lasso <- rep(0,TRIALS)
MSE.mcp   <- rep(0,TRIALS)
MSE.scad  <- rep(0,TRIALS)
MSE.emvs  <- rep(0,TRIALS)
MSE.bas   <- rep(0,TRIALS)
MSE.bms   <- rep(0,TRIALS)
MSE.vb.screening <- rep(0,TRIALS)
MSE.vb    <- rep(0,TRIALS)
MSE.varbvs <- rep(0,TRIALS)
MSE.bms    <- rep(0,TRIALS)
MSE.cva    <- rep(0,TRIALS)
 

bias.lasso <- rep(0,TRIALS)
bias.mcp   <- rep(0,TRIALS)
bias.scad  <- rep(0,TRIALS)
bias.emvs  <- rep(0,TRIALS)
bias.bas   <- rep(0,TRIALS)
bias.bms   <- rep(0,TRIALS)
bias.vb.screening <- rep(0,TRIALS)
bias.vb    <- rep(0,TRIALS)
bias.varbvs<- rep(0,TRIALS)
bias.bms    <- rep(0,TRIALS)
bias.cva    <- rep(0,TRIALS)

vnum <- c()

RESTART <- TRUE
if (RESTART) {
	start <- 1
} else {
	load("QLT.Rdata")
	start <- trials + 1
}

doBAS <- FALSE
doBMS <- TRUE
doEMVS <- TRUE
doVARBVS <- TRUE
doVB <- TRUE
doVBscreen <- TRUE
doCVA <- TRUE


for (trials in start:TRIALS) 
{
	#############################################################################
	
    set.seed(trials)
		
    true.inds <- which(vbeta!=0)
    cat("True predictors are:\n")
    cat(which(vbeta!=0),"\n")
	
	vf <- mX%*%matrix(vbeta)
	vy <- vf + rnorm(nrow(mX),0,sqrt(sigma2.true)) ##how about intercept??
	vy.cent <- vy - mean(vy)
	mX.std <- mX
	for (j in 1:ncol(mX.std)) {
		mX.std[,j] <- (mX.std[,j] - mean(mX.std[,j]))/sd(mX.std[,j])
	}

	
	vgamma     <- as.numeric(vbeta!=0)
	vbeta.til  <- c(0,vbeta)
	
	#############################################################################
	
	a3 <- proc.time()[3]
	res.lasso <- ebic.ncvreg(vy,mX,penalty="lasso")
	vy.hat <- mX.til%*%res.lasso$beta
	vgamma.hat <- as.numeric(res.lasso$beta[-1]!=0)
	scores.lasso <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
	b3 <- proc.time()[3]     
    print(b3-a3) 
    t.lasso[trials]    <- b3-a3
    MSE.lasso[trials]  <- sum((vf - vy.hat)^2)
    bias.lasso[trials] <- sum((res.lasso$beta-vbeta.til)^2)
    SCORES.lasso       <- cbind(SCORES.lasso,scores.lasso)
 
	#############################################################################
	
	a3 <- proc.time()[3]
	res.mcp <- ebic.ncvreg(vy,mX,penalty="MCP")
	vy.hat <- mX.til%*%res.mcp$beta
	vgamma.hat <- as.numeric(res.mcp$beta[-1]!=0)
	scores.mcp <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
	b3 <- proc.time()[3]     
	print(b3-a3)
	t.mcp[trials]    <- b3-a3    	
    MSE.mcp[trials]  <- sum((vf - vy.hat)^2)
    bias.mcp[trials] <- sum((res.mcp$beta-vbeta.til)^2)
    SCORES.mcp       <- cbind(SCORES.mcp,scores.mcp)

	#############################################################################
	
	a3 <- proc.time()[3]
	res.scad <- ebic.ncvreg(vy,mX,penalty="SCAD")
	vy.hat <- mX.til%*%res.scad$beta
	vgamma.hat <- as.numeric(res.scad$beta[-1]!=0)
	scores.scad <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
	b3 <- proc.time()[3]     
    print(b3-a3)
    t.scad[trials] <- b3-a3
    MSE.scad[trials]  <- sum((vf - vy.hat)^2)
    bias.scad[trials] <- sum((res.scad$beta-vbeta.til)^2)  
    SCORES.scad <- cbind(SCORES.scad,scores.scad)

    
    ###########################################################################    
    
 

	v0=exp(seq(log(0.1),log(1),,5))
	#v0=exp(seq(log(0.1),log(1),,20))
	v1=1000
	beta_init=rep(1,p-1)
	a=b=1
	epsilon=10E-4
	sigma_init=1  
	
	if (doEMVS) { 
	
		a3 <- proc.time()[3]
	    res.emvs <- EMVS(vy.cent, mX.til[,-1], v0=v0,v1=v1,type="betabinomial",beta_init=beta_init,sigma_init=1,epsilon=epsilon,a=a,b=b)
	    res.emvs.best <- EMVSbest.my(res.emvs)
	    #EMVSplot(res.emvs,plot_type="both",omit.zeroes=TRUE)
	    vgamma.hat       <- rep(0,p-1)
		vgamma.hat[res.emvs.best$indices] <- 1
		beta.hat.origin <- transfer(mX, res.emvs.best$beta.hat)
		#vy.hat <- mX.std%*%res.emvs.best$beta.hat + mean(y)
		vy.hat <- mX.til%*%beta.hat.origin + mean(vy)  
	    scores.emvs <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
		b3 <- proc.time()[3]     
	    print(b3-a3)	
	    t.emvs[trials] <- b3-a3
	    MSE.emvs[trials]  <- sum((vf - vy.hat)^2)
	    bias.emvs[trials] <- sum((beta.hat.origin - vbeta.til)^2)  
	    SCORES.emvs <- cbind(SCORES.emvs,scores.emvs)
    }
 
 
	#############################################################################

	if (doBMS) {
	
	
		
		res.init.lasso   <- ncvreg(mX,vy,penalty="lasso",dfmax=n)
		mbeta.lasso      <- res.init.lasso$beta[-1,]
		screening.lasso <- as.vector(which(mbeta.lasso[,ncol(mbeta.lasso)]>0))
		
		res.init.scad    <- ncvreg(mX,vy,penalty="SCAD",dfmax=n)
		mbeta.scad       <- res.init.scad$beta[-1,]
		screening.scad  <- as.vector(which(mbeta.scad[,ncol(mbeta.scad)]>0))
		
		res.init.mcp    <- ncvreg(mX,vy,penalty="MCP",dfmax=n)
		mbeta.mcp       <- res.init.mcp$beta[-1,]
		screening.mcp  <- as.vector(which(mbeta.mcp[,ncol(mbeta.mcp)]>0))
		
		
		start.value <- sort(unique(c(screening.mcp,screening.lasso,screening.scad))) 
		
		
	
		a3 <- proc.time()[3]
		res.bms <- bms(cbind(vy,mX), burn = 1000, iter = 1000, nmodel = 1000, mcmc = "bd", start.value = start.value,
		    g = "hyper=3", mprior = "random", mprior.size = NA, user.int = TRUE, g.stats = TRUE, logfile = TRUE, logstep = 100, 
		    force.full.ols = FALSE, fixed.reg=numeric(0))	
		fit        <- predict(res.bms, exact=TRUE)
		coef.bms   <- coef(res.bms,order.by.pip=FALSE )
	    postmean   <- coef.bms[,2]
	    prob0      <- round(coef.bms[,1])
	    ind        <- coef.bms[,5]
		vgamma.hat      <- rep(0,p-1)
		vgamma.hat[ind] <- prob0
		scores.bms <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
	 	b3 <- proc.time()[3]     
	    print(b3-a3)
	    t.bms[trials] <- b3-a3 
	    MSE.bms[trials]  <- sum((vf - fit)^2)
	    bias.bms[trials] <- sum((vbeta - postmean)^2)  
	    
	    SCORES.bms <- cbind(SCORES.bms,scores.bms)
	}
	
	#ans <- readline()
	
	###############################################################################
	
	if (doVARBVS) {
		a3 <- proc.time()[3]
		varbvs.res  <- varbvs (X=mX.til, Z=NULL, y=vy, family = "gaussian", sigma=1, sa=1, logodds=log(0.01/0.99))
		vgamma.hat  <- round(varbvs.res$alpha[-1])
		vbeta.hat   <- varbvs.res$mu
		vy.hat      <- mX.til %*%  vbeta.hat
		scores.varbvs <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
		b3 <- proc.time()[3]     
	    print(b3-a3) 
	    t.varbvs[trials] <- b3-a3 
	    MSE.varbvs[trials]  <- sum((vf - vy.hat)^2)
	    bias.varbvs[trials] <- sum((vbeta.hat - vbeta.til)^2)  
	    SCORES.varbvs <- cbind(SCORES.varbvs,scores.varbvs)
	}
	
	if (doCVA) {
		library(correlation)
		K <- 20
	  initial_gamma <- matrix(rbinom(K * ncol(mX.til), 1, 10/(K * ncol(mX.til))), K, ncol(mX.til))
		a4 <- proc.time()[3]
		n <- nrow(mX.til)
		p <- ncol(mX.til)

		# Normalise
		y.n <- (vy - mean(vy))/sd(vy)
		mX.n <- matrix(0,n,p)
		for (j in 1:p) {
			mX.n[,j] <- (mX.til[,j] - mean(mX.til[,j]))/sd(mX.til[,j])
		} 

		cva.res  		<- cva (initial_gamma, y.n, mX.n, K)
		vR2 <- rep(0, K)
		vp_gamma <- rep(0, K)
		for (model in 1:K) {
			inds <- which(cva.res$models[model, ] == 1)
			vR2[model] <- (t(y.n) %*% mX.n[, inds] %*% solve(t(mX.n[, inds]) %*% mX.n[, inds])  %*% t(mX.n[, inds]) %*% y.n)/(t(y.n) %*% y.n)
			vp_gamma[model] <- sum(cva.res$models[model, ])
		}
		log_p <- function(n, p, vR2, vp_gamma)
		{
		  a <- 1
		  b <- p
		  return(-n / 2 * log(1 - vR2) - vp_gamma / 2 * log(n) + lbeta(a + vp_gamma, b + p - vp_gamma))
		}
		vlog_p <- log_p(n, P, vR2, vp_gamma)

		# Get gamma with maximum likelihood
		vgamma.hat  <- cva.res$models[which.max(vlog_p), ]
		vgamma.hat  <- apply(vlog_p * cva.res$models, 2, sum) / sum(vlog_p)
		vgamma.hat.inds <- which(vgamma.hat != 0)
		vbeta.hat.inds   <- solve(t(mX.til[, vgamma.hat.inds]) %*% mX.til[, vgamma.hat.inds]) %*% t(mX.til[, vgamma.hat.inds]) %*% vy
		vbeta.hat <- rep(0, ncol(mX.til))
		vbeta.hat[vgamma.hat.inds] <- vbeta.hat.inds
		vy.hat      <- mX.til %*%  vbeta.hat
		scores.cva 	<- CalcSelectionScores(c(vgamma),c(vgamma.hat[2:length(vgamma.hat)])) 
		b4 <- proc.time()[3]     
		print(scores.cva)
	    print(b4-a4) 
	    t.cva[trials] <- b4-a4 
	    MSE.cva[trials]  <- sum((vf - vy.hat)^2)
	    bias.cva[trials] <- sum((vbeta.hat - vbeta.til)^2)  
	    SCORES.cva <- cbind(SCORES.cva,scores.cva)
	}

	dat <- cbind(
				as.numeric(SCORES.lasso[10,] ),
				as.numeric(SCORES.scad[10,] ),
				as.numeric(SCORES.mcp[10,] ),
				as.numeric(SCORES.emvs[10,] ),
				as.numeric(SCORES.bms[10,] ),
				as.numeric(SCORES.varbvs[10,] ),
				as.numeric(SCORES.cva[10,] )
				)
				
	dev.off()
	
	boxplot(dat,names=c("lasso","scad","mcp","emvs","bms","varbvs","cva"),ylim=c(0,1))
 
}

