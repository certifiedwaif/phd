doBAS <- TRUE
doBMS <- FALSE
doEMVS <- FALSE
doVARBVS <- TRUE

library(ncvreg)
if (doBAS) {
	library(BAS)
}
if (doBMS) {	
	library(BMS)
}
library(Matrix)
#library(EMVS)
library(varbvs)
library(ncvreg)
 

################################################################################


generate_data_QTL <- function()
{
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
	vf <- mX%*%matrix(vbeta)
	sigma2.true <- 20
	vy <- vf + rnorm(nrow(mX),0,sqrt(sigma2.true)) ##how about intercept??

	return(list(mX=mX, vbeta=vbeta, n=n, p=P, vf=vf, vy=vy))
}


generate_data_high_dimensional <- function(n=150, K=50)
{
  p <- 200
  sigma2.true <- 1.
  mX <- matrix(0, n, p)
  mSigma_block <- matrix(0.999, 10, 10)
  diag(mSigma_block) <- 1

  mSigma <- as.matrix(do.call(bdiag, rep(list(mSigma_block), 20)))
  chol_mSigma <- chol(mSigma)
  for (i in 1:n) {
    mX[i, 1:200] <- t(chol_mSigma) %*% rnorm(p)
  }
  mX <- scale(mX)
  nonzero_locations <- c(1, 11, 21, 31)
  nonzero_effects <- c(1.5, 2, 2.5, 3)
  vbeta <- rep(0, p)
  vbeta[nonzero_locations] <- nonzero_effects
  vf <- vector("double", n)
  for (i in 1:n) {
    vf[i] <- sum(nonzero_effects * mX[i, nonzero_locations])
  }
  vy <- vf + rnorm(nrow(mX),0,sqrt(sigma2.true))
  vy <- (vy - mean(vy))
  vy <- sqrt(n) * vy / sqrt(sum(vy ^ 2))
  initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
  return(list(vbeta=vbeta, n=n, p=p, vf=vf, vy=vy, mX=mX, K=K, initial_gamma=initial_gamma))
}


generate_data_Hitters <- function()
{
	library(ISLR)
	Hitters <- na.omit(Hitters)
	
	# Get y vector and X matrix
	y.t <- Hitters$Salary
	X.f <- model.matrix(Salary~.,Hitters)[,-1] 
	varnames <- colnames(X.f)

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
	colnames(mX) <- varnames
	
	return(list(vbeta=NULL, vy=vy, mX=mX, n=n, p=p))
} 


generate_data_Bodyfat <- function()
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
	colnames(mX) <- varnames
	
	return(list(vy=vy, mX=mX, n=n, p=p))
}


generate_data_Wage <- function()
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
	colnames(mX) <- varnames
	
	return(list(vy=vy, mX=mX, n=n, p=p))
}


generate_data_College <- function()
{
	library(ISLR)
	College <- na.omit(College)
	
	# Get y vector and X matrix
	y.t <- College$Grad.Rate
	X.f <- model.matrix(Grad.Rate~.,College)[,-1]
 
	varnames <- colnames(X.f)

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
	colnames(mX) <- varnames
	
	return(list(vy=vy, mX=mX, n=n, p=p))
}


generate_data_USCrime <- function()
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
	colnames(mX) <- varnames
	
	return(list(vy=vy, mX=mX, n=n, p=p))
}


generate_data_comData <- function()
{
	cat("Start of generate_data_comData\n")
	load("~/Dropbox/phd/code/comData.Rdata")
	sum.na <- function(x)
	{
	  sum(is.na(x))
	}
	inds = which(apply(X,2,sum.na)==0)
	mX = X[,inds]
	i <- 1
	vy = Y[,i]

	full_fit <- lm(vy~mX)
	summ <- summary(full_fit)

	vbeta <- coef(full_fit)
	vbeta[summ$coefficients[,4] >= .05] <- 0.
	vbeta <- vbeta[2:102]
	n <- nrow(mX)
	p <- ncol(mX)
	cat("vbeta ", length(vbeta), " p ", p, "\n")
	# initial_gamma <- matrix(rbinom(K * p, 1, .5), K, p)
	# cva_result <- cva(initial_gamma, vy, mX, K)
	vf <- mX%*%matrix(vbeta)

	return(list(vbeta=vbeta, vy=vy, mX=mX, vf=vf, n=n, p=p))
}

################################################################################

ebic.ncvreg <- function(vy,mX,penalty)
{
	n <- length(vy)
	p <- ncol(mX)
	res <- ncvreg(mX,vy,penalty=penalty)
	ebic <- c()
	for (i in 1:ncol(res$beta)) {
		dof  <- sum(res$beta[,i]!=0)
		rss  <- res$loss[i]
		ebic[i] <- log(rss/n) + (dof/n)*(log(n) + 2*log(p))
	}
	best <- which.min(ebic)
	return(list(beta=res$beta[,best],rss=res$loss[best],res=res))
}

################################################################################

EMVSbest.my <- function(result) 
{
    posts = result$log_g_function
    posts[is.finite(posts) == F] = NaN
    which <- which.max(posts)
    logpost <- posts[which]
    print("Best Model Found")
    gamma.hat <- result$prob_inclusion[which, ]
    return(list(log_g_function=logpost, 
        	indices = (1:dim(result$betas)[2])[result$prob_inclusion[which, ] > 0.5], 
        	beta.hat = result$betas[which, ], 
        	gamma.hat=gamma.hat))
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


log_p <- function(n, p, vR2, vp_gamma)
{
  a <- 1
  b <- p
  return(-n / 2 * log(1 - vR2) - vp_gamma / 2 * log(n) + lbeta(a + vp_gamma, b + p - vp_gamma))
}

model_likelihood <- function(models, y.n, mX.n)
{
	nmodels <- nrow(models)
	n <- nrow(mX.n)
	p <- ncol(mX.n)
	vR2 <- rep(0, nmodels)
	vp_gamma <- rep(0, nmodels)
	for (model in 1:nmodels) {
		inds <- which(models[model, ] == 1)
		vR2[model] <- (t(y.n) %*% mX.n[, inds] %*% solve(t(mX.n[, inds]) %*% mX.n[, inds])  %*% t(mX.n[, inds]) %*% y.n)/(t(y.n) %*% y.n)
		vp_gamma[model] <- sum(models[model, ])
	}
	vlog_p <- log_p(n, p, vR2, vp_gamma)

	return(vlog_p)
}
################################################################################


fill_in <- function(initial_gamma, lower, upper, proportion = 0.05)
{
	for (k in lower:upper) {
		valid <- FALSE
		while (!valid) {
			initial_gamma[k, ] <- rbinom(ncol(initial_gamma), 1, proportion)
			p_gamma <- sum(initial_gamma[k, ])
			if (p_gamma > 0 && p_gamma < ncol(initial_gamma)) {
				valid <- TRUE
			}
		}
	}
	return(initial_gamma)
}


initialise_gamma <- function(start, K, p, y.n, mX.n, models=NULL)
{
  if (start == "cold_start") {
    initial_gamma <- matrix(0, K, p)
    initial_gamma <- fill_in(initial_gamma, 1, K, 10/ncol(mX.n))
  } else {
  	initial_gamma <- matrix(0, K, ncol(mX.n))
    models_dedup <- models[!duplicated(models), ]
    # Remove null model
    models_dedup <- models_dedup[2:nrow(models_dedup), 2:ncol(models_dedup)]
    if (start == "warm_start_covariates") {
      cov_count <- apply(models_dedup, 1, sum)
      cov_order <- order(cov_count, decreasing = TRUE, na.last = NA)
      initial_gamma[1:min(K, length(cov_order)), ] <- models_dedup[1:min(K, length(cov_order)), ]
      if (length(cov_order) < K) {
      	initial_gamma <- fill_in(initial_gamma, length(cov_order) + 1, K)
    	}
    } else if (start == "warm_start_likelihood") {
      models_vlogp <- model_likelihood(models_dedup, y.n, mX.n)
      vlogp_order <- order(models_vlogp, decreasing = TRUE, na.last = NA)
      initial_gamma[1:min(K, length(vlogp_order)), ] <- models_dedup[1:min(K, length(vlogp_order)), ]
      if (length(vlogp_order) < K) {
      	initial_gamma <- fill_in(initial_gamma, length(vlogp_order) + 1, K)
    	}
    }
    if (K == 1) {
      initial_gamma <- matrix(initial_gamma, 1, p)
    }
  }
  return(initial_gamma)
}


QTL <- function(K, data_fn, start, prior, bUnique=TRUE, seed=1)
{
	cat(K, start, prior, "\n")
	# Check parameters
	if (!(is.numeric(K) && K > 0)) {
		stop("K must be a positive number")
	}

	cat("start", start, "\n")
	if (!(start %in% c("cold_start", "warm_start_covariates", "warm_start_likelihood"))) {
		stop("start must be one of cold_start, warm_start_covariates or warm_start_likelihood")
	}

	if (!(prior %in% c("maruyama", "BIC", "ZE", "liang_g1", "liang_g2", "liang_g3", "robust_bayarri1", "robust_bayarri2"))) {
		stop("prior must be one of BIC, ZE, liang_g1, liang_g2, liang_g3, robust_bayarri1 and robust_bayarri2")
	}

	TRIALS <- 20

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
	t.bas <- rep(0,TRIALS)

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

	doBAS <- TRUE
	doBMS <- FALSE
	doEMVS <- FALSE
	doVARBVS <- FALSE
	doVB <- TRUE
	doVBscreen <- TRUE
	doCVA <- TRUE
	
	vtimes = c()

	start_iter <- 1
	for (trials in start_iter:TRIALS) 
	{
		#############################################################################
	  set.seed(seed + trials)
		
		res.gen <- data_fn()
		mX <- res.gen$mX
		vbeta <- res.gen$vbeta
		n <- res.gen$n
		P <- res.gen$p
		vf <- res.gen$vf
		vy <- res.gen$vy

		mX.til <- cbind(1,mX)
		mX.intercept <- cbind(1,mX)
		p <- ncol(mX.intercept)

	  true.inds <- which(vbeta!=0)
	  cat("True predictors are:\n")
	  cat(which(vbeta!=0),"\n")
		
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
    
    if (doCVA) {
      cat("Start of CVA section\n")
      library(blma)
      
      # Normalise
      n <- nrow(mX)
      p <- ncol(mX)
      y.n <- (vy - mean(vy))/sd(vy)
      mX.n <- matrix(0,n,p)
      for (j in 1:p) {
        mX.n[,j] <- (mX[,j] - mean(mX[,j]))/sd(mX[,j])
      } 
      
      scad_models <- ifelse(t(res.scad$res$beta) != 0, 1, 0)
      initial_gamma <- initialise_gamma(start, K, p, y.n, mX.n, models=scad_models)

      # If K > nrow(warm_start), generate the rest of the rows by randomly selecting covariates from
      # t(res.scad$res$beta) != 0. Must ensure uniqueness.
      # FIXME: Choose model with the highest marginal likelihood
      # warm_start <- as.matrix(warm_start)
      # warm_start <- t(warm_start)
      
      # if (ncol(warm_start) == 1) {
      # 	warm_start <- t(warm_start)
      # }
      # if (ncol(warm_start) != ncol(mX)) {
      # 	warm_start <- matrix(warm_start[, 2:ncol(warm_start)], K, ncol(mX))
      # }
      
      # K <- 100
      # cat("Generating data ...")
      # initial_gamma <- matrix(0, K, ncol(mX.til))
      a4 <- proc.time()[3]
      cat(c(K, 1, prior, "\n"))
      modelprior <- "beta-binomial"
      modelprior_vec <- c(1, p)
      cva.res <- cva(y.n, mX.n, initial_gamma, prior, modelprior, modelprior_vec, cores=48)
      cat("covariates in models ", apply(cva.res$mGamma, 1, sum), "\n")
      
      #vlog_p <- model_likelihood(cva.res$mGamma, y.n, mX.n)
      vlog_p <- cva.res$posterior_model_probabilities
      cat("vlog_p ", vlog_p, "\n")
      
      # Get gamma with maximum likelihood
      vgamma.hat  <- cva.res$mGamma[which.max(vlog_p), ]
      vgamma.hat  <- apply(vlog_p * cva.res$mGamma, 2, sum) / sum(vlog_p)
      vgamma.hat.inds <- which(vgamma.hat > 0.5)
      cat("vgamma.hat.inds ", vgamma.hat.inds, "\n")
      if (length(vgamma.hat.inds) == 0) {
        vbeta.hat <- rep(0, ncol(mX.n))
        vy.hat    <- rep(mean(y.n), n)
      } else {
        vbeta.hat.inds   <- solve(t(mX.n[, vgamma.hat.inds]) %*% mX.n[, vgamma.hat.inds]) %*% t(mX.n[, vgamma.hat.inds]) %*% vy
        vbeta.hat <- rep(0, ncol(mX.n))
        vbeta.hat[vgamma.hat.inds] <- vbeta.hat.inds
        vy.hat      <- mX.n %*%  vbeta.hat
      }
      #cat("vgamma ", vgamma, " vgamma.hat ", vgamma.hat, "\n")
      scores.cva 	<- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
      b4 <- proc.time()[3]     
      cat("F1", scores.cva$F1, "\n")
      print(b4-a4) 
      t.cva[trials] <- b4-a4 
      MSE.cva[trials]  <- sum((vf - vy.hat)^2)
      bias.cva[trials] <- sum((vbeta.hat - vbeta.til)^2)  
      SCORES.cva <- cbind(SCORES.cva,scores.cva)
    }
    

		v0=exp(seq(log(0.1),log(1),,5))
		#v0=exp(seq(log(0.1),log(1),,20))
		v1=1000
		beta_init=rep(1,p-1)
		a=b=1
		epsilon=10E-4
		sigma_init=1  
		
		if (doEMVS) { 
			cat("Start of EMVS section\n")
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
			cat("Start of BMS section\n")
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
			res.bms <- bms(cbind(vy,mX), burn = 1000, iter = 1000000, nmodel = 1000, mcmc = "bd", start.value = start.value,
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
			cat("Start of VARBVS section\n")
			a3 <- proc.time()[3]
			varbvs.res  <- varbvs (X=mX.til, Z=NULL, y=vy, family = "gaussian", sigma=10.0, sa=10.0, logodds=log(0.001/0.999))
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

		if (doBAS) {
			cat("Start of BAS section\n")
		  a4 <- proc.time()[3]
			bas.res <- bas.lm(vy ~ mX, prior="g-prior", modelprior=beta.binomial(1, p), initprobs="uniform", method="MCMC", bestmodel=NULL, MCMC.iterations=1e5)
			vgamma.hat <- rep(0, length(vgamma))
			vgamma.hat[ bas.res$which[[which.max(bas.res$logmarg)]] ] <- 1
			vbeta.hat <- bas.res$mle[[which.max(bas.res$logmarg)]]
			vy.hat <- mX.til[, bas.res$which[[which.max(bas.res$logmarg)]] ] %*%  vbeta.hat[2:length(vbeta.hat)] + vbeta.hat[1]
			scores.bas <- CalcSelectionScores(c(vgamma),c(vgamma.hat)) 
			cat("F1", scores.bas$F1, "\n")
			b4 <- proc.time()[3]     
	    print(b4-a4) 
	    t.bas[trials] <- b4-a4 
	    MSE.bas[trials]  <- sum((vf - vy.hat)^2)
	    #bias.bas[trials] <- sum((vbeta.hat - vbeta.til)^2)  
	    SCORES.bas <- cbind(SCORES.bas,scores.bas)
		}

		dat <- cbind(
					as.numeric(SCORES.lasso[10,] ),
					as.numeric(SCORES.scad[10,] ),
					as.numeric(SCORES.mcp[10,] ),
					as.numeric(SCORES.emvs[10,] ),
					as.numeric(SCORES.bms[10,] ),
					as.numeric(SCORES.varbvs[10,] ),
					as.numeric(SCORES.bas[10,] ),
					as.numeric(SCORES.cva[10,] )
					)
		colnames(dat) <- c("lasso", "scad", "mcp", "bas", "cva")		
		#dev.off()
		
	 
	}
	
	mTimes = cbind(t.lasso,t.mcp,t.scad,t.emvs,t.varbvs,t.bas,t.cva)

	return(list(dat=dat,bas.res=bas.res,mTimes=mTimes,vgamma.hat=vgamma.hat,vgamma=vgamma))
}



 
