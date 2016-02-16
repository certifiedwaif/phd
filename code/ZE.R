library(Rcpp)
library(RcppArmadillo)

sourceCpp(file="ZE.cpp")

################################################################################

greycode <- function(p, standard_order=TRUE) 
{
	A <- matrix(c(0,1),2,1)
	if (p!=1) {
		for (i in 2:p) {
			P <- nrow(A)
			inds <- rev(1:P)
			A1 <- matrix(c(rep(0,P),rep(1,P)),2*P,1)
			A2 <- rbind(A,matrix(A[inds,],nrow(A),ncol(A)))
			A <- cbind(A1,A2)
		}
	}
	if (standard_order) {
		A_standard_order <- matrix(0, nrow(A), ncol(A))
		for (i in 1:nrow(A)) {
			for (j in 1:ncol(A)) {
				A_standard_order[i, (ncol(A) + 1) - j] <- A[i, j]
			}
		}
		return(A_standard_order)
	} else
	{
		return(A)
	}
}

################################################################################

ZE.constants <- function(n,pmax,LARGEP) 
{
	vcon <- c()
	vpen <- c()
	for (q in 0:pmax) 
	{
		con <- 0.5*(n-q) - 0.75
		dof <- q
		pen <- lbeta(0.5*q + 0.25,con) - lbeta(0.25,con)  
	if (LARGEP) {
		pen <- pen + lbeta(1+q,1+p-q)
	}
		vcon[q+1] <- con
		vpen[q+1] <- pen
	}
	return(list(vcon=vcon,vpen=vpen))
}

################################################################################

ZE.exact.slow <- function(vy,mX,LARGEP) 
{
	n <- length(vy)
	p <- ncol(mX)
	A <- greycode(p) 
	res.con <- ZE.constants(n,p,LARGEP) 
 
	vlog.ZE <- rep(0,nrow(A))
	vlog.ZE[1] <- res.con$vpen[1]
	vR2 <- rep(0, nrow(A))
	for (j in 2:nrow(A))
	{
		inds <- which(A[j,]==1)
		q <- length(inds)
		R2 <- summary(lm(vy~-1+mX[,inds]))$r.squared
		vR2[j] <- R2
		vlog.ZE[j] <-  -res.con$vcon[q+1]*log(1 - R2) + res.con$vpen[q+1]
	}
	return(list(A=A,vlog.ZE=vlog.ZE, vR2=vR2))
}

################################################################################

ZE.exact <- function(vy,mX,LARGEP,STORE.XTX=TRUE) 
{
	n <- length(vy)
	p <- ncol(mX)
	A <- greycode(p)
	vq <- A%*%matrix(1,p,1)
	linds <- apply(A==1, 1, which) 
	res.con <- ZE.constants(n,p,LARGEP) 
	if (STORE.XTX) {
		XTX <- t(mX)%*%mX
	}
	XTy <- t(mX)%*%vy
	yTy <- t(vy)%*%vy
	vR2 <- rep(0,nrow(A))
	for (j in 2:nrow(A))
	{
		inds <- linds[[j]]
		b <- XTy[inds]
		bbT <- b%*%t(b)
		if (STORE.XTX) {
			Q <- XTX[inds,inds]
		} else {
			X <- mX[,inds]
			Q <- t(X)%*%X
		}
		if (vq[j]==1) {
			Z <- 1/Q
		} else {
			Z <- solve(Q)
		}
		vR2[j] <- sum(bbT*Z) 
	}
	vR2 <- vR2/yTy
	vlog.ZE  <-  -res.con$vcon[vq+1]*log(1 - vR2) + res.con$vpen[vq+1]
	return(list(A=A,vlog.ZE=vlog.ZE, vR2=vR2))
}

################################################################################
 
ZE.exact.fast <- function(vy,mX,LARGEP) 
{
	n <- length(vy)
	p <- ncol(mX)
	q <- 0
	mA <- greycode(p)
	mD <- diff(mA)
	vonep <- matrix(1,p,1)
	vs <- (mD%*%vonep)==1
	vq <- mA%*%vonep
	vw <- abs(mD%*%matrix(1:p,p,1))
	res.con <- ZE.constants(n,p,LARGEP) 
	
	lmZ <- list()
	linds11 <- list()
	linds12 <- list()
	linds21 <- list()
	linds22 <- list()
	for (i in 1:p) {
		lmZ[[i]] <- matrix(0,i,i)
		inds <- 1:(i-1)
		mZ11 <- matrix(FALSE,i,i); mZ11[inds,inds] <- TRUE
		mZ12 <- matrix(FALSE,i,i); mZ12[inds,i]    <- TRUE
		mZ21 <- matrix(FALSE,i,i); mZ21[i,inds]    <- TRUE
		mZ22 <- matrix(FALSE,i,i); mZ22[i,i]       <- TRUE
		linds11[[i]] <- which(mZ11)
		linds12[[i]] <- which(mZ12)
		linds21[[i]] <- which(mZ21)
		linds22[[i]] <- which(mZ22)
	}
	
	XTX <- t(mX)%*%mX
	XTy <- t(mX)%*%vy
	yTy <- t(vy)%*%vy
	vR2  <- rep(0,nrow(mA))
	inds <- c()
	for (j in 2:nrow(mA)) 
	{
		cat("Iteration", j, "\n")
		if (vs[j-1]) {
			# Adding variable
			cat("Adding", vw[j-1], "\n")
			indsNew <- c(inds,vw[j-1])
			cat("indsNew", indsNew, "\n")
		} else {
			# Removing varable
			cat("Removing", vw[j-1], "\n")
			k <- which(inds==vw[j-1])
			indsNew <- inds[-k]
			cat("indsNew", indsNew, "\n")
		}
		b <- XTy[indsNew]
		q <- vq[j]
		if (q==1) {
			lmZ[[1]] <- 1/XTX[indsNew,indsNew] 
			Zb <- lmZ[[1]]*b
		} else {
			if (vs[j-1]) {
				v <- XTX[inds,vw[j-1]]
				Zv <- lmZ[[q-1]]%*%v
				d <- 1/(XTX[vw[j-1],vw[j-1]] - sum(v*Zv))
				cat("v", v, "\n")
				cat("Zv", Zv, "\n")
				cat("d", d, "\n")
				lmZ[[q]][linds11[[q]]] <- lmZ[[q-1]] + d*Zv%*%t(Zv)
				lmZ[[q]][linds12[[q]]] <- -d*Zv
				lmZ[[q]][linds21[[q]]] <- -d*Zv
				lmZ[[q]][linds22[[q]]] <- d
			} else {
				Z12 <- lmZ[[q+1]][-k,k]
				ZO <- Z12%*%t(Z12)
				lmZ[[q]] <- lmZ[[q+1]][-k,-k] -  ZO/lmZ[[q+1]][k,k]
			}
			Zb <- lmZ[[q]]%*%b
		} 
		vR2[j]  <- sum(b*Zb) 
		inds <- indsNew
	}
	vR2 <- vR2/yTy 
	vlog.ZE  <-  -res.con$vcon[vq+1]*log(1 - vR2) + res.con$vpen[vq+1]
	return(list(mA=mA,vlog.ZE=vlog.ZE, vR2=vR2))
}

################################################################################