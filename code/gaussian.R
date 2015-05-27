# John Ormerod's Poisson regression code, slightly modified to be used by my
# zero-inflated model code.

library(limma)
library(Matrix)
library(Rcpp)
library(RcppEigen)
library(digest)

source("CalculateB.R")
sourceCpp(file="fastdiag.cpp")

# Laplace's method of approxmation ----

f.lap <- function(vmu, vy, vr, mC, mSigma.inv, mLambda) 
{
  d <- length(vmu)
  veta <- mC %*% vmu
  mSigma <- solve(mSigma.inv)
  mDiag <- fastdiag(mC, mLambda)
  f <- sum(vy * vr * veta - vr * exp(veta + 0.5 * mDiag)) - 0.5 * t(vmu) %*% mSigma.inv %*% vmu - 
           0.5 * sum(diag(mLambda %*% mSigma))
  return(f)
}

norm <- function(v) sqrt(sum(v ^ 2))

vg.lap <- function(vmu, vy, vr, mC, mSigma.inv, mLambda) 
{       
  vg <- t(mC) %*% (vr * vy - vr * exp(mC %*% vmu)) - mSigma.inv %*% vmu
  
  return(vg)
}

mH.lap <- function(vmu, vy, vr, mC, mSigma.inv, mLambda) 
{
  vw <- exp(mC %*% vmu)
  dim(vw) <- NULL
  mH <- -t(mC * vw) %*% mC - mSigma.inv
  return(mH)
}

fit.Lap <- function(vmu, vy, vr, mC, mSigma.inv, mLambda) 
{
  MAXITER <- 100
  
  for (ITER in 1:MAXITER) {
    f  <- f.lap(vmu, vy, vr, mC, mSigma.inv, mLambda)
    vg <- vg.lap(vmu, vy, vr, mC, mSigma.inv, mLambda)
    mH <- mH.lap(vmu, vy, vr, mC, mSigma.inv, mLambda)

    # This algorithm can be numerically unstable. If the derivative becomes
    # infinite, we should return the last vmu and mLambda, rather than
    # updating them.
    if (any(is.nan(vg) || is.infinite(vg))) {
      break
    }

    old_mLambda <- mLambda
    old_vmu <- vmu
    mLambda <- solve(-mH, tol=1.0E-99)
    vmu <- vmu + mLambda %*% vg
    
    if (any(diag(mLambda < 0.0))) {
      # We've gone out of the allowable parameter space. Take the last known
      # good value.
      vmu <- old_vmu
      mLambda <- old_mLambda
      break
    }
    
    if (max(abs(vg)) < 1.0E-8) {
        break
    }
  }

  f <- f.lap(vmu, vy, vr, mC, mSigma.inv, mLambda) + 0.5*log(det(mLambda %*% mSigma.inv))
  return(list(vmu=vmu, mLambda=mLambda, f=f))
}

vtheta_enc <- function(vmu, mR)
{
  diag(mR) <- log(diag(mR))
  Rinds <- which(lower.tri(mR, diag=TRUE))
  c(vmu, mR[Rinds])
}

vtheta_dec <- function(vtheta, d)
{
  vmu <- vtheta[1:d]
  mR <- matrix(0, d, d)
  Rinds <- which(lower.tri(mR, diag=TRUE))
  mR[Rinds] <- vtheta[(d + 1):length(vtheta)]
  diag(mR) <- exp(diag(mR))
  diag(mR) <- min(c(1.0E3, diag(mR)))
  return(list(vmu=vmu, mR=mR, Rinds=Rinds))
}

# Gaussian variational approxmation, mLambda = mR mR^T ----
f.G <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, gh) 
{
  d <- length(vmu)
  
  vmu.til     <- mC %*% vmu
  vsigma2.til <- fastdiag(mC, mLambda)
  vB0 <- B0.fun("POISSON", vmu.til, vsigma2.til, gh) 
  f <- sum(vr * (vy * vmu.til - vB0)) - 0.5 * t(vmu) %*% mSigma.inv %*% vmu - 0.5 * tr(mSigma.inv %*% mLambda)
  f <- f - 0.5 * d * log(2 * pi) + 0.5 * log(det(mSigma.inv)) 
  return(f)
}

f.GVA <- function(vtheta, vy, vr, mC, mSigma.inv, gh)
{
  d <- ncol(mC)
  decode <- vtheta_dec(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR
  
  mLambda <- tcrossprod(mR)
  
  # f <- sum(log(diag(mR))) + f.G(vmu, mLambda, vy, vr, mC, mSigma.inv, gh)
  f <- -sum(log(diag(mR))) + f.G(vmu, mLambda, vy, vr, mC, mSigma.inv, gh)
  f <- f + 0.5 * d * log(2 * pi) + 0.5 * d
  
  if (!is.finite(f)) {
    f <- -1.0E16
  }
  #cat("f.GVA: f", round(f, 2), "vmu", round(vmu, 2), "diag(mR)", round(diag(mR), 2), "\n")
  return(f)
}

vg.G <- function(vmu, vy, vr, mC, mSigma.inv, vB1) 
{
  vg <- t(mC) %*% (vr * (vy - vB1)) - mSigma.inv %*% vmu     
  return(vg)
}

mH.G <- function(vmu, vy, vr, mC, mSigma.inv, vB2) 
{
  vw <-  vB2
  dim(vw) <- NULL
  mH <- -t(mC * (vr * vw)) %*% (mC) - mSigma.inv
  return(mH)    
}

vg.GVA.approx <- function(vmu, vy, vr, mC, mSigma.inv, gh)
{
  n <- length(vy)
  d <- ncol(mC)
	eps <- 1.0E-6
	f <- f.GVA(vmu, vy, vr, mC, mSigma.inv, gh)
	vg.approx <- matrix(0, dmLambda_check, 1)
	for (i in 1:d) {
		vmup <- vmu 
		vmup[i] <- vmu[i] + eps
		fp <- f.GVA(vmup, vy, vr, mC, mSigma.inv, gh)

		vmum <- vmu 
		vmum[i] <- vmu[i] - eps
		fm <- f.GVA(vmum, vy, vr, mC, mSigma.inv, gh)

		vg.approx[i] <- (fp - fm) / (2 * eps)
	}
	return(vg.approx)
}

vg.GVA.mat_approx <- function(vmu, mR, vy, vr, mC, mSigma.inv, gh)
{
  d <- ncol(mC)

  eps <- 1.0E-6
  vg.approx <- matrix(0, d, d)
  mR_copy <- mR
  for (i in 1:d) {
    for (j in 1:d) {
      mR_copy[i, j] <- mR_copy[i, j] + eps
      vtheta <- vtheta_enc(vmu, mR_copy)
      fp <- f.GVA(vtheta, vy, vr, mC, mSigma.inv, gh)

      mR_copy[i, j] <- mR_copy[i, j] - 2 * eps
      vtheta <- vtheta_enc(vmu, mR_copy)
      fm <- f.GVA(vtheta, vy, vr, mC, mSigma.inv, gh)
      mR_copy[i, j] <- mR_copy[i, j] + eps

      vg.approx[i, j] <- (fp - fm) / (2 * eps)
    }
  }
  return(vg.approx)
}

vg.GVA <- function(vtheta, vy, vr, mC, mSigma.inv, gh)
{
  d <- ncol(mC)
  decode <- vtheta_dec(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR
  Rinds <- decode$Rinds

  mLambda <- tcrossprod(mR)

  vmu.til     <- mC %*% vmu
  # Idea: Could multiply by mR and then square?
  vsigma2.til <- fastdiag(mC, mLambda)
  res.B12 <- B12.fun("POISSON", vmu.til, vsigma2.til, gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  vg <- rep(0, length(vtheta))
  vg[1:d] <- vg.G(vmu, vy, vr, mC, mSigma.inv, vB1) 

  mLambda.inv <- solve(mLambda, tol=1e-99)
  mH <- mH.G(vmu, vy, vr, mC, mSigma.inv, vB2)
  dmLambda <- (-mLambda.inv + mH) %*% mR
  # dmLambda <- (mLambda.inv + mH) %*% mR

  # Check derivative numerically
  dmLambda_check <- vg.GVA.mat_approx(vmu, mR, vy, vr, mC, mSigma.inv, gh)
  if (any(!is.finite(dmLambda) || is.nan(dmLambda))) {
   dmLambda <- matrix(0, d, d)
  }

  if (any(abs(dmLambda[Rinds] - dmLambda_check[Rinds]) > 1)) {
    cat("Analytic and numeric derivatives disagree.\n")
    # print(round(dmLambda, 2))
    # print(round(dmLambda_check, 2))
    print(round(dmLambda - dmLambda_check, 2))
    # print(round(diag(dmLambda - dmLambda_check), 2))
    # browser()
  } else {
    cat("Analytic and numeric derivatives agree.\n")
  }
  
  # # cat("vg.GVA: dmLambda", dmLambda, "\n")
  vg[(1 + d):length(vtheta)] <- dmLambda[Rinds]
  # vg[(1 + d):length(vtheta)] <- dmLambda_check[Rinds]
  #cat("vg.GVA: vg", round(vg, 2), "norm", sqrt(sum(vg^2)), "\n")
 
  return(vg)
}

fit.GVA <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, method, reltol=1.0e-12)
{
  gh2 <- NULL #list(x=gh$nodes, w=gh$weights, w.til=gh$weights*exp(gh$nodes^2))
		
  d <- length(vmu)
  mR <- t(chol(mLambda + diag(1.0E-8, d)))
	vtheta <- vtheta_enc(vmu, mR)
  lower_constraint <- rep(-Inf, length(vtheta))
  #lower_constraint[(d + 1):length(vmu)] <- -15
  
  if (method == "L-BFGS-B") {
    controls <- list(maxit=100, trace=6, fnscale=-1, REPORT=1, factr=1.0E-5, lmm=10, pgtol=reltol)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000, trace=0, fnscale=-1, REPORT=1, reltol=reltol) 
  } else {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, reltol=reltol) 
  }
  res <- optim(par=vtheta, fn=f.GVA, gr=vg.GVA, 
                method=method, lower=lower_constraint, upper=Inf, control=controls, 
                vy=vy, vr=vr, mC=mC, mSigma.inv=mSigma.inv, gh=gh2)        
      
  vtheta <- res$par 
  decode <- vtheta_dec(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR
  
  mLambda <- mR %*% t(mR)
  
  return(list(res=res, vmu=vmu, mLambda=mLambda))
}

# Gaussian variational approxmation, mLambda = (mR mR^T)^-1 ----
calc_cache <- function(hash, vtheta, mC)
{
  d <- ncol(mC)
  decode <- vtheta_dec(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR
  Rinds <- decode$Rinds
  
  mLambda.inv <- tcrossprod(mR)

  vmu.til <- mC %*% vmu
  va <- forwardsolve(mR, t(mC))
  vsigma2.til <- colSums(va^2)

  return(list(hash=hash, vmu=vmu,
              mR=mR, Rinds=Rinds, mLambda.inv=mLambda.inv,
              vmu.til=vmu.til, vsigma2.til=vsigma2.til))
}

f.G_new <- function(vmu, mR, vy, vr, mC, vmu.til, vsigma2.til, mSigma.inv, gh) 
{
  d <- length(vmu)
  vB0 <- B0.fun("POISSON", vmu.til, vsigma2.til, gh) 
  
  f <- sum(vr * (vy * vmu.til - vB0)) - 0.5 * t(vmu) %*% mSigma.inv %*% vmu
  f <- f - 0.5 * tr(forwardsolve(mR, backsolve(mR, mSigma.inv, transpose=TRUE)))
  f <- f - 0.5 * d * log(2 * pi) + 0.5 * log(det(mSigma.inv)) 
  return(f)
}

f.GVA_new <- function(vtheta, vy, vr, mC, mSigma.inv, gh, mR, Rinds, Dinds)
{
  hash <- digest(vtheta, algo="md5")
  if (is.null(cache$hash) || cache$hash != hash) {
    cache <<- calc_cache(hash, vtheta, mC)
  }
  vmu <- cache$vmu
  mR <- cache$mR
  mLambda.inv <- cache$mLambda.inv
  vmu.til <- cache$vmu.til
  vsigma2.til <- cache$vsigma2.til
  
  d <- length(vmu)
  f <- -sum(log(diag(mR))) + f.G_new(vmu, mR, vy, vr, mC, vmu.til, vsigma2.til, mSigma.inv, gh) 
  f <- f + 0.5 * d * log(2 * pi) + 0.5*d
  
  if (!is.finite(f)) {
    f <- -1.0E16
  }
  
  return(f)
}

vg.GVA.approx_new <- function(vtheta, vy, vr, mC, mSigma.inv, gh, mR, Rinds, Dinds)
{
  n <- length(vy)
  P <- length(vtheta)
  d <- ncol(mC)
  eps <- 1.0E-6
  f <- f.GVA_new(vtheta, vy, vr, mC, mSigma.inv, gh, mR, Rinds, Dinds)
  vg.approx <- matrix(0, P, 1)
  for (i in 1:P) {
    vthetap <- vtheta
    vthetap[i] <- vtheta[i] + eps
    fp <- f.GVA_new(vthetap, vy, vr, mC, mSigma.inv, gh, mR, Rinds, Dinds)
    
    vthetam <- vtheta
    vthetam[i] <- vtheta[i] - eps
    fm <- f.GVA_new(vthetam, vy, vr, mC, mSigma.inv, gh, mR, Rinds, Dinds)    
    
    vg.approx[i] <- (fp - fm) / (2 * eps)
  }
  return(vg.approx)
}

vg.GVA_new <- function(vtheta, vy, vr, mC, mSigma.inv, gh, mR, Rinds, Dinds)
{
  hash <- digest(vtheta, algo="md5")
  if (is.null(cache$hash) || cache$hash != hash) {
    cache <<- calc_cache(vtheta, mC)
  }
  vmu <- cache$vmu
  mR <- cache$mR
  Rinds <- cache$Rinds
  mLambda.inv <- cache$mLambda.inv
  vmu.til <- cache$vmu.til
  vsigma2.til <- cache$vsigma2.til
  
  res.B12 <- B12.fun("POISSON", vmu.til, vsigma2.til, gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  
  d <- length(vmu)
  vg <- 0 * vmu
  vg[1:d] <- vg.G(vmu, vy, vr, mC, mSigma.inv, vB1) 
  
  mH <- mH.G(vmu, vy, vr, mC, mSigma.inv, vB2)
  dmLambda <- (mLambda.inv + mH) %*% mR
  # dmLambda <- (0.5 * tr(mLambda.inv) + mH)
  #dmLambda_dmR <- -solve(mLambda.inv, solve(mLambda.inv, mR))
  #dmLambda_dmR <- -solve(mLambda.inv, solve(mLambda.inv, (t(mR) + mR)))
  dmLambda_dmR <- -forwardsolve(mR, backsolve(t(mR), forwardsolve(mR, backsolve(t(mR), (t(mR) + mR)))))
  #dmLambda_dmR <- -forwardsolve(mR, backsolve(t(mR), forwardsolve(mR, backsolve(t(mR), mR))))
  dmR <- dmLambda %*% dmLambda_dmR
  dmR[Dinds] <- dmR[Dinds] * mR[Dinds]

  # Check derivative numerically
  # dmR_check <- matrix(grad(func, as.vector(mR)), nrow(mR), ncol(mR))
  
  # cat("GVA2 vmu", vmu[9:10], "\ndvmu", vg[9:10], "\nmR", mR[9:10, 9:10], "\ndmR", 
		# 	dmR[9:10, 9:10]) #, "\ndmR_check", dmR_check[9:10, 9:10])
  
  vg[(1 + d):length(vtheta)] <- dmR[Rinds]    
  # vg[(1+d):length(vtheta)] <- dmR_check[Rinds]    
  
  return(vg)
}

# Swap fixed and random effects in mLambda so that inverse of mR is quick to
# calculate due to sparsity. If you do this, you have to re-order mSigma.inv, 
# vmu and mC as well.
swap_mX_mZ <- function(vmu, mLambda, mSigma.inv, mC, p, u_dim)
{
  beta_idx <- 1:p
  u_idx <- (p + 1):(p + u_dim)
  new_beta_idx <- (u_dim + 1):(u_dim + p)
  new_u_idx <- 1:u_dim
  mLambda_new <- blockDiag(mLambda[u_idx, u_idx], mLambda[beta_idx, beta_idx])
  mLambda_new[new_u_idx, new_beta_idx] <- mLambda[u_idx, beta_idx]
  mLambda_new[new_beta_idx, new_u_idx] <- t(mLambda[u_idx, beta_idx])
  mLambda <- mLambda_new
  mSigma.inv <- blockDiag(mSigma.inv[u_idx, u_idx], mSigma.inv[beta_idx, beta_idx])
  vmu <- c(vmu[u_idx], vmu[beta_idx])
  mC <- mC[, c(u_idx, beta_idx)]
  return(list(vmu=vmu, mLambda=mLambda, mSigma.inv=mSigma.inv, mC=mC))
}

# Swap everything back
swap_mX_mZ_back <- function(vmu, mLambda, mSigma.inv, mC, p, u_dim)
{
  beta_idx <- (u_dim + 1):(u_dim + p)
  u_idx <- 1:u_dim
  new_beta_idx <- 1:p
  new_u_idx <- (p + 1):(p + u_dim)
  mLambda_new <- blockDiag(mLambda[beta_idx, beta_idx], mLambda[u_idx, u_idx])
  mLambda_new[new_beta_idx, new_u_idx] <- mLambda[beta_idx, u_idx]
  mLambda_new[new_u_idx, new_beta_idx] <- t(mLambda[beta_idx, u_idx])
  mLambda <- mLambda_new
  mSigma.inv <- blockDiag(mSigma.inv[beta_idx, beta_idx], mSigma.inv[u_idx, u_idx])
  vmu <- c(vmu[beta_idx], vmu[u_idx])
  mC <- mC[, c(beta_idx, u_idx)]
  return(list(vmu=vmu, mLambda=mLambda, mSigma.inv=mSigma.inv, mC=mC))
}

fit.GVA_new <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, method, reltol=1.0e-12, p=NA, m=NA, 
                        blocksize=NA, spline_dim=NA)
{
  #gh  <- gauss.quad(N, kind="hermite")
  gh2 <- NULL 
  cache <<- list(hash=NULL)
  
  d <- length(vmu)
  Dinds <- d * ((1:d) - 1) + (1:d)
  
  u_dim <- (m - 1) * blocksize + spline_dim

  result <- swap_mX_mZ(vmu, mLambda, mSigma.inv, mC, p, u_dim)
  vmu <- result$vmu
  mLambda <- result$mLambda
  mSigma.inv <- result$mSigma.inv
  mC <- result$mC

  mR <- t(chol(solve(mLambda, tol=1.0E-99)))
  vtheta <- vtheta_enc(vmu, mR)

  if (method == "L-BFGS-B") {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, factr=1.0E-5, lmm=10)
  } else if (method == "Nelder-Mead") {
    controls <- list(maxit=100000000, trace=0, fnscale=-1, REPORT=1000, reltol=reltol) 
  } else {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, reltol=reltol) 
  }
  res <- optim(par=vtheta, fn=f.GVA_new, gr=vg.GVA_new, 
               method=method, upper=Inf, control=controls, 
               vy=vy, vr=vr, mC=mC, mSigma.inv=mSigma.inv,
               gh=gh2, mR=mR*0, Rinds=Rinds, Dinds=Dinds)        
  
  vtheta <- res$par 
  decode <- vtheta_decode(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR
  mLambda <- solve(tcrossprod(mR), tol=1.0E-99)

  result <- swap_mX_mZ_back(vmu, mLambda, mSigma.inv, mC, p, u_dim)
  vmu <- result$vmu
  mLambda <- result$mLambda
  mSigma.inv <- result$mSigma.inv
  mC <- result$mC
  
  return(list(res=res, vmu=vmu, mLambda=mLambda))
}

# Newton-Raphson Gaussian variational approximation ----
vg.G_nr <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, vB1) 
{
  vg <- t(mC) %*% (vr*(vy - vB1)) - mSigma.inv %*% vmu     
  return(vg)
}

mH.G_nr <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, vB2) 
{
  mH <- -t(mC*as.vector(vr*vB2)) %*% mC - mSigma.inv 
  return(mH)
}

fit.GVA_nr <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, method, reltol=1.0e-12, m=NA, p=NA, 
                        blocksize=NA, spline_dim=NA)
{
  MAXITER <- 1000
  TOL <- reltol
  
  for (ITER in 1:MAXITER) {
    vmu.old <- vmu
    
    # calculate B1
    # calculate B2
    vmu.til <- mC %*% vmu
    vsigma2.til <- fastdiag(mC, mLambda)
    res.B12 <- B12.fun("POISSON", vmu.til, vsigma2.til, gh)
    vB1 <- res.B12$vB1
    vB2 <- res.B12$vB2      

    vg <- vg.G_nr(vmu, mLambda, vy, vr, mC, mSigma.inv, vB1) 
    mH <- mH.G_nr(vmu, mLambda, vy, vr, mC, mSigma.inv, vB2) 
    
    # Use block inverse formula to speed computation
    # Let -mH = [A B]
    #           [B D]
    u_dim <- (m - 1) * blocksize + spline_dim
    A <- -mH[1:p, 1:p]
    B <- -mH[1:p, (p + 1):(p + u_dim)]
    D <- -mH[(p + 1):(p + u_dim), (p + 1):(p + u_dim)]
    # Then -mH^{-1} = [(A - B D^-1 B^T)^-1, -(A-B D^-1 B^T)^-1 B D^-1]
    #                 [-D^-1 B^T (A - B D^-1 B^T)^-1, D^-1 + D^-1 B^T (A - B D^-1 B^T)^-1 B D^-1]
    # D^-1 and (A - B D^-1 B^T)^-1 appear repeatedly, so we precalculate them
    D.inv <- solve(D)
    A_BDB.inv <- solve(A - B %*% D.inv %*% t(B))
    beta_idx <- 1:p
    u_idx <- (p + 1):(p + u_dim)
    mLambda[beta_idx, beta_idx] <- A_BDB.inv
    mLambda[beta_idx, u_idx] <- -A_BDB.inv %*% B %*% D.inv
    mLambda[u_idx, u_idx] <- D.inv + D.inv %*% t(B) %*% A_BDB.inv %*% B %*% D.inv
    mLambda[u_idx, beta_idx] <- t(mLambda[beta_idx, u_idx])
    
    vmu <- vmu + mLambda %*% vg
      
    err <- max(abs(vmu - vmu.old)) 
    if (err < TOL) {
      break;
    }
  } 
  
  return(list(vmu=vmu, mLambda=mLambda))
}
