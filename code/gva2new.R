f.G_new2 <- function(vmu,mR,vy,vr,mC,mSigma.inv,gh,p, m, blocksize, spline_dim, mC_sp)
{
  d <- length(vmu)
  
  vmu.til     <- mC%*%vmu
  #vsigma2.til <- diag(mC%*%mLambda%*%t(mC))
  #vsigma2.til = rep(0, nrow(mC))
  #for (i in 1:length(nrow(mC))) {
  #  vsigma2.til[i] = as.double(mC[i,]%*%mLambda%*%mC[i,])
  #}
  #vsigma2.til <- sapply(1:nrow(mC), function(i) {mC[i,]%*%mLambda%*%mC[i,]})
  #vsigma2.til <- apply(mC, 1, function(row) {row%*%mLambda%*%row})
  #vsigma2.til <- apply(mC, 1, function(row) {crossprod(crossprod(mLambda, row), row)})
  #vsigma2.til <- diag(mC %*% (mLambda %*%t(mC)))
  #vsigma2.til <- apply(mC, 1, function(row) {
  #  a <- forwardsolve(mR, row)
  #  sum(a^2)
  #})
  #mR <- Matrix(mR, sparse=TRUE)
  #mR_sp = sparse_R(mR, p, m, blocksize, spline_dim)
  #mR.inv = fastinv2(mR, p=p, m=m, blocksize=blocksize, spline_dim=spline_dim)
  #mR.inv2 = as.matrix(fastinv(mR, p, m, blocksize, spline_dim))
  #a <- solve(mR, t(mC))
  a <- forwardsolve(mR, t(mC))
  #a <- fastsolve(mR_sp, mC)
  vsigma2.til <- crossprod(a^2, rep(1, ncol(mC)))
  #vsigma2.til <- crossprod(a^2, rep(1, ncol(mC)))                           
  #mR.inv = fastinv(mR_sp)  
  #mLambda <- crossprod(mR.inv)
  mLambda <- chol2inv(mR)
  
  #vsigma2.til <- apply(mC, 1, function(row) {
  #  x <- row%*%t(mR.inv)
  #  return(t(x)%*%x)
  #)
  vB0 <- B0.fun("POISSON",vmu.til,vsigma2.til,gh) 
  f <- sum(vr*(vy*vmu.til - vB0)) - 0.5*t(vmu)%*%mSigma.inv%*%vmu - 0.5*tr(mSigma.inv%*%mLambda)
  f <- f - 0.5*d*log(2*pi) + 0.5*log(det(mSigma.inv)) 
  return(f)
}

###############################################################################

f.GVA_new2 <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds, p, m, blocksize, spline_dim, mC_sp)
{
  d <- ncol(mC)  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]
  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
  }   
  #mLambda.inv = mR%*%t(mR)
  #mLambda <- solve(mLambda.inv, tol=1.0E-99)
  #mR.inv = fastinv2(mR, p=p, m=m, blocksize=blocksize, spline_dim=spline_dim)
  #mLambda <- t(mR.inv) %*% mR.inv
  # TODO: You only need the _trace_ of this matrix?
  #mLambda = chol2inv(t(mR))
  
  f <- -sum(log(diag(mR))) + f.G_new2(vmu,mR, vy,vr,mC,mSigma.inv,gh, p, m, blocksize, spline_dim, mC_sp) 
  f <- f + 0.5*d*log(2*pi) + 0.5*d
  
  if (!is.finite(f)) {
    f <- -1.0E16
  }
  #cat("f.GVA_new2: vmu ", vmu, " f ", f, "\n")
  return(f)
}

###############################################################################

vg.G_new2 <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1) 
{
  vg <- t(mC)%*%(vr*(vy - vB1)) - mSigma.inv%*%vmu     
  #cat("vg.G", vg, "\n")
  return(vg)
}

###############################################################################

mH.G_new2 <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2) 
{
  vw <-  vB2; dim(vw) <- NULL
  #mH <- -t(mC*vw)%*%mC - mSigma.inv
  #mH <- -t(mC*(vr*vw))%*%(mC*vr) - mSigma.inv
  mH <- -t(mC*(vr*vw))%*%(mC) - mSigma.inv
  return(mH)    
}

###############################################################################

vg.GVA.approx_new2 <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
{
  n <- length(vy); P <- length(vtheta); d <- ncol(mC)
  eps <- 1.0E-6
  f <- f.GVA_new(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
  vg.approx <- matrix(0,P,1)
  for (i in 1:P) {
    vthetap <- vtheta
    vthetap[i] <- vtheta[i] + eps
    fp <- f.GVA_new(vthetap,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
    
    vthetam <- vtheta
    vthetam[i] <- vtheta[i] - eps
    fm <- f.GVA_new(vthetam,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)

    
    vg.approx[i] <- (fp - fm)/(2*eps)
  }
  return(vg.approx)
}

###############################################################################
fastinv2 = function(mR, p=NA, m=NA, blocksize=1, spline_dim=0)
{
  # Idea: This could be made faster by replacing with an equivalent
  # fastsolve function.
  
  u_dim = (m-1)*blocksize+spline_dim
  mR.inv.mZ = matrix(0, nrow=u_dim, ncol=u_dim)
  if (blocksize == 1 && m > 0) {
    # If the blocksize is 1, we can simply take the reciprocal of the
    # diagonal.
    diag(mR.inv.mZ) = 1/diag(mR[1:u_dim, 1:u_dim])
  } else {
    # Iterate through the blocks, taking the inverse of each
    if (m > 0)
      # FIXME: This section of code is _slow_.
      for (i in 1:m) {
        idx = ((i-1)*blocksize+1):(i*blocksize)
        mR.inv.mZ[idx, idx] = solve(mR[idx, idx])
        #mR.inv.mZ[idx, idx] = chol2inv(t(mR[idx, idx]))
      }
    
    if (spline_dim > 0) {
      spline_idx = (m*blocksize+1):(m*blocksize+spline_dim)
      spline_mR = mR[spline_idx, spline_idx]
      #spline_mR = tril(spline_mR, -3)
      spline_mR = band(spline_mR, -5, 0)
      mR.inv.mZ[spline_idx, spline_idx] = as.matrix(solve(spline_mR))
    }
  }
  
  # Last p lines can be solved as normal
  mR.inv.mX = solve(mR[(u_dim+1):(u_dim+p), (u_dim+1):(u_dim+p)], tol=1e-99)
  # Construct inverse of mR.
  mR.inv = blockDiag(mR.inv.mZ, mR.inv.mX)
  # Because of diagonal form of mR.inv.mZ, this could probably be optimised
  # further.
  mR.inv[(u_dim+1):(u_dim+p), 1:(u_dim)] = -mR.inv.mX %*% mR[(u_dim+1):(u_dim+p), 1:(u_dim)] %*% mR.inv.mZ
  return(mR.inv)
  #return(solve(mR, tol=1e-99))
  #return(Matrix(data=mR.inv, sparse=TRUE))
}
#require(compiler)
#fastinv = cmpfun(fastinv2)

###############################################################################
vg.GVA_new2 <- function(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds, p, m, blocksize, spline_dim, mC_sp)
{
  d <- ncol(mC)
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):length(vtheta)]

  mR[Dinds] <- exp(mR[Dinds]) 
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E3,mR[Dinds[i]]))
  }    
  
  # mR is lower triangular. Can you rewrite this using forward solves and
  # backsolves?
  # New
  #mR_sp = sparse_R(mR, p, m, blocksize, spline_dim)
  #mR.inv = fastinv2(mR, p=p, m=m, blocksize=blocksize, spline_dim=spline_dim)
  #mR.inv = fastinv(mR_sp)  
  #mLambda <- crossprod(mR.inv)
  mR.inv = solve(mR, tol=1.0E-99)
  #mLambda <- chol2inv(t(mR)) # Should this be chol2inv(t(mR))?
  mLambda.inv <- mR%*%t(mR)
  mLambda <- solve(mLambda.inv, tol=1.0E-99)
  
  #mR.inv = as.matrix(fastinv(mR, p, m, blocksize, spline_dim))
  #mR.inv = fastinv(mR_sp)
  #print(mR.inv - mR.inv2)
  # Old
  #mR.inv = solve(mR, tol=1.0E-99)
  
  # Old
  #mLambda.inv <- mR%*%t(mR)
  #mLambda <- solve(mLambda.inv, tol=1.0E-99)
  # New
  #mLambda <- t(mR.inv) %*% mR.inv
  #mLambda <- crossprod(mR.inv)
  #mLambda <- chol2inv(mR_sp)
  vmu.til     <- mC%*%vmu
  # TODO: Fast solve
  a <- forwardsolve(mR, t(mC))
  vsigma2.til <- crossprod(a^2, rep(1, ncol(mC)))
  #vsigma2.til <- fastdiag(mC, mLambda)

  res.B12 <- B12.fun("POISSON",vmu.til,vsigma2.til,gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  vg <- 0*vmu
  vg[1:d] <- vg.G_new2(vmu,mLambda,vy,vr,mC,mSigma.inv,vB1)
  mH <- mH.G_new2(vmu,mLambda,vy,vr,mC,mSigma.inv,vB2)
  dmLambda <- -solve(mR, tol=1.0E-99)%*%(mLambda.inv + mH)%*%mLambda  
  #dmLambda <- -solve(mR, tol=1.0E-99)%*%(diag(rep(1, ncol(mC))) + mH%*%mLambda)
  dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  #dmLambda <- -t(mR.inv%*%(mLambda.inv +  mH)%*%mLambda)
  #dmLambda <- -mR.inv%*%(diag(1, ncol(mC)) + mH %*% mLambda)
  #dmLambda <- -solve(mR, diag(1, ncol(mC)) + mH %*% mLambda)
  #dmLambda <- -solve(mR, (mLambda.inv + mH)%*%mLambda)
  
  #dmLambda[Rinds] <- dmLambda[Rinds]*mR[Rinds]
  #dmLambda[Rinds] <- t(dmLambda)[Rinds]
  #dmLambda[Dinds] <- dmLambda[Dinds]*mR[Dinds]
  #cat("vg.GVA_new2: vg ", vg, " dmLambda ", as.matrix(dmLambda), "\n")
  
  #require(gridExtra)
  #require(lattice)
  #plot1 <- image(Matrix(mR))
  #plot2 <- image(Matrix(mLambda))
  #plot3 <- image(Matrix(dmLambda))
  #plot4 <- image(Matrix(mH))
  #grob <- grid.arrange(plot1, plot2, plot3, plot4, ncol=2, nrow=2)
  #ans <- readline()
  #if (ans == "Q") stop("Time to go!")
  #print(sum(eigen(dmLambda)$values^2))
  #print("\n")
  #print(det(dmLambda))
  
  #res <- vg.GVA.approx_new(vtheta,vy,vr,mC,mSigma.inv,gh,mR,Rinds,Dinds)
 
  #print("john test")

  #print("mLambda.inv")
  #print(mLambda.inv)
  
  #print("mH")
  #print(mH)
  
  #print("mLambda")
  #print(mLambda)

  #print("mSigma.inv")
  #print(mSigma.inv)

  #print("diff")
  #print(mLambda - solve(mLambda.inv,tol=1.0E-99))

  #print("mR")
  #print(mR)

  #print(res)
  #bit <- res[(1+d):length(vtheta)]
  
  #print(cbind(bit,dmLambda[Rinds],bit-dmLambda[Rinds]))    

  #ans <- readline()

  # Check with numeric derivative
  #h = 1e-8
  #dmLambda2 = matrix(0, nrow(mLambda), ncol(mLambda))
  #for (i in 1:length(Rinds)) {
  #  mR_high = mR
  #  mR_high[Rinds[i]] = mR_high[Rinds[i]] + h
  #  mR_low = mR
  #  mR_low[Rinds[i]] = mR_low[Rinds[i]] - h
  #mLambda_high = solve(mR_high%*%t(mR_high), tol=1e-99)
  #mLambda_low = solve(mR_low%*%t(mR_low), tol=1e-99)
  #dmLambda2[Rinds[i]] = (-0.5*log(det(mLambda_high)) + f.G_new(vmu,mLambda_high,vy,vr,mC,mSigma.inv,gh)  + 0.5*log(det(mLambda_low)) - f.G_new(vmu,mLambda_low,vy,vr,mC,mSigma.inv,gh))/(2*h)
  #vtheta_high = c(vmu, mR_high[Rinds])
  #vtheta_low = c(vmu, mR_low[Rinds])
  #dmLambda2[Rinds[i]] = (f.GVA_new(vtheta_high,vy,vr,mC,mSigma.inv,gh,mR_high,Rinds,Dinds) - f.GVA_new(vtheta_low,vy,vr,mC,mSigma.inv,gh,mR_low,Rinds,Dinds))/(2*h)
  #  mLambda_high = solve(mR_high%*%t(mR_high), tol=1e-99)
  #  mLambda_low = solve(mR_low%*%t(mR_low), tol=1e-99)
  #  dmLambda2[Rinds[i]] = (f.G_new(vmu,mLambda_high,vy,vr,mC,mSigma.inv,gh) - f.G_new(vmu,mLambda_low,vy,vr,mC,mSigma.inv,gh))/(2*h)
  #}  
  
  #cat("dmLambda", dmLambda, "\n")
  #cat("dmLambda2", dmLambda2, "\n")
  #cat("dmLambda diff", dmLambda - dmLambda2, "\n")
  #for (i in 1:ncol(dmLambda))
  #  for (j in 1:nrow(dmLambda))
  #    cat(j, i, dmLambda[i, j], dmLambda2[i, j], dmLambda[i, j] - dmLambda2[i, j], "\n")
  #cat("\n")
  
  #printMatrix(dmLambda)
  #printMatrix(dmLambda2)
  #printMatrix(dmLambda-dmLambda2)
  
  vg[(1+d):length(vtheta)] <- dmLambda[Rinds]    
  #vg[(1+d):length(vtheta)] <- t(chol(dmLambda))[Rinds]    
  #vg[(1+d):length(vtheta)] <- dmLambda2[Rinds]    
  
  return(vg)
}

###############################################################################

fit.GVA_new2 <- function(vmu,mLambda,vy,vr,mC,mSigma.inv,method,reltol=1.0e-12, p=NA, m=NA, blocksize=NA, spline_dim=NA)
{
  #library(statmod)
  #N <- 15
  #gh  <- gauss.quad(N,kind="hermite")
  gh2 <- NULL #list(x=gh$nodes,w=gh$weights,w.til=gh$weights*exp(gh$nodes^2))
  
  d <- length(vmu)
  Dinds <- d*((1:d)-1)+(1:d)
  
  u_dim = (m-1)*blocksize+spline_dim
  # Swap fixed and random effects in mLambda so that inverse of mR is quick to
  # calculate due to sparsity. If you do this, you have to re-order mSigma.inv,
  # vmu and mC as well.
  beta_idx = 1:p
  u_idx = (p+1):(p+u_dim)
  new_beta_idx = (u_dim+1):(u_dim+p)
  new_u_idx = 1:u_dim
  mLambda_new = blockDiag(mLambda[u_idx, u_idx], mLambda[beta_idx, beta_idx])
  mLambda_new[new_u_idx, new_beta_idx] = mLambda[u_idx, beta_idx]
  mLambda_new[new_beta_idx, new_u_idx] = t(mLambda[u_idx, beta_idx])
  mLambda = mLambda_new
  mSigma.inv = blockDiag(mSigma.inv[u_idx, u_idx], mSigma.inv[beta_idx, beta_idx])
  vmu = c(vmu[u_idx], vmu[beta_idx])
  mC = mC[,c(u_idx, beta_idx)]
  #mC = Matrix(mC, sparse=TRUE)
  
  mR <- t(chol(solve(mLambda, tol=1.0E-99) + diag(1.0E-8,d)))

  #cat("mR", round(mR, 3), "\n")
  #printMatrix(round(mR, 3))
  #cat("mSigma.inv", mSigma.inv, "\n")
  mR[Dinds] <- log(mR[Dinds])
  # Idea: If Rinds included only the diagonal for the mZ part, and the
  # lower triangular mX part, we could optimise over a space of lower dimension.
  low.tri <- lower.tri(mR,diag=TRUE)
  low.tri[1:m,] <- FALSE
  diag(low.tri) <- TRUE
  Rinds <- which(low.tri)
  #Rinds <- lower.tri(mR,diag=TRUE)
  vmu <- c(vmu,mR[Rinds])
  P <- length(vmu)
  lower_constraint <- rep(-Inf, length(vmu))
  
  if (method=="L-BFGS-B") {
    #controls <- list(maxit=1000,trace=1,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
    controls <- list(maxit=10,trace=1,fnscale=-1,REPORT=1,factr=1.0E-5,lmm=10)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000,trace=0,fnscale=-1,REPORT=1000,reltol=reltol) 
  } else {
    controls <- list(maxit=1000,trace=0,fnscale=-1,REPORT=1,reltol=reltol) 
  }
  res <- optim(par=vmu, fn=f.GVA_new2, gr=vg.GVA_new2,
               method=method,lower=lower_constraint, upper=Inf, control=controls,
               vy=vy,vr=vr,mC=mC,mSigma.inv=mSigma.inv,gh=gh2,mR=mR*0,Rinds=Rinds,Dinds=Dinds, p=p, m=m, blocksize=blocksize, spline_dim=spline_dim, mC_sp=mC_sp)        
  
  vtheta <- res$par 
  
  vmu <- vtheta[1:d]
  mR[Rinds] <- vtheta[(1+d):P]
  mR[Dinds] <- exp(mR[Dinds])  
  mLambda <- solve(mR%*%t(mR), tol=1.0E-99)
  
  browser()

  # Swap everything back
  beta_idx = (u_dim+1):(u_dim+p)
  u_idx = 1:u_dim
  new_beta_idx = 1:p
  new_u_idx = (p+1):(p+u_dim)
  mLambda_new = blockDiag(mLambda[beta_idx, beta_idx], mLambda[u_idx, u_idx])
  mLambda_new[new_beta_idx, new_u_idx] = mLambda[beta_idx, u_idx]
  mLambda_new[new_u_idx, new_beta_idx] = t(mLambda[beta_idx, u_idx])
  mLambda = mLambda_new
  mSigma.inv = blockDiag(mSigma.inv[beta_idx, beta_idx], mSigma.inv[u_idx, u_idx])
  vmu = c(vmu[beta_idx], vmu[u_idx])
  mC = mC[,c(beta_idx, u_idx)]
  
  return(list(res=res,vmu=vmu,mLambda=mLambda))
}


