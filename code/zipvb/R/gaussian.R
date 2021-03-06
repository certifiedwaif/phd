# John Ormerod's Poisson regression code, slightly modified to be used by my
# zero-inflated model code.

# Laplace's method of approxmation ----

f.lap <- function(vmu, vy, vr, mC, mSigma.inv, mLambda)
{
  d <- length(vmu)
  veta <- mC %*% vmu
  mSigma <- solve(mSigma.inv)
  mDiag <- fastdiag(mLambda, mC)
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
    # FIXME: Need to handle the case where solve can't successfully invert
    mLambda <- tryCatch(solve(-mH + diag(1e-8, length(vmu)), tol=1.0E-99), error = function(e) {NULL})
    # If we can't invert -mH, keep old mLambda and return
    if (is.null(mLambda)) {
      mLambda <- old_mLambda
      break
    }
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

safeexp <- function(x)
{
  t <- getOption("threshold")
  if (is.null(t)) {
    t <- 2
  }
  a <- exp(t) / 2
  b <- (1 - t) * exp(t)
  c <- (1 - (t^2) / 2 - (1 - t) * t) * exp(t)

  ifelse(x < t, exp(x), a * x^2 + b * x + c)
}

safeexp_diff <- function(x)
{
  t <- getOption("threshold")
  if (is.null(t)) {
    t <- 2
  }
  a <- exp(t) / 2
  b <- (1 - t) * exp(t)
  c <- (1 - (t^2) / 2 - (1 - t) * t) * exp(t)

  ifelse(x < t, exp(x), 2 * a * x + b)
}

safeexp_diff_diff <- function(x)
{
  t <- getOption("threshold")
  if (is.null(t)) {
    t <- 2
  }
  a <- exp(t) / 2
  b <- (1 - t) * exp(t)
  c <- (1 - (t^2) / 2 - (1 - t) * t) * exp(t)

  ifelse(x < t, exp(x), 2 * a)
}

safeexp_inv <- function(x)
{
  t <- getOption("threshold")
  if (is.null(t)) {
    t <- 2
  }
  a <- exp(t) / 2
  b <- (1 - t) * exp(t)
  c <- (1 - (t^2) / 2 - (1 - t) * t) * exp(t)

  # Complete the square
  h <- -b / (2 * a)
  k <- c - a * h^2
  ifelse(x < exp(t), log(x), sqrt((x - k) / a) + h)
}

vtheta_enc <- function(vmu, mR)
{
  safe_exp <- !is.null(getOption("safe_exp"))
  if (safe_exp)
    diag(mR) <- safeexp_inv(diag(mR))
  else
    diag(mR) <- log(diag(mR))
  Rinds <- which(lower.tri(mR, diag=TRUE))
  c(vmu, mR[Rinds])
}

vtheta_dec <- function(vtheta, d)
{
  safe_exp <- !is.null(getOption("safe_exp"))
  vmu <- vtheta[1:d]
  mR <- matrix(0, d, d)
  Dinds <- d*((1:d)-1)+(1:d)
  Rinds <- which(lower.tri(mR, diag=TRUE))
  mR[Rinds] <- vtheta[(d + 1):length(vtheta)]
  if (safe_exp)
    diag(mR) <- safeexp(diag(mR))
  else
    diag(mR) <- exp(diag(mR))
  # Threshold entries in the diagonal of mR at 100,000
  for (i in 1:length(Dinds)) {
      mR[Dinds[i]] <- min(c(1.0E5,mR[Dinds[i]]))
  }
  return(list(vmu=vmu, mR=mR, Dinds=Dinds, Rinds=Rinds))
}

# Gaussian variational approxmation, mLambda = mR mR^T ----
f.G <- function(vmu, mR, vy, vr, mC, mSigma.inv, gh)
{
  d <- length(vmu)

  # vmu.til     <- min(mC %*% vmu, 1000)
  vmu.til     <- mC %*% vmu
  vsigma2.til <- fastdiag2(mR, mC)
  mLambda <- tcrossprod(mR)
  vB0 <- B0.fun("POISSON", vmu.til, vsigma2.til, gh)
  f <- sum(vr * (vy * vmu.til - vB0))
  f <- f - 0.5 * t(vmu) %*% mSigma.inv %*% vmu - 0.5 * tr(mSigma.inv %*% mLambda)
  f <- f - 0.5 * d * log(2 * pi) + 0.5 * log(det(mSigma.inv))
  return(f)
}

f.GVA <- function(vtheta, vy, vr, mC, mSigma.inv, gh)
{
  d <- ncol(mC)
  decode <- vtheta_dec(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR
  Dinds <- decode$Dinds

  f <- sum(log(diag(mR))) + f.G(vmu, mR, vy, vr, mC, mSigma.inv, gh)
  f <- f + 0.5 * d * log(2 * pi) + 0.5 * d

  if (is.nan(f) || !is.finite(f)) {
    f <- -1.0E16
  }

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
  Dinds <- decode$Dinds
  Rinds <- decode$Rinds

  # mLambda <- tcrossprod(mR)

  vmu.til     <- mC %*% vmu
  # Idea: Could multiply by mR and then square?
  vsigma2.til <- fastdiag2(mR, mC)
  res.B12 <- B12.fun("POISSON", vmu.til, vsigma2.til, gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2

  vg <- rep(0, length(vtheta))
  vg[1:d] <- vg.G(vmu, vy, vr, mC, mSigma.inv, vB1)

  # mLambda.inv <- solve(mLambda + diag(1e-8, d), tol=1e-99)
  mLambda.inv <- chol2inv(t(mR))
  mH <- mH.G(vmu, vy, vr, mC, mSigma.inv, vB2)
  dmLambda <- (mLambda.inv + mH) %*% mR
  dmLambda[Dinds] <- dmLambda[Dinds] * mR[Dinds]
  # dmLambda <- (mLambda.inv + mH) %*% mR

  # Check derivative numerically
  # dmLambda_check <- vg.GVA.mat_approx(vmu, mR, vy, vr, mC, mSigma.inv, gh)
  if (any(!is.finite(dmLambda) || is.nan(dmLambda))) {
   dmLambda <- matrix(0, d, d)
  }

  # if (any(abs(dmLambda[Rinds] - dmLambda_check[Rinds]) > 1)) {
  #   cat("Analytic and numeric derivatives disagree.\n")
  #   # print(round(dmLambda, 2))
  #   # print(round(dmLambda_check, 2))
  #   print(round(dmLambda - dmLambda_check, 2))
  #   # print(round(diag(dmLambda - dmLambda_check), 2))
  #   # browser()
  # } else {
  #   cat("Analytic and numeric derivatives agree.\n")
  # }

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
  # Dinds <- d*((1:d)-1)+(1:d)
  # lower_constraint[d + Dinds] <- -8

  if (method == "L-BFGS-B") {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, factr=1.0E-5, lmm=10, pgtol=reltol)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000, trace=0, fnscale=-1, REPORT=1, reltol=reltol)
  } else {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, reltol=reltol)
  }
  res <- optim(par=vtheta, fn=f.GVA, gr=vg.GVA,
                method=method, lower=-Inf, upper=Inf, control=controls,
                vy=vy, vr=vr, mC=mC, mSigma.inv=mSigma.inv, gh=gh2)

  vtheta <- res$par
  decode <- vtheta_dec(vtheta, d)
  vmu <- decode$vmu
  mR <- decode$mR

  mLambda <- mR %*% t(mR)

  return(list(res=res, vmu=vmu, mLambda=mLambda))
}

# Gaussian variational approxmation, mLambda = (mR mR^T)^-1 ----
vtheta_enc_new <- function(vmu, mR, Rinds)
{
  # cat("vtheta_enc_new: vmu ", vmu, "\n")
  safe_exp <- !is.null(getOption("safe_exp"))
  if (safe_exp) 
    diag(mR) <- -safeexp_inv(diag(mR))
  else
    diag(mR) <- -log(diag(mR))
  return(c(vmu, mR[Rinds]))
}

vtheta_dec_new <- function(vtheta, d, Rinds)
{
  # browser()
  safe_exp <- !is.null(getOption("safe_exp"))
  vmu <- vtheta[1:d]
  mR <- matrix(0, d, d)
  Dinds <- d*((1:d)-1)+(1:d)
  mR[Rinds] <- vtheta[(d + 1):length(vtheta)]
  if (safe_exp)
    diag(mR) <- safeexp(-diag(mR))
  else
    diag(mR) <- exp(-diag(mR))
  return(list(vmu=vmu, mR=mR, Dinds=Dinds))
}

f.G_new <- function(vmu, mR, vy, vr, mC, mSigma.inv, gh)
{
  d <- length(vmu)
  # cat("vmu", vmu, "\n")
  vmu.til     <- mC %*% vmu
  # cat("vmu.til", vmu.til, "\n")
  # va <- forwardsolve(mR, t(mC))
  # vsigma2.til <- colSums(va^2)
  vsigma2.til <- fastsolve(mR, mC)
  # cat("vsigma2.til", vsigma2.til, "\n")
  # for (i in 1:length(vsigma2.til)) {
  #   if (vsigma2.til[i] > 1) {
  #     vsigma2.til[i] <- 1
  #   }
  # }
  mLambda <- chol2inv(t(mR) + diag(1e-8, d))
  # vsigma2.til <- fastdiag(mLambda, mC)
  vB0 <- B0.fun("POISSON", vmu.til, vsigma2.til, gh)
  # cat("vB0", vB0, "\n")
  # cat("sum(vy * vmu.til - vB0)", sum(vy * vmu.til - vB0), "\n")
  # cat("sum(vr * (vy * vmu.til - vB0))", sum(vr * (vy * vmu.til - vB0)), "\n")
  if (any(!is.finite(vB0))) {
    f <- -1.0E16
  } else {
    f <- sum(vr * (vy * vmu.til - vB0))
    f <- f - 0.5 * t(vmu) %*% mSigma.inv %*% vmu
    f <- f - 0.5 * tr(mSigma.inv %*% mLambda)
    f <- f - 0.5 * d * log(2 * pi) + 0.5 * log(det(mSigma.inv))
  }
  return(f)
}

f.GVA_new <- function(vtheta, vy, vr, mC, mSigma.inv, gh, Rinds)
{
  d <- ncol(mC)
  decode <- vtheta_dec_new(vtheta, d, Rinds)
  vmu <- decode$vmu
  mR <- decode$mR
  Dinds <- decode$Dinds
  for (i in 1:length(Dinds)) {
    mR[Dinds[i]] <- min(c(1.0E5, mR[Dinds[i]]))
    mR[Dinds[i]] <- max(c(1.0E-5, mR[Dinds[i]]))
  }

  # Threshold entries of mR at 1e-2
  #cat("mR[Dinds]", mR[Dinds], "\n")

  # mLambda.inv <- tcrossprod(mR)

  # cat("sum(log(diag(mR)))", sum(log(diag(mR))), "\n")
  f <- -sum(log(diag(mR))) + f.G_new(vmu, mR, vy, vr, mC, mSigma.inv, gh)
  f <- f + 0.5 * d + 0.5 * d * log(2 * pi)

  if (!is.finite(f) || abs(f) > 1e36) {
    f <- -1.0E16
  }
  # cat("f.GVA_new: f", round(f, 2), "vmu", round(vmu, 2), "diag(mR)", round(diag(mR), 2), "\n")
  return(f)
}

vg.G_new <- function(vmu, vy, vr, mC, mSigma.inv, vB1)
{
  vg <- t(mC) %*% (vr * (vy - vB1)) - mSigma.inv %*% vmu
  return(vg)
}

mH.G_new <- function(vmu, vy, vr, mC, mSigma.inv, vB2)
{
  vw <-  vB2
  dim(vw) <- NULL
  mH <- -t(mC * (vr * vw)) %*% (mC) - mSigma.inv
  return(mH)
}

vg.GVA_new.approx <- function(vmu, vy, vr, mC, mSigma.inv, gh, Rinds)
{
  n <- length(vy)
  d <- ncol(mC)
  eps <- 1.0E-6
  f <- f.GVA_new(vmu, vy, vr, mC, mSigma.inv, gh, Rinds)
  vg.approx <- matrix(0, dmLambda_check, 1)
  for (i in 1:d) {
    vmup <- vmu
    vmup[i] <- vmu[i] + eps
    fp <- f.GVA_new(vmup, vy, vr, mC, mSigma.inv, gh, Rinds)

    vmum <- vmu
    vmum[i] <- vmu[i] - eps
    fm <- f.GVA_new(vmum, vy, vr, mC, mSigma.inv, gh, Rinds)

    vg.approx[i] <- (fp - fm) / (2 * eps)
  }
  return(vg.approx)
}

vg.GVA_new.mat_approx <- function(vmu, mR, vy, vr, mC, mSigma.inv, gh, Rinds)
{
  d <- ncol(mC)

  eps <- 1.0E-6
  vg.approx <- matrix(0, d, d)
  mR_copy <- mR
  for (i in 1:d) {
    for (j in 1:d) {
      mR_copy[i, j] <- mR_copy[i, j] + eps
      vtheta <- vtheta_enc_new(vmu, mR_copy, Rinds)
      fp <- f.GVA_new(vtheta, vy, vr, mC, mSigma.inv, gh, Rinds)

      mR_copy[i, j] <- mR_copy[i, j] - 2 * eps
      vtheta <- vtheta_enc_new(vmu, mR_copy, Rinds)
      fm <- f.GVA_new(vtheta, vy, vr, mC, mSigma.inv, gh, Rinds)
      mR_copy[i, j] <- mR_copy[i, j] + eps

      vg.approx[i, j] <- (fp - fm) / (2 * eps)
    }
  }
  return(vg.approx)
}

vg.GVA_new <- function(vtheta, vy, vr, mC, mSigma.inv, gh, Rinds)
{
  d <- ncol(mC)
  decode <- vtheta_dec_new(vtheta, d, Rinds)
  vmu <- decode$vmu
  mR <- decode$mR
  # mR <- new("dtrMatrix", uplo="L", diag="N", x=as.vector(decode$mR), Dim=as.integer(c(d, d)))
  Dinds <- decode$Dinds

  for (i in 1:length(Dinds)) {
      mR[Dinds[i]] <- min(c(1.0E5, mR[Dinds[i]]))
      mR[Dinds[i]] <- max(c(1.0E-5, mR[Dinds[i]]))
  }
  # cat("mR[Dinds]", mR[Dinds], "\n")

  # mLambda.inv <- tcrossprod(mR)
  # mLambda <- solve(mLambda.inv, tol=1e-99)
  mLambda <- chol2inv(t(mR) + diag(1e-8, d))
  # cat("diag(mLambda)", diag(mLambda), "\n")

  vmu.til     <- mC %*% vmu
  # cat("vmu.til", vmu.til, "\n")
  # Idea: Could multiply by mR and then square?
  # va <- forwardsolve(mR, t(mC))
  # vsigma2.til <- colSums(va^2)
  # vsigma2.til <- fastdiag(mLambda, mC)
  vsigma2.til <- fastsolve(mR, mC)
  # cat("vsigma2.til", vsigma2.til, "\n")
  res.B12 <- B12.fun("POISSON", vmu.til, vsigma2.til, gh)
  vB1 <- res.B12$vB1
  vB2 <- res.B12$vB2
  if (any(!is.finite(vB1)) || any(!is.finite(vB2)) || any(abs(vB1) > 1e100)) {
    # cat("vB1", vB1, "\n")
    # cat("vB2", vB2, "\n")
    vB1 <- 1
    vB2 <- 1
  }

  vg <- rep(0, length(vtheta))
  vg[1:d] <- vg.G(vmu, vy, vr, mC, mSigma.inv, vB1)

  mH <- mH.G(vmu, vy, vr, mC, mSigma.inv, vB2)
  # cat("rcond(mLambda.inv)", rcond(mLambda.inv), "\n")
  # cat("rcond(mH)", rcond(mH), "\n")
  # cat("rcond(mLambda.inv + mH)", rcond(mLambda.inv + mH), "\n")
  # cat("rcond(mLambda)", rcond(mLambda), "\n")
  # cat("rcond(mR)", rcond(mR), "\n")
  # dmLambda <- (mLambda.inv + mH) %*% (-mLambda %*% mR %*% mLambda)
  mR_mLambda <- mR %*% mLambda
  # mR_mLambda <- t(backsolve(mR, forwardsolve(mR, t(mR)), transpose=TRUE))
  dmLambda <- -(mR_mLambda + mH %*% mLambda %*% mR_mLambda)
  # 97% for the fixed slope, but 89.6% for the fixed intercept
  # dmLambda <- -(mR_mLambda + mH %*% forwardsolve(tcrossprod(mR) + diag(1e-8, d), mR_mLambda))
  # 90.28% on the fixed intercept, 92.65% on the slope
  # dmLambda <- -(mR_mLambda + mH %*% solve(tcrossprod(mR) + diag(1e-8, d), mR_mLambda, tol=1e-99))
  dmLambda[Dinds] <- dmLambda[Dinds] * -mR[Dinds]
  # This is very clever, but seems to mess up the accuracy
  # dmLambda <- -(mR_mLambda + mH %*% backsolve(mR, forwardsolve(mR, mR_mLambda), transpose=TRUE))
  # Accuracy of beta_1 decreased
  # dmLambda <- -mLambda %*% (mLambda.inv + mH) %*% mR %*% mLambda
  # dmLambda <- -(diag(1, d) + mLambda %*% mH) %*% mR %*% mLambda
  # dmLambda <- forwardsolve(mR, -(mR + mH / mR[Dinds]) / mR[Dinds]^2)

  # Check derivative numerically
  # dmLambda_check <- vg.GVA_new.mat_approx(vmu, mR, vy, vr, mC, mSigma.inv, gh)
  # if (any(!is.finite(dmLambda) || is.nan(dmLambda))) {
  #  dmLambda <- matrix(0, d, d)
  # }

  # if (any(abs(dmLambda[Rinds] - dmLambda_check[Rinds]) > 1)) {
  #   cat("Analytic and numeric derivatives disagree.\n")
  #   # print(round(dmLambda, 2))
  #   # print(round(dmLambda_check, 2))
  #   print(round(dmLambda - dmLambda_check, 2))
  #   # print(round(diag(dmLambda - dmLambda_check), 2))
  #   # browser()
  # } else {
  #   cat("Analytic and numeric derivatives agree.\n")
  # }

  # # cat("vg.GVA_new: dmLambda", dmLambda, "\n")
  vg[(1 + d):length(vtheta)] <- dmLambda[Rinds]
  # vg[(1 + d):length(vtheta)] <- dmLambda_check[Rinds]
  if (any(!is.finite(vg)) || any(is.nan(vg))) {
    cat("vmu", vmu, "\n")
    cat("diag(mR)", diag(mR), "\n")
    cat("vg.GVA_new: vg", round(vg, 2), "norm", sqrt(sum(vg^2)), "\n")
  }
  return(vg)
}

# Swap fixed and random effects in mLambda so that inverse of mR is quick to
# calculate due to sparsity. If you do this, you have to re-order mSigma.inv,
# vmu and mC as well.
swap_mX_mZ <- function(vmu, mLambda, mSigma.inv, mC, p, u_dim)
{
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
  return(list(vmu=vmu, mLambda=mLambda, mSigma.inv=mSigma.inv, mC=mC))
}

# Swap everything back
swap_mX_mZ_back <- function(vmu, mLambda, mSigma.inv, mC, p, u_dim)
{
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
  return(list(vmu=vmu, mLambda=mLambda, mSigma.inv=mSigma.inv, mC=mC))
}

fit.GVA_new <- function(vmu, mLambda, vy, vr, mC, mSigma.inv, method, reltol=1.0e-12, p=NA, m=NA, blocksize=NA, spline_dim=NA)
{

  gh2 <- NULL #list(x=gh$nodes, w=gh$weights, w.til=gh$weights*exp(gh$nodes^2))

  d <- length(vmu)

  # Swap order of parameters so that random effects are before fixed effects
  u_dim = (m-1)*blocksize+spline_dim
  if (m > 0) {
    result <- swap_mX_mZ(vmu, mLambda, mSigma.inv, mC, p, u_dim)
    vmu <- result$vmu
    mLambda <- result$mLambda
    mSigma.inv <- result$mSigma.inv
    mC <- result$mC
  }

  mR <- t(chol(solve(mLambda + diag(1.0E-8, d))))
  mC2 <- mC
  mC2[mC != 0] <- 1
  ones <- crossprod(mC2)# + solve(mSigma.inv)
  ones[ones != 0] <- 1
  Rinds <- intersect(which(lower.tri(mR, diag=TRUE)), which(ones == 1))
  # Rinds <- which(lower.tri(mR, diag=TRUE))
  # print(length(Rinds))
  vtheta <- vtheta_enc_new(vmu, mR, Rinds)

  if (method == "L-BFGS-B") {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, factr=1.0E-5, lmm=10, pgtol=reltol)
  } else if (method=="Nelder-Mead") {
    controls <- list(maxit=100000000, trace=0, fnscale=-1, REPORT=1, reltol=reltol)
  } else {
    controls <- list(maxit=1000, trace=0, fnscale=-1, REPORT=1, reltol=reltol)
  }
  res <- optim(par=vtheta, fn=f.GVA_new, gr=vg.GVA_new,
                method=method, lower=-Inf, upper=Inf, control=controls,
                vy=vy, vr=vr, mC=mC, mSigma.inv=mSigma.inv, gh=gh2, Rinds=Rinds)

  vtheta <- res$par
  decode <- vtheta_dec_new(vtheta, d, Rinds)
  vmu <- decode$vmu
  mR <- decode$mR
  # browser()

  mLambda <- solve(crossprod(t(mR)) + diag(1e-8, d), tol=1e-99)

  # Swap back
  if (m > 0) {
    result <- swap_mX_mZ_back(vmu, mLambda, mSigma.inv, mC, p, u_dim)
    vmu <- result$vmu
    mLambda <- result$mLambda
    mSigma.inv <- result$mSigma.inv
    mC <- result$mC
  }

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

  d <- ncol(mC)
  for (ITER in 1:MAXITER) {
    vmu.old <- vmu
    mLambda.old <- mLambda

    # calculate B1
    # calculate B2
    vmu.til <- mC %*% vmu
    vsigma2.til <- fastdiag(mLambda, mC)
    res.B12 <- B12.fun("POISSON", vmu.til, vsigma2.til, gh)
    vB1 <- res.B12$vB1
    vB2 <- res.B12$vB2

    vg <- vg.G_nr(vmu, mLambda, vy, vr, mC, mSigma.inv, vB1)
    mH <- mH.G_nr(vmu, mLambda, vy, vr, mC, mSigma.inv, vB2)

    # Use block inverse formula to speed computation
    # Let -mH = [A B]
    #           [B D]
    if (m > 0) {
      u_dim <- (m - 1) * blocksize + spline_dim
      A <- -mH[1:p, 1:p]
      B <- -mH[1:p, (p + 1):(p + u_dim)]
      D <- -mH[(p + 1):(p + u_dim), (p + 1):(p + u_dim)]
      # Then -mH^{-1} = [(A - B D^-1 B^T)^-1, -(A-B D^-1 B^T)^-1 B D^-1]
      #                 [-D^-1 B^T (A - B D^-1 B^T)^-1, D^-1 + D^-1 B^T (A - B D^-1 B^T)^-1 B D^-1]
      # D^-1 and (A - B D^-1 B^T)^-1 appear repeatedly, so we precalculate them
      D.inv <- tryCatch(solve(D), error=function(e) NULL)
      # D.inv <- solve(D)
      # We can't invert D, so bail out with the last known good vmu and mLambda.
      # This probably means that the answer we've come up with to our optimisation
      # problem isn't very good, because we're in a part of the parameter space that
      # we shouldn't be in!
      if (is.null(D.inv)) {
        vmu <- vmu.old
        mLambda <- mLambda.old
        break
      }
      A_BDB.inv <- solve(A - B %*% D.inv %*% t(B))
      beta_idx <- 1:p
      u_idx <- (p + 1):(p + u_dim)
      mLambda[beta_idx, beta_idx] <- A_BDB.inv
      mLambda[beta_idx, u_idx] <- -A_BDB.inv %*% B %*% D.inv
      mLambda[u_idx, u_idx] <- D.inv + D.inv %*% t(B) %*% A_BDB.inv %*% B %*% D.inv
      mLambda[u_idx, beta_idx] <- t(mLambda[beta_idx, u_idx])
    } else {
      mLambda <- solve(-mH)
    }

    vmu <- vmu + mLambda %*% vg

    err <- max(abs(vmu - vmu.old))
    if (is.nan(err)) {
      # vmu and mLambda are probably full of NaNs as well, so return last good
      # vmu and mLambda
      vmu <- vmu.old
      mLambda <- mLambda.old
      break
    }
    if (err < TOL) {
      break;
    }
  }

  return(list(vmu=vmu, mLambda=mLambda))
}
