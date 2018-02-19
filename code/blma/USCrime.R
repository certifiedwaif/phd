dat <- read.table(file="uscrime.txt",header=TRUE)
dat <- na.omit(dat)

# Get y vector and X matrix
y.t <- dat$Crime
X.f <- model.matrix(Crime~.,dat)[,-1]

# Get dinemsins
n <- nrow(X.f)
p <- ncol(X.f)

mu.y      <- mean( y.t )
sigma2.y  <- (n-1)*var(y.t)/n
y.n <- (y.t - mu.y)/sqrt(sigma2.y)

# Nornalise covariates
X.n <- matrix(0,n,p)
mu.x <- c()
sigma2.x <- c()
for (j in 1:p)
{
  mu.x[j]      <- mean( X.f[,j] )
  sigma2.x[j]  <- (n-1)*var(X.f[,j])/n
  X.n[,j] <- (X.f[,j] - mu.x[j])/sqrt(sigma2.x[j])
}


library(correlation)

vy = y.n
mX = X.n

vR2 <- correlation::all_correlations_mX(vy, mX)$vR2

mGamma = graycode(p,0)

vq = mGamma%*%matrix(1,p,1)

M = nrow(mGamma)

XTy = t(mX)%*%vy
XTX = t(mX)%*%mX

# Check Graycode order
verr = c()
for (i in 1:M) {
  inds = which(mGamma[i,]==1)
  if (length(inds)>2) {
    R2 = t(XTy[inds])%*%solve(XTX[inds,inds])%*%XTy[inds]/n
    verr <- c(verr,abs(R2-vR2[i]))
    #cat(i,R2,vR2[i],abs(R2-vR2[i]),"\n")
  }
}

print("The maximum error is:")
print(max(verr))


vBIC = n*log(1 - vR2) + vq*log(n) # + 2*vq*log(vq) - 2*vq

a = -0.75
vb = 0.5*(n - vq - 5) - a
c = 0.5*(n - 1)
vd = 0.5*vq + a

log.vp.ZE = -(vb+1)*log(1 - vR2) + lbeta(vd+1,vb+1) - lbeta(a+1,vb+1)
vZE <- -2*log.vp.ZE

model.post.ZE <- exp(log.vp.ZE - max(log.vp.ZE))/sum(exp(log.vp.ZE - max(log.vp.ZE)))


ord = order(log.vp.ZE,decreasing=TRUE)

TOPMODELS = 20

best = which.min(vZE)

tab = cbind( mGamma, round(100*model.post.ZE,2), round(vR2,3), vq,  round(vZE - vZE[best],3), round(vBIC - vBIC[best],3)  )[ord[1:TOPMODELS],]


var.post.ZE <- t(mGamma)%*%model.post.ZE

print(tab)
print(var.post.ZE)

models <- apply(mGamma[ord[1:TOPMODELS],] == 1, 1, which)
stats <- cbind(round(100*model.post.ZE, 2), round(vR2, 2),   round(vZE - vZE[best], 2), round(vBIC - vBIC[best], 2)  )[ord[1:TOPMODELS],]
for (i in 1:TOPMODELS) {
  # Print model indices
  for (idx in 1:(length(models[[i]])-1)) {
    cat(models[[i]][idx], ", ")
  }
  if (length(models[[i]]) != 1) {
    cat(models[[i]][length(models[[i]])])
  }
  # Print stats
  for (j in 1:ncol(stats)) {
    cat("& ", stats[i, j])
  }
  cat("\\\\\n")
}
