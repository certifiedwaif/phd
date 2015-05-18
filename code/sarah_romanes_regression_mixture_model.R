library(rje)

#--------------------------------------------------------------------------------------------------------------------#
  
#Setting up the distribution and priors 

n=500
xs=rnorm(n)
X=cbind(1,xs)
X=matrix(X,n,2)

sigma1s<- 1
alpha <- c(8,4)

sigma0s <- 2
beta <-c(1,-9)

rho=0.3
z=rbinom(n,1,rho)

y=matrix(NA,n,1)

for (i in 1:n) {
  
  y[i]=(X[i,]%*%alpha + rnorm(1,sd=sqrt(sigma1s)))*z[i] + (X[i,]%*%beta + rnorm(1,sd=sqrt(sigma0s)))*(1-z[i])
}

plot(xs, y)

I=diag(2)

#Beta (need DIFFUSE priors). Note Beta ~ N(0, SIGMA2_B)
Sigma2_B = 1.0E8*I
Sigma2_B.inv <- 1.0E-8*I

#Alpha (need DIFFUSE priors). Note Alpha ~ N(0, SIGMA2_A)
Sigma2_A = 1.0E8*I
Sigma2_A.inv <- 1.0E-8*I


#sigma                       Note sigma ~ IG(A,B)
A1=0.01
B1=0.01

A0=0.01
B0=0.01

#--------------------------------------------------------------------------------------------------------------------#

#necessary functions

tr <- function(A) {
  return(sum(diag(A)))
  }
#--------------------------------------------------------------------------------------------------------------------#

#initiating variables 

sa=matrix(1,2,2)
ma=matrix(1,2,1)

sb=matrix(2,2,2)
mb=matrix(2,2,1)

eta=matrix(1,n,1)
p=rep(0.5,n)
P=diag(p)

vec1=rep(1,n)

r1=1
r0=1
ph1=1
ph0=1


#--------------------------------------------------------------------------------------------------------------------#

#implementing the algorithm

TOL <- 1.0E-5
MAXITER <- 100

theta.old <- c(sa,ma,sb,mb,r1,ph1,r0,ph0) 
for(i in 1:MAXITER)

  {
  
  XTPX <- t(X)%*%P%*%X
  XTPy <- t(X)%*%P%*%y
  XTIPX <- t(X)%*%(diag(n)-P)%*%X
  XTIPY <- t(X)%*%(diag(n)-P)%*%y
  
 
  sa <- solve((r1/ph1)*XTPX + Sigma2_A.inv)
  ma <- sa%*%XTPy*(r1/ph1)
  
  sb <- solve((r0/ph0)* XTIPX + Sigma2_B.inv)
  mb <- sb%*%XTIPY *(r0/ph0)
  
  
  r1= as.numeric(A1+0.5*(t(vec1)%*%p))
  
  ph1=as.numeric(B1 +0.5*(t(p)%*%((y-X%*%ma)^2) + tr(sa%*%XTPX)))
  
  r0=as.numeric(A0 + 0.5*t(vec1)%*%(vec1-p))
  
  ph0=as.numeric(B0 +0.5*(t(vec1-p)%*%((y-X%*%mb)^2) + tr(sb%*%XTPX)))
  
  for(j in 1:500){
  
     eta[j]=logit(rho)+ (-0.5*log(ph1)+0.5*digamma(r1)-0.5*(r1/ph1)*(((y[j])-t(X[j,])%*%ma)^2 + t(X[j,])%*%sa%*%t(t(X[j, ])))) - (-0.5*log(ph0)+0.5*digamma(r0)-0.5*(r0/ph0)*(((y[j])-t(X[j,])%*%mb)^2 + t(X[j,])%*%sb%*%t(t(X[j, ]))))
  
  p[j]=expit(eta[j])
  P=diag(p)
  }

  
  theta <- c(sa,ma,sb,mb,r1,ph1,r0,ph0)
  err <- max(abs(theta - theta.old))
  theta.old <- theta
  if (err < TOL) {  
      break;
  }
}


param=c(sa,ma,sb,mb,r1,ph1,r0,ph0)
param
