
dyn.load("Gassianintegrals.dll")

B0.fun <- function(vmu,vsigma2) {
    n <- length(vmu)
    res.B0 <- .C("R_vB0",mu=vmu,sigma=sqrt(vsigma2),N=as.integer(50),val=rep(0,n),n=as.integer(n))
    return(res.B0$val)
}

B12.fun <- function(vmu,vsigma2) {
    n <- length(vmu)
    res.B12 <- .C("R_vB12",mu=vmu,sigma=sqrt(vsigma2),N=as.integer(50),val=rep(0,n),n=as.integer(n))
    return(list(vB1=res.B12$val1,vB2=res.B12$val2))
}

trapint <- function(xgrid, fgrid) 
{
    ng <- length(xgrid)
    xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
    fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
    integ <- sum(xvec * fvec)/2
    return(integ)
}

mu = 0
sigma = 1
xgrid = seq(from = -3.9527, to = 5.1807, length.out = 50)
fgrid = log(1+exp(xgrid)) * dnorm(xgrid, mu, sigma)
trapint(xgrid, fgrid)