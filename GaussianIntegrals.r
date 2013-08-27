
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

