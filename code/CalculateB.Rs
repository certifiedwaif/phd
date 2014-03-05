###############################################################################

b0 <- function(x) { log(1+exp(x)); }
b1 <- function(x) { 1/(1+exp(-x)); }
b2 <- function(x) { mu <- b1(x); mu*(1-mu); }
b3 <- function(x) { mu <- b1(x); (mu-2*(mu^2))*(1 - mu); }
b4 <- function(x) { mu <- b1(x); (mu - 6*(mu^2)  + 6*(mu^3))*(1 - mu); }
b5 <- function(x) { mu <- b1(x); (mu - 14*(mu^2) + 36*(mu^3)  - 24*(mu^4))*(1 - mu); }
b6 <- function(x) { mu <- b1(x); (mu - 30*(mu^2) + 150*(mu^3) - 240*(mu^4) + 120*(mu^5))*(1 - vmu); }

##############################################################################

B0.fun <- function(family,mu,sigma2,gh) 
{ 
    if (family=="POISSON") {
        vB0 <- exp(mu+0.5*sigma2)
    } else {
		    vB0 <- B0.aghq(mu,sigma2,gh)$vB0
		}
    return(vB0)
}

##############################################################################

B0.sn.fun <- function(family,mu,sigma2,delta,gh) 
{ 
    if (family=="POISSON") {
        vB0 <- 2*exp(mu+0.5*sigma2)*pnorm(delta)
    } else {
		    vB0 <- B0.sn.aghq(mu,sigma2,delta,gh)$vB0
		}
    return(vB0)
}


###############################################################################

B12.fun <- function(family,mu,sigma2,gh) 
{
    if (family=="POISSON") {
        vB1 <- exp(matrix(mu)+0.5*matrix(sigma2))
        vB2 <- vB1
    } else {
				res <- B12.aghq(mu,sigma2,gh)
				vB1 <- res$vB1
				vB2 <- res$vB2	
		}		
		return(list(vB1=vB1,vB2=vB2))
}


###############################################################################

B12.sn.fun <- function(family,mu,sigma2,delta,gh) 
{
    if (family=="POISSON") {
        vB1 <- 2*exp(mu+0.5*sigma2)*pnorm(delta)
        vB2 <- vB1
    } else {
				res <- B12.aghq(mu,sigma2,gh)
				vB1 <- res$vB1
				vB2 <- res$vB2	
		}		
		return(list(vB1=vB1,vB2=vB2))
}

###############################################################################

B1234.fun <- function(family,mu,sigma2,gh) 
{
    if (family=="POISSON") {
        vB1 <- exp(mu+0.5*sigma2)
        vB2 <- vB1
        vB3 <- vB1
        vB4 <- vB1                
    } else {
				res <- B1234.aghq(mu,sigma2,gh)
				vB1 <- res$vB1
				vB2 <- res$vB2	
				vB3 <- res$vB3
				vB4 <- res$vB4					
		}		
		return(list(vB1=vB1,vB2=vB2,vB3=vB3,vB4=vB4))
}

###############################################################################

solve.cubic <- function(a,b,c) {

    Q <- (a^2 - 3*b)/9
    R <- (2*(a^3) - 9*a*b + 27*c)/54
 
    if ((R^2)<=(Q^3)) {
        theta <- acos(R/sqrt(Q^3))
        x1 <- -2*sqrt(Q)*cos(theta/3) - a/3
        x2 <- -2*sqrt(Q)*cos((theta+2*pi)/3) - a/3
        x3 <- -2*sqrt(Q)*cos((theta-2*pi)/3) - a/3
        return(c(x1,x2,x3))
    } else {
        A <- - sign(R)*(abs(R) + sqrt(R^2 - Q^3))^(1/3)
        if (A!=0) {
            B <- Q/A
        } else {
            B <- 0
        }
        return( A + B - a/3 )
    }        
}

###############################################################################

logistic.normal.approx.mode <- function(vmu,vsigma2)
{    
    a1 <- 0.0128663
    p0 <- 0.353417
    p1 <- 0.918588
    p2 <- 0.0179221
    
    vc2 <-  a1*vsigma2 + p1/p2 - vmu
    vc1 <-  a1*p1*vsigma2/p2 - 2*vsigma2 - vmu*p1/p2 + p0/p2
    vc0 <- -p1*vsigma2/p2 + a1*p0*vsigma2/p2 - vmu*p0/p2
        
    vx <- solve.cubic(vc2,vc1,vc0)
    
    return(vx)
}

###############################################################################

faster.logistic.normal.approx.mode <- function(vmu,vsigma2)
{    
    a1 <- 0.0128663
    p0 <- 0.353417
    p1 <- 0.918588
    p2 <- 0.0179221
    
    vc2 <-  a1*vsigma2 + p1/p2 - vmu
    vc1 <-  a1*p1*vsigma2/p2 - 2*vsigma2 - vmu*p1/p2 + p0/p2
    vc0 <- -p1*vsigma2/p2 + a1*p0*vsigma2/p2 - vmu*p0/p2
    
    # Almost always the middle root of the cubic solver.    
    vQ <- (vc2^2 - 3*vc1)/9
    vR <- (2*(vc2^3) - 9*vc2*vc1 + 27*vc0)/54    
    vtheta <- acos(vR/sqrt(vQ^3))
    vx <- -2*sqrt(vQ)*cos((vtheta+2*pi)/3) - vc2/3
    
    return(vx)
}

###############################################################################

find.modes <- function(vmu,vsigma2,vx.hat) 
{
    MAXITER  <- 100
    EPS.TERM <- 1.0e-5
    vsigma <- sqrt(vsigma2)
		for (ITER in 1:MAXITER) {
				vx.til <- vmu + vsigma*vx.hat
				vb0 <- b0(vx.til)
				vb1 <- b1(vx.til)
				vb2 <- b2(vx.til)		   
			  vg  <- vsigma*(vb1/vb0)  - vx.hat
				mH  <- vsigma2*(vb2/vb0 - (vb1/vb0)^2) - 1
				dx  <- vg/mH				   
				vx.hat <- vx.hat - dx
				if (all(abs(vg)<EPS.TERM)) { break; }
		}		
		return(list(vx.hat=vx.hat,mH=mH))
}

###############################################################################

aghq.grid <- function(vmu,vsigma2,gh) 
{
    n <- length(vmu)
    N <- length(gh$x)
    
    vsigma <- sqrt(vsigma2)
    
    # Calculate Initial Guess  
    #vx.0 <- matrix(0,n,1)        
    #for (i in 1:n) {    
    #   # Calculate intitial guesses for modes
    #   vx.hats <- logistic.normal.approx.mode(vmu[i],vsigma2[i])  
    #   vf.hats <- b0(vx.hats)*dnorm(vx.hats,vmu[i],vsigma2[i])  
    #   # Select the intitial guess with the highest value
    #   print(c(vmu[i],vsigma2[i],vf.hats))
    #   vx.0[i] <- vx.hats[which.max(vf.hats)]
    #}
    
    
    vx.0 <- faster.logistic.normal.approx.mode(vmu,vsigma2)      
    vx.0 <- (vx.0 - vmu)/vsigma
       
		# Polish the Initial Guess using Newton's method
		res <- find.modes(vmu,vsigma2,vx.0) 
		vmu.star     <- as.vector(res$vx.hat)
		vsigma2.star <- as.vector(-1/res$mH)
    
		sqrt.2vsigma2.star <- as.vector(sqrt(2*vsigma2.star))
		vmu.til <- vmu + vsigma*vmu.star
		
		# Calculate the weights and abscissae
		mX     <- matrix(vmu.star,n,N) + sqrt.2vsigma2.star%o%gh$x
		mX.til <- matrix(vmu.til,n,N)  + (sqrt.2vsigma2.star*vsigma)%o%gh$x
		mW.til <- (sqrt.2vsigma2.star%o%(gh$w.til/sqrt(2*pi)))*exp(-0.5*mX^2)		
		
		return(list(mX=mX,mX.til=mX.til,mW=mW.til))
}

###############################################################################

find.modes.sn <- function(vmu,vsigma2,vd,vx.hat) 
{
    MAXITER  <- 100
    EPS.TERM <- 1.0e-5
    vsigma <- sqrt(vsigma2)
    vd.til <- vd*vsigma
		for (ITER in 1:MAXITER) {
				vx.til <- vmu + vsigma*vx.hat
				vb0 <- b0(vx.til)
				vb1 <- b1(vx.til)
				vb2 <- b2(vx.til)
				
				vx.bar <- vd.til*vx.hat
				
        phi <- dnorm(vx.bar)
        Phi <- pnorm(vx.bar)
        R0 <- zeta(1,vx.bar)
        #R0 <- phi/Phi
        R1 <- -R0*(vx.bar + R0)				
				
			  vg  <- vsigma*(vb1/vb0)  - vx.hat  +  vd.til*R0
				mH  <- vsigma2*(vb2/vb0 - (vb1/vb0)^2) - 1    +  (vd.til^2)*R1
				dx  <- vg/mH				   
				vx.hat <- vx.hat - dx
				if (all(abs(vg)<EPS.TERM)) { break; }
		}		
		return(list(vx.hat=vx.hat,mH=mH))
}


###############################################################################

aghq.sn.grid <- function(mu,sigma,d,b.fun) 
{
		LG <-  -200
		RG <-  200
		NG <-  500
				
	  x.g <- seq(LG,RG,,NG) 
		f.g <- b.fun(x.g)*dnorm(x.g,mu,sigma)*2*pnorm(d*(x.g-mu))
		ind <- ((f.g/max(f.g)) > 1.0E-12)		           
				
		LG <-  min(x.g[ind])  
		RG <-  max(x.g[ind]) 	    
		NG <- 500 
		x.g <- seq(LG,RG,,NG) 
    
    return(x.g)
}

###############################################################################

B0.aghq <- function(vmu,vsigma2,gh) 
{
    negind <- (vmu<0)
		vmu    <- abs(vmu)
		
		# Calculate the weights and abscissae
		res <- aghq.grid(as.vector(vmu),as.vector(vsigma2),gh) 
		
		# Evalute the integrals
		vB0 <- (res$mW*b0(res$mX.til))%*%matrix(1,length(gh$x),1)
		
		# Fix up values corresponding to originally negative vmu
		vB0[negind] <- vB0[negind] - vmu[negind]
	
		return(list(vB0=vB0))    
}

###############################################################################


B0.sn.aghq <- function(vmu,vsigma2,vdelta,gh)
{
    n <- length(vmu)
		vsigma <- sqrt(vsigma2)
		vd <- (vdelta/vsigma2)/sqrt(1 - vdelta^2/vsigma2)
 	
		vB0 <- c()
		for (i in 1:n) {
				x.g <- aghq.sn.grid(vmu[i],vsigma[i],vd[i],b0) 
				f.g <- b0(x.g)*dnorm(x.g,vmu[i],vsigma[i])*2*pnorm(vd[i]*(x.g-vmu[i]))
				#print("B0")
				#print(i)
				#plot(x.g,f.g,type="l")
				#ans <- readline()
				vB0[i] <- trapint(x.g,f.g)		
		} 	
 	
		return(list(vB0=vB0))    
}


###############################################################################

B12.aghq <- function(vmu,vsigma2,gh) 
{
    negind <- (vmu<0)
    vmu    <- abs(vmu)
			
		# Calculate the weights and abscissae
		res <- aghq.grid(as.vector(vmu),as.vector(vsigma2),gh) 	
		
		# Evalute the integrals
		vB1 <- (res$mW*b1(res$mX.til))%*%matrix(1,length(gh$x),1)
		vB2 <- (res$mW*b2(res$mX.til))%*%matrix(1,length(gh$x),1)
		
		# Fix up values corresponding to originally negative vmu
		vB1[negind] <- 1 - vB1[negind]	
			
		return(list(vB1=vB1,vB2=vB2))    
}

###############################################################################

B12.sn.aghq <- function(vmu,vsigma2,vdelta,gh) 
{
    n <- length(vmu)
		vsigma <- sqrt(vsigma2)
		vd <- (vdelta/vsigma2)/sqrt(1 - (vdelta^2)/vsigma2)
 	
		vB1 <- c()
		vB2 <- c()
		for (i in 1:n) {
				x.g <- aghq.sn.grid(vmu[i],vsigma[i],vd[i],b1) 
				f.g <- b1(x.g)*dnorm(x.g,vmu[i],vsigma[i])*2*pnorm(vd[i]*(x.g-vmu[i]))
				vB1[i] <- trapint(x.g,f.g)		
				#print("B1")
				#print(i)
				#plot(x.g,f.g,type="l")
				#ans <- readline()
				
				x.g <- aghq.sn.grid(vmu[i],vsigma[i],vd[i],b2) 
				f.g <- b2(x.g)*dnorm(x.g,vmu[i],vsigma[i])*2*pnorm(vd[i]*(x.g-vmu[i]))
				vB2[i] <- trapint(x.g,f.g)					
				#print("B2")
				#print(i)
				#plot(x.g,f.g,type="l")
				#ans <- readline()				
		} 		
			
		return(list(vB1=vB1,vB2=vB2))    
}


###############################################################################

B1234.aghq <- function(vmu,vsigma2,gh) 
{
    negind <- (vmu<0)
    vmu    <- abs(vmu)
		
		# Calculate the weights and abscissae
		res <- aghq.grid(as.vector(vmu),as.vector(vsigma2),gh) 	
							  						
		# Evalute the integrals
		vB1 <- (res$mW*b1(res$mX.til))%*%matrix(1,length(gh$x),1)
		vB2 <- (res$mW*b2(res$mX.til))%*%matrix(1,length(gh$x),1)	
		vB3 <- (res$mW*b3(res$mX.til))%*%matrix(1,length(gh$x),1)	
		vB4 <- (res$mW*b4(res$mX.til))%*%matrix(1,length(gh$x),1)	
		
		# Fix up values corresponding to originally negative vmu
		vB1[negind] <- 1 - vB1[negind]	
		vB3[negind] <-   - vB3[negind]	
			
		return(list(vB1=vB1,vB2=vB2,vB3=vB3,vB4=vB4))    
}

###############################################################################
