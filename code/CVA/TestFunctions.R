####################################################################################################

set.seed(1)

source("functions.R")


n = 30

ln = vector("list", 6)
ln[[1]] = n
ln[[2]] = n
ln[[3]] = n
ln[[4]] = n
ln[[5]] = n
ln[[6]] = n

lp = vector("list", 6)
lp[[1]] = 3 
lp[[2]] = 10
lp[[3]] = 25
lp[[4]] = 3
lp[[5]] = 10
lp[[6]] = 25

lR2 = vector("list", 6)
lR2[[1]] = 0.1
lR2[[2]] = 0.1
lR2[[3]] = 0.1
lR2[[4]] = 0.95
lR2[[5]] = 0.95
lR2[[6]] = 0.95

PLOT.G = TRUE
PLOT.U = TRUE
PLOT.SIGMA = TRUE
PLOT.ALPHA = TRUE

PDF = TRUE
CHECK = FALSE 

####################################################################################################

if (PLOT.G) 
{
	if (PDF) {
		pdf("gGivenY.pdf",width=12,height=8)
	}
	
	par(mfrow=c(2,3))
	
	for (i in 1:6) {
	
		n = ln[[i]]
		p = lp[[i]]
		R2 = lR2[[i]]
		
		plot.gGivenY(g=NULL,n,p,R2,N=1000,LEGEND=(i==1),CHECK=TRUE)
	}
	
	if (PDF) {
		dev.off()
	}
}
	
####################################################################################################

if (PLOT.U) 
{
	if (PDF) {
		pdf("uGivenY.pdf",width=12,height=8)
	}
	
	par(mfrow=c(2,3))
	
	for (i in 1:6) {
	
		n = ln[[i]]
		p = lp[[i]]
		R2 = lR2[[i]]
		
		plot.uGivenY(u=NULL,n,p,R2,N=1000,LEGEND=(i==1),CHECK=TRUE)
	}
	
	if (PDF) {
		dev.off()
	}
}

####################################################################################################

if (PLOT.SIGMA) {
	if (PDF) {
		pdf("sigma2GivenY.pdf",width=12,height=8)
	}
	
	par(mfrow=c(2,3))
	
	for (i in 1:6) {
	
		n = ln[[i]]
		p = lp[[i]]
		R2 = lR2[[i]]

		plot.sigma2GivenY(sigma2=NULL,n,p,R2,N=1000,LEGEND=(i==1),CHECK=TRUE)
	}
	
	if (PDF) {
		dev.off()
	}
}

####################################################################################################

if (PLOT.ALPHA) {

	if (PDF) {
		pdf("alphaGivenY.pdf",width=12,height=8)
	}
	
	par(mfrow=c(2,3))
	
	for (i in 1:6) {
	
		n = ln[[i]]
		p = lp[[i]]
		R2 = lR2[[i]]
		
		plot.alphaGivenY(alpha=NULL,n,p,R2,N=1000,LEGEND=(i==1),CHECK=TRUE)
	}
	
	if (PDF) {
		dev.off()
	}
}
	
	
	 
