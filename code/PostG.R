library(gsl)

R2 = 0.5
n = 100
p = 10

g = seq(0,150,,10000)

a = -0.75
b = (n - p - 5)/2 - a
c = (n - 1)/2
d = p/2 + a

log.f = (b + 1)*log(1 - R2) + b*log(g) - (d + b + 2)*log(1 + g*(1 - R2)) - lbeta(d+1,b+1)

plot(g,exp(log.f),type="l")


trapint <- function(xgrid, fgrid) 
{
	ng <- length(xgrid)
	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)]
	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng]
	integ <- sum(xvec * fvec)/2
	return(integ)
}

print(trapint(g,exp(log.f)))

N = 100
vB = rbeta(N,b+1,d+1)
vg.tilde = vB/(1 - vB)
vg = vg.tilde/(1 - R2)

lines(density(vg),col="red")






library(gsl)


vsigma2 = seq(0.1,1.5,,1000)

log.f = c*log(n/2)
log.f = log.f + (b + 1)*log(1 - R2)
log.f = log.f - lgamma(c)
log.f = log.f - (c + 1)*log(vsigma2)
log.f = log.f - 0.5*n/vsigma2
log.f = log.f + log( hyperg_1F1(b + 1, c, 0.5*n*R2/vsigma2, give=FALSE, strict=TRUE) )


plot(vsigma2,exp(log.f),type="l")

print(trapint(vsigma2,exp(log.f)))


vsigma2.post = 1/rgamma(N,c,0.5*n*(1 - R2*vg/(1 + vg)))
lines(density(vsigma2.post),col="red")

dinvgamma <- function(x,a,b) {
	log.f = a*log(b) - lgamma(a) - (a + 1)*log(x) - b/x
	return(exp(log.f))
}

vf = 0*exp(log.f)
for (i in 1:N) {
	vf = vf + (1/N)*dinvgamma(vsigma2,c,0.5*n*(1 - R2*vg[i]/(1 + vg[i])))
}

lines(vsigma2,vf,col="blue")

# Rao-Blackwellised alpha|y
valpha <- seq(-0.3, 0.3,, 1e5)
vf_alpha = 0 * exp(log.f)
for (i in 1:N) {
	vf_alpha = vf_alpha + (1 / N) * dnorm(valpha, 0, sqrt(vsigma2.post[i]) / sqrt(n))
}
plot(valpha,vf_alpha,col="blue",type="l")

analytic <- exp(lgamma(c + 0.5) + (b + 1) * log(1 - R2) - lgamma(c) - 0.5 * log(pi) -n/2 * log(1 + valpha^2)) * hyperg_2F1(b + 1, c + 0.5, c, R2 / (1 + valpha^2))
lines(valpha, analytic, col="red")
legend("topright", lty=1, legend = c("Rao-Blackwellised estimator", "Analytic expression"), col = c("blue", "red"))

# Rao-Blackwellised beta|y
vbeta <- seq(-10, 10,, 1e5)
vf_beta = 0 * exp(log.f)
beta_hat <- 1.
for (i in 1:N) {
	vf_beta = vf_beta + (1 / N) * dnorm(vbeta, g[i]/(1 + g[i]) * beta_hat, sqrt(g[i]/(1 + g[i]) * vsigma2.post[i]))
}
plot(vbeta,vf_beta,col="blue",type="l")


log_K <- (b + 1) * log(1-R2) - lbeta(p/2 + a + 1, b + 1)
analytic <- exp(log_K + lbeta(p / 2 + a + 1, b + 2)) * hyperg_2F1(n/2, b + 2, n / 2 + 1, R2)
# analytic <- exp(log_K + lbeta(p / 2 + a + 1, b + 2)) * hyperg_2F1(p/2 + a + 1, 1, n / 2 + 1, R2)
lines(vbeta, analytic, col="red")
legend("topright", lty=1, legend = c("Rao-Blackwellised estimator", "Analytic expression"), col = c("blue", "red"))

trapint(vbeta, vbeta * vf_beta)
