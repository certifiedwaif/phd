// ExampleUScrime_main.cpp

#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include <gsl/gsl_sf_hyperg.h>

using Eigen::VectorXd;
using namespace std;

const double PI  =3.141592653589793238463;

// # Trapezoidal integration over a potentially irregular grid

// trapint <- function(xgrid, fgrid) 
// { 
// 	ng <- length(xgrid) 
// 	xvec <- xgrid[2:ng] - xgrid[1:(ng - 1)] 
// 	fvec <- fgrid[1:(ng - 1)] + fgrid[2:ng] 
// 	integ <- sum(xvec * fvec)/2 
// 	return(integ) 
// }   

double trapint(VectorXd xgrid, VectorXd fgrid)
{
	auto sum = 0.0;

	for (auto i = 0; i < xgrid.size(); i++) {
		sum += 0.5 * (xgrid(i + 1) - xgrid(i)) * (fgrid(i + 1) - fgrid(i));
	}

	return sum;
}

// # Solve a*x*x + b*x + c = 0
// # Probably a numerically bad way of doing this see Num Rec in C

// solveQuad <- function(a,b,c) 
// {
// 	disc <- b*b - 4*a*c
// 	#cat("disc=",disc,"\n")
// 	if (disc<0) {
// 		val <- c(NA)
// 	} else {
// 		val1 <-(-b + sqrt(disc))/(2*a)
// 		val2 <-(-b - sqrt(disc))/(2*a) 
// 		val <- c(val1,val2)
// 	}
// 	return( val )
// }


VectorXd solveQuad(double a, double b, double c)
{
	VectorXd val(2);
	auto disc = b*b - 4*a*c;
	if (disc < 0) {
		auto NEGATIVE_INFINITY = -numeric_limits<double>::infinity();
		val << NEGATIVE_INFINITY, NEGATIVE_INFINITY;
	} else {
		val << -(-b + sqrt(disc))/(2*a), -(-b - sqrt(disc))/(2*a);
	}
	return val;
}

// ################################################################################

// # The log the q-density for g

// log.qg <- function(x,A,B,C) {
// 	return( A*log(x) + B*log(1 + x) - C/x )	
// }

double log_qg(double x, double A, double B, double C)
{
	return A*log(x) + B*log(1 + x) - C/x;
}


VectorXd log_qg(VectorXd x, double A, double B, double C)
{
	VectorXd result(x.size());
	for (auto i = 0; i < x.size(); i++)
		result(i) = A*log(x(i)) + B*log(1 + x(i)) - C/x(i);
	return result;
}

// ################################################################################

// # The log the q-density for h (where h is the inverse of g)

// log.qh <- function(x,U,B,C) {
// 	return( U*log(x) + B*log(1 + x) - C*x )	
// }

double log_qh(double x, double U, double B, double C)
{
	return U*log(x) + B*log(1 + x) - C*x;
}

VectorXd log_qh(VectorXd x, double U, double B, double C)
{
	VectorXd result(x.size());
	for (auto i = 0; i < x.size(); i++)
		result(i) = U*log(x(i)) + B*log(1 + x(i)) - C*x(i);
	return result;
}


// # Calculate the mode of the q-density for g

// mode.qg <- function(A,B,C) {
// 	res <- solveQuad(A+B,A+C,C) 
// 	if (length(res)>1) {
// 		res <- res[res>0]
// 	}
// 	return(res)
// }

double mode_qg(double A, double B, double C)
{
	auto result = solveQuad(A+B, A+C, C);
	if (result(0) > 0.0) return result(0);
	if (result(1) > 0.0) return result(1);
	return result(1);
}

// ################################################################################

// # Calculate the mode of the q-density for h (where h is the inverse of g)

// mode.qh <- function(U,B,C) {
// 	res <- solveQuad(C,C-B-U,-U) 
// 	if (length(res)>1) {
// 		res <- res[res>0]
// 	}
// 	return(res)
// }

double mode_qh(double U, double B, double C)
{
	auto result = solveQuad(C, C-B-U, -U);
	if (result(0) > 0.0) return result(0);
	if (result(1) > 0.0) return result(1);
	return result(1);
}

// ################################################################################

// # Calculate the Hessian of the q-density for g evaluated at x

// hessian.qg <- function(x,A,B,C) 
// {
// 	x2 <- x*x
// 	x3 <- x2*x
// 	xp12 <- (1+x)*(1+x)
// 	val <-  -A/x2 - B/xp12 - 2*C/x3
// 	return(val)
// }

double hessian_qg(double x, double A, double B, double C)
{
	auto x2 = x*x;
	auto x3 = x2*x;
	auto xp12 = (1.0 + x) * (1.0 + x);
	return -A/x2 - B/xp12 - 2.0*C/x3;
}

// ################################################################################

// # Calculate the Hessian of the q-density for h evaluated at x (where h is the inverse of g)

// hessian.qh <- function(x,U,B,C) 
// {	
// 	x2 <- x*x
// 	xp12 <- (1+x)*(1+x)
// 	return( -U/x2 - B/xp12 )	
// }

double hessian_qh(double x, double U, double B, double C) 
{	
	auto x2 = x*x;
	auto xp12 = (1+x)*(1+x);
	return -U/x2 - B/xp12;
}

// ################################################################################

// # Calculate the Laplace approximation of the normalizing constant for the
// # q-density for g

// Z.g.Laplace <- function(A,B,C)
// {
// 	res <- mode.qg(A,B,C)
// 	g.hat <- res[1]
// 	qg.hat <- exp(log.qg(g.hat,A,B,C))
// 	sigma2.inv <- -hessian.qg(g.hat,A,B,C) 
// 	sigma2 <- 1/sigma2.inv
// 	return( sqrt(2*pi*sigma2)*qg.hat )
// }

double Z_g_Laplace(double A,double B,double C)
{
	auto res = mode_qg(A,B,C);
	auto g_hat = res;
	auto qg_hat = exp(log_qg(g_hat,A,B,C));
	auto sigma2_inv = -hessian_qg(g_hat,A,B,C) ;
	auto sigma2 = 1/sigma2_inv;

	return sqrt(2*PI*sigma2)*qg_hat;
}

// ################################################################################

// # Calculate the Laplace approximation of the normalizing constant for the
// # q-density for h (where h is the inverse of g)

// Z.h.Laplace <- function(U,B,C)
// {
// 	res <- mode.qh(U,B,C)
// 	h.hat <- res[1]
// 	qh.hat <- exp(log.qh(h.hat,U,B,C))
// 	sigma2.inv <- -hessian.qh(h.hat,U,B,C) 
// 	sigma2 <- 1/sigma2.inv
// 	return( sqrt(2*pi*sigma2)*qh.hat )
// }

double Z_h_Laplace(double U,double B,double C)
{
	auto res = mode_qh(U,B,C);
	auto h_hat = res;
	auto qh_hat = exp(log_qh(h_hat,U,B,C));
	auto sigma2_inv = -hessian_qh(h_hat,U,B,C) ;
	auto sigma2 = 1/sigma2_inv;

	return sqrt(2*PI*sigma2)*qh_hat;
}

// ################################################################################

// # Approximate the expected value of g with respect to q(g) using the mode

// E.g.Laplace <- function(A,B,C) {
// 	res <- mode.qg(A,B,C)
// 	g.hat <- res[1]
// 	return(g.hat)
// }

double E_g_Laplace(double A,double B,double C) {
	auto res = mode_qg(A,B,C);
	auto g_hat = res;
	return(g_hat);
}

// ################################################################################

// # Approximate the expected value of h with respect to q(h) using the mode

// E.h.Laplace <- function(U,B,C) {
// 	res <- mode.qh(U,B,C)
// 	h.hat <- res[1]
// 	return(h.hat)
// }

double E_h_Laplace(double U,double B,double C) {
	auto res = mode_qh(U,B,C);
	auto h_hat = res;
	return h_hat;
}

// ################################################################################

// # Approximate the mth moment of g using the Fully-Exponential-Laplace method

// moment.g.FullyExponentialLaplace <- function(m,A,B,C) {
// 	val1 <- Z.g.Laplace(A+m,B,C)
// 	val2 <- Z.g.Laplace(A,B,C)
// 	return(val1/val2)
// }

// ################################################################################

// # Approximate the mth moment of h using the Fully-Exponential-Laplace method

// moment.h.FullyExponentialLaplace <- function(m,U,B,C) {
// 	val1 <- Z.h.Laplace(U+m,B,C)
// 	val2 <- Z.h.Laplace(U,B,C)
// 	return(val1/val2)
// }

double moment_g_FullyExponentialLaplace(double m,double A,double B,double C) {
	auto val1 = Z_g_Laplace(A+m,B,C);
	auto val2 = Z_g_Laplace(A,B,C);
	return(val1/val2);
}

double moment_h_FullyExponentialLaplace(double m,double U,double B,double C) {
	auto val1 = Z_h_Laplace(U+m,B,C);
	auto val2 = Z_h_Laplace(U,B,C);
	return(val1/val2);
}

// ################################################################################

// # Approximate the normalizing constant for the q-density for g using trapezoidal
// # integration

// Z.g.trapint <- function(A,B,C,N=10000,PLOT=FALSE)
// {
// 	TOL <- 1.0E-14
// 	TOR <- 1.0E-5
	
// 	# Find the mode of g
// 	res <- mode.qg(A,B,C)
// 	g.hat <- res[1]
	
// 	# Evaluate q(g) at the mode
// 	log.qg.max <- log.qg(g.hat,A,B,C)
// 	qg.hat <- exp(log.qg.max)
	
// 	# Use a normal approximation at the mode
// 	sigma2.inv <- -hessian.qg(g.hat,A,B,C) 
// 	sigma2 <- 1/sigma2.inv
// 	delta   <- 0.5*sqrt(sigma2)
	
// 	# Loop using steps of size delta to the right until the difference in the 
// 	# log-densities between the maximum and the current location is smaller 
// 	# than TOR
// 	g.curr  <- g.hat + delta
// 	log.qg.curr <- log.qg(g.curr,A,B,C)
// 	#count <- 0	
// 	#cat(count,g.hat,g.curr,log.qg.curr,log.qg.max,log.qg.curr - log.qg.max,log(TOL),"\n")
// 	while ( (log.qg.curr - log.qg.max) > log(TOR) ) {
// 		#cat(count,g.hat,g.curr,log.qg.curr,log.qg.max,log.qg.curr - log.qg.max,log(TOL),"\n")
// 		#count <- count + 1
// 		g.curr <- g.curr + delta
// 		log.qg.curr <- log.qg(g.curr,A,B,C)
// 	}
// 	R <- g.curr
	
// 	# Loop using steps of size delta to the left until the difference in the 
// 	# log-densities between the maximum and the current location is smaller 
// 	# than TOL
// 	delta <- g.hat/4
// 	g.curr  <- g.hat/2
// 	log.qg.curr <- log.qg(g.curr,A,B,C)
// 	while ( (log.qg.curr - log.qg.max) > log(TOL) ) {
// 		while ((g.curr - delta)<=0) {
// 			# Reduce the step size if necessary
// 			delta <- delta/5
// 		}
// 		g.curr <- g.curr - delta
// 		log.qg.curr <- log.qg(g.curr,A,B,C)
// 	}
// 	L <- g.curr
	
// 	#print(L)
	
// 	# Calculate a grid between L and R with N points
// 	vg <- seq(L,R,,N)
// 	log.qg.vg <- log.qg(vg,A,B,C)
	
// 	# Use trapezoidal integration
// 	intVal <- trapint(vg,exp(log.qg.vg))
	
// 	# Plot the points for diagnostic purposes
// 	if (PLOT) {
// 		plot(vg,exp(log.qg.vg),type="l")
// 		points(g.hat,exp(log.qg.max))
// 	}
	
// 	return(list(intVal=intVal,vg=vg,log.qg.vg=log.qg.vg))
// }

//################################################################################

// # Approximate the normalizing constant for the q-density for g using trapezoidal
// # integration

struct Trapint {
	double intVal;
	VectorXd vg;
	VectorXd log_qg_vg;
};

VectorXd exp(VectorXd x) {
	VectorXd result(x.size());
	for (auto i = 0; i < x.size(); i++)
		result(i) = exp(x(i));
	return result;
}

VectorXd seq(double L, double R, size_t N)
{
	VectorXd result(N);
	auto delta = (R - L) / N;
	for (auto i = 0; i < N; i++) {
		result(i) = L + delta * i;
	}
	return result;
}

Trapint Z_g_trapint(double A, double B, double C, size_t N=10000)
{
	auto TOL = 1.0E-14;
	auto TOR = 1.0E-5;
	
	// Find the mode of g
	auto res = mode_qg(A,B,C);
	auto g_hat = res;
	
	// Evaluate q(g) at the mode
	auto log_qg_max = log_qg(g_hat,A,B,C);
	auto qg_hat = exp(log_qg_max);
	
	// Use a normal approximation at the mode
	auto sigma2_inv = -hessian_qg(g_hat,A,B,C);
	auto sigma2 = 1/sigma2_inv;
	auto delta_1   = 0.5*sqrt(sigma2);
	
	// Loop using steps of size delta to the right until the difference in the 
	// log-densities between the maximum and the current location is smaller 
	// than TOR
	auto g_curr  = g_hat + delta_1;
	auto log_qg_curr = log_qg(g_curr,A,B,C);
	//count = 0	
	//cat(count,g_hat,g_curr,log_qg_curr,log_qg_max,log_qg_curr - log_qg_max,log(TOL),"\n")
	while ( (log_qg_curr - log_qg_max) > log(TOR) ) {
		//cat(count,g_hat,g_curr,log_qg_curr,log_qg_max,log_qg_curr - log_qg_max,log(TOL),"\n")
		//count = count + 1
		g_curr = g_curr + delta_1;
		log_qg_curr = log_qg(g_curr,A,B,C);
	}
	auto R = g_curr;
	
	// Loop using steps of size delta to the left until the difference in the 
	// log-densities between the maximum and the current location is smaller 
	// than TOL
	auto delta_2 = g_hat/4;
	g_curr  = g_hat/2;
	log_qg_curr = log_qg(g_curr,A,B,C);
	while ( (log_qg_curr - log_qg_max) > log(TOL) ) {
		while ((g_curr - delta_2)<=0) {
			// Reduce the step size if necessary
			delta_2 = delta_2/5;
		}
		g_curr = g_curr - delta_2;
		log_qg_curr = log_qg(g_curr,A,B,C);
	}
	auto L = g_curr;
	
	//print(L)
	
	// Calculate a grid between L and R with N points
	// vg = seq(L,R,,N)
	VectorXd vg = seq(L, R, N);
	auto log_qg_vg = log_qg(vg,A,B,C);
	
	// Use trapezoidal integration
	double intVal = trapint(vg,exp(log_qg_vg));
	
	// return(list(intVal=intVal,vg=vg,log_qg_vg=log_qg_vg))
	Trapint result = {intVal, vg, log_qg_vg};
	return result;
}

// ################################################################################

// Z.h.trapint <- function(U,B,C,N=10000,PLOT=FALSE)
// {
// 	TOL <- 1.0E-5
// 	L <- 1.0E-3
// 	res <- mode.qh(U,B,C)
// 	h.hat <- res[1]
// 	log.qh.max <- log.qh(h.hat,U,B,C)
// 	log.qh.L <- log.qh(L,U,B,C)
// 	delta   <- h.hat
// 	h.curr  <- h.hat
// 	log.qh.curr <- log.qh(h.curr,U,B,C)
// 	while ( (log.qh.curr - log.qh.max)>log(TOL) ) {
// 		h.curr <- h.curr + delta
// 		log.qh.curr <- log.qh(h.curr,U,B,C)
// 	}
// 	vh <- seq(L,h.curr,,N)
// 	log.qh.vh <- log.qh(vh,U,B,C)
// 	intVal <- trapint(vh,exp(log.qh.vh))
	
// 	if (PLOT) {
// 		plot(vh,exp(log.qh.vh),type="l")
// 		points(h.hat,exp(log.qh.max))
// 	}

// 	return(intVal)
// }


double Z_h_trapint(double U, double B, double C, size_t N=10000)
{
	auto TOL = 1.0E-5;
	auto L = 1.0E-3;
	auto res = mode_qh(U,B,C);
	auto h_hat = res;
	auto log_qh_max = log_qh(h_hat,U,B,C);
	auto log_qh_L = log_qh(L,U,B,C);
	auto delta   = h_hat;
	auto h_curr  = h_hat;
	auto log_qh_curr = log_qh(h_curr,U,B,C);

	while ( (log_qh_curr - log_qh_max)>log(TOL) ) {
		h_curr = h_curr + delta;
		log_qh_curr = log_qh(h_curr,U,B,C);
	}

	// vh = seq(L,h_curr,,N)
	VectorXd vh = seq(L, h_curr, N);
	auto log_qh_vh = log_qh(vh,U,B,C);
	auto intVal = trapint(vh,exp(log_qh_vh));

	return intVal;
}

// ################################################################################

// # Calculate the mth moment of q(g) using trapezoidal integration

// moment.g.trapint <- function(m,A,B,C,N=1000000)
// {
// 	val1 <- Z.g.trapint(A+m,B,C,N)$intVal
// 	val2 <- Z.g.trapint(A,B,C,N)$intVal
// 	return(val1/val2)
// }

// ################################################################################

// # Calculate the mth moment of q(h) using trapezoidal integration

// moment.h.trapint <- function(m,U,B,C,N=1000000)
// {
// 	val1 <- Z.h.trapint(U+m,B,C,N)
// 	val2 <- Z.h.trapint(U,B,C,N)
// 	return(val1/val2)
// }

// Calculate the mth moment of q(g) using trapezoidal integration

double moment_g_trapint(double m, double A, double B, double C, size_t N=1000000)
{
	auto val1 = Z_g_trapint(A+m,B,C,N).intVal;
	auto val2 = Z_g_trapint(A,B,C,N).intVal;
	return val1/val2;
}

// Calculate the mth moment of q(h) using trapezoidal integration

double moment_h_trapint(double m, double U, double B, double C, size_t N=1000000)
{
	auto val1 = Z_h_trapint(U+m,B,C,N);
	auto val2 = Z_h_trapint(U,B,C,N);
	return val1/val2;
}

double whittakerW(double z, double kappa, double mu)
{
	return exp(-.5*z) * pow(z, mu + .5) * gsl_sf_hyperg_U(mu - kappa + .5, 1 + 2 * mu, z);
}

// # Calculate the normalizing constant using the exact result

// Z.g.exact <- function(A,B,C) 
// {
// 	nu <- A + 1
// 	mu <- B + 1
// 	beta <- C
	
// 	# Check the conditions under which the result holds
// 	if ( ((1-mu)>nu)&(nu>0) ) {
// 		#print("fine")
// 	} else {
// 		#print("not fine")
// 	}
	
// 	val1 <- 0.5*(nu-1)*log(beta) + lgamma(1 - mu - nu) + 0.5*beta
// 	val2 <- Re( whittakerW(beta, 0.5*(nu-1) + mu, -0.5*nu) )
	
// 	return( list(val=exp(val1)*val2,val1=val1,val2=val2) )
// }

// Calculate the normalizing constant using the exact result

struct Z_g_exact_result {
	double val;
	double val1;
	double val2;
};

Z_g_exact_result Z_g_exact(double A, double B, double C) 
{
	auto nu = A + 1;
	auto mu = B + 1;
	auto beta = C;
	
	// Check the conditions under which the result holds
	if ( ((1-mu)>nu)&(nu>0) ) {
		cout << "fine" << endl;
	} else {
		cout << "not fine" << endl;
	}
	
	auto val1 = 0.5*(nu-1)*log(beta) + lgamma(1 - mu - nu) + 0.5*beta;
	auto val2 = whittakerW(beta, 0.5*(nu-1) + mu, -0.5*nu);
	
	Z_g_exact_result result = {exp(val1)*val2,val1,val2};
	return result;
}

// ################################################################################

// # Calculate the mth moment of q(g) using the exact result

// moment.g.exact <- function(m,A,B,C)
// {
// 	res1 <- Z.g.exact(A+m,B,C)
// 	res2 <- Z.g.exact(A,B,C)
// 	return(list(val=res1$val/res2$val,res1=res1,res2=res2))
// }


// Calculate the mth moment of q(g) using the exact result

struct moment_g_exact_result {
	double val;
	Z_g_exact_result res1;
	Z_g_exact_result res2;
};

moment_g_exact_result moment_g_exact(double m, double A, double B, double C)
{
	auto res1 = Z_g_exact(A+m,B,C);
	auto res2 = Z_g_exact(A,B,C);
	moment_g_exact_result result = {res1.val/res2.val,res1,res2};
	return result;
}


// ################################################################################

// # Use the plug in approximation for tau

// tau.plugin  <- function(a,n,p,R2)
// {
// 	b <- (n-p)/2 - a - 2
// 	A <- b - p/2
// 	B <- -(n-p)/2
	
// 	tau <- (1 - R2)*(1 + (0.5*p + a + 1)/b) 
// 	return(tau)
// }

// Use the plug in approximation for tau

double tau_plugin(uint a, uint n, uint p, double R2)
{
	auto b = (n-p)/2 - a - 2;
	auto A = b - p/2;
	auto B = -(n-p)/2;
	
	auto tau = (1 - R2)*(1 + (0.5*p + a + 1)/b) ;
	return tau;
}

// ################################################################################

// # Apply the update for tau_g using trapezoidal integration for ITER iterations

// tau.g.trapint <- function(a,n,p,R2,ITER,N=1000)
// {
// 	tau <- tau.plugin(a,n,p,R2)
	
// 	b <- (n-p)/2 - a - 2
// 	A <- b - p/2
// 	B <- -(n-p)/2
	
// 	for (i in 1:ITER) {
// 		C <- 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau)
// 		tau <- moment.g.trapint(m=-1,A,B,C,N)
// 	}
	
// 	return(tau)		
// }


// Apply the update for tau_g using trapezoidal integration for ITER iterations

double tau_g_trapint(double a, uint n, uint p, double R2, uint ITER, uint N=1000)
{
	auto tau = tau_plugin(a,n,p,R2);
	
	auto b = (n-p)/2 - a - 2;
	auto A = b - p/2;
	auto B = -(n-p)/2;
	
	for (auto i = 0; i < ITER; i++) {
		auto C = 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau);
		tau = moment_g_trapint(-1,A,B,C,N);
	}
	
	return tau;
}

// ################################################################################

// # Use the Laplace approximation for tau_g

// tau.g.Laplace <- function(a,n,p,R2)
// {
// 	A <- 2*a + p
// 	res <- A*(1 - R2)/(n - A) 
// 	return(res)
// }

// Use the Laplace approximation for tau_g

double tau_g_Laplace(double a, uint n, uint p, double R2)
{
	auto A = 2*a + p;
	auto res = A*(1 - R2)/(n - A);
	return res;
}

// ################################################################################

// # Apply the update for tau_g using FEL for ITER iterations

// tau.g.FullyExponentialLaplace <- function(a,n,p,R2,ITER)
// {
// 	tau <- tau.g.Laplace(a,n,p,R2)
	
// 	b <- (n-p)/2 - a - 2
// 	A <- b - p/2
// 	B <- -(n-p)/2
	
// 	U <-  -(A+B+2)
	
// 	for (i in 1:ITER) {
// 		C <- 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau)
// 		#cat(i,tau,C,"\n")
// 		tau <- moment.h.FullyExponentialLaplace(m=1,U,B,C)
// 		#cat(i,tau,C,"\n")
// 		if (is.na(tau)) {
// 			tau <- tau.plugin(a,n,p,R2)
// 			break;
// 		}
// 	}
	
// 	return(tau)		
// }

// Apply the update for tau_g using FEL for ITER iterations

double tau_g_FullyExponentialLaplace(double a, uint n, uint p, double R2, uint ITER)
{
	auto tau = tau_g_Laplace(a,n,p,R2);
	
	auto b = (n-p)/2 - a - 2;
	auto A = b - p/2;
	auto B = -(n-p)/2;
	
	auto U =  -(A+B+2);
	
	for (auto i = 0; i < ITER; i++) {
		auto C = 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau);
		// #cat(i,tau,C,"\n")
		cout << i << tau << C << endl;
		tau = moment_h_FullyExponentialLaplace(1,U,B,C);
		//  #cat(i,tau,C,"\n")
		cout << i << tau << C << endl;
		// TODO: How to signal and cope with failure here?
		#ifdef FALSE
		if (is_na(tau)) {
			tau = tau_plugin(a,n,p,R2);
			break;
		}
		#endif
	}
	
	return tau;
}

// ################################################################################

// # Apply the update for tau_g using Laplace's method for ITER iterations

// tau.g.IterativeLaplace  <- function(a,n,p,R2,ITER) {

// 	tau <- tau.plugin(a,n,p,R2)
	
// 	b <- (n-p)/2 - a - 2
// 	A <- b - p/2
// 	B <- -(n-p)/2
	
// 	U <-  -(A+B+2)
	
// 	for (i in 1:ITER) {
// 	C <- 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau)
// 	tau <- E.h.Laplace(U,B,C)
// 	}
	
// 	return(tau)	
// }	

// Apply the update for tau_g using Laplace's method for ITER iterations

double tau_g_IterativeLaplace(double a, uint n, uint p, double R2, uint ITER) {

	auto tau = tau_plugin(a,n,p,R2);
	
	auto b = (n-p)/2 - a - 2;
	auto A = b - p/2;
	auto B = -(n-p)/2;
	
	auto U =  -(A+B+2);
	
	for (auto i = 0; i < ITER; i++) {
		auto C = 0.5*n*R2/((1 + tau)*(1 - R2 + tau)) + 0.5*p/(1 + tau);
		tau = E_h_Laplace(U,B,C);
	}
	
	return tau;
}	


int main()
{
	cout << "Hello world\n";

	return 0;
}