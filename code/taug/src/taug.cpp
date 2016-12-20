// taug.cpp

#include "taug.h"
#include <ctgmath>
#include <functional>
#include <limits>
#include <Eigen/Dense>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>

using namespace std;
using namespace Rcpp;

const auto NEGATIVE_INFINITY = -numeric_limits<double>::infinity();

// Trapezoidal integration over a potentially irregular grid
template <typename Derived1, typename Derived2>
double trapint(const DenseBase<Derived1>& xgrid, const DenseBase<Derived2>& fgrid)
{
	auto sum = 0.0;

	#pragma omp parallel for simd reduction(+:sum)
	for (auto i = 0; i < xgrid.size() - 1; i++) {
		sum += 0.5 * (xgrid(i + 1) - xgrid(i)) * (fgrid(i) + fgrid(i + 1));
	}

	return sum;
}


double trapint(std::function<double(int)> x, std::function<double(double)> f)
{
	auto sum = 0.0;

  #pragma omp parallel for simd reduction(+:sum)
	for (auto i = 0; i < GRID_POINTS - 1; i++) {
		sum += 0.5 * (x(i + 1) - x(i)) * (f(x(i + 1)) + f(x(i)));
    // #ifndef _OPENMP
	    // Rcout << sum << std::endl;
    // #endif
	}

	return sum;
}

// Solve a*x*x + b*x + c = 0
// Probably a numerically bad way of doing this see Num Rec in C
VectorXd solveQuad(double a, double b, double c)
{
	VectorXd val(2);
	auto disc = b*b - 4.*a*c;
	if (disc < 0) {
		val << NEGATIVE_INFINITY, NEGATIVE_INFINITY;
	}
	else {
		val << (-b + sqrt(disc))/(2.*a), (-b - sqrt(disc))/(2.*a);
	}
	return val;
}


// The log the q-density for g
double log_qg(double x, double A, double B, double C)
{
	return A*log(x) + B*log(1. + x) - C/x;
}


template <typename Derived>
VectorXd log_qg(const DenseBase<Derived>& x, double A, double B, double C)
{
	VectorXd result(x.size());
	#pragma omp parallel for simd
	for (auto i = 0; i < x.size(); i++)
		result(i) = A*log(x(i)) + B*log(1. + x(i)) - C/x(i);
	return result;
}


// The log the q-density for h (where h is the inverse of g)
double log_qh(double x, double U, double B, double C)
{
	auto result = U*log(x) + B*log(1. + x) - C*x;
	// Rcout << "log_qh: x " << x << " U " << U << " B " << B << " C " << C << " result " << result << std::endl;
	return result;
}


template <typename Derived>
VectorXd log_qh(const DenseBase<Derived>& x, double U, double B, double C)
{
	VectorXd result(x.size());
	#pragma omp parallel for simd
	for (auto i = 0; i < x.size(); i++)
		result(i) = U*log(x(i)) + B*log(1. + x(i)) - C*x(i);
	return result;
}


// Calculate the mode of the q-density for g
double mode_qg(double A, double B, double C)
{
	auto result = solveQuad(A+B, A+C, C);
	if (result(0) > 0.0) return result(0);
	else return result(1);
}


// Calculate the mode of the q-density for h (where h is the inverse of g)
double mode_qh(double U, double B, double C)
{
	auto result = solveQuad(C, C-B-U, -U);
	if (result(0) > 0.0) return result(0);
	else return result(1);
}


// Calculate the Hessian of the q-density for g evaluated at x
double hessian_qg(double x, double A, double B, double C)
{
	auto x2 = x*x;
	auto x3 = x2*x;
	auto xp12 = (1.0 + x) * (1.0 + x);
	return -A/x2 - B/xp12 - 2.0*C/x3;
}


// Calculate the Hessian of the q-density for h evaluated at x (where h is the inverse of g)
double hessian_qh(double x, double U, double B, double C)
{
	auto x2 = x*x;
	auto xp12 = (1+x)*(1+x);
	auto result = -U/x2 - B/xp12;
	// Rcout << "hessian_qh: x " << x << " U " << U << " B " << B << " C " << C << " result " << result << std::endl;
	return result;
}


// Calculate the Laplace approximation of the normalizing constant for the
// q-density for g
double Z_g_Laplace(double A,double B,double C)
{
	auto res = mode_qg(A,B,C);
	auto g_hat = res;
	auto qg_hat = exp(log_qg(g_hat,A,B,C));
	auto sigma2_inv = -hessian_qg(g_hat,A,B,C) ;
	auto sigma2 = 1/sigma2_inv;

	return sqrt(2*PI*sigma2)*qg_hat;
}


// Calculate the Laplace approximation of the normalizing constant for the
// q-density for h (where h is the inverse of g)
double log_Z_h_Laplace(double U,double B,double C)
{
	auto res = mode_qh(U,B,C);
	auto h_hat = res;
	auto log_qh_hat = log_qh(h_hat,U,B,C);
	auto sigma2_inv = -hessian_qh(h_hat,U,B,C) ;
	auto sigma2 = 1/sigma2_inv;
	// Rcout << "Z_h_Laplace: res " << res;
	// Rcout << " h_hat " << h_hat;
	// Rcout << " log_qh_hat " << log_qh_hat;
	// Rcout << " sigma2_inv " << sigma2_inv;
	// Rcout << " sigma2 " << sigma2 << std::endl;

	return log(sqrt(2*PI)) + log(sigma2) + log_qh_hat;
}


// Approximate the expected value of g with respect to q(g) using the mode
double E_g_Laplace(double A,double B,double C)
{
	auto res = mode_qg(A,B,C);
	auto g_hat = res;
	return(g_hat);
}


// Approximate the expected value of h with respect to q(h) using the mode
double E_h_Laplace(double U,double B,double C)
{
	auto res = mode_qh(U,B,C);
	auto h_hat = res;
	return h_hat;
}


// Approximate the mth moment of g using the Fully-Exponential-Laplace method
double moment_g_FullyExponentialLaplace(double m,double A,double B,double C)
{
	auto val1 = Z_g_Laplace(A+m,B,C);
	auto val2 = Z_g_Laplace(A,B,C);
	return(val1/val2);
}


// Approximate the mth moment of h using the Fully-Exponential-Laplace method
double moment_h_FullyExponentialLaplace(double m1, double m2, double U,double B,double C)
{
	auto val1 = log_Z_h_Laplace(U+m1,B+m2,C);
	auto val2 = log_Z_h_Laplace(U,B,C);
	// Rcout << "val1 " << val1 << " val2 " << val2 << std::endl;
	return exp(val1 - val2);
}


// Approximate the normalizing constant for the q-density for g using trapezoidal
// integration

struct Trapint
{
	double intVal;
	double log_intVal;
	VectorXd vg;
	VectorXd log_qg_vg;
};

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
	auto delta_2 = g_hat/4.;
	g_curr  = g_hat/2.;
	log_qg_curr = log_qg(g_curr,A,B,C);
	while ( (log_qg_curr - log_qg_max) > log(TOL) ) {
		while ((g_curr - delta_2)<=0.) {
			// Reduce the step size if necessary
			delta_2 = delta_2/5.;
		}
		g_curr = g_curr - delta_2;
		log_qg_curr = log_qg(g_curr,A,B,C);
	}
	auto L = g_curr;

	//print(L)

	// Calculate a grid between L and R with N points
	// vg = seq(L,R,,N)
	const VectorXd vg = VectorXd::LinSpaced(N, L, R);
	auto log_qg_vg = log_qg(vg,A,B,C);
	auto min_log_qg = log_qg_vg.minCoeff();

	// Use trapezoidal integration
	double intVal = trapint(vg,(log_qg_vg.array() - min_log_qg).exp());

	// return(list(intVal=intVal,vg=vg,log_qg_vg=log_qg_vg))
	Trapint result;
	result.intVal = exp(min_log_qg) * intVal;
	result.log_intVal = log(intVal) + min_log_qg;
	result.vg = vg;
	result.log_qg_vg = log_qg_vg;
	return result;
}


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
	const VectorXd vh = VectorXd::LinSpaced(N, L, h_curr);
	auto log_qh_vh = log_qh(vh,U,B,C);
	auto intVal = trapint(vh,log_qh_vh.array().exp());

	return intVal;
}


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
	gsl_sf_result_e10 gsl_result;
	gsl_sf_hyperg_U_e10_e(mu - kappa + .5, 1. + 2. * mu, z, &gsl_result);

	auto U = gsl_result.val * pow(10, gsl_result.e10);
	auto result = exp(-.5*z) * pow(z, mu + .5) * U;

	// Rcout << "whittakerW: z " << z;
	// Rcout << " kappa " << kappa;
	// Rcout << " mu " << mu;
	// Rcout << " gsl_result.val " << gsl_result.val;
	// Rcout << " gsl_result.err " << gsl_result.err;
	// Rcout << " gsl_result.e10 " << gsl_result.e10;
	// Rcout << " U " << U;
	// Rcout << " result " << result << endl;

	return result;
}


// Calculate the normalizing constant using the exact result

struct Z_g_exact_result
{
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
	if ( ((1.-mu)>nu)&&(nu>0.) ) {
		Rcpp::Rcout << "fine" << endl;
	}
	else {
		Rcpp::Rcout << "not fine" << endl;
	}

	auto val1 = 0.5*(nu-1)*log(beta) + lgamma(1. - mu - nu) + 0.5*beta;
	auto val2 = whittakerW(beta, 0.5*(nu-1) + mu, -0.5*nu);

	Z_g_exact_result result = {exp(val1)*val2,val1,val2};
	return result;
}


// Calculate the mth moment of q(g) using the exact result

struct moment_g_exact_result
{
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


// Use the plug in approximation for tau

double tau_plugin(double a, double n, double p, double R2)
{
	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	double tau = (1. - R2)*(1. + (0.5*p + a + 1.)/b) ;
	return tau;
}


// Apply the update for tau_g using trapezoidal integration for ITER iterations

double tau_g_trapint(double a, double n, double p, double R2, uint ITER, uint N=1000)
{
	double tau = tau_plugin(a,n,p,R2);

	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	for (uint i = 0; i < ITER; i++) {
		double C = 0.5*n*R2/((1. + tau)*(1. - R2 + tau)) + 0.5*p/(1. + tau);
		tau = moment_g_trapint(-1.,A,B,C,N);
	}

	return tau;
}


// Use the Laplace approximation for tau_g

double tau_g_Laplace(double a, double n, double p, double R2)
{
	double A = 2.*a + p;
	double res = A*(1. - R2)/(n - A);
	return res;
}


// Apply the update for tau_g using FEL for ITER iterations

double tau_g_FullyExponentialLaplace(double a, double n, double p, double R2, uint ITER)
{
	double tau = tau_g_Laplace(a,n,p,R2);
	// Rcout << "Initial value of tau " << tau << std::endl;

	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	double U =  -(A+B+2.);

	for (uint i = 0; i < ITER; i++) {
		double C = 0.5*n*R2/((1. + tau)*(1. - R2 + tau)) + 0.5*p/(1. + tau);
		// #cat(i,tau,C,"\n")
		tau = moment_h_FullyExponentialLaplace(1., 0.,U,B,C);
		//  #cat(i,tau,C,"\n")
		// Rcout << "Iteration " << i << " C " << C << " tau " << tau << std::endl;
		if (isnan(tau) || fabs(tau) == 0.0) {
			// Rcout << "In NAN branch, using plug-in estimate: tau ";
			tau = tau_plugin(a,n,p,R2);
			// Rcout << tau << std::endl;
			break;
		}
	}

	return tau;
}


//' Variance of g over one plug g under q
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation co-efficient, squared
//' @return The variance of g/(1+g) under q
//' @export
// [[Rcpp::export]]
double var_g_over_one_plus_g(double n, double p, double R2)
{
	auto ITER = 20;
	auto a = -0.75;
	double tau = tau_g_Laplace(a,n,p,R2);
	// Rcout << "Initial value of tau " << tau << std::endl;

	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	double U =  -(A+B+2.);

	double E_g_over_one_plus_g;
	double E_g_over_one_plus_g_2;

	for (uint i = 0; i < ITER; i++) {
		double C = 0.5*n*R2/((1. + tau)*(1. - R2 + tau)) + 0.5*p/(1. + tau);
		tau = moment_h_FullyExponentialLaplace(1.,0.,U,B,C);
		E_g_over_one_plus_g = moment_h_FullyExponentialLaplace(0.,1.,U,B,C);
		E_g_over_one_plus_g_2 = moment_h_FullyExponentialLaplace(0.,2.,U,B,C);
		// Rcout << "var_g_over_one_plus_g:";
		// Rcout << " tau " << tau;
		// Rcout << " E_g_over_one_plus_g " << E_g_over_one_plus_g;
		// Rcout << " E_g_over_one_plus_g_2 " << E_g_over_one_plus_g_2 << std::endl;

		// Rcout << "Iteration " << i << " C " << C << " tau " << tau << std::endl;
		if (isnan(tau) || fabs(tau) == 0.0) {
			// Rcout << "In NAN branch, using plug-in estimate: tau ";
			tau = tau_plugin(a,n,p,R2);
			// Rcout << tau << std::endl;
			break;
		}
	}

	return (E_g_over_one_plus_g_2 - E_g_over_one_plus_g * E_g_over_one_plus_g);
}


// Apply the update for tau_g using Laplace's method for ITER iterations

double tau_g_IterativeLaplace(double a, double n, double p, double R2, uint ITER)
{

	double tau = tau_plugin(a,n,p,R2);

	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	double U =  -(A+B+2.);

	for (uint i = 0; i < ITER; i++) {
		double C = 0.5*n*R2/((1. + tau)*(1. - R2 + tau)) + 0.5*p/(1. + tau);
		tau = E_h_Laplace(U,B,C);
	}

	return tau;
}


struct ZE_constants_result
{
	VectorXd vcon;
	VectorXd vpen;
};


ZE_constants_result ZE_constants(double n, int pmax, bool LARGEP = false)
{
	VectorXd vcon(pmax + 1);
	VectorXd vpen(pmax + 1);
	ZE_constants_result result;

	#pragma omp parallel for
	for (auto q = 0; q < pmax + 1; q++) {
		double con = .5 * (n - q) - 0.75;
		auto dof = q;
		double pen = gsl_sf_lnbeta(0.5*q + 0.25,con) - gsl_sf_lnbeta(0.25,con);
		vcon(q) = con;
		vpen(q) = pen;
	}

	result.vcon = vcon;
	result.vpen = vpen;
	return result;
}


VectorXd ZE_exact(VectorXd vn, VectorXd vp, VectorXd vR2)
{
	VectorXd vlog_ZE(vR2.size());
	#pragma omp parallel for
	for (auto i = 0; i < vlog_ZE.size(); i++) {
		auto n = vn(i);
		auto p = vp(i);
		ZE_constants_result res_con = ZE_constants(n, p);
		// Rcpp::Rcout << "i " << i << " p_gamma " << p_gamma << endl;
		vlog_ZE(i) = -res_con.vcon(p) * log(1.0 - vR2(i)) + res_con.vpen(p);
	}
	return vlog_ZE;
}


//' Calculate tau_g from n, p and R2
//'
//' @param n The number of observations
//' @param p The number of parameters
//' @param R2 The correlation co-efficient
//' @return The log-lower bound of tau_g
//'
//' @export
// [[Rcpp::export]]
double tau_g(double n, double p, double R2)
{
	auto a = -0.75;

	// Note res$mA contains the gray-code
	double b = (n - p)/2. - a - 2. ;

	// Calculate
	double A =  n/2. - p - a - 2.;
	double B = -static_cast<double>(n-p)/2.;
	double U = -(A+B+2.);

	// Calculate tau_g
	double tau_g = tau_g_FullyExponentialLaplace(a, n, p, R2, 20);
	// Rcpp::Rcout << "n " << n << " p " << p << " tau_g " << tau_g << std::endl;
	return tau_g;
}


void tau_g(int n, const MatrixXd& mGraycode, const VectorXd& vR2, const VectorXd& vlog_ZE,
						VectorXd& vp, VectorXd& vq)
{
	auto logpy = vlog_ZE;
	auto logpy_til = logpy.array() - logpy.maxCoeff();

	// Calculate the marginal variable inclusion probability
	vp = logpy_til.array().exp() / logpy_til.array().exp().sum();
	auto a = -0.75;
	VectorXd velbo(vR2.size());

	#pragma omp parallel for
	for (auto i = 0; i < vR2.size(); i++) {
		// Note res$mA contains the gray-code
		uint p = mGraycode.row(i).sum();
		double b = (n - p)/2. - a - 2. ;

		// Calculate
		double A =  n/2. - p - a - 2.;
		double B = -static_cast<double>(n-p)/2.;
		double U = -(A+B+2.);

		// Calculate tau_g
		double tau_g = tau_g_FullyExponentialLaplace(a,n,p,vR2(i),20);

		// Calculate the constant C (needed to calculate the normalizing constant for q(g)
		double C = 0.5*n*vR2(i)/((1 + tau_g)*(1 - vR2(i) + tau_g)) + 0.5*p/(1 + tau_g);

		// Calculate the
							 // Z.h.Laplace(U,B,C)
		double log_Z = Z_g_trapint(A,B,C,1000).log_intVal;

		// Calculate the lower bound on the log-likelihood
		velbo(i) = 0.5*p - 0.5*n*log(2*PI) - gsl_sf_lnbeta(a+1,b+1)  - 0.5*n*log(1 + tau_g - vR2(i)) ;
		velbo(i) = velbo(i) - 0.5*(n+p)*log(0.5*(n+p))+ gsl_sf_lngamma(0.5*(n+p)) + C*tau_g + log_Z + 0.5*(n-p)*log(1 + tau_g);

		// cout << " " << i << " " << velbo[i] << " " << tau_g << " " << C << " " << Z << "\n";

		// How are failures from Z_g_trapint signalled?
		// If there is an error stop here and have a look
		// if (Z.failed) {
		// 	print("error! press escape and have a look")
		// 	string ans;
		// 	cin >> ans;
		// }

	}

	auto logqy_til = velbo.array() - velbo.maxCoeff();
	vq = logqy_til.array().exp() / logqy_til.array().exp().sum();

}


//' Calculate marginal log-likelihood log p(y)
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation co-efficient, squared
//' @return The value of the marginal log-likelihood log p(y)
//' @export
// [[Rcpp::export]]
double log_p(double n, double p, double R2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;

  auto result = gsl_sf_lngamma(p/2. + a + 1.) - n/2. * log(n * PI) + gsl_sf_lngamma((n - p) / 2.);
  result = result - gsl_sf_lngamma(a + 1.) - ((n - p) / 2. - a - 1.) * log(1 - R2);

  return result;
}


//' Calculate AIC
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation co-efficient, squared
//' @return The value of the marginal log-likelihood log p(y)
//' @export
// [[Rcpp::export]]
double aic(double n, double p, double R2)
{
	return -2 * log_p(n, p, R2) + 2 * p;
}


//' Calculate BIC
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation co-efficient, squared
//' @return The value of the marginal log-likelihood log p(y)
//' @export
// [[Rcpp::export]]
double bic(double n, double p, double R2)
{
	return -2 * log_p(n, p, R2) + p * log(n);
}

//' Calculate the variational lower bound
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param c The variational parameter c
//' @param s The variational parameter s
//' @param tau_g The variational parameter tau_g
//' @param log_det_XTX The log of the determinant of X^T X
//' @param log_det_mSigma The log of the determinant of mSigma
//' @return The variational lower bound
//' @export
// [[Rcpp::export]]
double elbo(double n, double p, double c, double s, double tau_g, double log_det_XTX, double log_det_mSigma)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
  auto beta = c;
  auto nu = n / 2. - p - a - 1.;
  auto mu = -(n - p/2.) + 1.;
  auto old_Z = pow(beta, (nu - 1.) / 2.) * gsl_sf_gamma(1. - mu - nu) * exp(- beta/(2.*nu)) * whittakerW(beta, (nu - 1.)/2 + mu, -nu / 2.);

	// Calculate
	double A =  n/2. - p - a - 2.;
	double B = -static_cast<double>(n-p)/2.;
	double U = -(A+B+2.);

	// Calculate the
						 // Z.h.Laplace(U,B,C)
	double log_new_Z = Z_g_trapint(A,B,c,1000).log_intVal;

  auto result =  p / 2. - n / 2. * log(2 * PI) + 0.5 * log_det_XTX - gsl_sf_lnbeta(a + 1., b + 1.) \
          - (n + p) / 2. * log(s) + gsl_sf_lngamma((n + p) / 2.) + c * tau_g + log_new_Z + 0.5 * log_det_mSigma;
  // Rcout << "elbo: n " << n;
  // Rcout << " p " << p;
  // Rcout << " c " << c;
  // Rcout << " s " << s;
  // Rcout << " tau_g " << tau_g;
  // Rcout << " log_det_XTX " << log_det_XTX;
  // Rcout << " log_det_mSigma " << log_det_mSigma;
  // Rcout << " a " << a;
  // Rcout << " b " << b;
  // Rcout << " beta " << beta;
  // Rcout << " nu " << nu;
  // Rcout << " mu " << mu;
  // Rcout << " old_Z " << old_Z;
  // Rcout << " log_new_Z " << log_new_Z;
  // Rcout << " result " << result << std::endl;

  return result;
}


double fall_fact(double x, double k)
{
	return gsl_sf_gamma(x + 1) / gsl_sf_gamma(x - k + 1);
}


double hyperg_1F2(int a1, int b1, int b2, double z)
{
	double sum = 0.;
	for (auto k = 0; k < 10; k++) {
		sum += (fall_fact(a1, k) * pow(z, k)) / (fall_fact(b1, k) * fall_fact(b2, k) * gsl_sf_gamma(k + 1));
	}
	return sum;
}


//' Calculate the exact precision
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation coefficient squared
//' @return The exact precision
//' @export
// [[Rcpp::export]]
double exact_precision(double n, double p, double R2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;

	return pow((1-R2), b + 1.) * gsl_sf_hyperg_2F1(n / 2. + 1., b + 1., n/2., R2);
	// return pow((1-R2), b + 1.) * hyperg_1F2(n / 2. + 1, b + 1., n/2., R2);
}


//' Calculate the exact posterior expectation of g
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation coefficient squared
//' @return The exact posterior expectation of g
//' @export
// [[Rcpp::export]]
double E_g_y(double n, double p, double R2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
	auto A =  n/2. - p - a - 2.;
	auto B = -static_cast<double>(n-p)/2.;

	return (b + 2.) / ((1. - R2) * (p / 2. + a + 1.));
	// const auto GRID_POINTS = 10000;
	// VectorXd x(GRID_POINTS), f(GRID_POINTS);
	// for (auto i = 1; i < GRID_POINTS; i++) {
	// 	double g = i / 100.;
	// 	x(i) = g;
	// 	// f(i) = pow((g/(1 + g)), (n - p) / 2.) * pow(g, -p/2. - a - 1.) * exp(-c/g);
	// 	// f(i) = exp((n - p) / 2. * (log(g) - log(1. + g)) + (-p/2. - a - 1.) * log(g) + -c/g);
	// 	// f(i) = pow(g, b + 1.) * pow(1 + g * (1 - R2), -n / 2.);
	// 	f(i) = exp((b + 1.) * log(g) - (n / 2.) * log(1 + g * (1 - R2)));
	// 	// Rcout << x(i) << " " <<  f(i) << std::endl;
	// }
	// return pow(1. - R2, b + 1.) * trapint(x, f) / gsl_sf_beta(p / 2. + a + 1., b + 1.);
}


//' Calculate the approximate posterior expectation of g
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation coefficient squared
//' @param c The co-efficient c of the exponential function exp(-c/g)
//' @return The approximate posterior expectation of g
//' @export
// [[Rcpp::export]]
double E_q_g_c(double n, double p, double R2, double c)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
	auto A =  n/2. - p - a - 2.;
	auto B = -static_cast<double>(n-p)/2.;
	auto log_new_Z = Z_g_trapint(A,B,c,1000).log_intVal;

	auto result = exp(-log_new_Z +  Z_g_trapint(A + 1., B, c, 1000).log_intVal);
	return result;

	// Rcout << "n " << n << " p " << p << " R2 " << R2 << " c " << c <<  std::endl;
	// auto result = exp(-log_new_Z + c) * pow(c, (n/2. - p - a - 1) / 2.) * gsl_sf_gamma((n - p) / 2. - n / 2. + p + a) \
	// 			 * whittakerW(-3. * n / 4. - a / 2. + 1. / 2., -n / 4. + p / 2. + a / 2.,  c);;
	// auto result = whittakerW(-3. * n / 4. - a / 2. + 1. / 2., -n / 4. + p / 2. + a / 2.,  c);
	// const auto GRID_POINTS = 10000;
	// VectorXd x(GRID_POINTS), f(GRID_POINTS);
	// for (auto i = 1; i < GRID_POINTS; i++) {
	// 	double g = i / 100.;
	// 	x(i) = g;
	// 	// f(i) = pow((g/(1 + g)), (n - p) / 2.) * pow(g, -p/2. - a - 1.) * exp(-c/g);
	// 	// f(i) = exp((n - p) / 2. * (log(g) - log(1. + g)) + (-p/2. - a - 1.) * log(g) + -c/g);
	// 	f(i) = exp((n / 2. - p - a - 1.) * log(g) - ((n - p) / 2.) * log(1. + g) - c/g);
	// 	// Rcout << x(i) << " " <<  f(i) << std::endl;
	// }
	// auto result = trapint(x, f);
	// Rcout << " = " << result << std::endl;
	// return moment_g_FullyExponentialLaplace(1.0, A,  B, c);
}


double p_sigma2_y(double n, double p, double R2, double sigma2)
{
  auto a = -3./4.;
  auto b = (n - p) / 2. - 2. - a;
	auto result =  pow(1. - R2, b + 1.) * pow(n / 2., n / 2.) * pow(sigma2, -(n / 2. + 1.)) / gsl_sf_gamma(n / 2.) \
				* exp(- n / (2. * sigma2)) * gsl_sf_hyperg_1F1(b + 1., n / 2., (n * R2) / (2 * sigma2));

	if (isnan(result)) {
	  result = exp((b + 1.) * log(1. - R2) + n / 2. * log(n / 2.) - (n / 2. + 1.) * log(sigma2) \
	                - gsl_sf_lngamma(n / 2.) - n / (2. * sigma2) \
	                + log(gsl_sf_hyperg_1F1(b + 1., n / 2., (n * R2) / (2 * sigma2))));
	}

	if (isnan(result)) {
		Rcout << "p_sigma2_y: n " << n;
		Rcout << " p " << p;
		Rcout << " R2 " << R2;
		Rcout << " sigma2 " << sigma2;
		Rcout << " = " << result << std::endl;
	}

	return result;
}


double q_sigma2(double r, double s, double sigma2)
{
	auto result = pow(s, r) / gsl_sf_gamma(r) * pow(sigma2, -r - 1.) * exp(- s / sigma2);

	if (isnan(result)) {
		result = exp(r * log(s) - gsl_sf_lngamma(r) - (r + 1.) * log(sigma2) - s / sigma2);
	}

	if (isnan(result)) {
		Rcout << "q_sigma2: r " << r;
		Rcout << " s " << s;
		Rcout << " sigma2 " << sigma2;
		Rcout << " = " << result << std::endl;
	}

	return result;
}


//' Calculate the accuracy of the approximation of sigma2
//'
//' @param n The number of observations
//' @param p The number of covariates
//' @param R2 The correlation co-efficient, squared
//' @param r The first parameter of the q(sigma2) distribution
//' @param s The second parameter of the q(sigma2) distribution
//' @return The accuracy of the approximation of sigma2
//' @export
// [[Rcpp::export]]
double accuracy_sigma2(double n, double p, double R2, double r, double s)
{
	const auto GRID_POINTS = 1000;
	VectorXd x(GRID_POINTS), f(GRID_POINTS);
	for (auto i = 1; i < GRID_POINTS; i++) {
		double sigma2 = i / 100.;
		x(i) = sigma2;
		f(i) = abs(p_sigma2_y(n, p, R2, sigma2) - q_sigma2(r, s, sigma2));
		// Rcout << x(i) << " " <<  f(i) << std::endl;
	}

	return 1. - 0.5 * trapint(x, f);
	// return trapint(x, f);
}