// ExampleUScrime_main.cpp

#include <ctgmath>
#include <iostream>
#include <limits>
#include <Eigen/Dense>
#include <gsl/gsl_sf_gamma.h>
#include <gsl/gsl_sf_hyperg.h>
#include "graycode.hpp"
#include "correlation.hpp"

using namespace std;
using Eigen::DenseBase;

const auto PI  =3.141592653589793238463;
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
	return U*log(x) + B*log(1. + x) - C*x;
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
	return -U/x2 - B/xp12;
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
double Z_h_Laplace(double U,double B,double C)
{
	auto res = mode_qh(U,B,C);
	auto h_hat = res;
	auto qh_hat = exp(log_qh(h_hat,U,B,C));
	auto sigma2_inv = -hessian_qh(h_hat,U,B,C) ;
	auto sigma2 = 1/sigma2_inv;

	return sqrt(2*PI*sigma2)*qh_hat;
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
double moment_h_FullyExponentialLaplace(double m,double U,double B,double C)
{
	auto val1 = Z_h_Laplace(U+m,B,C);
	auto val2 = Z_h_Laplace(U,B,C);
	return(val1/val2);
}


// Approximate the normalizing constant for the q-density for g using trapezoidal
// integration

struct Trapint
{
	double intVal;
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

	// Use trapezoidal integration
	double intVal = trapint(vg,log_qg_vg.array().exp());

	// return(list(intVal=intVal,vg=vg,log_qg_vg=log_qg_vg))
	Trapint result;
	result.intVal = intVal;
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
	return exp(-.5*z) * pow(z, mu + .5) * gsl_sf_hyperg_U(mu - kappa + .5, 1. + 2. * mu, z);
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
	if ( ((1.-mu)>nu)&(nu>0.) ) {
		cout << "fine" << endl;
	}
	else {
		cout << "not fine" << endl;
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

double tau_plugin(int a, int n, int p, double R2)
{
	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	double tau = (1. - R2)*(1. + (0.5*p + a + 1.)/b) ;
	return tau;
}


// Apply the update for tau_g using trapezoidal integration for ITER iterations

double tau_g_trapint(double a, int n, int p, double R2, uint ITER, uint N=1000)
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

double tau_g_Laplace(double a, int n, int p, double R2)
{
	double A = 2.*a + p;
	double res = A*(1. - R2)/(n - A);
	return res;
}


// Apply the update for tau_g using FEL for ITER iterations

double tau_g_FullyExponentialLaplace(double a, int n, int p, double R2, uint ITER)
{
	double tau = tau_g_Laplace(a,n,p,R2);

	double b = (n-p)/2. - a - 2.;
	double A = b - p/2.;
	double B = -(n-p)/2.;

	double U =  -(A+B+2.);

	for (uint i = 0; i < ITER; i++) {
		double C = 0.5*n*R2/((1. + tau)*(1. - R2 + tau)) + 0.5*p/(1. + tau);
		// #cat(i,tau,C,"\n")
		// cout << i << " " << tau << " " << C << endl;
		tau = moment_h_FullyExponentialLaplace(1,U,B,C);
		//  #cat(i,tau,C,"\n")
		// cout << i << tau << C << endl;
		if (isnan(tau) || fabs(tau) == 0.0) {
			tau = tau_plugin(a,n,p,R2);
			break;
		}
	}

	return tau;
}


// Apply the update for tau_g using Laplace's method for ITER iterations

double tau_g_IterativeLaplace(double a, int n, int p, double R2, uint ITER)
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

ZE_constants_result ZE_constants(int n, int pmax, bool LARGEP = false)
{
	VectorXd vcon(pmax + 1);
	VectorXd vpen(pmax + 1);
	ZE_constants_result result;

	#pragma omp parallel for simd
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


struct ZE_exact_result
{
	MatrixXi mGraycode;
	VectorXd vR2;
	VectorXd vlog_ZE;
};

ZE_exact_result ZE_exact(VectorXd vy, MatrixXd mX)
{
	auto n = mX.rows();
	auto p = mX.cols();
	ZE_constants_result res_con = ZE_constants(n, p);
	Graycode graycode(p);
	auto mGraycode = graycode.to_MatrixXi();
	auto vR2 = all_correlations_mX_cpp(vy, mX, 0, true, true);
	VectorXi vq = mGraycode * MatrixXi::Ones(p, 1);
	VectorXd vlog_ZE(vq.size());
	#pragma omp parallel for simd
	for (auto i = 0; i < vlog_ZE.size(); i++) {
		auto p_gamma = vq(i);
		// cout << "i " << i << " p_gamma " << p_gamma << endl;
		vlog_ZE(i) = -res_con.vcon(p_gamma) * log(1.0 - vR2(i)) + res_con.vpen(p_gamma);
	}
	ZE_exact_result result;
	result.mGraycode = mGraycode;
	result.vR2 = vR2;
	result.vlog_ZE = vlog_ZE;
	return result;
}


int main()
{
	try
	{
		string mX_fname, vy_fname, vR2_fname, vlog_ZE_fname, mGraycode_fname;

		const auto SETTING = 2;
		if (SETTING == 1) {
			mX_fname = "Hitters_mX.csv";
			vy_fname = "Hitters_vy.csv";
			vR2_fname = "Hitters_vR2.csv";;
			vlog_ZE_fname = "Hitters_vlog_ZE.csv";
			mGraycode_fname = "Hitters_mGraycode.csv";
		}

		if (SETTING == 2) {
			mX_fname = "bodyfat_mX.csv";
			vy_fname = "bodyfat_vy.csv";
			vR2_fname = "bodyfat_vR2.csv";
			vlog_ZE_fname = "bodyfat_vlog_ZE.csv";
			mGraycode_fname = "bodyfat_mGraycode.csv";
		}

		if (SETTING == 3) {
			mX_fname = "Wage_mX.csv";
			vy_fname = "Wage_vy.csv";
			vR2_fname = "Wage_vR2.csv";
			vlog_ZE_fname = "Wage_vlog_ZE.csv";
			mGraycode_fname = "Wage_mGraycode.csv";
		}

		if (SETTING == 4) {
			mX_fname = "GradRate_mX.csv";
			vy_fname = "GradRate_vy.csv";
			vR2_fname = "GradRate_vR2.csv";
			vlog_ZE_fname = "GradRate_vlog_ZE.csv";
			mGraycode_fname = "GradRate_mGraycode.csv";
		}
		MatrixXd mX = parseCSVfile_double(mX_fname);
		VectorXd vy = parseCSVfile_double(vy_fname);
		VectorXd vR2 = parseCSVfile_double(vR2_fname);
		VectorXd vlog_ZE = parseCSVfile_double(vlog_ZE_fname);
		MatrixXd mGraycode = parseCSVfile_double(mGraycode_fname);

		auto n = mX.rows();
		auto p = mX.cols();
		auto logpy = vlog_ZE;
		auto logpy_til = logpy.array() - logpy.maxCoeff();

		// Calculate the marginal variable inclusion probability
		VectorXd vp = logpy_til.array().exp() / logpy_til.array().exp().sum();
		auto a = -0.75;
		VectorXd velbo(vR2.size());

		#pragma omp parallel for
		for (auto i = 0; i < vR2.size(); i++) {
			// Note res$mA contains the gray-code
			uint p = mGraycode.row(i).sum();
			double b = (n - p)/2. - a - 2. ;

			// Calculate
			double A =  n/2. - p - a - 2.;
			double B = -(n-p)/2.;
			double U = -(A+B+2.);

			// Calculate tau_g
			double tau_g = tau_g_FullyExponentialLaplace(a,n,p,vR2[i],20);

			// Calculate the constant C (needed to calculate the normalizing constant for q(g)
			double C = 0.5*n*vR2[i]/((1 + tau_g)*(1 - vR2[i] + tau_g)) + 0.5*p/(1 + tau_g);

			// Calculate the
								 // Z.h.Laplace(U,B,C)
			double Z = Z_g_trapint(A,B,C,1000).intVal;

			// Calculate the lower bound on the log-likelihood
			velbo[i] = 0.5*p - 0.5*n*log(2*PI) - gsl_sf_lnbeta(a+1,b+1)  - 0.5*n*log(1 + tau_g - vR2[i]) ;
			velbo[i] = velbo[i] - 0.5*(n+p)*log(0.5*(n+p))+ gsl_sf_lngamma(0.5*(n+p)) + C*tau_g + log(Z)  + 0.5*(n-p)*log(1 + tau_g);

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
		VectorXd vq = logqy_til.array().exp() / logqy_til.array().exp().sum();

		// Calculate the variable inclusion probabilities
		auto vw1 = vp.transpose() * mGraycode.cast<double>();
		auto vw2 = vq.transpose() * mGraycode.cast<double>();

		cout << vw1 << endl;
		cout << vw2 << endl;
	}
	catch (std::runtime_error& e) {
		cerr << "There was a problem loading the matrix files." << e.what() << endl;
		return -1;
	}

	return 0;
}
