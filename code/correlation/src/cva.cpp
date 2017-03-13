#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <vector>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include "graycode.h"
#include "correlation.h"

using namespace std;
using namespace Rcpp;


double prob(const uint n, const uint p, const double sigma2, const double a, const double b)
{
	double log_p = -n / 2. * log(sigma2) - p / 2. * log(n) + gsl_sf_lnbeta(a, b);

	return exp(log_p);
}


//' Run a Collapsed Variational Approximation to find the K best linear models
//'
//' @param vy Vector of responses
//' @param mX Matrix of covariates
//' @param K The number of particles in the population
//' @return Good question!
//' @export
// [[Rcpp::export]]
double cva(NumericVector vy_in, NumericMatrix mX_in, const uint K)
{
	VectorXd vy(vy_in.length()); // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
	for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
	MatrixXd mX(mX_in.nrow(), mX_in.ncol()); // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
	for (auto i = 0; i < mX_in.nrow(); i++)
		for (auto j = 0; j < mX_in.ncol(); j++)
			mX(i, j) = mX_in(i, j);
	const auto a = 1.;
	const auto b = 1.;
	const auto n = mX.rows();
	const auto p = mX.cols();
	const MatrixXd mXTX = mX.transpose() * mX;
	const MatrixXd mXTy = mX.transpose() * vy;
	VectorXd probs(K);
	VectorXd w(K);
	vector< dbitset > gamma(K);
	vector< MatrixXd > mXTX_inv(K);
	VectorXd sigma2(K);
	const gsl_rng_type *T;
	gsl_rng *r;
	const auto RHO = 0.1;
	auto f_lambda_prev = 0.;
	const auto EPSILON = 1e-8;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// Initialise population of K particles randomly
	for (auto k = 0; k < K; k++) {
		gamma[k].resize(p);
		// Empty models are pointless, so we keep regenerating gammas until we get a non-zero one
		while (gamma[k].count() == 0) {
			for (auto j = 0; j < p; j++) {
				if (gsl_rng_uniform(r) <= RHO)
					gamma[k][j] = true;
				else
					gamma[k][j] = false;
			}
		}
	}

	// Initialise mXTX_inv
	// Initialise sigma2
	for (auto k = 0; k < K; k++) {
		auto p_gamma = gamma[k].count();
		MatrixXd mX_gamma_Ty(p_gamma, 1);
		get_rows(mXTy, gamma[k], mX_gamma_Ty);
		sigma2[k] = 1. - (mX_gamma_Ty.transpose() * mXTX_inv[k] * mX_gamma_Ty).value() / n;
	}

	// Loop until convergence
	auto converged = false;
	while (!converged) {
		for (auto k = 0; k < K; k++) {
			for (auto j = 0; j < p; j++) {
				dbitset gamma_0 = gamma[k];
				gamma_0[j] = false;
				auto p_gamma_0 = gamma_0.count();
				MatrixXd mX_gamma_0(n, p_gamma_0);
				get_cols(mX, gamma_0, mX_gamma_0);
				MatrixXd mX_gamma_0_Ty = mX_gamma_0.transpose() * vy;

				dbitset gamma_1 = gamma[k];
				gamma_1[j] = true;
				auto p_gamma_1 = gamma_1.count();
				MatrixXd mX_gamma_1(n, p_gamma_1);
				get_cols(mX, gamma_1, mX_gamma_1);
				MatrixXd mX_gamma_1_Ty = mX_gamma_1.transpose() * vy;

				uint p_gamma;
				bool bUpdate = !gamma[k][j];
				MatrixXd mXTX_inv_prime;
				Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
				if (bUpdate) {
					p_gamma = p_gamma_1;
					Rcpp::Rcout << "p_gamma_1 " << p_gamma_1 << std::endl;
						// Update mXTX_inv
						mXTX_inv_prime.resize(p_gamma_1, p_gamma_1);
						bool bLow;
						rank_one_update(gamma_1, j, gamma_1.find_first(),
														0,
														mXTX,
														mXTX_inv[k],
														mXTX_inv_prime,
														bLow);
						mXTX_inv[k] = mXTX_inv_prime;
				} else {
					p_gamma = p_gamma_0;
					// Downdate mXTX_inv
					mXTX_inv_prime.resize(p_gamma_0, p_gamma_0);
					rank_one_downdate(j, gamma_0.find_first(), 0,
														mXTX_inv[k], mXTX_inv_prime);
				}

				MatrixXd mX_gamma_0TX_gamma_0_inv(p_gamma_0, p_gamma_0);
				MatrixXd mX_gamma_1TX_gamma_1_inv(p_gamma_1, p_gamma_1);

				MatrixXd mXgammaTXgamma_inv(p_gamma, p_gamma);
				MatrixXd mXgamma_Ty = get_rows(mXTy, gamma_1, mXgamma_Ty);
				double sigma2_0;
				double sigma2_1;

 				if (bUpdate){
 					sigma2_0 = sigma2[k];
					sigma2_1 = 1. - (mX_gamma_1_Ty.transpose() * mX_gamma_1TX_gamma_1_inv * mX_gamma_1_Ty).value() / n;
 				} else {
 					sigma2_0 = 1. - (mX_gamma_0_Ty.transpose() * mX_gamma_0TX_gamma_0_inv * mX_gamma_0_Ty).value() / n;
 					sigma2_1 = sigma2[k];
 				}

				double p_0 = prob(n, p, sigma2_0, a, b);
				double p_1 = prob(n, p, sigma2_1, a, b);
				if (p_0 > p_1) {
					gamma[k][j] = true;
					sigma2[k] = sigma2_0;
					if (!bUpdate) mXTX_inv[k] = mXTX_inv_prime;
				}	else {
					gamma[k][j] = false;
					sigma2[k] = sigma2_1;
					mXTX_inv[k] = mXTX_inv_prime;
					if (bUpdate) mXTX_inv[k] = mXTX_inv_prime;
				}
			}
		}

		for (auto k = 0; k < K; k++) {
			probs[k] = prob(n, p, sigma2[k], a, b);
		}

		for (auto k = 0; k < K; k++) {
			w[k] = probs[k] / probs.sum();
		}

		// Check for convergence - is f_lambda changed from the last iteration?
		const auto lambda = 1.;
		VectorXd p_l(K);
		p_l = VectorXd::Zero(K);
		for (auto l = 0; l < K; l++) {
			for (auto k = 0; k < K; k++) {
				p_l[l] += w[k] * (gamma[l][k] ? 1. : 0.);
			}
		}
		auto H = -(p_l.array() * p_l.array().log()).sum();
		double f_lambda = w.dot(probs) + lambda * H;
		if ((f_lambda - f_lambda_prev) > EPSILON) {
			f_lambda_prev = f_lambda;
		} else {
			converged = true;
		}
	}
}