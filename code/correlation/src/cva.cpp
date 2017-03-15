#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
// [[Rcpp::depends(RcppGSL)]]
#include <RcppGSL.h>
#include <vector>
#include <cmath>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_sf.h>
#include "graycode.h"
#include "correlation.h"

using namespace std;
using namespace Rcpp;


double prob(const int n, const int p, const double sigma2, const double a, const double b)
{
	double log_sigma2 = std::log(sigma2);
	double log_n = std::log(n);
	// Rcpp::Rcout << "log_sigma2 " << log_sigma2 << std::endl;
	// Rcpp::Rcout << "log_n " << log_n << std::endl;
	double log_p = -n / 2. * log_sigma2 - p / 2. * log_n + gsl_sf_lnbeta(a, b);
	// Rcpp::Rcout << "log_p " << log_p << std::endl;

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
void cva(NumericVector vy_in, NumericMatrix mX_in, const int K)
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

	// // Centre vy
	// centre(vy);

	// // Centre non-intercept columns of mX
	// for (uint i = 0; i < mX.cols(); i++) {
	// 	VectorXd vcol = mX.col(i);
	// 	centre(vcol);
	// 	mX.col(i) = vcol;
	// }

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
		Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
	}

	// Initialise mXTX_inv
	// Initialise sigma2
	for (auto k = 0; k < K; k++) {
		auto p_gamma = gamma[k].count();
		MatrixXd mX_gamma(n, p_gamma);
		get_cols(mX, gamma[k], mX_gamma);
		MatrixXd mX_gamma_Ty(p_gamma, 1);
		get_rows(mXTy, gamma[k], mX_gamma_Ty);
		mXTX_inv[k] = (mX_gamma.transpose() * mX_gamma).inverse();
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

				dbitset gamma_1 = gamma[k];
				gamma_1[j] = true;
				auto p_gamma_1 = gamma_1.count();

				bool bUpdate = !gamma[k][j];
				MatrixXd mXTX_inv_prime;
				Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
				if (bUpdate) {
					Rcpp::Rcout << "Updating " << j << std::endl;
					Rcpp::Rcout << "p_gamma_1 " << p_gamma_1 << std::endl;
					// Update mXTX_inv
					mXTX_inv_prime.resize(p_gamma_1, p_gamma_1);
					bool bLow; // gamma_1 or gamma[k]?
					rank_one_update(gamma[k], j, gamma_1.find_first(),
													0,
													mXTX,
													mXTX_inv[k],
													mXTX_inv_prime,
													bLow);
					Rcpp::Rcout << "mXTX_inv[k]\n " << mXTX_inv[k] << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime\n " << mXTX_inv_prime << std::endl;
					Rcpp::Rcout << "bLow " << bLow << std::endl;
					// mXTX_inv[k] = mXTX_inv_prime;
				} else {
					// If only one bit is set and downdating would lead to the null model, then we should skip
					// checking the downdate
					if (p_gamma_0 == 0)
						continue;
					Rcpp::Rcout << "Downdating " << j << std::endl;
					Rcpp::Rcout << "p_gamma_0 " << p_gamma_0 << std::endl;
					// Downdate mXTX_inv
					mXTX_inv_prime.resize(p_gamma_0, p_gamma_0);
					rank_one_downdate(j, gamma_0.find_first(), 0,
														mXTX_inv[k], mXTX_inv_prime);
					Rcpp::Rcout << "mXTX_inv[k]\n " << mXTX_inv[k] << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime\n " << mXTX_inv_prime << std::endl;
				}

				double sigma2_0;
				double sigma2_1;

 				if (bUpdate){
					MatrixXd mX_gamma_1(n, p_gamma_1);
					get_cols(mX, gamma_1, mX_gamma_1);
					MatrixXd mX_gamma_1_Ty = mX_gamma_1.transpose() * vy;
					// MatrixXd mX_gamma_1_Ty(n, p_gamma_1);
					// get_rows(mXTy, gamma_1, mX_gamma_1_Ty);
 					sigma2_0 = sigma2[k];
					Rcpp::Rcout << "mX_gamma_1.transpose() * mX_gamma_1\n" << mX_gamma_1.transpose() * mX_gamma_1 << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_1.transpose() * mX_gamma_1\n" << mXTX_inv_prime * mX_gamma_1.transpose() * mX_gamma_1 << std::endl;
 					Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_1_Ty " << mXTX_inv_prime * mX_gamma_1_Ty << std::endl;
 					Rcpp::Rcout << "(mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value() " << (mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value() << std::endl;
					sigma2_1 = 1. - (mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value() / n;
 				} else {
					MatrixXd mX_gamma_0(n, p_gamma_0);
					get_cols(mX, gamma_0, mX_gamma_0);
					MatrixXd mX_gamma_0_Ty = mX_gamma_0.transpose() * vy;
					// MatrixXd mX_gamma_0_Ty(n, p_gamma_0);
					// get_rows(mXTy, gamma_0, mX_gamma_0_Ty);
					Rcpp::Rcout << "mX_gamma_0.transpose() * mX_gamma_0\n" << mX_gamma_0.transpose() * mX_gamma_0 << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_0.transpose() * mX_gamma_0\n" << mXTX_inv_prime * mX_gamma_0.transpose() * mX_gamma_0 << std::endl;
 					Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_0_Ty " << mXTX_inv_prime * mX_gamma_0_Ty << std::endl;
 					Rcpp::Rcout << "(mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value() " << (mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value() << std::endl;
 					sigma2_0 = 1. - (mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value() / n;
 					sigma2_1 = sigma2[k];
 				}

				Rcpp::Rcout << "sigma2_0 " << sigma2_0 << std::endl;
				Rcpp::Rcout << "sigma2_1 " << sigma2_1 << std::endl;
				double p_0 = prob(n, p_gamma_0, sigma2_0, a, b);
				double p_1 = prob(n, p_gamma_1, sigma2_1, a, b);
				Rcpp::Rcout << "p_0 " << p_0 << std::endl;
				Rcpp::Rcout << "p_1 " << p_1 << std::endl;
				if (p_0 > p_1) {
					if (!bUpdate) {
						gamma[k][j] = false;
						Rcpp::Rcout << "Keep downdate" << std::endl;
						sigma2[k] = sigma2_0;
						mXTX_inv[k] = mXTX_inv_prime;
					}
				}	else {
					if (bUpdate) {
						gamma[k][j] = true;
						Rcpp::Rcout << "Keep update" << std::endl;
						sigma2[k] = sigma2_1;
						mXTX_inv[k] = mXTX_inv_prime;
					}
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