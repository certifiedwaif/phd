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

void cva(VectorXd& vy, MatrixXd& mX, const uint K, const uint p)
{
	const double a = 1.;
	const double b = 1.;
	const uint n = mX.rows();
	const MatrixXd mXTy = mX.transpose() * vy;
	VectorXd w(K);
	vector< dbitset > gamma(K);
	vector< MatrixXd > mXTX_inv;
	VectorXd sigma2(K);
	const gsl_rng_type *T;
	gsl_rng *r;
	const auto RHO = 0.1;

	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);

	// Initialise population of K particles randomly
	for (auto k = 0; k < K; k++) {
		gamma[k].resize(p);
		for (auto j = 0; j < p; j++) {
			if (gsl_rng_uniform(r) <= RHO)
				gamma[k][j] = true;
			else
				gamma[k][j] = false;
		}
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
				bool bUpdate = gamma[k][j];
				if (bUpdate) {
					p_gamma = p_gamma_1;
				} else {
					p_gamma = p_gamma_0;
				}

				MatrixXd mX_gamma_0TX_gamma_0_inv(p_gamma_0, p_gamma_0);
				MatrixXd mX_gamma_1TX_gamma_1_inv(p_gamma_1, p_gamma_1);

				MatrixXd mXgammaTXgamma_inv(p_gamma, p_gamma);
				MatrixXd mXgamma_Ty = get_rows(mXTy, gamma_1, mXgamma_Ty);
				double sigma2_0 = 1. - (mX_gamma_0_Ty.transpose() * mX_gamma_0TX_gamma_0_inv * mX_gamma_0_Ty).value() / n;
				double sigma2_1 = 1. - (mX_gamma_1_Ty.transpose() * mX_gamma_1TX_gamma_1_inv * mX_gamma_1_Ty).value() / n;

				double p_0 = prob(n, p, sigma2_0, a, b);
				double p_1 = prob(n, p, sigma2_1, a, b);
				if (p_0 > p_1) {
					gamma[k][j] = true;
					sigma2[k] = sigma2_0;
				}	else {
					gamma[k][j] = false;
					sigma2[k] = sigma2_1;
				}
			}
		}

		double sum_p = 0.;
		for (auto k = 0; k < K; k++) sum_p += prob(n, p, sigma2[k], a, b);
		for (auto k = 0; k < K; k++) {
			w[k] = prob(n, p, sigma2[k], a, b) / sum_p;
		}

		// Check for convergence
	}
}