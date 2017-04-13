#include <Rcpp.h>
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_map>
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/functional/hash.hpp>
#include "graycode.h"
#include "correlation.h"

#define DEBUG

using namespace std;
using namespace Rcpp;

namespace boost {
template <typename B, typename A>
std::size_t hash_value(const boost::dynamic_bitset<B, A>& bs) {
  return boost::hash_value(bs.m_bits);
}
}

double log_prob(const int n, const int p, int p_gamma, const double sigma2, const double a, const double b)
{
	double log_sigma2 = std::log(sigma2);
	double log_n = std::log(n);
	#ifdef DEBUG
	// Rcpp::Rcout << "log_sigma2 " << log_sigma2 << std::endl;
	// Rcpp::Rcout << "log_n " << log_n << std::endl;
	// Rcpp::Rcout << "Rf_lbeta(" << a + p_gamma << ", " << b + p - p_gamma << ") " << Rf_lbeta(a + p_gamma, b + p - p_gamma) << std::endl;
	#endif
	double log_p = -n / 2. * log_sigma2 - p / 2. * log_n + Rf_lbeta(a + p_gamma, b + p - p_gamma);
	#ifdef DEBUG
	// Rcpp::Rcout << "log_p " << log_p << std::endl;
	#endif

	return log_p;
}


void gamma_to_NumericMatrix(const vector< dbitset >& gamma, NumericMatrix& nm)
{
	auto K = gamma.size();
	auto p = gamma[0].size();
	for (auto k = 0; k < K; k++) {
		for (auto j = 0; j < p; j++) {
			nm(k, j) = gamma[k][j] ? 1. : 0.;
		}
	}
}


//' Run a Collapsed Variational Approximation to find the K best linear models
//'
//' @param gamma_initial Matrix of initial models, a K by p logical matrix
//' @param vy Vector of responses
//' @param mX Matrix of covariates
//' @param K The number of particles in the population
//' @param lambda The weighting factor for the entropy in f_lambda. Defaults to 1.
//' @return A list containing the named element models, which is a K by p matrix of the models
//'					selected by the algorithm, and the named element trajectory, which includes a list
//'					of the populations of models for each iteration of the algorithm until it converged
//' @export
// [[Rcpp::export]]
List cva(NumericMatrix gamma_initial, NumericVector vy_in, NumericMatrix mX_in, const int K, const double lambda = 1.)
{
	VectorXd vy(vy_in.length()); // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
	for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
	MatrixXd mX(mX_in.nrow(), mX_in.ncol()); // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
	for (auto i = 0; i < mX_in.nrow(); i++)
		for (auto j = 0; j < mX_in.ncol(); j++)
			mX(i, j) = mX_in(i, j);
	const auto n = mX.rows();
	const auto p = mX.cols();
	auto a = 1.;
	auto b = p;
	const MatrixXd mXTX = mX.transpose() * mX;
	const MatrixXd mXTy = mX.transpose() * vy;
	VectorXd probs(K);
	VectorXd w(K);
	vector< dbitset > gamma(K);
	vector< vector< dbitset > > trajectory;
	vector< MatrixXd > mXTX_inv(K);
	VectorXd sigma2(K);
	std::unordered_map< std::size_t, bool > hash;
	// const gsl_rng_type *T;
	// gsl_rng *r;
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

	// gsl_rng_env_setup();

	// T = gsl_rng_default;
	// r = gsl_rng_alloc(T);

	// Initialise population of K particles randomly
	// Rcpp::Rcout << "Generated" << std::endl;
	// for (auto k = 0; k < K; k++) {
	// 	gamma[k].resize(p);
	// 	// Empty models are pointless, so we keep regenerating gammas until we get a non-zero one
	// 	while (gamma[k].count() == 0) {
	// 		for (auto j = 0; j < p; j++) {
	// 			if (gsl_rng_uniform(r) <= RHO)
	// 				gamma[k][j] = true;
	// 			else
	// 				gamma[k][j] = false;
	// 		}
	// 	}
	// 	Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
	// }

	if (gamma_initial.nrow() != K || gamma_initial.ncol() != p) {
		stringstream ss;
		ss << "initial_gamma is " << to_string(gamma_initial.nrow()) << " by " << to_string(gamma_initial.ncol()) << ", expected " << K << " by " << p;
		stop(ss.str());
	}

	Rcpp::Rcout << "initial_gamma" << std::endl;
	for (auto k = 0; k < K; k++) {
		gamma[k].resize(p);
		for (auto j = 0; j < p; j++) {
			if (gamma_initial(k, j) >= EPSILON) {
				gamma[k][j] = true;
			} else {
				gamma[k][j] = false;
			}
		}
		Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
		hash.insert({boost::hash_value(gamma[k]), true});
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
	auto iteration = 1;
	while (!converged) {
		Rcpp::Rcout << "Iteration " << iteration << std::endl;
		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
			#endif
			for (auto j = 0; j < p; j++) {
				dbitset gamma_0 = gamma[k];
				gamma_0[j] = false;
				auto p_gamma_0 = gamma_0.count();

				dbitset gamma_1 = gamma[k];
				gamma_1[j] = true;
				auto p_gamma_1 = gamma_1.count();

				bool bUpdate = !gamma[k][j];
				MatrixXd mXTX_inv_prime;
				if (bUpdate) {
					// If we've seen this bitstring before, don't do the update
					auto h = boost::hash_value(gamma_1);
					#ifdef DEBUG
					// Rcpp::Rcout << "h " << h << std::endl;
					#endif
					auto search = hash.find(h);
					// Rcpp::Rcout << "search->first " << search->first << std::endl;
					if (search != hash.end()) {
						#ifdef DEBUG
						Rcpp::Rcout << "Skipping" << std::endl;
						#endif
						continue;
					} else {
						#ifdef DEBUG
						Rcpp::Rcout << "Inserting" << std::endl;
						#endif
						hash.insert({h, true});
					}
					#ifdef DEBUG
					Rcpp::Rcout << "Updating " << j << std::endl;
					// Rcpp::Rcout << "p_gamma_1 " << p_gamma_1 << std::endl;
					#endif
					// Update mXTX_inv
					mXTX_inv_prime.resize(p_gamma_1, p_gamma_1);
					bool bLow; // gamma_1 or gamma[k]?
					rank_one_update(gamma[k], j, gamma_1.find_first(),
													0,
													mXTX,
													mXTX_inv[k],
													mXTX_inv_prime,
													bLow);
					#ifdef DEBUG
					// Rcpp::Rcout << "mXTX_inv[k]\n " << mXTX_inv[k] << std::endl;
					// Rcpp::Rcout << "mXTX_inv_prime\n " << mXTX_inv_prime << std::endl;
					Rcpp::Rcout << "bLow " << bLow << std::endl;
					#endif
				} else {
					// If only one bit is set and downdating would lead to the null model, then we should skip
					// checking the downdate
					if (p_gamma_0 == 0)
						continue;
					// If we've seen this bitstring before, don't do the update
					auto h = boost::hash_value(gamma_0);
					// Rcpp::Rcout << "h " << h << std::endl;
					auto search = hash.find(h);
					// Rcpp::Rcout << "search->first " << search->first << std::endl;
					if (search != hash.end()) {
						Rcpp::Rcout << "Skipping" << std::endl;
						continue;
					} else {
						Rcpp::Rcout << "Inserting" << std::endl;
						hash.insert({h, true});
					}
					#ifdef DEBUG
					Rcpp::Rcout << "Downdating " << j << std::endl;
					// Rcpp::Rcout << "p_gamma_0 " << p_gamma_0 << std::endl;
					#endif
					// Downdate mXTX_inv
					mXTX_inv_prime.resize(p_gamma_0, p_gamma_0);
					rank_one_downdate(j, gamma_0.find_first(), 0,
														mXTX_inv[k], mXTX_inv_prime);
					#ifdef DEBUG
					// Rcpp::Rcout << "mXTX_inv[k]\n " << mXTX_inv[k] << std::endl;
					// Rcpp::Rcout << "mXTX_inv_prime\n " << mXTX_inv_prime << std::endl;
					#endif
				}

				double sigma2_0;
				double sigma2_1;

 				if (bUpdate){
					MatrixXd mX_gamma_1(n, p_gamma_1);
					get_cols(mX, gamma_1, mX_gamma_1);
					MatrixXd mX_gamma_1_Ty = mX_gamma_1.transpose() * vy;
 					sigma2_0 = sigma2[k];
					#ifdef DEBUG
					// Rcpp::Rcout << "mX_gamma_1.transpose() * mX_gamma_1\n" << mX_gamma_1.transpose() * mX_gamma_1 << std::endl;
					// Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_1.transpose() * mX_gamma_1\n" << mXTX_inv_prime * mX_gamma_1.transpose() * mX_gamma_1 << std::endl;
 				// 	Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_1_Ty\n" << mXTX_inv_prime * mX_gamma_1_Ty << std::endl;
 				// 	Rcpp::Rcout << "(mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value() " << (mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value() << std::endl;
					#endif
 					double nR2 = (mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value();
 					if (nR2 <= n) {
						sigma2_1 = 1. - nR2 / n;
					} else {
						// It's mathematically impossible that nR2 can be that high. Therefore, our inverse must be bad.
						// Recalculate it from scratch and try again.
						#ifdef DEBUG
						// Rcpp::Rcout << "Re-calculating mXTX_inv_prime" << std::endl;
						#endif
						mXTX_inv_prime = (mX_gamma_1.transpose() * mX_gamma_1).inverse();
						#ifdef DEBUG
						// Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_1.transpose() * mX_gamma_1\n" << mXTX_inv_prime * mX_gamma_1.transpose() * mX_gamma_1 << std::endl;
						#endif
						nR2 = (mX_gamma_1_Ty.transpose() * mXTX_inv_prime * mX_gamma_1_Ty).value();
						sigma2_1 = 1. - nR2 / n;
					}
 				} else {
					MatrixXd mX_gamma_0(n, p_gamma_0);
					get_cols(mX, gamma_0, mX_gamma_0);
					MatrixXd mX_gamma_0_Ty = mX_gamma_0.transpose() * vy;
					#ifdef DEBUG
					// Rcpp::Rcout << "mX_gamma_0.transpose() * mX_gamma_0\n" << mX_gamma_0.transpose() * mX_gamma_0 << std::endl;
					// Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_0.transpose() * mX_gamma_0\n" << mXTX_inv_prime * mX_gamma_0.transpose() * mX_gamma_0 << std::endl;
 				// 	Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_0_Ty\n" << mXTX_inv_prime * mX_gamma_0_Ty << std::endl;
 				// 	Rcpp::Rcout << "(mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value() " << (mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value() << std::endl;
					#endif
 					double nR2 = (mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value();
 					sigma2_0 = 1. - (mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value() / n;
 					if (nR2 <= n) {
						sigma2_0 = 1. - nR2 / n;
					} else {
						// It's mathematically impossible that nR2 can be that high. Therefore, our inverse must be bad.
						// Recalculate it from scratch and try again.
						#ifdef DEBUG
						// Rcpp::Rcout << "Re-calculating mXTX_inv_prime" << std::endl;
						#endif
						mXTX_inv_prime = (mX_gamma_0.transpose() * mX_gamma_0).inverse();
						#ifdef DEBUG
						// Rcpp::Rcout << "mXTX_inv_prime * mX_gamma_0.transpose() * mX_gamma_0\n" << mXTX_inv_prime * mX_gamma_0.transpose() * mX_gamma_0 << std::endl;
						#endif
						nR2 = (mX_gamma_0_Ty.transpose() * mXTX_inv_prime * mX_gamma_0_Ty).value();
						sigma2_0 = 1. - nR2 / n;
					}
 					sigma2_1 = sigma2[k];
 				}

				#ifdef DEBUG
				// Rcpp::Rcout << "sigma2_0 " << sigma2_0 << std::endl;
				// Rcpp::Rcout << "sigma2_1 " << sigma2_1 << std::endl;
				#endif
				double log_p_0 = log_prob(n, p, p_gamma_0, sigma2_0, a, b);
				double log_p_1 = log_prob(n, p, p_gamma_1, sigma2_1, a, b);
				#ifdef DEBUG
				Rcpp::Rcout << "log_p_0 " << log_p_0;
				Rcpp::Rcout << " log_p_1 " << log_p_1 << std::endl;
				#endif
				if (log_p_0 > log_p_1) {
					if (!bUpdate) {
						hash.erase(boost::hash_value(gamma[k]));
						gamma[k][j] = false;
						#ifdef DEBUG
						Rcpp::Rcout << "Keep downdate" << std::endl;
						#endif
						sigma2[k] = sigma2_0;
						mXTX_inv[k] = mXTX_inv_prime;
					}
				}	else {
					if (bUpdate) {
						hash.erase(boost::hash_value(gamma[k]));
						gamma[k][j] = true;
						#ifdef DEBUG
						Rcpp::Rcout << "Keep update" << std::endl;
						#endif
						sigma2[k] = sigma2_1;
						mXTX_inv[k] = mXTX_inv_prime;
					}
				}
			}
		}

		for (auto k = 0; k < K; k++) {
			auto p_gamma = gamma[k].count();
			auto b = p;
			probs[k] = exp(log_prob(n, p, p_gamma, sigma2[k], a, b));
			#ifdef DEBUG
			// Rcpp::Rcout << "probs[" << k << "] " << probs[k] << std::endl;
			#endif
		}

		for (auto k = 0; k < K; k++) {
			w[k] = probs[k] / probs.sum();
			#ifdef DEBUG
			// Rcpp::Rcout << "w[" << k << "] " << w[k] << std::endl;
			#endif
		}

		// Check for convergence - is f_lambda changed from the last iteration?
		auto H = -(w.array() * w.array().log()).sum();
		Rcpp::Rcout << "H " << H << std::endl;
		double f_lambda = w.dot(probs) + lambda * H;
		Rcpp::Rcout << "f_lambda_prev " << f_lambda_prev << " f_lambda " << f_lambda << std::endl;
		if ((f_lambda - f_lambda_prev) > EPSILON) {
			f_lambda_prev = f_lambda;
		} else {
			converged = true;
		}
		for (auto k = 0; k < K; k++) {
			Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
		}
		iteration++;
		trajectory.push_back(gamma);
	}
	Rcpp::Rcout << "Converged" << std::endl;
	NumericMatrix bitstrings(K, p);
	gamma_to_NumericMatrix(gamma, bitstrings);
	
	List trajectory_bitstrings;
	for (auto i = 0; i < trajectory.size(); i++) {
		NumericMatrix bitstrings2(K, p);
		gamma_to_NumericMatrix(trajectory[i], bitstrings2);
		trajectory_bitstrings.push_back(bitstrings2);
	}

	List result = List::create(Named("models") = bitstrings,
															Named("trajectory") = trajectory_bitstrings);
	return result;
}
