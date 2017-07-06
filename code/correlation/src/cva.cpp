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

// #define DEBUG

using namespace std;
using namespace Rcpp;

namespace boost
{
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
	double log_p;
	if (sigma2 == 0.) {
		log_p = -INFINITY;
	}
	log_p = -n / 2. * log_sigma2 - p_gamma / 2. * log_n + Rf_lbeta(a + p_gamma, b + p - p_gamma);
	#ifdef DEBUG
	// Rcpp::Rcout << "log_p " << log_p << std::endl;
	#endif

	return log_p;
}


void calculate_log_probabilities(const vector< dbitset >& gamma, const VectorXd& sigma2, const int n,
																	VectorXd& log_probs)
{
	auto K = gamma.size();
	auto p = gamma[0].size();
	auto a = 1.;

	for (auto k = 0; k < K; k++) {
		if (sigma2[k] == 0.) {
			log_probs[k] = -INFINITY;
		} else {
			auto p_gamma = gamma[k].count();
			auto b = p;
			log_probs[k] = log_prob(n, p, p_gamma, sigma2[k], a, b);
		}
		#ifdef DEBUG
		// Rcpp::Rcout << "log_probs[" << k << "] " << log_probs[k] << std::endl;
		#endif
	}
	// Re-normalise
	log_probs = log_probs.array() - log_probs.maxCoeff();
}


void calculate_weights(const VectorXd& sigma2, const VectorXd& log_probs, VectorXd& w, const int K)
{
	for (auto k = 0; k < K; k++) {
		// Need to use log-sum-exp trick here
		// Have to skip 0 probabilities somehow.
		if (sigma2[k] == 0.) {
			w[k] = 0.;
		} else {
			w[k] = exp(log_probs[k]) / log_probs.array().exp().sum();
		}
		#ifdef DEBUG
		Rcpp::Rcout << "w[" << k << "] " << w[k] << std::endl;
		#endif
	}
}


double calculate_entropy(const VectorXd& w, const int K)
{
	// Don't add 0 weights
	auto H = 0.;
	for (auto k = 0; k < K; k++) {
		if (w[k] == 0.) {
			continue;
		}	else {
			H += -w[k] * log(w[k]);
		}
	}

	return H;
}


double calculate_w_dot_prob(const VectorXd& w, const VectorXd& log_probs, const int K)
{
	auto w_dot_prob = 0.;
	for (auto k = 0; k < K; k++) {
		if (w[k] == 0.) {
			continue;
		} else {
			w_dot_prob += w[k] * exp(log_probs[k]);
		}
		#ifdef DEBUG
		Rcpp::Rcout << "w[k] " << w[k] << " log_probs[k] " << log_probs[k] << " w_dot_prob " << k << " " << w_dot_prob << " " << std::endl;
		#endif
	}
	return w_dot_prob;
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
	VectorXd vy(vy_in.length());   // = Rcpp::as<Eigen::Map<Eigen::VectorXd>>(vy_in);
	for (auto i = 0; i < vy_in.length(); i++) vy[i] = vy_in[i];
																 // = Rcpp::as<Eigen::Map<Eigen::MatrixXd>>(mX_in);
	MatrixXd mX(mX_in.nrow(), mX_in.ncol());
	for (auto i = 0; i < mX_in.nrow(); i++)
		for (auto j = 0; j < mX_in.ncol(); j++)
			mX(i, j) = mX_in(i, j);
	const auto n = mX.rows();
	const auto p = mX.cols();
	auto a = 1.;
	auto b = p;
	const MatrixXd mXTX = mX.transpose() * mX;
	const MatrixXd mXTy = mX.transpose() * vy;
	VectorXd log_probs(K);
	VectorXd w(K);
	vector< dbitset > gamma(K);
	vector< vector< dbitset > > trajectory;
	vector< VectorXd > trajectory_probs;
	vector< MatrixXd > mXTX_inv(K);
	VectorXd sigma2(K);
	std::unordered_map< std::size_t, bool > hash;
	// const gsl_rng_type *T;
	// gsl_rng *r;
	const auto RHO = 0.1;
	auto f_lambda_prev = 0.;
	const auto EPSILON = 1e-8;

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

	#ifdef DEBUG
	Rcpp::Rcout << "initial_gamma" << std::endl;
	#endif
	for (auto k = 0; k < K; k++) {
		gamma[k].resize(p);
		for (auto j = 0; j < p; j++) {
			if (gamma_initial(k, j) >= EPSILON) {
				gamma[k][j] = true;
			}
			else {
				gamma[k][j] = false;
			}
		}
		#ifdef DEBUG
		Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
		#endif
		hash.insert({boost::hash_value(gamma[k]), true}
		);
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
		// Rcpp::Rcout << "sigma2[" << k << "] " << sigma2[k] << std::endl;
		calculate_log_probabilities(gamma, sigma2, n, log_probs);
	}
	trajectory.push_back(gamma);
	trajectory_probs.push_back(log_probs.array().exp());

	// Loop until convergence
	auto converged = false;
	auto iteration = 1;
	while (!converged) {
		#ifdef DEBUG
		Rcpp::Rcout << "Iteration " << iteration << std::endl;
		#endif

		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
			#endif
			for (auto j = 0; j < p; j++) {
				dbitset gamma_prime = gamma[k];
				gamma_prime[j] = !gamma_prime[j];
				auto p_gamma = gamma[k].count();
				auto p_gamma_prime = gamma_prime.count();
				if (p_gamma_prime == 0)
					continue;
				bool bUpdate = !gamma[k][j];

				// If we've seen this bitstring before, don't do the update
				auto h = boost::hash_value(gamma_prime);
				auto search = hash.find(h);
				if (search != hash.end()) {
					continue;
				}	else {
					hash.insert({h, true});
				}
				#ifdef DEBUG
				if (bUpdate)
					Rcpp::Rcout << "Updating " << j << std::endl;
				else
					Rcpp::Rcout << "Downdating " << j << std::endl;
				#endif

				// Update or downdate mXTX_inv
				MatrixXd mXTX_inv_prime;
				mXTX_inv_prime.resize(p_gamma_prime, p_gamma_prime);
				bool bLow;             // gamma_1 or gamma[k]?
				uint min_idx = std::min(gamma[k].find_first(), gamma_prime.find_first());
				#ifdef DEBUG
				Rcpp::Rcout << "min_idx " << min_idx << std::endl;
				#endif

				// This is totally evil
				// Explanation: mXTX is addressed absolutely, but mXTX_inv is addressed relatively. To account for
				// this, we abuse fixed to adjust for the gaps in gamma_prime
				int fixed = 0;
				for (auto idx = min_idx; idx < j; idx++) {
					if (!gamma_prime[idx])
						fixed--;
				}

				if (bUpdate) {
					rank_one_update(gamma[k], j, min_idx,
						fixed,
						mXTX,
						mXTX_inv[k],
						mXTX_inv_prime,
						bLow);
						#ifdef DEBUG
						Rcpp::Rcout << "bLow " << bLow << std::endl;
						#endif
				} else {
					rank_one_downdate(j, min_idx, fixed,
						mXTX_inv[k], mXTX_inv_prime);
				}

				// Calculate sigma2_prime
				double sigma2_prime;

				MatrixXd mX_gamma_prime(n, p_gamma_prime);
				get_cols(mX, gamma_prime, mX_gamma_prime);
				MatrixXd mX_gamma_prime_Ty = mX_gamma_prime.transpose() * vy;
				double nR2 = (mX_gamma_prime_Ty.transpose() * mXTX_inv_prime * mX_gamma_prime_Ty).value();
				sigma2_prime = 1. - nR2 / n;
				// #ifdef DEBUG
				// It's mathematically impossible that nR2 can be that high. Therefore, our inverse must be bad.
				// Recalculate it from scratch and try again.
				
				#ifdef DEBUG
				MatrixXd mXTX_inv_prime_check = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
				if (mXTX_inv_prime != mXTX_inv_prime_check) {
					Rcpp::Rcout << "gamma[k]    " << gamma[k] << std::endl;
					Rcpp::Rcout << "gamma_prime " << gamma_prime << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime " << std::endl << mXTX_inv_prime << std::endl;
					Rcpp::Rcout << "mXTX_inv_prime_check " << std::endl << mXTX_inv_prime_check << std::endl;
					Rcpp::Rcout << "Difference " << std::endl << mXTX_inv_prime - mXTX_inv_prime_check << std::endl;
					// throw std::range_error("Rank one update failed. I wonder why?");
				}

				double nR2_check = (mX_gamma_prime_Ty.transpose() * mXTX_inv_prime * mX_gamma_prime_Ty).value();
				double sigma2_prime_check = 1. - nR2_check / n;
				if (abs(sigma2_prime - sigma2_prime_check) > EPSILON) {
					Rcpp::Rcout << "sigma2_prime " << sigma2_prime << " sigma2_prime_check " << sigma2_prime_check << std::endl;
				}
				mXTX_inv_prime = mXTX_inv_prime_check;
				#endif
				// Rcpp::Rcout << "sigma2_1 " << sigma2_1 << std::endl;
				// #endif

				#ifdef DEBUG
				Rcpp::Rcout << "sigma2[k] " << sigma2[k] << std::endl;
				Rcpp::Rcout << "sigma2_prime " << sigma2_prime << std::endl;
				#endif
				double log_p_0;
				double log_p_1;
				if (bUpdate) {
					log_p_0 = log_prob(n, p, p_gamma, sigma2[k], a, b);
					log_p_1 = log_prob(n, p, p_gamma_prime, sigma2_prime, a, b);
				} else {
					log_p_0 = log_prob(n, p, p_gamma_prime, sigma2_prime, a, b);
					log_p_1 = log_prob(n, p, p_gamma, sigma2[k], a, b);
				}
				#ifdef DEBUG
				Rcpp::Rcout << "log_p_0 " << log_p_0;
				Rcpp::Rcout << " log_p_1 " << log_p_1 << std::endl;
				#endif
				if ((log_p_0 > log_p_1 && !bUpdate) || (log_p_1 > log_p_0 && bUpdate)) {
					hash.erase(boost::hash_value(gamma[k]));
					gamma[k][j] = bUpdate;
					hash.insert({boost::hash_value(gamma[k]), true});
					#ifdef DEBUG
					if (bUpdate)
						Rcpp::Rcout << "Keep update" << std::endl;
					else
						Rcpp::Rcout << "Keep downdate" << std::endl;
					#endif
					sigma2[k] = sigma2_prime;
					mXTX_inv[k] = mXTX_inv_prime;
				}
			}
		}

		calculate_log_probabilities(gamma, sigma2, n, log_probs);

		// Calculate weights
		calculate_weights(sigma2, log_probs, w, K);

		// Calculate entropy
		auto H = calculate_entropy(w, K);

		// Rcpp::Rcout << "w.dot(probs) " << w.dot(probs) << std::endl;
		auto w_dot_prob = calculate_w_dot_prob(w, log_probs, K);
		#ifdef  DEBUG
		Rcpp::Rcout << "w_dot_prob " << w_dot_prob << std::endl;
		Rcpp::Rcout << "H " << H << std::endl;
		#endif

		// Calculate f_lambda
		double f_lambda = w_dot_prob + lambda * H;
		#ifdef DEBUG
		Rcpp::Rcout << "f_lambda_prev " << f_lambda_prev << " f_lambda " << f_lambda << std::endl;
		#endif

		// Check for convergence - is f_lambda changed from the last iteration?
		if ((f_lambda - f_lambda_prev) > EPSILON) {
			f_lambda_prev = f_lambda;
		}
		else {
			converged = true;
		}
		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k + 1 << "] " << gamma[k] << std::endl;
			#endif
		}
		iteration++;
		trajectory.push_back(gamma);
		trajectory_probs.push_back(log_probs.array().exp());
	}
	#ifdef DEBUG
	Rcpp::Rcout << "Converged" << std::endl;
	#endif
	NumericMatrix bitstrings(K, p);
	gamma_to_NumericMatrix(gamma, bitstrings);

	List trajectory_bitstrings;
	NumericMatrix trajectory_probabilities(K, trajectory.size());
	for (auto i = 0; i < trajectory.size(); i++) {
		NumericMatrix bitstrings2(K, p);
		gamma_to_NumericMatrix(trajectory[i], bitstrings2);
		trajectory_bitstrings.push_back(bitstrings2);
		for (auto k = 0; k < K; k++) {
			trajectory_probabilities(k, i) = trajectory_probs[i](k);
		}
	}

	List result = List::create(Named("models") = bitstrings,
															Named("trajectory") = trajectory_bitstrings,
															Named("trajectory_probs") = trajectory_probabilities);
	return result;
}
