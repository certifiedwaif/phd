#include <Rcpp.h>
// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <vector>
#include <cmath>
#include <sstream>
#include <unordered_map>
#define BOOST_DYNAMIC_BITSET_DONT_USE_FRIENDS
#include <boost/functional/hash.hpp>

// [[Rcpp::depends(RcppGSL)]]
#include <gsl/gsl_sf_hyperg.h>
#include <gsl/gsl_errno.h>


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


// Which one was this one again? Maruyama?
double maruyama(const int n, const int p, const double R2, int p_gamma)
{
	const auto sigma2 = 1. - R2;
	const auto a = 1.;
	const auto b = p;
	#ifdef DEBUG
	Rcpp::Rcout << "sigma2 " << sigma2 << std::endl;
	Rcpp::Rcout << "n " << n << std::endl;
	#endif
	const auto log_sigma2 = std::log(sigma2);
	const auto log_n = std::log(n);
	#ifdef DEBUG
	Rcpp::Rcout << "log_sigma2 " << log_sigma2 << std::endl;
	Rcpp::Rcout << "p_gamma " << p_gamma << std::endl;
	Rcpp::Rcout << "log_n " << log_n << std::endl;
	Rcpp::Rcout << "Rf_lbeta(" << a + p_gamma << ", " << b + p - p_gamma << ") " << Rf_lbeta(a + p_gamma, b + p - p_gamma) << std::endl;
	#endif
	double log_p;
	log_p = -n / 2. * log_sigma2 - p_gamma / 2. * log_n + Rf_lbeta(a + p_gamma, b + p - p_gamma);
	if (sigma2 == 0. || std::isnan(sigma2)) {
		log_p = -INFINITY;
		// throw std::range_error("Invalid sigma2");
	}
	#ifdef DEBUG
	Rcpp::Rcout << "log_p " << log_p << std::endl;
	#endif
	// if (isnan(log_p)) {
	// // If a floating point number is not equal to itself, then it must be NaN
	// // if (log_p != log_p) {
	// 	log_p = -INFINITY;
	// }
	#ifdef DEBUG
	Rcpp::Rcout << "log_p " << log_p << std::endl;
	#endif

	return log_p;
}


double BIC(const int n, const int p, double R2, int vp_gamma)
{
	const auto sigma2 = 1. - R2;
	#ifdef DEBUG
	Rcpp::Rcout << "n " << n << " p " << p << " R2 " << R2 << " vp_gamma " << vp_gamma << std::endl;
	#endif
	auto BIC = n * log(sigma2) + vp_gamma * log(n);
	if (sigma2 == 0. || std::isnan(sigma2)) {
		BIC = -INFINITY;
		// throw std::range_error("Invalid sigma2");
	}
	#ifdef DEBUG
	Rcpp::Rcout << "n * log(1 - R2) " << n * log(1 - R2) << std::endl;
	Rcpp::Rcout << "vp_gamma * log(n) " << vp_gamma * log(n) << std::endl;
	#endif
	return -0.5 * BIC;
}


double ZE(const int n, const int p, double R2, int p_gamma)
{
	auto a = -0.75;
	auto b = 0.5 * (n - p_gamma - 5) - a;
	auto c = 0.5 * (n - 1);
	auto d = 0.5 * p_gamma + a;

	auto log_p = -(b+1)*log(1 - R2) + Rf_lbeta(d+1,b+1) - Rf_lbeta(a+1,b+1);
	auto ZE = -2*log_p;
	return log_p;
}


double log_hyperg_2F1(double b, double c, double x)
{
	if (x == 0.)
		return 0.;
	auto val = 0.;
	val += log(c-1);
	val += (1-c)*log(x);
	val += (c-b-1)*log(1-x);
	val += Rf_lbeta(c-1,b-c+1);
	val += Rf_pbeta(x, (c-1), (b-c+1), true, true);
	return val;
}


double log_hyperg_2F1_naive(double b, double c, double x)
{
	auto val = log(gsl_sf_hyperg_2F1( b, 1, c, x));
	return val;
}


// Liang's hyper g-prior
double liang_g1(const int n, const int p, double R2, int p_gamma)
{
	auto a = 3.;
	double log_p_g;
	log_p_g = log(a - 2) - log(p_gamma + a - 2) + log(gsl_sf_hyperg_2F1(0.5*(n-1), 1, 0.5*(p_gamma + a), R2));
	return log_p_g;
}


// Liang's g prior
double liang_g2(const int n, const int p, double R2, int p_gamma)
{
	auto a = 3.;
	auto log_vp_g2 = log(a - 2) - log(p_gamma + a - 2) + log_hyperg_2F1( 0.5*(n-1), 0.5*(p_gamma+a), R2);
	return log_vp_g2;
}


// Liang's g prior
double liang_g3(const int n, const int p, double R2, int p_gamma)
{
	double log_vp_gprior5;
	if (p_gamma == 0)
		return 0.;
	log_vp_gprior5 -= 0.5*p_gamma*log(n+1);
	log_vp_gprior5 += 0.5*p_gamma*log(p_gamma+1);
	log_vp_gprior5 -= 0.5*(n - 1)*log(R2);
	log_vp_gprior5 -= log(p_gamma+1);
	// Check for errors
	gsl_sf_result result;
	int error_code = gsl_sf_hyperg_2F1_e( 0.5*(p_gamma+1), 0.5*(n-1), 0.5*(p_gamma+3), (1-1/R2)*(p_gamma+1)/(n+1), &result);
	if (error_code == GSL_EMAXITER) {
		log_vp_gprior5 = -INFINITY;
	}	else {
		log_vp_gprior5 += log(result.val);
	}
	return log_vp_gprior5;
}


// Robust Bayarri
double robust_bayarri1(const int n, const int p, double R2, int p_gamma)
{
	#ifdef DEBUG
	Rcpp::Rcout << "n " << n;
	Rcpp::Rcout << " p " << p;
	Rcpp::Rcout << " R2 " << R2;
	Rcpp::Rcout << " p_gamma " << p_gamma;
	#endif
	auto L = (1. + n)/(1. + p_gamma) - 1.;
	auto sigma2 = 1. - R2;
	auto z = R2/(1. + L*sigma2);

	double log_vp_gprior6;

	if (p_gamma == 0)
		return 0.;
	log_vp_gprior6 += 0.5*(n - p_gamma - 1)*log( n + 1 );
	log_vp_gprior6 -= 0.5*(n - p_gamma - 1)*log( p_gamma + 1);
	log_vp_gprior6 -= 0.5*(n - 1)*log(1 + L*sigma2);
	log_vp_gprior6 -= log (p_gamma + 1);
	log_vp_gprior6 += log(gsl_sf_hyperg_1F1( 0.5*(n-1), 0.5*(p_gamma+3), z ));
	#ifdef DEBUG
	Rcpp::Rcout << " log_vp_gprior6 " << log_vp_gprior6 << std::endl;
	#endif
	return log_vp_gprior6;
}


double robust_bayarri2(const int n, const int p, double R2, int p_gamma)
{
	auto sigma2 = 1. - R2;
	auto L = (1. + n)/(1. + p_gamma) - 1.;
	auto z = R2/(1. + L*sigma2);

	if (p_gamma == 0)
		return 0.;
	double log_vp_gprior7;
	log_vp_gprior7 += 0.5*(n - p_gamma - 1)*log( n + 1 );
	log_vp_gprior7 -= 0.5*(n - p_gamma - 1)*log( p_gamma + 1);
	log_vp_gprior7 -= 0.5*(n - 1)*log(1 + L*sigma2);
	log_vp_gprior7 -= log (p_gamma + 1);
	log_vp_gprior7 += log(gsl_sf_hyperg_2F1( 0.5*(n-1), 1, 0.5*(p_gamma + 3), z));
	return log_vp_gprior7;
}


void calculate_log_probabilities(const vector< dbitset >& gamma, const VectorXd& sigma2, const int n,
																	VectorXd& log_probs, std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob)
{
	auto K = gamma.size();
	auto p = gamma[0].size();
	auto a = 1.;

	for (auto k = 0; k < K; k++) {
		auto p_gamma = gamma[k].count();
		auto b = p;
		log_probs[k] = log_prob(n, p, 1. - sigma2[k], p_gamma);
		#ifdef DEBUG
		// Rcpp::Rcout << "log_probs[" << k << "] " << log_probs[k] << std::endl;
		#endif
	}
	// Re-normalise
	log_probs = log_probs.array() - log_probs.maxCoeff();
}


void calculate_weights(const VectorXd& sigma2, const VectorXd& log_probs, VectorXd& w)
{
	const auto K = log_probs.size();
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


double calculate_entropy(const VectorXd& w)
{
	const auto K = w.size();
	auto H = 0.;
	for (auto k = 0; k < K; k++) {
		// Don't add 0 weights
		if (w[k] == 0.) {
			continue;
		}	else {
			H += -w[k] * log(w[k]);
		}
	}
	return H;
}


double calculate_w_dot_prob(const VectorXd& w, const VectorXd& log_probs)
{
	const auto K = w.size();
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
//' @param log_lik The log likelihood function to use in ranking models. One of
//' 							 "log_prob1" "BIC" "ZE" "3" "4" "5" "6" "7".
//' @param bUnique Whether to ensure uniqueness in the population of particles or not. Defaults to true.
//' @return A list containing the named element models, which is a K by p matrix of the models
//'					selected by the algorithm, and the named element trajectory, which includes a list
//'					of the populations of models for each iteration of the algorithm until it converged
//' @export
// [[Rcpp::export]]
List cva(NumericMatrix gamma_initial, NumericVector vy_in, NumericMatrix mX_in, const int K,
				 const double lambda = 1., std::string log_lik = "maruyama", const bool bUnique = true)
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

	std::function<double (const int n, const int p, double vR2, int vp_gamma)> log_prob;
	if (log_lik == "maruyama") {
		log_prob = maruyama;
	} else if (log_lik == "BIC") {
		log_prob = BIC;
	} else if (log_lik == "ZE") {
		log_prob = ZE;
	} else if (log_lik == "liang_g1") {
		log_prob = liang_g1;
	} else if (log_lik == "liang_g2") {
		log_prob = liang_g2;
	} else if (log_lik == "liang_g3") {
		log_prob = liang_g3;
	} else if (log_lik == "robust_bayarri1") {
		log_prob = robust_bayarri1;
	} else if (log_lik == "robust_bayarri2") {
		log_prob = robust_bayarri2;
	} else {
		log_prob = maruyama;
	}
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
			if (gamma_initial(k, j) == 1.) {
				gamma[k][j] = true;
			}
			else {
				gamma[k][j] = false;
			}
		}
		#ifdef DEBUG
		Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
		#endif
		if (bUnique) {
			hash.insert({boost::hash_value(gamma[k]), true});
		}
	}

	// Initialise mXTX_inv
	// Initialise sigma2
	for (auto k = 0; k < K; k++) {
		auto p_gamma = gamma[k].count();
		if (p_gamma == 0) {
			stringstream ss;
			ss << "gamma[" << k + 1 << "] has no bits set" << endl;
			throw domain_error(ss.str());
		}
		MatrixXd mX_gamma(n, p_gamma);
		get_cols(mX, gamma[k], mX_gamma);
		MatrixXd mX_gamma_Ty(p_gamma, 1);
		get_rows(mXTy, gamma[k], mX_gamma_Ty);
		mXTX_inv[k] = (mX_gamma.transpose() * mX_gamma).inverse();
		sigma2[k] = 1. - (mX_gamma_Ty.transpose() * mXTX_inv[k] * mX_gamma_Ty).value() / n;
		#ifdef DEBUG
		Rcpp::Rcout << "sigma2[" << k << "] " << sigma2[k] << std::endl;
		#endif
	}
	calculate_log_probabilities(gamma, sigma2, n, log_probs, log_prob);
	trajectory.push_back(gamma);
	trajectory_probs.push_back(log_probs.array().exp());

	// Loop until convergence
	auto converged = false;
	auto iteration = 1;
	while (!converged) {
		#ifdef DEBUG
		Rcpp::Rcout << "Iteration " << iteration << std::endl;
		#endif

		// #pragma omp parallel for
		for (auto k = 0; k < K; k++) {
			#ifdef DEBUG
			Rcpp::Rcout << "gamma[" << k << "] " << gamma[k] << std::endl;
			#endif
			for (auto j = 0; j < p; j++) {
				dbitset gamma_prime = gamma[k];
				gamma_prime[j] = !gamma_prime[j];
				auto p_gamma = gamma[k].count();
				auto p_gamma_prime = gamma_prime.count();
				if ((p_gamma_prime == 0) || (p_gamma_prime >= n - 1))
					continue;
				bool bUpdate = !gamma[k][j];

				// If we've seen this bitstring before, don't do the update
				if (bUnique) {
					auto h = boost::hash_value(gamma_prime);
					auto search = hash.find(h);
					if (search != hash.end()) {
						continue;
					}	else {
						hash.insert({h, true});
					}
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

				// Rcpp::Rcout << "gamma[" << k << "]    " << gamma[k] << std::endl;
				// Rcpp::Rcout << "gamma_prime " << gamma_prime << std::endl;
				#ifdef DEBUG
				MatrixXd mXTX_inv_prime_check = (mX_gamma_prime.transpose() * mX_gamma_prime).inverse();
				if (!mXTX_inv_prime.isApprox(mXTX_inv_prime_check, 1e-8)) {
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
					log_p_0 = log_prob(n, p, 1. - sigma2[k], p_gamma);
					log_p_1 = log_prob(n, p, 1. - sigma2_prime, p_gamma_prime);
				} else {
					log_p_0 = log_prob(n, p, 1. - sigma2_prime, p_gamma_prime);
					log_p_1 = log_prob(n, p, 1. - sigma2[k], p_gamma);
				}
				#ifdef DEBUG
				Rcpp::Rcout << "log_p_0 " << log_p_0;
				Rcpp::Rcout << " log_p_1 " << log_p_1 << std::endl;
				#endif
				if ((log_p_0 > log_p_1 && !bUpdate) || (log_p_1 > log_p_0 && bUpdate)) {
					if (bUnique) {
						hash.erase(boost::hash_value(gamma[k]));
					}
					gamma[k][j] = bUpdate;
					if (bUnique) {
						hash.insert({boost::hash_value(gamma[k]), true});
					}
					#ifdef DEBUG
					if (bUpdate)
						Rcpp::Rcout << "Keep update" << std::endl;
					else
						Rcpp::Rcout << "Keep downdate" << std::endl;
					#endif
					sigma2[k] = sigma2_prime;
					mXTX_inv[k] = mXTX_inv_prime;
				} else {
					#ifdef DEBUG
					if (bUpdate)
						Rcpp::Rcout << "Don't keep update" << std::endl;
					else
						Rcpp::Rcout << "Don't keep downdate" << std::endl;
					#endif
				}
			}
		}

		calculate_log_probabilities(gamma, sigma2, n, log_probs, log_prob);

		// Calculate weights
		calculate_weights(sigma2, log_probs, w);

		// Calculate entropy
		auto H = calculate_entropy(w);

		// Rcpp::Rcout << "w.dot(probs) " << w.dot(probs) << std::endl;
		auto w_dot_prob = calculate_w_dot_prob(w, log_probs);
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
