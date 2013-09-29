//============================================================================
// Name        : armadillo_variational_approximation_to_finite_gaussian_mixture.cpp
// Author      : Mark Greenaway
// Version     :
// Copyright   : Your copyright notice
// Description : Implement finite Gaussian mixture model in C++ using the Armadillo
//               library.
//============================================================================

#include <iostream>
#include <armadillo>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include "faithful.hpp"

using namespace std;
using namespace arma;

double digamma(const double x)
{
	return gsl_sf_psi(x);
}

double square(const double x)
{
	return x*x;
}

vec lgamma(vec x)
{
	vec result(x.n_elem);

	for (unsigned int i = 0; i < x.n_elem; i++) {
		result(i) = lgamma(x(i));
	}

	return result;
}

#define n 272
#define K 2

struct FiniteGaussianMixture {
	vec mu_mu, sigma2_mu, prior_A, prior_B;
	const double prior_alpha;
	vec mu, alpha, sigma2, A, B, x;
	mat v, w;

	FiniteGaussianMixture(): mu_mu(K), sigma2_mu(K), prior_A(K), prior_B(K), prior_alpha(.1),
			mu(K), alpha(K), sigma2(K),
			A(K), B(K), x(n), v(n, K), w(n, K)
	{
	}

	void initialise()
	{
		initialise_priors();
		initialise_data();
		initialise_parameters();
	}

	void initialise_priors()
	{
		// Initialise priors
		mu_mu << 0.0 << 1.0 << endr;
		sigma2_mu << 10000 << 10000 << endr;
		prior_A.fill(1.0);
		prior_B = prior_A;
	}

	void initialise_data()
	{
		// Initialise x from faithful, then standardise.
		for (int i=0;i<n;i++) {
			x(i) = faithful[i].eruptions;
		}
		x = (x - mean(x))/stddev(x);
	}

	void initialise_parameters()
	{
		int i, k;
		// Initialise x from faithful, then standardise.
		for (i=0;i<n;i++) {
			x(i) = faithful[i].eruptions;
		}
		x = (x - mean(x))/stddev(x);

		// Initialise parameters
		for (k=0; k<K; k++) {
			mu(k) = mean(x) + k;
			sigma2(k) = var(x);
			alpha(k) = .1;
			A(k) = .1;
			B(k) = .1;
		}

		// Start with equal weighting of all components. We do this because otherwise the initialise
		// evaluation of the log-likelihood produces NaNs.
		// This seems to break things even worse!
		for (i = 0; i < n; i++) {
			for (k = 0; k < K; k++) {
				w(i, k) = 1.0/K;
			}
		}
	}

	void cycle()
	{
		int i, k;

		for (i = 0; i<n; i++) {
			for (k=0; k < K; k++) {
				v(i, k) = digamma(alpha(k)) + .5*digamma(A(k)) - .5*log(B(k)) - .5*A(k)*(square(x(i) - mu(k)) + sigma2(k))/B(k);
			}
		}

		for (i = 0; i<n; i++) {
			for (k=0; k < K; k++) {
				// TODO: There's the potential for overflow here
				// Remove the maximum using the trick that John showed you.
				w(i, k) = exp(v(i,k)) / sum(exp(v.row(i)));
			}
		}

		for (k=0;k<K;k++) {
			sigma2(k) = 1.0/(1.0/sigma2_mu(k) + A(k)*sum(w.col(k))/B(k));
			mu(k) = sigma2(k) * (mu_mu(k)/sigma2_mu(k) + A(k)*sum(w.col(k).t() * x)/B(k));
			alpha(k) = prior_alpha + sum(w.col(k));
			A(k) = prior_A(k) + .5*sum(w.col(k));

			B(k) = prior_B(k) + .5*sum(w.col(k).t()*square(x - mu(k)) + sigma2(k));
		}
	}

	const double log_likelihood(void) {
		double log_lik;

		log_lik = .5*K*(1 - n*log(2* M_PI)) + lgamma(K*prior_alpha);
		log_lik += -K * lgamma(prior_alpha) - lgamma(n + K * prior_alpha);

		vec log_lik_v(K);

		for (int k = 0; k < K; k++) {
			log_lik_v(k) = prior_A(k) * log(prior_B(k)) - A(k) * log(B(k));
			log_lik_v(k) += lgamma(alpha(k)) + .5 * log(sigma2(k)/sigma2_mu(k));
			log_lik_v(k) += -.5*(square(mu(k) - mu_mu(k)) + sigma2(k))/sigma2_mu(k);
			log_lik_v(k) += - sum(w.col(k).t() * log(w.col(k)));
		}

		log_lik += sum(log_lik_v);

		return log_lik;
	}

};

int main(int argc, char **argv) {
	FiniteGaussianMixture m;

	cout << "Variational approximation to finite Gaussian mixture" << endl;
	m.initialise();

	// Cycle until convergence
	int iter = 1;
	double log_likelihood = m.log_likelihood(), last_log_likelihood = -INFINITY;
	while (log_likelihood > last_log_likelihood) {
		last_log_likelihood = log_likelihood;
		cout << "Iteration: " << iter << " , log. lik. " << log_likelihood << endl;
		m.cycle();
		log_likelihood = m.log_likelihood();
		iter++;
	}

	return 0;
}
