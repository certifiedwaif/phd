//============================================================================
// Name        : armadillo_demo.cpp
// Author      : Mark Greenaway
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <armadillo>
#include <math.h>
#include <gsl/gsl_sf.h>
#include <gsl/gsl_statistics.h>
#include "faithful.hpp"

double digamma(const double x)
{
	return gsl_sf_psi(x);
}

double square(const double x)
{
	return x*x;
}

#define n 272
#define K 2

using namespace std;
using namespace arma;

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
	}

	void cycle()
	{
		int i, k;

		for (i = 0; i<n; i++) {
			for (k=0; k < K; k++) {
				v(i, k) = digamma(alpha(k)) + .5*digamma(A(k)) - .5*log(B(k)) - .5*A(k)*(square(x(i) - mu(k)) + sigma2(k))/B(k);
			}
		}
		cout << v << endl;

		for (i = 0; i<n; i++) {
			for (k=0; k < K; k++) {
				// TODO: There's the potential for overflow here
				// Remove the maximum using the trick that John showed you.
				w(i, k) = exp(v(i,k)) / sum(exp(v.row(i)));
			}
		}
		cout << w << endl;

		for (k=0;k<K;k++) {
			sigma2(k) = 1.0/(1.0/sigma2_mu(k) + A(k)*sum(w.col(k))/B(k));
			mu(k) = sigma2(k) * (mu_mu(k)/sigma2_mu(k) + A(k)*sum(w.col(k).t() * x)/B(k));
			alpha(k) = prior_alpha + sum(w.col(k));
			A(k) = prior_A(k) + .5*sum(w.col(k));

			B(k) = prior_B(k) + .5*sum(w.col(k).t()*square(x - mu(k)) + sigma2(k));

			cout << mu(k) << " " << sigma2(k) << " " << alpha(k) << " " << A(k) << " " << B(k) << endl;
		}
	}

	const double log_likelihood(void) {
		return 0.0;
	}

};

int main(int argc, char **argv) {
	FiniteGaussianMixture m;

	cout << "Variational approximation to finite Gaussian mixture" << endl;
	m.initialise();

	// Cycle until convergence
	for (int z=0; z<100; z++) {
		cout << "Iteration: " << z << endl;
		m.cycle();
	}

	return 0;
}
