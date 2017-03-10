
double p(const uint n, const uint p, const double sigma2, const double a, const double b)
{
	double log_p = -n / 2. * log(sigma2) - p / 2. * log(n) + gsl_sf_lnbeta(a, b);

	return exp(log_p);
}

void cva(VectorXd vy, const uint K, const uint p)
{
	dbitset gamma(p);
	VectorXd w(K);
	// Loop until convergence
	for (auto k = 0; k < K; k++) {
		for (auto j = 0; j < p; j++) {
		}
	}

	double sum_p = 0.0;
	for (k = 0; k < K; k++) sum_p += p(vy, gamma[k];)
	for (auto k = 0; k < K; k++) {
		w[k] = p(vy, gamma[k]) / sum_p;
	}
}