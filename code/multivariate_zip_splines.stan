data {
  int<lower=1> N; // Number of observations
  int<lower=1> P; // Number of fixed effects covariates
  int<lower=0> B; // Block size
  int<lower=1> M; // Number of subjects
  int<lower=0> y[N]; // Estimated treatment effects
  vector[P] X[N]; // Fixed effects covariate matrix
  vector[B*(M-1)] Z[N]; // Random effects covariate matrix
}

transformed data {
  vector[P] zeros_beta;
  vector[B] zeros_u;

  zeros_beta <- rep_vector(0, P);
  zeros_u <- rep_vector(0, B);
}

parameters {
  vector[P] vbeta; 
  vector[B] vu[M];
  cov_matrix[B] psi;
  cov_matrix[P] BetaPrior;
  real<lower=0, upper=1> rho;
  cov_matrix[B] sigma_u;
}

model {
  real eta;
  matrix[P,P] chol_BetaPrior;
  vector[B*(M-1)] u;

  rho ~ beta(1.0, 1.0);
  sigma_u ~ inv_wishart(B + 1, psi);

  chol_BetaPrior <- cholesky_decompose(BetaPrior);
  vbeta ~ multi_normal_cholesky(zeros_beta, chol_BetaPrior);
  
  // This definitely works, but it's slow.
  for (m in 1:M) {
    vu[m] ~ multi_normal(zeros_u, sigma_u);
  }
  // This is faster. But does it work? Seems to, but distributions are weird.
  //vu ~ multi_normal(zeros_u, sigma_u);
  
  for (n in 1:N) {
    // to_vector?
    // i.e. u <- to_vector(vu);
    for (m in 1:M) {
      for (b in 1:B) {
        u[(m-1)*B+b] <- vu[m][b];
      }
    }

    eta <- dot_product(X[n], vbeta) + dot_product(Z[n], u);

    if (y[n] == 0)
      increment_log_prob(log_sum_exp(bernoulli_log(0, rho),
        bernoulli_log(1, rho) + poisson_log_log(y[n], eta)));
    else
      increment_log_prob(bernoulli_log(1, rho) + poisson_log_log(y[n], eta));
  };
}
