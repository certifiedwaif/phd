data {
  int<lower=1> N; // Number of observations
  # int<lower=1> P; // Number of fixed effects covariates
  int<lower=1> spline_dim; // Number of dimensions for spline random effects
  int<lower=0> y[N]; // Estimated treatment effects
  # vector[P] X[N]; // Fixed effects covariate matrix
  vector[spline_dim] Z[N]; // Random effects covariate matrix
}

transformed data {
  # vector[P] zeros_beta;
  vector[spline_dim] zeros_u;

  # zeros_beta <- rep_vector(0, P);
  zeros_u <- rep_vector(0, spline_dim);
}

parameters {
  # vector[P] vbeta; 
  vector[spline_dim] vu;
  cov_matrix[spline_dim] psi;
  # cov_matrix[P] BetaPrior;
  real<lower=0, upper=1> rho;
  cov_matrix[spline_dim] sigma_u;
}

model {
  real eta;
  # matrix[P, P] chol_BetaPrior;
  matrix[spline_dim, spline_dim] chol_sigma_u;

  rho ~ beta(1.0, 1.0);
  sigma_u ~ inv_wishart(spline_dim + 1, psi);

  # chol_BetaPrior <- cholesky_decompose(BetaPrior);
  # vbeta ~ multi_normal_cholesky(zeros_beta, chol_BetaPrior);
  
  chol_sigma_u <- cholesky_decompose(sigma_u);
  vu ~ multi_normal_cholesky(zeros_u, chol_sigma_u);
  
  for (n in 1:N) {
    # eta <- dot_product(X[n], vbeta) + dot_product(Z[n], vu);
    eta <- dot_product(Z[n], vu);

    if (y[n] == 0)
      increment_log_prob(log_sum_exp(bernoulli_log(0, rho), bernoulli_log(1, rho) + poisson_log_log(y[n], eta)));
    else
      increment_log_prob(bernoulli_log(1, rho) + poisson_log_log(y[n], eta));
  };
}
