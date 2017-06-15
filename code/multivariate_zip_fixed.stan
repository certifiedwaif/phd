data {
  int<lower=1> N; // Number of observations
  int<lower=1> P; // Number of fixed effects covariates
  int<lower=0> y[N]; // Estimated treatment effects
  vector[P] X[N]; // Fixed effects covariate matrix
}

transformed data {
  vector[P] zeros_beta;

  zeros_beta = rep_vector(0, P);
}

parameters {
  vector[P] vbeta; 
  cov_matrix[P] BetaPrior;
  real<lower=0, upper=1> rho;
}

model {
  real eta;
  matrix[P, P] chol_BetaPrior;

  rho ~ beta(1.0, 1.0);

  chol_BetaPrior = cholesky_decompose(BetaPrior);
  vbeta ~ multi_normal_cholesky(zeros_beta, chol_BetaPrior);
  
  for (n in 1:N) {
    eta = dot_product(X[n], vbeta);

    if (y[n] == 0)
      # target += log_sum_exp(bernoulli_lpmf(0|rho),
      #                       bernoulli_lpmf(1|rho) + poisson_log_lpmf(y[n]|eta));
      target += bernoulli_lpmf(0|rho);
    else
      target += bernoulli_lpmf(1|rho) + poisson_log_lpmf(y[n]|eta);
  };
}
