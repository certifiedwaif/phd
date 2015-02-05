data {
  int<lower=1> N; // Number of observations
  int<lower=1> P; // Number of fixed effects covariates
  int<lower=1> B; // Block size
  int<lower=1> M; // Number of subjects
  int<lower=0> y[N]; // Estimated treatment effects
  vector[P] X[N]; // Fixed effects covariate matrix
  vector[M*B] Z[N]; // Random effects covariate matrix
}

parameters {
  vector[P] vbeta; 
  vector[B] u;
  cov_matrix[B] S;
  real<lower=0, upper=1> rho;
  real<lower=0> sigma_u;
}

model {
  real eta;

  rho ~ beta(.01, .01);
  sigma_u ~ inv_wishart(B+1, S);

  vbeta ~ normal(0, 100.0);
  
  u ~ normal(0, sigma_u);
  
  for (n in 1:N) {
    eta <- dot_product(X[n], vbeta) + dot_product(Z[n], u);
    if (y[n] == 0)
      increment_log_prob(log_sum_exp(bernoulli_log(0, rho),
        bernoulli_log(1, rho) + poisson_log(y[n], exp(eta))));
    else
      increment_log_prob(bernoulli_log(1, rho) + poisson_log(y[n], exp(eta)));
  };
}