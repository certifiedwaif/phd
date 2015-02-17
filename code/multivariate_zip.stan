data {
  int<lower=1> N; // Number of observations
  int<lower=1> P; // Number of fixed effects covariates
  int<lower=1> B; // Block size
  int<lower=1> M; // Number of subjects
  int<lower=0> y[N]; // Estimated treatment effects
  vector[P] X[N]; // Fixed effects covariate matrix
  vector[B*M] Z[N]; // Random effects covariate matrix
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
  cov_matrix[B] S;
  cov_matrix[P] BetaPrior;
  real<lower=0, upper=1> rho;
  cov_matrix[B] sigma_u;
  //real<lower=0> v;
}

model {
  real eta;
  vector[B*M] u;
  matrix[P,P] chol_BetaPrior;

  rho ~ beta(1.0, 1.0);
  // sigma_u ~ inv_wishart(v, S);
  sigma_u ~ inv_wishart(B+1, S);


  chol_BetaPrior <- cholesky_decompose(BetaPrior);
  vbeta ~ multi_normal_cholesky(zeros_beta, chol_BetaPrior);
  
  vu ~ multi_normal(zeros_u, sigma_u);
  
  for (n in 1:N) {
    for (m in 1:M)
      for (b in 1:B)
        u[(m-1)*B+b] <- vu[m][b];

    eta <- dot_product(X[n], vbeta) + dot_product(Z[n], u);
    //for (m in 1:M)
    //  for (b in 1:B)
    //    eta <- eta + Z[n,(m-1)*B+b]*u[m][b];

    if (y[n] == 0)
      increment_log_prob(log_sum_exp(bernoulli_log(0, rho),
        bernoulli_log(1, rho) + poisson_log_log(y[n], eta)));
    else
      increment_log_prob(bernoulli_log(1, rho) + poisson_log_log(y[n], eta));
  };
}