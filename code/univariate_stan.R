# univariate_stan.R
require(rstan)

zip_code <- '
  data {
    int<lower=1> N; // number of schools 
    int<lower=0> y[N]; // estimated treatment effects
  }
  parameters {
    real lambda; 
    real<lower=0, upper=1> rho;
  }
  model {
    lambda ~ gamma(.01, .01);
    rho ~ beta(.01, .01);

    for (i in 1:N) {
      if (y[i] == 0)
        increment_log_prob(log_sum_exp(bernoulli_log(1, rho),
                                        bernoulli_log(0, rho) + poisson_log(y[i], lambda)));
      else
        increment_log_prob(bernoulli_log(0, rho) + poisson_log(y[i], lambda));
    };
  }
'

zip_dat <- list(N = 50, 
                    y = c(0,7,3,4,5,3,2,6,5,0,0,1,
                          0,0,5,0,2,3,6,4,0,5,4,0,
                          7,0,0,0,7,0,6,6,0,3,0,5,
                          0,4,0,0,0,2,3,0,3,4,5,0,
                          8,0))

fit <- stan(model_code = zip_code, data = zip_dat, 
            iter = 1e6, chains = 4)
plot(fit)
