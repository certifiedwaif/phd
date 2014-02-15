# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC
source("zero_inflated_model.R")
n = 100
rho = .5
lambda = 100

x = generate_test_data(n, rho, lambda)

a = 1000
b = 10

# Variational approximation ----
start = Sys.time()
result_var = zero_infl_var(x, a, b)
var_runtime = Sys.time() - start
# Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
# or 200 iterations a second.
var_lambda = result_var$a_lambda / result_var$b_lambda
var_rho = result_var$a_rho / (result_var$a_rho + result_var$b_lambda)
