# Test variational approximation. Unit tests? ----
# Simulate data
# Variational approximation
# Check that lower bounds are monotonically increasing
# Compare accuracy against MCMC
setwd("~/phd/code")
source("zero_inflated_model.R")
n = 1000
rho = .5
lambda = 1

x = generate_test_data(n, rho, lambda)

a = 0.01
b = 0.01

# Variational approximation ----
start = Sys.time()
#pdf("lower_bound_convergence.pdf")
result_var = zero_infl_var(x, a, b)
#dev.off()
var_runtime = Sys.time() - start
# Variational approximation takes .05 seconds to run 10 iterations. So 5ms per iteration,
# or 200 iterations a second.
var_lambda = result_var$a_lambda / result_var$b_lambda
var_rho = result_var$a_rho / (result_var$a_rho + result_var$b_rho)

main_check_accuracy()
