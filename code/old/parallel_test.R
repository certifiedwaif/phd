# parallel_test.R
require(multicore)

sum_rnorm = function (x)
{
  sum_r = sum(rnorm(1e3))
  return(sum_r)
}

rand_norm = lapply(1:1e3, sum_rnorm)
rand_norm = simplify2array(rand_norm)
mean(rand_norm)
