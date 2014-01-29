# Statistical inference Assignment 2 Question 1 ----
# Simulate the data, assuming Pi is known ----
set.seed(2014)
pi = .5
true_lambda = 5
n = 1000
x = rep(NA, n)
for (i in 1:n) {
  x[i] = rexp(1, true_lambda)
  if (rbinom(1, 1, .5) == 1) {
    x[i] = x[i] + 1
  }
}
hist(x)

# EM algorithm ----
# Initialise ----
iterations = 10
lambda = rep(NA, iterations)
lambda[1] = 1/mean(x)
k = 1
# E step ----
e = rep(NA, n)
for (i in 1:n) {
  num = pi * dexp(x[i], lambda[k])
  if ((x[i] - 1) <= 0) {
    denom = pi * dexp(x[i], lambda[k])
  } else {
    denom = pi * dexp(x[i], lambda[k]) + (1 - pi) * dexp(x[i] - 1, lambda[k])
  }
  e[i] = num/denom
  print(sprintf("i %f x[i] %f num %f denom %f e[i] %f", i, x[i], num, denom, e[i]))
}
# M step ----
sum(e*x)
sum((1-e)*(x-1))
lambda[k + 1] = n/(sum(e*x) + sum((1-e)*(x-1))) # Sum part is wrong?
k = k + 1
