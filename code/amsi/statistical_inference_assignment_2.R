# Statistical inference Assignment 2 Question 1 ----
# Simulate the data, assuming Pi is known ----
set.seed(2012)
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
hist(x, breaks=100)

# EM algorithm ----
em = function(estimate_pi=FALSE)
{
  # Initialise ----
  iterations = 10
  lambda = rep(NA, iterations)
  lambda[1] = 1/mean(x)
  for (k in 1:iterations) {    
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
      #print(sprintf("i %f x[i] %f num %f denom %f e[i] %f", i, x[i], num, denom, e[i]))
    }
    # M step ----
    #sum(e*x)
    #sum((1-e)*(x-1))
    lambda[k + 1] = n/(sum(e*x) + sum((1-e)*(x-1))) # Sum part is wrong?
    print(sprintf("Iteration %f: lambda %f", k, lambda[k + 1]))
    if (estimate_pi==TRUE) {
      pi = sum(e)/n
      print(sprintf("Iteration %f: pi %f", k, pi))
    }
  }
  return(list(lambda=lambda[k], pi=pi))
}
em(estimate_pi=FALSE)
em(estimate_pi=TRUE)

# Question 2 ----
# Visualising influence functions
# (e) ----
IF = function(y)
{
  multiplicand = 2*sqrt(pi)
  return(multiplicand*(pnorm(y) - .5))
}
pdf("influence_function1.pdf")
curve(IF, xlim=c(0, 4))
dev.off()

IF2 = function(y)
{
  multiplicand = 4*sqrt(pi)
  return(multiplicand*y*dnorm(y))
}
pdf("influence_function2.pdf")
curve(IF2, xlim=c(0, 4))
dev.off()
