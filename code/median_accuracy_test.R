library(optparse)
library(parallel)
library(zipvb)
source("generate.R")
source("accuracy.R")
source("median_accuracy.R")

i <- 4
test <- "slope"
mult <- generate_test_case(i, test)
mult$vy[67] <- 20
stan_fit <- mcmc(mult, iterations = 30000, warmup = 5000, mc.cores = 1)
