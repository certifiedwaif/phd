library(optparse)
library(parallel)
source("generate.R")
library(zipvb)
source("zero_inflated_model.R")
source("accuracy.R")
source("median_accuracy.R")

i <- 4
test <- "slope"
mult <- generate_test_case(i, test)
mult$vy[67] <- 20
stan_fit <- mcmc(mult, iterations=3e4, warmup=5e3, mc.cores = 1)
