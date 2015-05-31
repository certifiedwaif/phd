# Contents:
# -- linear mixed-effects models---random intercepts only, independent random intercepts and slopes, correlated random intercepts and slopes (simulated data, R analysis, JAGS analysis)
# -- modeling 2 correlated random effects with a scaled Inverse Wishart instead (simulated data, R analysis, JAGS analysis)
# -- general correlated random effects with a scaled Inverse Wishart (simulated data, R analysis, JAGS analysis)
# -- the material in this presentation is based on Kery (2010), Gelman & Hill (2007), and notes here http://www.unc.edu/courses/2010fall/ecol/563/001/docs/lectures/lecture28.htm and here http://people.tamu.edu/~daveamp/Bugs.ppt


# Linear mixed-effects models

# Mixed-effects or mixed models contain factors, or more generally covariates, with both fixed and random effects.
# -- we constrain the values for at least one set of effects (intercepts and/or slopes) to come from a normal distribution
# -- this is what the random-effects assumption means (usually; in general, the random effects can come from any distribution)

# There are at least three sets of assumptions that one may make about the random effects for the intercept and/or the slope of regression lines that are fitted to grouped data:
### 1. only intercepts are random, but slopes are identical for all groups (e.g., subjects, items etc.)
### 2. both intercepts and slopes are random, but they are independent
### 3. both intercepts and slopes are random and there is a correlation between them

### (an additional case, where slopes are random and intercepts are fixed, is not a sensible model in most circumstances)


# Model No. 1 is often called a random-intercepts model.
# Both models No. 2 and 3 are called random-coefficients models.
# Model No. 3 is the default in R's function lmer() in package lme4 when fitting a random-coefficients model.


# The plan:
# -- [Model 1 and Model 2] first, we generate a random-coefficients dataset under model No. 2, where both intercepts and slopes are uncorrelated random effects; we then fit both a random-intercepts (No. 1) and a random-coefficients model without correlation (No. 2) to this dataset
# -- [Model 3] next, we generate a second data set that includes a correlation between random intercepts and random slopes and adopt the random-coefficients model with correlation between intercepts and slopes (No. 3) to analyze it
# -- [Model 3, ctd.] we estimate this model with JAGS in 2 different ways: (i) by separately modeling the sd.s of the 2 random-effect distributions (the intercepts and the slopes) and the correlation between the 2 distributions; (ii) by modeling all three parameters simultaneously with a scaled inverse Wishart prior; this second way has much better mixing and it also generalizes immediately to more complex random effects structures involving 3 or more correlated random effects
# -- [Model 4] finally, we generate a third data set that includes 3 correlated random effects (intercepts and 2 distinct slope vectors) and show how the random-coefficients model with correlation between intercepts and slopes generalizes to account for this kind of random effect structures; we only estimate this with a scaled inverse Wishart prior


# A close examination of how such a dataset can be assembled (i.e., simulated) will help us better understand how analogous datasets are broken down (i.e., analyzed) using mixed models.
# -- as Kery (2010) puts it: very few strategies can be more effective to understand this type of mixed model than the combination of simulating data sets and describing the models fitted in BUGS syntax



# Data generation

# The factor for subjects:
n_subjects <- 56 # Number of subjects
n_sample <- 10 # Number of observations for each subject
(n_obs <- n_subjects*n_sample) # Total number of data points
(subjects <- gl(n=n_subjects, k=n_sample)) # Indicator for subjects

# Continuous predictor x
original_x <- runif(n_obs, 45, 70)
summary(original_x)

# We standardize it (always good to do with continuous predictors when using JAGS/BUGS)
(mean_orig_x <- mean(original_x))
(sd_orig_x <- sd(original_x))
x <- (original_x-mean_orig_x)/sd_orig_x
round(x, 2)
summary(x)

hist(x, col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(x), col="lightsteelblue3", lwd=2)


# This is the model matrix for the means parametrization of the interaction model between the subjects and the continuous covariate x:
Xmat <- model.matrix(~subjects*x-1-x)
dim(Xmat)

# Q: where do these dimensions come from?
head(Xmat)

# -- there are 560 observations (rows) and 112 regression terms / variables (columns)
dimnames(Xmat)[[1]]
dimnames(Xmat)[[2]]

# -- there are 56 terms for subjects, the coefficients of which will provide the subject-specific intercepts (i.e., the intercept random effects, or intercept effects for short)
# -- there are 56 terms for interactions between each subject and the continuous covariate x, the coefficients of which will provide the subject-specific slopes (i.e., the slope random effects, or slope effects for short)

round(Xmat[1, ], 2) 		# Print the top row for each column
Xmat[, 1]               # Print all rows for column 1 (group 1)
round(Xmat[, 57], 2)    # Print all rows for column 57 (group 1:x)


# Parameters for the distributions of the random coefficients / random effects (note that the intercepts and slopes comes from two independent Gaussian distributions):

intercept_mean <- 230 # mu_alpha
intercept_sd <- 20 # sigma_alpha

slope_mean <- 60 # mu_beta
slope_sd <- 30 # sigma_beta


# Generate the random coefficients:

intercept_effects <- rnorm(n=n_subjects, mean=intercept_mean, sd=intercept_sd)
slope_effects <- rnorm(n=n_subjects, mean=slope_mean, sd=slope_sd)

par(mfrow=c(1, 2))
hist(intercept_effects, col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(intercept_effects), col="lightsteelblue3", lwd=2)
hist(slope_effects, col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(slope_effects), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))

all_effects <- c(intercept_effects, slope_effects) # Put them all together
round(all_effects, 2)

# -- thus, we have two stochastic components in our model IN ADDITION TO the usual stochastic component for the individual-level responses, to which we now turn


# Generating the continuous response variable:

# -- the deterministic part
lin_pred <- Xmat %*% all_effects # Value of lin_predictor
str(lin_pred)

# -- the stochastic part
sigma_res <- 30
normal_error <- rnorm(n=n_obs, mean=0, sd=sigma_res) # residuals
str(normal_error)

# -- put the two together
y <- lin_pred+normal_error
str(y)
# or, alternatively
y <- rnorm(n=n_obs, mean=lin_pred, sd=sigma_res)
str(y)

# We take a look at the response variable
hist(y, col="lightsteelblue1", border="white", breaks=30, freq=FALSE)
lines(density(y), col="lightsteelblue3", lwd=2)
summary(y)

library("lattice")
xyplot(y~x|subjects)


# Analysis under a random-intercepts model

# REML analysis using R
library("lme4")
lme_fit1 <- lmer(y~x+(1|subjects))
print(lme_fit1, cor=FALSE)

fixef(lme_fit1)
ranef(lme_fit1)
coef(lme_fit1)

# Compare with true values:
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(lme_fit1, cor=FALSE)


# Bayesian analysis using JAGS

# Write model
cat("model {
# Priors
mu_int~dnorm(0, 0.0001) # Mean hyperparameter for random intercepts
sigma_int~dunif(0, 100) # SD hyperparameter for random intercepts
tau_int <- 1/(sigma_int*sigma_int)
for (i in 1:n_subj) {
    alpha[i]~dnorm(mu_int, tau_int) # Random intercepts
}
beta~dnorm(0, 0.0001) # Common slope
sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res)
# Likelihood
for (i in 1:n_obs) {
    mu[i] <- alpha[subjects[i]]+beta*x[i] # Expectation
    y[i]~dnorm(mu[i], tau_res) # The actual (random) responses
}
}", fill=TRUE, file="lme_model1.txt")

# Bundle data
jags_data <- list(y=as.numeric(y), subjects=as.numeric(subjects), x=as.numeric(x), n_subj=max(as.numeric(subjects)), n_obs=as.numeric(n_obs))
# use as.numeric across the board for the data passed to JAGS; it might work w/o it, but this is often needed for other BUGS packages

# Inits function
inits <- function() {
    list(alpha=rnorm(n_subjects, 0, 2), beta=rnorm(1, 1, 1), mu_int=rnorm(1, 0, 1), sigma_int=rlnorm(1), sigma_res=rlnorm(1))
}

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "sigma_res")

# MCMC settings
ni <- 11000; nb <- 1000; nt <- 20; nc <- 3

# Start Gibbs sampling
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model1.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

traceplot(outj)

# Inspect results
print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(out$mean, dig=4)
print(lme_fit1, cor=FALSE)

# The same model can be parametrized by making the mean of the intercept ranefs a fixed effect and modeling the subject ranefs as coming from a normal distribution centered at 0. This is in fact how the lmer function models ranefs.

# alternative model
cat("model {
# Priors
sigma_int~dunif(0, 100) # SD hyperparameter for random intercepts
tau_int <- 1/(sigma_int*sigma_int)
for (i in 1:n_subj) {
    alpha[i]~dnorm(0, tau_int) # Random by-subject deflections to the intercept
}
mu_int~dnorm(0, 0.0001) # The mean intercept
beta~dnorm(0, 0.0001) # Common slope
sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res)
# Likelihood
for (i in 1:n_obs) {
    mu[i] <- mu_int + alpha[subjects[i]] + beta*x[i] # Expectation
    y[i]~dnorm(mu[i], tau_res) # The actual (random) responses
}
}", fill=TRUE, file="lme_model1_1.txt")

# Bundle data
jags_data <- list(y=as.numeric(y), subjects=as.numeric(subjects), x=as.numeric(x), n_subj=max(as.numeric(subjects)), n_obs=as.numeric(n_obs))
# use as.numeric across the board for the data passed to JAGS; it might work w/o it, but this is often needed for other BUGS packages

# Inits function
inits <- function() {
    list(alpha=rnorm(n_subjects, 0, 2), beta=rnorm(1, 1, 1), mu_int=rnorm(1, 0, 1), sigma_int=rlnorm(1), sigma_res=rlnorm(1))
}

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "sigma_res")

# MCMC settings
ni <- 11000; nb <- 1000; nt <- 20; nc <- 3

# Start Gibbs sampling
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model1_1.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

traceplot(outj)

# Inspect results
print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(out$mean, dig=4)
print(lme_fit1, cor=FALSE)


# Analysis under a random-coefficients model without correlation between intercept and slope

# REML analysis using R
library("lme4")
lme_fit2 <- lmer(y~x+(1|subjects)+(0+x|subjects))
print(lme_fit2, cor=F)
fixef(lme_fit2)
ranef(lme_fit2)
coef(lme_fit2)


# Compare with true values:
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(lme_fit2, cor=FALSE)


# Bayesian analysis using JAGS

# Define model
cat("model {
# Priors
mu_int~dnorm(0, 0.001) # Mean hyperparameter for random intercepts
sigma_int~dunif(0, 100) # SD hyperparameter for random intercepts
tau_int <- 1/(sigma_int*sigma_int)
mu_slope~dnorm(0, 0.001) # Mean hyperparameter for random slopes
sigma_slope~dunif(0, 100) # SD hyperparameter for slopes
tau_slope <- 1/(sigma_slope*sigma_slope)
for (i in 1:n_subj) {
    alpha[i]~dnorm(mu_int, tau_int) # Random intercepts
    beta[i]~dnorm(mu_slope, tau_slope) # Random slopes
}
sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res) # Residual precision
# Likelihood
for (i in 1:n_obs) {
    mu[i] <- alpha[subjects[i]]+beta[subjects[i]]*x[i]
    y[i]~dnorm(mu[i], tau_res)
}
}", fill=TRUE, file="lme_model2.txt")

# Bundle data
jags_data <- list(y=as.numeric(y), subjects=as.numeric(subjects), x=as.numeric(x), n_subj=max(as.numeric(subjects)), n_obs=as.numeric(n_obs))

# Inits function
inits <- function() {
    list(alpha=rnorm(n_subjects, 0, 2), beta=rnorm(n_subjects, 10, 2), mu_int=rnorm(1, 0, 1), sigma_int=rlnorm(1), mu_slope=rnorm(1, 0, 1), sigma_slope=rlnorm(1), sigma_res=rlnorm(1))
}

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "mu_slope", "sigma_slope", "sigma_res")

# MCMC settings
ni <- 11000; nb <- 1000; nt <- 10; nc <- 3

# Start Gibbs sampling
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model2.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

traceplot(outj)

# -- you can see that the chains are mixing much better when we use the appropriate model for the data
# -- see Gelman's `folk theorem of statistical computing' (http://andrewgelman.com/2008/05/13/the_folk_theore/):
# when you have computational problems, often there's a problem with your model

print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, sigma_res=sigma_res)
print(out$mean, dig=4)
print(lme_fit2, cor=FALSE)

# Using simulated data and successfully recovering the input values makes us fairly confident that the JAGS analysis has been correctly specified.

# This is very helpful for more complex models b/c it's easy to make mistakes:
# -- a good way to check the JAGS analysis for a custom model that is needed for a particular phenomenon is to simulate the data and run the JAGS model on that data


# [Model 3] The random-coefficients model with correlation between intercept and slope

# Data generation

n_subjects <- 56
n_sample <- 10
(n_obs <- n_subjects*n_sample)
(subjects <- gl(n=n_subjects, k=n_sample))

# Standardized continuous covariate:
original_x <- runif(n_obs, 45, 70)
(mean_orig_x <- mean(original_x))
(sd_orig_x <- sd(original_x))
x <- (original_x-mean_orig_x)/sd_orig_x

hist(x, col="lightsteelblue1", border="white", breaks=20, freq=FALSE)
lines(density(x), col="lightsteelblue3", lwd=2)

# Design matrix:
Xmat <- model.matrix(~subjects*x-1-x)

# -- there are 560 observations (rows) and 112 regression terms / variables (columns), just as before
dimnames(Xmat)[[1]]
dimnames(Xmat)[[2]]

round(Xmat[1, ], 2) # Print the top row for each column

# Generate the correlated random effects for intercept and slope:


# Assembling the parameters for the multivariate normal distribution

intercept_mean <- 230 # Values for five hyperparameters
intercept_sd <- 20
slope_mean <- 60
slope_sd <- 30
intercept_slope_covariance <- 10
intercept_slope_correlation <- intercept_slope_covariance/(intercept_sd*slope_sd)

(mu_vector <- c(intercept_mean, slope_mean))
(var_covar_matrix <- matrix(c(intercept_sd^2, intercept_slope_covariance, intercept_slope_covariance, slope_sd^2), 2, 2))

# Generating the correlated random effects for intercepts and slopes:

library("MASS") # Load MASS to sample from a multivariate normal
effects <- mvrnorm(n=n_subjects, mu=mu_vector, Sigma=var_covar_matrix)
round(effects, 2)

par(mfrow=c(1, 2))
hist(effects[, 1], col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(effects[, 1]), col="lightsteelblue3", lwd=2)

hist(effects[, 2], col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(effects[, 2]), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))

# Plotting the bivariate distribution:
effects_kde <- kde2d(effects[, 1], effects[, 2], n=50) # kernel density estimate
par(mfrow=c(1, 3))
contour(effects_kde)
image(effects_kde)
persp(effects_kde, phi=45, theta=30)

# even better:
par(mfrow=c(1, 2))
image(effects_kde); contour(effects_kde, add=T)
persp(effects_kde, phi=45, theta=-30, shade=.1, border=NULL, col="lightsteelblue1", ticktype="detailed", xlab="", ylab="", zlab="")
par(mfrow=c(1, 1))

apply(effects, 2, mean)
apply(effects, 2, sd)
apply(effects, 2, var)
cov(effects[, 1], effects[, 2])
var(effects)

# Population parameters
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, intercept_slope_covariance=intercept_slope_covariance, sigma_res=sigma_res)

# Sampling error for intercept-slope covariance (200 samples of 50, 500, 5000 and 50000 group/random effects each) -- just to illustrate how difficult it is to estimate measures of association, even from fairly large data sets with 5000 observations:

par(mfrow=c(1, 4))
cov_temp1 <- numeric()
for (i in 1:200) {
    temp1 <- mvrnorm(50, mu=mu_vector, Sigma=var_covar_matrix)
    cov_temp1[i] <- var(temp1)[1, 2]
}
hist(cov_temp1, col="lightsteelblue1", border="white", freq=FALSE, main="n_obs=50")
lines(density(cov_temp1), col="lightsteelblue3", lwd=2)

cov_temp2 <- numeric()
for (i in 1:200) {
    temp2 <- mvrnorm(500, mu=mu_vector, Sigma=var_covar_matrix)
    cov_temp2[i] <- var(temp2)[1, 2]
}
hist(cov_temp2, col="lightsteelblue1", border="white", freq=FALSE, main="n_obs=500")
lines(density(cov_temp2), col="lightsteelblue3", lwd=2)

cov_temp3 <- numeric()
for (i in 1:200) {
    temp3 <- mvrnorm(5000, mu=mu_vector, Sigma=var_covar_matrix)
    cov_temp3[i] <- var(temp3)[1, 2]
}
hist(cov_temp3, col="lightsteelblue1", border="white", freq=FALSE, main="n_obs=5000")
lines(density(cov_temp3), col="lightsteelblue3", lwd=2)

cov_temp4 <- numeric()
for (i in 1:200) {
    temp4 <- mvrnorm(50000, mu=mu_vector, Sigma=var_covar_matrix)
    cov_temp4[i] <- var(temp4)[1, 2]
}
hist(cov_temp4, col="lightsteelblue1", border="white", freq=FALSE, main="n_obs=50000")
lines(density(cov_temp4), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))


intercept_effects <- effects[, 1]
round(intercept_effects, 2)
slope_effects <- effects[, 2]
round(slope_effects, 2)
all_effects <- c(intercept_effects, slope_effects) # Put them all together
round(all_effects, 2)


# Generate the response variable:
# -- the deterministic part
lin_pred <- Xmat %*% all_effects
round(as.vector(lin_pred), 2)

# -- the stochastic part
sigma_res <- 30
(normal_error <- rnorm(n=n_obs, mean=0, sd=sigma_res))	# residuals

# -- add them together
y <- lin_pred+normal_error
# or, in one go:
y <- rnorm(n=n_obs, mean=lin_pred, sd=sigma_res)

hist(y, col="lightsteelblue1", border="white", breaks=15, freq=FALSE)
lines(density(y), col="lightsteelblue3", lwd=2)

library("lattice")
xyplot(y~x|subjects, pch=20)


# REML analysis using R
library("lme4")
lme_fit3 <- lmer(y~x+(x|subjects))
print(lme_fit3, cor=FALSE)

# Compare with the true values:
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, intercept_slope_correlation=intercept_slope_covariance/(intercept_sd*slope_sd), sigma_res=sigma_res)



# Bayesian analysis using JAGS

# This is one way in which we can specify a Bayesian analysis of the random-coefficients model with correlation.
# -- it is more intuitive but does not generalize well to more than 2 correlated random effects.

# We will introduce a different and more general way to allow for correlation among two or more sets of random effects in a model after this.

# Define model
cat("model {
# Priors
mu_int~dnorm(0, 0.0001) # mean for random intercepts
mu_slope~dnorm(0, 0.0001) # mean for random slopes
sigma_int~dunif(0, 100) # SD of intercepts
sigma_slope~dunif(0, 100) # SD of slopes
rho~dunif(-1, 1) # correlation between intercepts and slopes
Sigma_B[1, 1] <- pow(sigma_int, 2) # We start assembling the var-covar matrix for the random effects
Sigma_B[2, 2] <- pow(sigma_slope, 2)
Sigma_B[1, 2] <- rho*sigma_int*sigma_slope
Sigma_B[2, 1] <- Sigma_B[1, 2]
covariance <- Sigma_B[1, 2]
Tau_B[1:2, 1:2] <- inverse(Sigma_B[,])
for (i in 1:n_subj) {
    B_hat[i, 1] <- mu_int
    B_hat[i, 2] <- mu_slope
    B[i, 1:2]~dmnorm(B_hat[i, ], Tau_B[,]) # the pairs of correlated random effects
    alpha[i] <- B[i, 1] # random intercept
    beta[i] <- B[i, 2] # random slope
}
sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res)
# Likelihood
for (i in 1:n_obs) {
    mu[i] <- alpha[subjects[i]]+beta[subjects[i]]*x[i]
    y[i]~dnorm(mu[i], tau_res)
}
}", fill=TRUE, file="lme_model3.txt")

# Bundle data
jags_data <- list(y=as.numeric(y), subjects=as.numeric(subjects), x=as.numeric(x), n_subj=max(as.numeric(subjects)), n_obs=as.numeric(n_obs))

# Inits function
inits <- function() {
    list(mu_int=rnorm(1, 0, 1), sigma_int=rlnorm(1), mu_slope=rnorm(1, 0, 1), sigma_slope=rlnorm(1), rho=runif(1, -1, 1), sigma_res=rlnorm(1))
}

# Parameters to estimate
params <- c("alpha", "beta", "mu_int", "sigma_int", "mu_slope", "sigma_slope", "rho", "covariance", "sigma_res")

# MCMC settings
ni <- 3200; nb <- 200; nt <- 6; nc <- 3 # more than this probably needed for a good approx. of the posterior distribution

# Start Gibbs sampler
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model3.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
traceplot(outj)

# -- this type of model does not converge very fast in addition to the fact that it does not generalize very well; we increase the number of iterations

ni <- 25000; nb <- 5000; nt <- 40; nc <- 3
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model3.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)
traceplot(outj)

print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, intercept_slope_correlation=intercept_slope_covariance/(intercept_sd*slope_sd), sigma_res=sigma_res)
print(out$mean, dig=4)
print(lme_fit3, cor=FALSE)

# Note the very large SD for the posterior distribution of the covariance (relative to the mean):
print(out$mean$covariance, dig=2)
print(out$sd$covariance, dig=2)

print(out$mean$rho, dig=2)
print(out$sd$rho, dig=2)

# -- R does not even provide an SE for the covariance estimator (equivalently, for the correlation of random effects)
# -- covariances are even harder to reliably estimate than variances, which are harder than mean estimators (it's easy to estimate measures of center / location, harder to estimate measures of dispersion and even harder to estimate measures of association)


# MODELING CORRELATED RANEFS WITH A SCALED INVERSE WISHART: the 2 correlated ranefs case first.

# -- we will now introduce an alternative way of placing priors over correlated random effects that both converges faster and generalizes to structures with more than 2 correlated random effects.
# -- see the relevant chapter of Gelman & Hill (2007) for more introductory discussion and references

# Bayesian analysis using a scaled Inverse Wishart

# Define model
cat("model {

# Set up the means for the multivariate ranef distribution
for (i in 1:2) {
    xi[i]~dunif(0, 100) # scaling for the multivariate ranef distribution (for means, sds, and the ranefs themselves)
    mu_raw[i]~dnorm(0, .0001) # unscaled means for the multivariate ranef distribution
    mu[i] <- xi[i]*mu_raw[i] # scaled means for the multivariate ranef distribution
}
mu_int <- mu[1] # mean for random intercepts
mu_slope <- mu[2] # mean for random slopes

# Set up the var-covar matrix for the multivariate ranef distribution
Tau_B_raw[1:2, 1:2] ~ dwish(W[,], 3) # W is the identity matrix, provided as data; we have 3 dofs, i.e., 2 ranefs + 1, to ensure a uniform (-1, 1) prior for the correlation between ranefs
Sigma_B_raw[1:2, 1:2] <- inverse(Tau_B_raw[,])
for (i in 1:2) {
    sigma[i] <- xi[i]*sqrt(Sigma_B_raw[i, i])
}
sigma_int <- sigma[1] # SD of intercepts
sigma_slope <- sigma[2] # SD of slopes
for (i in 1:2) { for (j in 1:2) {
    rho[i, j] <- Sigma_B_raw[i, j]/sqrt(Sigma_B_raw[i, i]*Sigma_B_raw[j, j])
} }
rho_int_slope <- rho[1, 2]
covariance <- rho_int_slope*sigma_int*sigma_slope

# The multivariate ranef distribution, i.e., modeling the correlated ranefs
for (j in 1:n_subj) {
	B_raw_hat[j, 1] <- mu_raw[1]
	B_raw_hat[j, 2] <- mu_raw[2]
	B_raw[j, 1:2] ~ dmnorm(B_raw_hat[j, ], Tau_B_raw[, ]) # the pairs of unscaled (raw) correlated random effects
	alpha[j] <- xi[1]*B_raw[j, 1] # random intercept
    beta[j] <- xi[2]*B_raw[j, 2] # random slope
}

# Model the resid. sd independently
sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res)

# Likelihood
for (i in 1:n_obs) {
    mu_obs[i] <- alpha[subjects[i]]+beta[subjects[i]]*x[i]
    y[i]~dnorm(mu_obs[i], tau_res)
}

# Sampling from the prior: given that we do not place hyperpriors directly on the means, sds and correlation(s) of the multivariate ranef distribution, we want to sample from the prior to make sure we didn't accidentally make it more informed than we wanted (and we want it very vague)
for (i in 1:2) {
    xi_prior[i]~dunif(0, 100)
    mu_raw_prior[i]~dnorm(0, .0001)
    mu_prior[i] <- xi_prior[i]*mu_raw_prior[i]
}
mu_int_prior <- mu_prior[1]
mu_slope_prior <- mu_prior[2]
Tau_B_raw_prior[1:2, 1:2] ~ dwish(W[,], 3)
Sigma_B_raw_prior[1:2, 1:2] <- inverse(Tau_B_raw_prior[,])
for (i in 1:2) {
    sigma_prior[i] <- xi_prior[i]*sqrt(Sigma_B_raw_prior[i, i])
}
sigma_int_prior <- sigma_prior[1]
sigma_slope_prior <- sigma_prior[2]
for (i in 1:2) { for (j in 1:2) {
    rho_prior[i, j] <- Sigma_B_raw_prior[i, j]/sqrt(Sigma_B_raw_prior[i, i]*Sigma_B_raw_prior[j, j])
} }
rho_int_slope_prior <- rho_prior[1, 2]
}", fill=TRUE, file="lme_model4.txt")

# Bundle data
jags_data <- list(y=as.numeric(y), subjects=as.numeric(subjects), x=as.numeric(x), n_subj=max(as.numeric(subjects)), n_obs=as.numeric(n_obs), W=diag(2))

#install.packages("bayesm")
library("bayesm")
var_vec <- apply(coef(lme_fit3)$subjects, 2, var)

# Inits function
inits <- function() {
    list(xi=rlnorm(2), mu_raw=rnorm(2), Tau_B_raw=rwishart(3, diag(2)*var_vec)$W, sigma_res=rlnorm(1), xi_prior=rlnorm(2), mu_raw_prior=rnorm(2), Tau_B_raw_prior=rwishart(3, diag(2)*var_vec)$W)
}

# Parameters to estimate
params <- c("mu", "mu_int", "mu_slope", "sigma", "sigma_int", "sigma_slope", "rho", "rho_int_slope", "covariance", "alpha", "beta", "sigma_res", "mu_int_prior", "mu_slope_prior", "sigma_int_prior", "sigma_slope_prior", "rho_int_slope_prior")

# MCMC settings
ni <- 7000; nb <- 1000; nt <- 6; nc <- 3 # more than this probably needed for a good approx. of the posterior distribution

# Start Gibbs sampler
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model4.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

traceplot(outj)

print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
data.frame(intercept_mean=intercept_mean, slope_mean=slope_mean, intercept_sd=intercept_sd, slope_sd=slope_sd, intercept_slope_correlation=intercept_slope_covariance/(intercept_sd*slope_sd), sigma_res=sigma_res)
print(out$mean, dig=4)
print(lme_fit3, cor=FALSE)

# Once again, note the very large SD for the posterior distribution of the covariance (relative to the mean):
print(out$mean$covariance, dig=2)
print(out$sd$covariance, dig=2)

print(out$mean$rho, dig=2)
print(out$sd$rho, dig=2)

# Finally, we compare the prior and posterior distributions for the (derived) ranef distribution parameters:

par(mfrow=c(3, 4))
hist(out$sims.list$mu_int_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_int_prior")
lines(density(out$sims.list$mu_int_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$mu_int, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_int")
lines(density(out$sims.list$mu_int), col="lightsteelblue3", lwd=2)

hist(out$sims.list$mu_slope_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_slope_prior")
lines(density(out$sims.list$mu_slope_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$mu_slope, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_slope")
lines(density(out$sims.list$mu_slope), col="lightsteelblue3", lwd=2)

hist(out$sims.list$sigma_int_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_int_prior")
lines(density(out$sims.list$sigma_int_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$sigma_int, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_int")
lines(density(out$sims.list$sigma_int), col="lightsteelblue3", lwd=2)

hist(out$sims.list$sigma_slope_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_slope_prior")
lines(density(out$sims.list$sigma_slope_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$sigma_slope, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_slope")
lines(density(out$sims.list$sigma_slope), col="lightsteelblue3", lwd=2)

hist(out$sims.list$rho_int_slope_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_int_slope_prior")
lines(density(out$sims.list$rho_int_slope_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$rho_int_slope, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_int_slope")
lines(density(out$sims.list$rho_int_slope), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))


# [Model 4] Finally, we simulate data from a model with 3 correlated ranefs and then provide the R and JAGS analyses of this data

# Data generation

n_subjects <- 56
n_sample <- 10
(n_obs <- n_subjects*n_sample)
(subjects <- gl(n=n_subjects, k=n_sample))

# Two standardized continuous covariates:
original_x1 <- runif(n_obs, 45, 70)
(mean_orig_x1 <- mean(original_x1))
(sd_orig_x1 <- sd(original_x1))
x1 <- (original_x1-mean_orig_x1)/sd_orig_x1

original_x2 <- rgamma(n_obs, 10, 1)
(mean_orig_x2 <- mean(original_x2))
(sd_orig_x2 <- sd(original_x2))
x2 <- (original_x2-mean_orig_x2)/sd_orig_x2

par(mfrow=c(1, 2))
hist(x1, col="lightsteelblue1", border="white", breaks=15, freq=FALSE)
lines(density(x1), col="lightsteelblue3", lwd=2)

hist(x2, col="lightsteelblue1", border="white", breaks=15, freq=FALSE)
lines(density(x2), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))

# Design matrix:
Xmat <- model.matrix(~subjects*x1+subjects*x2-1-x1-x2)

# -- there are 560 observations (rows) and 112 regression terms / variables (columns), just as before
dimnames(Xmat)[[1]]
dimnames(Xmat)[[2]]

round(Xmat[1, ], 2) # Print the top row for each column

# Generate the correlated random effects for intercept and slope:


# Assembling the parameters for the multivariate normal distribution

intercept_mean <- 230
intercept_sd <- 20
slope1_mean <- 60
slope1_sd <- 30
slope2_mean <- 40
slope2_sd <- 25
intercept_slope1_covariance <- 10
(intercept_slope1_correlation <- intercept_slope1_covariance/(intercept_sd*slope1_sd))
intercept_slope2_covariance <- 300
(intercept_slope2_correlation <- intercept_slope2_covariance/(intercept_sd*slope2_sd))
slope1_slope2_covariance <- -200
(slope1_slope2_correlation <- slope1_slope2_covariance/(slope1_sd*slope2_sd))

(mu_vector <- c(intercept_mean, slope1_mean, slope2_mean))
(var_covar_matrix <- matrix(c(intercept_sd^2, intercept_slope1_covariance, intercept_slope2_covariance, intercept_slope1_covariance, slope1_sd^2, slope1_slope2_covariance, intercept_slope2_covariance, slope1_slope2_covariance, slope2_sd^2), 3, 3))
eigen(var_covar_matrix)$values

# Generating the correlated random effects for intercepts and slopes:

library("MASS")
effects <- mvrnorm(n=n_subjects, mu=mu_vector, Sigma=var_covar_matrix)
str(effects)
round(effects, 2)

par(mfrow=c(1, 3))
hist(effects[, 1], col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(effects[, 1]), col="lightsteelblue3", lwd=2)

hist(effects[, 2], col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(effects[, 2]), col="lightsteelblue3", lwd=2)

hist(effects[, 3], col="lightsteelblue1", border="white", breaks=10, freq=FALSE)
lines(density(effects[, 3]), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))

# Plotting the bivariate distributions:
effects_kde12 <- kde2d(effects[, 1], effects[, 2], n=50)
effects_kde13 <- kde2d(effects[, 1], effects[, 3], n=50)
effects_kde23 <- kde2d(effects[, 2], effects[, 3], n=50)

par(mfrow=c(3, 2))
image(effects_kde12); contour(effects_kde12, add=T)
persp(effects_kde12, phi=45, theta=-30, shade=.1, border=NULL, col="lightsteelblue1", ticktype="detailed", xlab="", ylab="", zlab="")
image(effects_kde13); contour(effects_kde13, add=T)
persp(effects_kde13, phi=45, theta=-30, shade=.1, border=NULL, col="lightsteelblue1", ticktype="detailed", xlab="", ylab="", zlab="")
image(effects_kde23); contour(effects_kde23, add=T)
persp(effects_kde23, phi=45, theta=-30, shade=.1, border=NULL, col="lightsteelblue1", ticktype="detailed", xlab="", ylab="", zlab="")
par(mfrow=c(1, 1))

apply(effects, 2, mean)
apply(effects, 2, sd)
apply(effects, 2, var)
cov(effects[, 1], effects[, 2])
cov(effects[, 1], effects[, 3])
cov(effects[, 2], effects[, 3])
var(effects)

sigma_res <- 30

# Population parameters
data.frame(intercept_mean=intercept_mean, slope1_mean=slope1_mean, slope2_mean=slope2_mean, intercept_sd=intercept_sd, slope1_sd=slope1_sd, slope2_sd=slope2_sd, intercept_slope1_covariance=intercept_slope1_covariance, intercept_slope2_covariance=intercept_slope2_covariance, slope1_slope2_covariance=slope1_slope2_covariance, sigma_res=sigma_res)

intercept_effects <- effects[, 1]
round(intercept_effects, 2)
slope1_effects <- effects[, 2]
round(slope1_effects, 2)
slope2_effects <- effects[, 3]
round(slope2_effects, 2)
all_effects <- c(intercept_effects, slope1_effects, slope2_effects)


# Generate the response variable:
# -- the deterministic part
lin_pred <- Xmat %*% all_effects
round(as.vector(lin_pred), 2)

# -- the stochastic part
sigma_res <- 30
(normal_error <- rnorm(n=n_obs, mean=0, sd=sigma_res))	# residuals

# -- add them together
y <- lin_pred+normal_error
# or, in one go:
y <- rnorm(n=n_obs, mean=lin_pred, sd=sigma_res)

hist(y, col="lightsteelblue1", border="white", breaks=15, freq=FALSE)
lines(density(y), col="lightsteelblue3", lwd=2)

library("lattice")
xyplot(y~x1+x2|subjects, pch=20)

data.frame(y, x1, x2, subjects)
# REML analysis using R
library("lme4")
lme_fit4 <- lmer(y~x1+x2+(x1+x2|subjects))
print(lme_fit4, cor=FALSE)

# Compare with the true values:
data.frame(intercept_mean=intercept_mean, slope1_mean=slope1_mean, slope2_mean=slope2_mean, intercept_sd=intercept_sd, slope1_sd=slope1_sd, slope2_sd=slope2_sd, intercept_slope1_covariance=intercept_slope1_covariance, intercept_slope2_covariance=intercept_slope2_covariance, slope1_slope2_covariance=slope1_slope2_covariance, sigma_res=sigma_res)

# Bayesian analysis using a scaled Inverse Wishart

# Define model
cat("model {

# Set up the means for the multivariate ranef distribution
for (i in 1:3) {
    xi[i]~dunif(0, 100)
    mu_raw[i]~dnorm(0, .0001)
    mu[i] <- xi[i]*mu_raw[i]
}
mu_int <- mu[1]
mu_slope1 <- mu[2]
mu_slope2 <- mu[3]

Tau_B_raw[1:3, 1:3] ~ dwish(W[,], 4)
Sigma_B_raw[1:3, 1:3] <- inverse(Tau_B_raw[,])
for (i in 1:3) {
    sigma[i] <- xi[i]*sqrt(Sigma_B_raw[i, i])
}
sigma_int <- sigma[1]
sigma_slope1 <- sigma[2]
sigma_slope2 <- sigma[3]
for (i in 1:3) { for (j in 1:3) {
    rho[i, j] <- Sigma_B_raw[i, j]/sqrt(Sigma_B_raw[i, i]*Sigma_B_raw[j, j])
} }
rho_int_slope1 <- rho[1, 2]
rho_int_slope2 <- rho[1, 3]
rho_slope1_slope2 <- rho[2, 3]

for (j in 1:n_subj) {
	B_raw_hat[j, 1] <- mu_raw[1]
	B_raw_hat[j, 2] <- mu_raw[2]
	B_raw_hat[j, 3] <- mu_raw[3]
	B_raw[j, 1:3] ~ dmnorm(B_raw_hat[j, ], Tau_B_raw[, ])
	alpha[j] <- xi[1]*B_raw[j, 1]
    beta1[j] <- xi[2]*B_raw[j, 2]
    beta2[j] <- xi[3]*B_raw[j, 3]
}

sigma_res~dunif(0, 100) # Residual standard deviation
tau_res <- 1/(sigma_res*sigma_res)

for (i in 1:n_obs) {
    mu_obs[i] <- alpha[subjects[i]]+beta1[subjects[i]]*x1[i]+beta2[subjects[i]]*x2[i]
    y[i]~dnorm(mu_obs[i], tau_res)
}

for (i in 1:3) {
    xi_prior[i]~dunif(0, 100)
    mu_raw_prior[i]~dnorm(0, .0001)
    mu_prior[i] <- xi_prior[i]*mu_raw_prior[i]
}
mu_int_prior <- mu_prior[1]
mu_slope1_prior <- mu_prior[2]
mu_slope2_prior <- mu_prior[3]
Tau_B_raw_prior[1:3, 1:3] ~ dwish(W[,], 4)
Sigma_B_raw_prior[1:3, 1:3] <- inverse(Tau_B_raw_prior[,])
for (i in 1:3) {
    sigma_prior[i] <- xi_prior[i]*sqrt(Sigma_B_raw_prior[i, i])
}
sigma_int_prior <- sigma_prior[1]
sigma_slope1_prior <- sigma_prior[2]
sigma_slope2_prior <- sigma_prior[3]
for (i in 1:3) { for (j in 1:3) {
    rho_prior[i, j] <- Sigma_B_raw_prior[i, j]/sqrt(Sigma_B_raw_prior[i, i]*Sigma_B_raw_prior[j, j])
} }
rho_int_slope1_prior <- rho_prior[1, 2]
rho_int_slope2_prior <- rho_prior[1, 3]
rho_slope1_slope2_prior <- rho_prior[2, 3]
}", fill=TRUE, file="lme_model5.txt")

# Bundle data
jags_data <- list(y=as.numeric(y), subjects=as.numeric(subjects), x1=as.numeric(x1), x2=as.numeric(x2), n_subj=max(as.numeric(subjects)), n_obs=as.numeric(n_obs), W=diag(3))

#install.packages("bayesm")
library("bayesm")
var_vec <- apply(coef(lme_fit4)$subjects, 2, var)

# Inits function
inits <- function() {
    list(xi=rlnorm(3), mu_raw=rnorm(3), Tau_B_raw=rwishart(4, diag(3)*var_vec)$W, sigma_res=rlnorm(1), xi_prior=rlnorm(3), mu_raw_prior=rnorm(3), Tau_B_raw_prior=rwishart(4, diag(3)*var_vec)$W)
}

# Parameters to estimate
params <- c("mu", "mu_int", "mu_slope1", "mu_slope2", "sigma", "sigma_int", "sigma_slope1", "sigma_slope2", "rho", "rho_int_slope1", "rho_int_slope2", "rho_slope1_slope2", "alpha", "beta1", "beta2", "sigma_res", "mu_int_prior", "mu_slope1_prior", "mu_slope2_prior", "sigma_int_prior", "sigma_slope1_prior", "sigma_slope2_prior", "rho_prior", "rho_int_slope1_prior", "rho_int_slope2_prior", "rho_slope1_slope2_prior")

# MCMC settings
ni <- 7000; nb <- 1000; nt <- 6; nc <- 3 # more than this probably needed for a good approx. of the posterior distribution

# Start Gibbs sampler
library("R2jags")
outj <- jags(jags_data, inits=inits, parameters.to.save=params, model.file="lme_model5.txt", n.thin=nt, n.chains=nc, n.burnin=nb, n.iter=ni)

traceplot(outj)

print(outj, dig=3)

out <- outj$BUGSoutput

# Compare with MLEs and true values
print(out$mean, dig=4)
data.frame(intercept_mean=intercept_mean, slope1_mean=slope1_mean, slope2_mean=slope2_mean, intercept_sd=intercept_sd, slope1_sd=slope1_sd, slope2_sd=slope2_sd, intercept_slope1_covariance=intercept_slope1_covariance, intercept_slope2_covariance=intercept_slope2_covariance, slope1_slope2_covariance=slope1_slope2_covariance, sigma_res=sigma_res)
print(lme_fit4, cor=FALSE)
data.frame(intercept_mean=intercept_mean, slope1_mean=slope1_mean, slope2_mean=slope2_mean, intercept_sd=intercept_sd, slope1_sd=slope1_sd, slope2_sd=slope2_sd, intercept_slope1_covariance=intercept_slope1_covariance, intercept_slope2_covariance=intercept_slope2_covariance, slope1_slope2_covariance=slope1_slope2_covariance, sigma_res=sigma_res)

# Note the large sd.s for the measures of association
print(out$mean$rho, dig=2)
print(out$sd$rho, dig=2)

# Finally, we compare the prior and posterior distributions for the (derived) ranef distribution parameters:

par(mfrow=c(3, 6))
hist(out$sims.list$mu_int_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_int_prior")
lines(density(out$sims.list$mu_int_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$mu_int, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_int")
lines(density(out$sims.list$mu_int), col="lightsteelblue3", lwd=2)

hist(out$sims.list$mu_slope1_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_slope1_prior")
lines(density(out$sims.list$mu_slope1_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$mu_slope1, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_slope1")
lines(density(out$sims.list$mu_slope1), col="lightsteelblue3", lwd=2)

hist(out$sims.list$mu_slope2_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_slope2_prior")
lines(density(out$sims.list$mu_slope2_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$mu_slope2, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="mu_slope2")
lines(density(out$sims.list$mu_slope2), col="lightsteelblue3", lwd=2)

hist(out$sims.list$sigma_int_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_int_prior")
lines(density(out$sims.list$sigma_int_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$sigma_int, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_int")
lines(density(out$sims.list$sigma_int), col="lightsteelblue3", lwd=2)

hist(out$sims.list$sigma_slope1_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_slope1_prior")
lines(density(out$sims.list$sigma_slope1_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$sigma_slope1, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_slope1")
lines(density(out$sims.list$sigma_slope1), col="lightsteelblue3", lwd=2)

hist(out$sims.list$sigma_slope2_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_slope2_prior")
lines(density(out$sims.list$sigma_slope2_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$sigma_slope2, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="sigma_slope2")
lines(density(out$sims.list$sigma_slope2), col="lightsteelblue3", lwd=2)

hist(out$sims.list$rho_int_slope1_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_int_slope1_prior")
lines(density(out$sims.list$rho_int_slope1_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$rho_int_slope1, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_int_slope1")
lines(density(out$sims.list$rho_int_slope1), col="lightsteelblue3", lwd=2)

hist(out$sims.list$rho_int_slope2_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_int_slope2_prior")
lines(density(out$sims.list$rho_int_slope2_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$rho_int_slope2, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_int_slope2")
lines(density(out$sims.list$rho_int_slope2), col="lightsteelblue3", lwd=2)

hist(out$sims.list$rho_slope1_slope2_prior, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_slope1_slope2_prior")
lines(density(out$sims.list$rho_slope1_slope2_prior), col="lightsteelblue3", lwd=2)
hist(out$sims.list$rho_slope1_slope2, col="lightsteelblue1", border="white", breaks=10, freq=FALSE, main="rho_slope1_slope2")
lines(density(out$sims.list$rho_slope1_slope2), col="lightsteelblue3", lwd=2)
par(mfrow=c(1, 1))
