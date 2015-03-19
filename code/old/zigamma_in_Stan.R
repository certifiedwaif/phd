######################################
# zero-inflated gamma
library(MASS)
library(rstan)

# just my inv_logit function
logistic <- function (x) 
{
    p <- exp(x)/(1 + exp(x))
    p <- ifelse(x == Inf, 1, p)
    p
}

# simulate data
N <- 1000 # number of cases
J <- 20 # number of clusters
NperJ <- N/J # cases per cluster
theta <- 1 # scale
mu <- c( 0 , log(7) ) # means of varying intercepts
id <- rep( 1:J , each=NperJ )
Sigma <- matrix( 0 , nrow=2 , ncol=2 ) # var-cov matrix for varying intercepts
Sigma[1,1] <- 0.15
Sigma[2,2] <- 0.5
Sigma[1,2] <- Sigma[2,1] <- -0.8 * sqrt( Sigma[1,1] * Sigma[2,2] )
alpha <- mvrnorm( J , mu=mu , Sigma=Sigma )
y1 <- rbinom( N , size=1 , prob=logistic( alpha[id,1] ) ) # sim zeros
y2 <- (1-y1) * rgamma( N , shape=exp(alpha[id,2])*theta , scale=theta ) # sim observed

# prep data for Stan
dat <- list(
    y = ifelse( y2==0 , 20 , y2 ), # 20 to prevent gamma density gacking in vectorized code
    iszero = y1,
    N = N,
    J = J,
    Omega = diag(2),
    id = id
)

model_zigamma <- '
    data {
        int<lower=0> N;                 // number of cases
        int<lower=1> J;                 // number of clusters
        real<lower=0> y[N];             // observed outcome
        int<lower=0,upper=1> iszero[N]; // indicates a zero outcome
        int<lower=1> id[N];             // cluster number for each case
        cov_matrix[2] Omega;            // diagonal prior
    }
    parameters {
        vector[2] mu;                   // means of varying effects
        cov_matrix[2] Sigma;            // var-cov matrix for varying effects
        vector[2] alpha[J];             // varying effects for each cluster
        real theta;                     // log scale
    }
    model {
        real pi;                        // probability of zero GLM
        real mugamma;                   // mean of gamma GLM
        Sigma ~ inv_wishart( 3 , Omega );
        mu[1] ~ normal( 0 , 100 );
        mu[2] ~ normal( 0 , 100 );
        theta ~ normal( 0 , 100 );
        for (j in 1:J) alpha[j] ~ multi_normal( mu , Sigma );
        for (n in 1:N) {
            pi <- inv_logit( alpha[ id[n] , 1 ] );
            mugamma <- exp( alpha[ id[n] , 2 ] );
            lp__ <- lp__ + if_else( iszero[n] , log(pi) , log1m(pi) + gamma_log( y[n] , mugamma*exp(theta) , exp(theta) ) );
        }
    }
    generated quantities {
        real dev;                       // deviance of each set of samples
        real pi;                        // temp for computing dev
        real mugamma;                   // temp for computing dev
        real rho;                       // correlation btw varying intercepts
        real sd[2];                     // sd of varying intercepts
        dev <- 0;                       // not sure need to init to zero
        for ( n in 1:N ) {
            pi <- inv_logit( alpha[ id[n] , 1 ] );
            mugamma <- exp( alpha[ id[n] , 2 ] );
            dev <- dev + (-2) * if_else( iszero[n] , log(pi) , log1m(pi) + gamma_log( y[n] , mugamma*exp(theta) , exp(theta) ) );
        }
        for( k in 1:2 ) sd[k] <- sqrt( Sigma[k,k] );
        rho <- Sigma[1,2] / sqrt( Sigma[1,1] * Sigma[2,2] );
    }
'

# optional init list
init_alpha = matrix( 0 , nrow=J , ncol=2 )
init_alpha[ , 2 ] <- rep( log(7) , J )
initlist <- list(list(
    theta = 1 ,
    mu = c( 0 , log(7) ) ,
    Sigma = diag(2) ,
    alpha = init_alpha
))

fit <- stan( model_code = model_zigamma , data = dat , iter = 6000 , warmup = 1000 , init="random" , chains = 1 )

rstan::traceplot( fit , c("mu","theta","sd","rho") )

print( fit , digits=2 , pars=c("mu","theta","sd","rho") , probs=c(0.025,0.5,0.975) )
