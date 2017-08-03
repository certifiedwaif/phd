# cmd_QLT.R

source("QLT.R")

args <- commandArgs(TRUE)
K <- as.integer(args[1])

allowable_data_fn <- c("generate_data_QTL", "generate_data_high_dimensional")
data_fn <- args[2]
if (!(data_fn %in% allowable_data_fn)) {
	stop(paste0("Second argument must be one of ", allowable_data_fn))
}

allowable_start <- c("warm_start_covariates", "warm_start_likelihood", "cold_start")
start <- args[3]
if (!(start %in% allowable_start)) {
	stop(paste0("Second argument must be one of ", allowable_start))
}

allowable_prior <- c("log_prob1", "BIC", "ZE", "3", "4", "5", "6", "7")
prior <- args[4]
if (!(prior %in% allowable_prior)) {
	stop(paste0("Second argument must be one of ", allowable_prior))
}

if (data_fn == "generate_data_QTL") {
	dat <- QLT(K, generate_data_QTL, start, prior)
} else {
	dat <- QLT(K, generate_data_high_dimensional, start, prior)
}

pdf(sprintf("results/%s_%s_%s_%s.pdf", K, data_fn, start, prior))
boxplot(dat,names=c("lasso","scad","mcp",
# "emvs",
"bms","varbvs","cva"),ylim=c(0,1))
dev.off()
