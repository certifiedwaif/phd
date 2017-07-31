# test_QLT.R

source("QLT.R")

K <- c(1, 10, 20, 50, 100)
data_fn <- c("generate_data_QTL", "generate_data_high_dimensional")
start <- c("cold_start", "warm_start_covariates", "warm_start_likelihood")
prior <- c("BIC", "ZE", "3", "4", "5", "6", "7")
grid_df <- expand.grid(K, data_fn, start, prior)
colnames(grid_df) <- c("K", "data_fn", "start", "prior")

for (i in 1:nrow(grid_df)) {
	K <- grid_df[i, ]$K
	data_fn <- grid_df[i, ]$data_fn
	start <- grid_df[i, ]$start
	prior <- grid_df[i, ]$prior

	if (data_fn == "generate_data_QTL") {
		dat <- QLT(K, generate_data_QTL, start, prior)
	} else {
		dat <- QLT(K, generate_data_high_dimensional, start, prior)
	}
	pdf(sprintf("%s_%s_%s_%s.pdf", K, data_fn, start, prior))
	boxplot(dat,names=c("lasso","scad","mcp",
	# "emvs",
	"bms","varbvs","cva"),ylim=c(0,1))
	dev.off()
}
