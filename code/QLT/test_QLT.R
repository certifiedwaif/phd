# test_QLT.R

library(optparse)
source("QLT.R")


generate_F_1_data <- function(K, data_fn, start, prior)
{
  if (data_fn == "generate_data_QTL") {
    dat <- QLT(K, generate_data_QTL, start, prior)
  } else {
    dat <- QLT(K, generate_data_high_dimensional, start, prior)
  }
  write(dat, file = sprintf("results/%s_%s_%s_%s.dat", K, data_fn, start, prior))
}


generate_boxplots <- function(dat)
{
	# Base
  pdf(sprintf("results/%s_%s_%s_%s.pdf", K, data_fn, start, prior))
  boxplot(dat,names=c("lasso","scad","mcp",
                      # "emvs",
                      "bms",
                      #"varbvs",
                      "cva"),ylim=c(0,1))
  dev.off()

  # Ggplot
	library(tidyverse)
	library(forcats)
	dat2 <- data.frame(dat)
	colnames(dat2) <- c("lasso", "scad", "mcp", "bms", "varbvs", "cva")
	dat3 <- dat2 %>%
	  gather(method, F_1, lasso:cva) %>%
	  mutate(method = parse_factor(method, c("lasso", "scad", "mcp", "bms", "varbvs", "cva")))
	ggplot(dat3, aes(x=method, y=F_1)) + geom_boxplot() + ggtitle("Something")
}


main <- function()
{
  option_list <- list(make_option("--k", type="integer"),
  										make_option("--data"),
  										make_option("--prior"),
  										make_option("--start"))
  opt <- parse_args(OptionParser(option_list=option_list))
  print(str(opt))
  generate_F_1_data(opt$k, opt$data, opt$prior, opt$start)
}
main()
