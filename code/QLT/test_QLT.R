# test_QLT.R

library(optparse)
source("QLT.R")


generate_F_1_data <- function(K, data_fn, start, prior)
{
  if (data_fn == "QLT") {
    dat <- QLT(K, generate_data_QLT, start, prior)
  } else if (data_fn == "high_dimensional") {
    dat <- QLT(K, generate_data_high_dimensional, start, prior)
  } else {
  	stop("data_fn unknown")
  }
  # save(dat, file = sprintf("results/%s_%s_%s_%s.dat", K, data_fn, start, prior))
  return(dat)
}


variable_inclusion <- function(K, data_fn, start, prior)
{
  if (data_fn == "Hitters") {
    dat <- generate_data_Hitters()
  } else if (data_fn == "Bodyfat") {
    dat <- generate_data_Bodyfat()
  } else if (data_fn == "Wage") {
    dat <- generate_data_Wage()
  } else if (data_fn == "College") {
    dat <- generate_data_College()
  } else if (data_fn == "USCrime") {
    dat <- generate_data_USCrime()
  }

  p <- dat$p
  vy <- dat$vy
  mX <- dat$mX

  initial_gamma <- initialise_gamma(start, K, p, vy, mX, models=NULL)
  library(correlation)
  cva.res <- cva(initial_gamma, vy, mX, K, 1, prior)
  variable_inclusion <- apply(cva.res$models, 2, sum) / K

  return(variable_inclusion)
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


build_db <- function()
{
	library(tidyverse)

	orig_fnames = list.files(pattern = "*.dat$")
	fnames <- orig_fnames %>%
							str_replace("generate_data_", "") %>%
							str_replace("high_dimensional", "highdimensional") %>%
							str_replace("warm_start_", "") %>%
							str_replace("cold_start", "cold") %>%
							str_replace("log_prob1", "logprob1")
	pattern <- "^(.*)_(.*)_(.*)_(.*).dat$"
	db <- str_match(fnames, pattern)
	tbl <- tibble(fname = orig_fnames, k = as.numeric(db[, 2]), data = db[, 3], start = db[, 4], prior = db[, 5], f_1 = NA)

	path <- "~/Dropbox/phd/code/QLT/results/"
	for (i in 1:nrow(tbl)) {
		load(str_c(path, tbl$fname[i]))
		tbl$f_1[i] <- list(dat)
	}
	tbl <- tbl %>%
					mutate(f_1_mean = map_dbl(tbl$f_1, mean),
									f_1_median = map_dbl(tbl$f_1, median),
									f_1_sd = map_dbl(tbl$f_1, sd))

	return(tbl)
}

# tbl %>% filter(f_1_median < .5) %>% group_by(k, prior) %>% arrange(k, prior) %>% summarise(count = n())
# tbl %>% filter(f_1_median < .5) %>% group_by(data, prior) %>% arrange(k, prior) %>% summarise(count = n())
# tbl %>% filter(f_1_median >= .5) %>% group_by(k, prior) %>% arrange(k, prior) %>% summarise(count = n())
# tbl %>% filter(f_1_median >= .5) %>% group_by(data, prior) %>% arrange(k, prior) %>% summarise(count = n())

main <- function()
{
  option_list <- list(make_option("--k", type="integer"),
  										make_option("--data"),
  										make_option("--prior"),
  										make_option("--start"))
  opt <- parse_args(OptionParser(option_list=option_list))
  print(str(opt))
  generate_F_1_data(opt$k, opt$data, opt$start, opt$prior)
}
main()
