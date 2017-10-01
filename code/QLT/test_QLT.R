# test_QLT.R

library(tidyverse)
library(optparse)
library(latex2exp)
source("QLT.R")


generate_F_1_data <- function(K, data_fn, start, prior, bUnique)
{
  if (data_fn == "QLT") {
    dat <- QLT(K, generate_data_QLT, start, prior, bUnique=bUnique)
  } else if (data_fn == "high_dimensional") {
    dat <- QLT(K, generate_data_high_dimensional, start, prior, bUnique=bUnique)
  } else {
  	stop("data_fn unknown")
  }
  # save(dat, file = sprintf("results/%s_%s_%s_%s.dat", K, data_fn, start, prior))
  return(dat)
}


calc_var_probs <- function(vy, mX, mGamma, vR2)
{
  library(gsl)
  calcVarProbs <- function(score,mGamma) {
    model.prob = exp(score - max(score))/sum(exp(score - max(score)))
    var.prob = t(mGamma)%*%model.prob
    return(var.prob)
  }
  
  n <- length(vy)
  p <- ncol(mX)
  M = length(vR2)
  vq = mGamma%*%matrix(1,p,1)
  
  vBIC = n*log(1 - vR2) + vq*log(n) 
  var.prob1 = calcVarProbs(-0.5*vBIC,mGamma) 
  
  a  = -0.75
  vb = 0.5*(n - vq - 5) - a
  c  = 0.5*(n - 1)
  vd = 0.5*vq + a
  
  log.vp = -(vb+1)*log(1 - vR2) + lbeta(vd+1,vb+1) - lbeta(a+1,vb+1)
  vZE <- -2*log.vp
  
  var.prob2 = calcVarProbs(log.vp,mGamma) 
  
  log.hyperg_2F1 = function(b,c,x) {
    val = 0 
    val = val + log(c-1) 
    val = val + (1-c)*log(x) 
    val = val + (c-b-1)*log(1-x) 
    val = val + lbeta(c-1,b-c+1) 
    val = val + pbeta(x,shape1=(c-1),shape2=(b-c+1),log=TRUE)
    return(val)
  }
  
  log.hyperg_2F1.naive = function(b,c,x) {
    val = log( hyperg_2F1( b, 1, c, x, give=FALSE, strict=TRUE) )
    return(val)
  } 
  
  a = 3
  log.vp.g = log(a - 2) - log(vq + a - 2) + log( hyperg_2F1(0.5*(n-1), 1, 0.5*(vq+a), vR2, give=FALSE, strict=TRUE) )
  var.prob3 = calcVarProbs(log.vp.g,mGamma) 
  
  a = 3
  log.vp.g2 = log(a - 2) - log(vq + a - 2) + log.hyperg_2F1( 0.5*(n-1), 0.5*(vq+a), vR2 )
  var.prob4 = calcVarProbs(log.vp.g2,mGamma) 
  
  log.vp.gprior5 = rep(0,M)
  log.vp.gprior5[-1] = log.vp.gprior5[-1] - 0.5*vq[-1]*log(n+1)
  log.vp.gprior5[-1] = log.vp.gprior5[-1] + 0.5*vq[-1]*log(vq[-1]+1)
  log.vp.gprior5[-1] = log.vp.gprior5[-1] - 0.5*(n - 1)*log(vR2[-1])
  log.vp.gprior5[-1] = log.vp.gprior5[-1] - log(vq[-1]+1)
  log.vp.gprior5[-1] = log.vp.gprior5[-1] + log( hyperg_2F1( 0.5*(vq[-1]+1), 0.5*(n-1), 0.5*(vq[-1]+3), (1-1/vR2[-1])*(vq[-1]+1)/(n+1), give=FALSE, strict=TRUE) )
  
  model.prob5 = exp(log.vp.gprior5 - max(log.vp.gprior5))/sum(exp(log.vp.gprior5 - max(log.vp.gprior5)))
  var.prob5 = t(mGamma)%*%model.prob5
  
  vL = (1 + n)/(1 + vq) - 1
  vsigma2 = 1 - vR2
  vz = vR2/(1 + vL*vsigma2)
  
  log.vp.gprior6 = rep(0,M)
  log.vp.gprior6[-1] = log.vp.gprior6[-1] + 0.5*(n - vq[-1] - 1)*log( n + 1 )
  log.vp.gprior6[-1] = log.vp.gprior6[-1] - 0.5*(n - vq[-1] - 1)*log( vq[-1] + 1)
  log.vp.gprior6[-1] = log.vp.gprior6[-1] - 0.5*(n - 1)*log(1 + vL[-1]*vsigma2[-1])
  log.vp.gprior6[-1] = log.vp.gprior6[-1] - log (vq[-1] + 1)
  log.vp.gprior6[-1] = log.vp.gprior6[-1] + log.hyperg_2F1( 0.5*(n-1), 0.5*(vq[-1]+3), vz[-1] )
  
  model.prob6 = exp(log.vp.gprior6 - max(log.vp.gprior6))/sum(exp(log.vp.gprior6 - max(log.vp.gprior6)))
  var.prob6 = t(mGamma)%*%model.prob6
  
  log.vp.gprior7 = rep(0,M)
  log.vp.gprior7[-1] = log.vp.gprior7[-1] + 0.5*(n - vq[-1] - 1)*log( n + 1 )
  log.vp.gprior7[-1] = log.vp.gprior7[-1] - 0.5*(n - vq[-1] - 1)*log( vq[-1] + 1)
  log.vp.gprior7[-1] = log.vp.gprior7[-1] - 0.5*(n - 1)*log(1 + vL[-1]*vsigma2[-1])
  log.vp.gprior7[-1] = log.vp.gprior7[-1] - log (vq[-1] + 1)
  log.vp.gprior7[-1] = log.vp.gprior7[-1] + log( hyperg_2F1( 0.5*(n-1), 1, 0.5*(vq[-1] + 3), vz[-1], give=FALSE, strict=TRUE) )
  
  model.prob7 = exp(log.vp.gprior7 - max(log.vp.gprior7))/sum(exp(log.vp.gprior7 - max(log.vp.gprior7)))
  var.prob7 = t(mGamma)%*%model.prob7
  
  tab4 <- cbind(var.prob1,var.prob2,var.prob3,var.prob4,var.prob5,var.prob6,var.prob7)
  colnames(tab4) = c("BIC","ZE","g.naive","g.safe","Robust.naive","Robust","Robust.safe")
  return(tab4)
}


variable_inclusion_tbl <- function(K, data_fn, start, prior, bUnique=TRUE)
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
  n <- length(vy)

  initial_gamma <- initialise_gamma(start, K, p, vy, mX, models=NULL)
  library(correlation)
  cva.res <- cva(initial_gamma, vy, mX, K, 1, prior, bUnique)

	vR21 = correlation::all_correlations_mX(vy, mX)$vR2
	mGamma1 = graycode(p,0)
	tab1 <- calc_var_probs(vy, mX, mGamma1, vR21)
	tab_to_tbl <- function(tab, mX)
	{
    # Rows are covariates, columns are priors
  	varnames <- colnames(mX)
  	priors <- colnames(tab)
  	# I'm sure there's some clever tidyverse way of doing this
    var <- rep(NA, length(varnames) * length(priors))
    prior <- rep(NA, length(varnames) * length(priors))
    prob <- rep(NA, length(varnames) * length(priors))
    counter <- 1
    for (row in 1:nrow(tab)) {
  	  for (col in 1:ncol(tab)) {
  	    var[counter] <- varnames[row]
  	    prior[counter] <- priors[col]
  	    prob[counter] <- tab[row, col]
  	    counter <- counter + 1
  	  }
  	}
  	tibble(covariate=var, var_prior=prior, inclusion_probability=prob)
	}
	tbl1 <- tab_to_tbl(tab1, mX)
	
	mGamma2 = cva.res$models
	vR22 = apply(mGamma2, 1, function(vGamma) {
	  mX <- mX[, which(vGamma == 1)]
	  (t(vy) %*% mX %*% solve(t(mX) %*% mX) %*% t(mX) %*% vy) / n
	})
	tab2 <- calc_var_probs(vy, mX, mGamma2, vR22)
	tbl2 <- tab_to_tbl(tab2, mX)
	tbl <- rbind(tbl1 %>% mutate(source="Exact"),
	             tbl2 %>% mutate(source="CVA"))
	
	return(tbl %>% mutate(K=K, start=start, cva_prior=prior, data_set=data_fn))
}


variable_inclusion_graph <- function(tbl)
{
	ggplot(tbl, aes(x=covariate, y=inclusion_probability)) + facet_grid(var_prior~source) + geom_col() +
		xlab("Covariate") + ylab("Inclusion probability") + ggtitle(sprintf("%s, K=%s", tbl$data_set, tbl$K))
}


variable_inclusion_graphs <- function()
{
  pdf("/tmp/variable_inclusion.pdf")
  
  a1 <- proc.time()[3]
  variable_inclusion(20, "Hitters", "cold_start", "log_prob1", bUnique=TRUE) %>% variable_inclusion_graph
  b1 <- proc.time()[3]
  print(b1 - a1)
  
  a2 <- proc.time()[3]
  variable_inclusion(20, "Hitters", "cold_start", "log_prob1", bUnique=FALSE) %>% variable_inclusion_graph
  b2 <- proc.time()[3]
  print(b2 - a2)

  expand.grid(K = c(20, 50, 100),
  						data_set = c("Hitters", "Bodyfat", "Wage", "College", "USCrime")) %>%
  						as_tibble %>%
  						map2(K, data_set,
  									~variable_inclusion_tbl(.x, .y, "cold_start", "log_prob1")) %>%
  						variable_inclusion_graph
  variable_inclusion_tbl(20, "Hitters", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(20, "Bodyfat", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(20, "Wage", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(20, "College", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(20, "USCrime", "cold_start", "log_prob1") %>% variable_inclusion_graph
  
  variable_inclusion_tbl(50, "Hitters", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Bodyfat", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Wage", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "College", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "USCrime", "cold_start", "log_prob1") %>% variable_inclusion_graph
  
  variable_inclusion_tbl(100, "Hitters", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(100, "Bodyfat", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(100, "Wage", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(100, "College", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(100, "USCrime", "cold_start", "log_prob1") %>% variable_inclusion_graph
  dev.off()
  
  pdf("/tmp/variable_inclusion_model_selection.pdf")
  variable_inclusion_tbl(50, "Hitters", "cold_start", "log_prob1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "BIC") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "ZE") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "liang_g1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "liang_g2") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "liang_g3") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "robust_bayarri1") %>% variable_inclusion_graph
  variable_inclusion_tbl(50, "Hitters", "cold_start", "robust_bayarri2") %>% variable_inclusion_graph
  dev.off()
}


relative_error_tbl <- function(data_set)
{
  tbl <- bind_rows(map(.x = c(20, 50, 100),
  											.f = ~ variable_inclusion_tbl(.x, data_set, "cold_start", "log_prob1")))
  exact_tbl <- tbl %>% 
  							filter(source == "Exact") %>%
  							rename(exact_prob = inclusion_probability) %>%
  							dplyr::select(-source)
  cva_tbl <- tbl %>%
  						filter(source == "CVA") %>%
  						rename(est_prob = inclusion_probability) %>%
  						dplyr::select(-source)
  rel_error_tbl <- exact_tbl %>%
								  	full_join(cva_tbl) %>%
								  	mutate(relative_error = round(abs((exact_prob - est_prob)) / exact_prob, 2)) %>%
								  	dplyr::select(data_set, var_prior, covariate, K, exact_prob, relative_error) %>%
								  	arrange(data_set, var_prior, covariate, K)
	
	return(rel_error_tbl)
}


relative_error_tables <- function()
{
	hitters_re_tbl <- relative_error_tbl("Hitters")
	bodyfat_re_tbl <- relative_error_tbl("Bodyfat")
	wage_re_tbl <- relative_error_tbl("Wage")
	college_re_tbl <- relative_error_tbl("College")
	uscrime_re_tbl <- relative_error_tbl("USCrime")

	bind_rows(hitters_re_tbl,
						bodyfat_re_tbl,
						wage_re_tbl,
						college_re_tbl,
						uscrime_re_tbl) %>%
				mutate(greater_than_half = ifelse(exact_prob >= 0.5, ">0.5", "<=0.5")) %>%
				filter(var_prior != "Robust.naive" & var_prior != "g.naive" & var_prior != "Robust.safe") %>%
				mutate(K=sprintf("%03d", K)) %>%
				group_by(data_set, K, var_prior, greater_than_half) %>%
				summarise(relative_error = round(mean(relative_error), 2)) %>%
				unite(greater_than_half_K, greater_than_half, K) %>%
		 		spread(greater_than_half_K, relative_error) %>%
		 		write_csv("/tmp/summarised_relative_error_table.csv")
}


relative_error_plot <- function(tbl)
{
	 tbl %>% ggplot(aes(x=K, y=relative_error, col=var_prior)) +
	 					ggtitle(tbl$data_set) + geom_line() + facet_grid(~covariate) + ylab("Relative error")
}


relative_error_plots <- function()
{
	pdf("relative_error_plots.pdf")
	relative_error_tbl("Hitters") %>% relative_error_plot
	relative_error_tbl("Bodyfat") %>% relative_error_plot
	relative_error_tbl("Wage") %>% relative_error_plot
	relative_error_tbl("College") %>% relative_error_plot
	relative_error_tbl("USCrime") %>% relative_error_plot
	dev.off()
}


#' Build a database of the F_1 scores produced by various experiments
build_f_1_db <- function()
{
	library(stringr)
	path <- "~/Dropbox/phd/code/QLT/results/"

	# Load CVA F_1 data
	orig_fnames = list.files(path = path, pattern = "*.dat$")
	fnames <- orig_fnames %>%
							str_replace("generate_data_", "") %>%
							str_replace("high_dimensional", "highdimensional") %>%
							str_replace("warm_start_", "") %>%
							str_replace("cold_start", "cold") %>%
							str_replace("log_prob1", "logprob1")
	pattern <- "^(.*)_(.*)_(.*)_(.*).dat$"
	db <- str_match(fnames, pattern)
	cva_tbl <- tibble(fname = orig_fnames, k = as.numeric(db[, 2]), data = db[, 3], start = db[, 4], prior = db[, 5], f_1 = NA)
	cva_tbl <- cva_tbl %>% filter(!is.na(k)) # Remove comparison/non-matching data sets
	for (i in 1:nrow(cva_tbl)) {
		load(str_c(path, cva_tbl$fname[i]))
		cva_tbl$f_1[i] <- list(dat)
	}
	cva_tbl$f_1 <- map(cva_tbl$f_1, function(x) {
	  if (is.vector(x)) {
	    x
	  } else if (is.matrix(x)) {
	    x[, ncol(x)]
	  } else {
	    x
	  }
	})
	
	# Load acomparison data
	load("~/Dropbox/phd/code/QLT/results/QLT_comparison.dat")
	QLT_dat <- dat[, 1:5]
	loaded <- load("~/Dropbox/phd/code/QLT/results/high_dimensional_comparison.dat")
	high_dimensional_dat <- dat[, 1:5]
	columns <- c("Lasso", "Scad", "Mcp", "Bms", "Varbvs")
	colnames(QLT_dat) <- columns
	QLT_list <- map(.x=columns, ~QLT_dat[,.x])
	colnames(high_dimensional_dat) <- columns
	high_dim_list <- map(.x=columns, ~high_dimensional_dat[,.x])
	QLT_tbl <- tibble(data="QLT", method=columns, f_1=QLT_list)
	high_dim_tbl <- tibble(data="highdimensional", method=columns, f_1=high_dim_list)
	comparison_tbl <- QLT_tbl %>% bind_rows(high_dim_tbl)

	# Some CVA data got mixed in with this, so seperate that out and combine with the rest of the CVA data
	#cva_k1000_tbl <- comparison_tbl %>% filter(method == "CVA_K1000") %>% dplyr::select(-method)
	#comparison_no_cva_tbl <- comparison_tbl %>% filter(method != "CVA_K1000")
	#cva_k1000_tbl <- cva_k1000_tbl %>%
	#									mutate(k=1000, 
	#												 fname=str_c(data, "_comparison.dat"),
	#												 start="cold",
	#												 prior="BIC")
	#cva_tbl <- cva_tbl %>% bind_rows(cva_k1000_tbl)

	return(list(comparison_tbl=comparison_tbl, cva_tbl=cva_tbl))
}


generate_boxplot <- function(data_set, prior_, start_, cva_tbl, comparison_tbl)
{
	# Comparison F_1 scores are always the same
	# Compare against CVA for different values of K, start, prior
  cva_plot_tbl <- cva_tbl %>%
										filter(data == data_set & prior == prior_ & start == start_) %>%
										mutate(method = sprintf("CVA: K=%04d", k))
  # I don't know how this happened, but some of the f_1 data files contain 50 by c
  # matrices. We only want the last column.
  cva_plot_tbl$f_1 <- map(cva_plot_tbl$f_1, function(x) {
    if (is.vector(x)) {
      x
    } else {
      x[, 4]
    }
  })
  method_levels <- c(unique(comparison_tbl$method), sort(unique(cva_plot_tbl$method)))
  bind_rows(comparison_tbl %>% filter(data == data_set) %>% select(method, f_1),
						cva_plot_tbl %>% select(method, f_1)) %>%
    unnest %>%
		ggplot(aes(x=factor(method, levels=method_levels, ordered=TRUE), y=f_1)) + 
    geom_boxplot() + ylim(low=0, high=1) + xlab("Method") + ylab(TeX("F_1"))
}


generate_boxplots <- function()
{
  dbs <- build_f_1_db()
  
  setwd("~/Dropbox/phd/code/QLT")
  pdf("results/QTL_cold_maruyama.pdf")
  generate_boxplot("QLT", "maruyama", "cold", dbs$cva_tbl, dbs$comparison_tbl)
  dev.off()
  
  pdf("results/QTL_likelihood_maruyama.pdf")
  generate_boxplot("QLT", "maruyama", "likelihood", dbs$cva_tbl, dbs$comparison_tbl)
  dev.off()
  
  pdf("results/QTL_covariates_maruyama.pdf")
  generate_boxplot("QLT", "maruyama", "covariates", dbs$cva_tbl, dbs$comparison_tbl)
  dev.off()
  
  pdf("results/highdimensional_cold_logprob1.pdf")
  generate_boxplot("highdimensional", "logprob1", "cold", dbs$cva_tbl, dbs$comparison_tbl)
  dev.off()
  
  pdf("results/highdimensional_likelihood_logprob1.pdf")
  generate_boxplot("highdimensional", "logprob1", "likelihood", dbs$cva_tbl, dbs$comparison_tbl)
  dev.off()
  
  pdf("results/highdimensional_covariates_logprob1.pdf")
  generate_boxplot("highdimensional", "logprob1", "covariates", dbs$cva_tbl, dbs$comparison_tbl)
  dev.off()
}

# # Base
 #  pdf(sprintf("results/%s_%s_%s_%s.pdf", K, data_fn, start, prior))
 #  boxplot(dat,names=c("lasso","scad","mcp",
 #                      # "emvs",
 #                      "bms",
 #                      #"varbvs",
 #                      "cva"),ylim=c(0,1))
 #  dev.off()

	# tbl <- tbl %>%
	# 				mutate(f_1_mean = map_dbl(tbl$f_1, mean),
	# 								f_1_median = map_dbl(tbl$f_1, median),
	# 								f_1_sd = map_dbl(tbl$f_1, sd))

 #  # Ggplot
	# library(tidyverse)
	# library(forcats)
	# dat2 <- data.frame(dat)
	# colnames(dat2) <- c("lasso", "scad", "mcp", "bms", "varbvs", "cva")
	# dat3 <- dat2 %>%
	#   gather(method, F_1, lasso:cva) %>%
	#   mutate(method = parse_factor(method, c("lasso", "scad", "mcp", "bms", "varbvs", "cva")))
	# ggplot(dat3, aes(x=method, y=F_1)) + geom_boxplot() + ggtitle("Something")
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
