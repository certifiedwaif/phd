# plot_local_solutions.R
library(optparse)

display_error_locations <- function(error_locations)
{
	plot(error_locations[, 1], error_locations[, 2], xlab=expression(beta[0]), ylab=expression(beta[1]))
}

main <- function()
{
  option_list <- list(make_option(c("-i", "--input")),
  										make_option(c("-o", "--output")))
  opt <- parse_args(OptionParser(option_list=option_list))
  input_file <- opt$input
  error_locations <- read.csv(file=input_file)
  pdf(opt$output)
  display_error_locations(error_locations[, 2:3])
  dev.off()
}
main()
