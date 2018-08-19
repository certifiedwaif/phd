# plot_local_solutions.R
library(optparse)
library(tidyverse)
library(latex2exp)

plot_error_locations <- function(error_locations) {
    # plot(error_locations[, 1], error_locations[, 2], xlab = expression(beta[0]),
    # ylab = expression(beta[1]))
    ggplot(error_locations, aes(x = x, y = y)) + geom_point() + xlab(TeX("$\\beta_1$")) + 
        ylab(TeX("$\\beta_2$"))
}

main <- function() {
    option_list <- list(make_option(c("-i", "--input")), make_option(c("-o", "--output")))
    opt <- parse_args(OptionParser(option_list = option_list))
    input_file <- opt$input
    error_locations <- read_csv(file = input_file)
    pdf(opt$output)
    p <- plot_error_locations(error_locations)
    print(p)
    dev.off()
}
main()
