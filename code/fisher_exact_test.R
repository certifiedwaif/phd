# fisher_exact_test.R
pvalue <- function(n) {
	p_non_effect <- 0.05
	p_effect <- 0.1
	p_group <- .5
	mat <- matrix(c(n * p_non_effect * p_group, n * (1 - p_non_effect)  * (1 - p_group),
									n * p_effect * p_group, n * (1 - p_effect) * p_group), ncol=2)
	fisher.test(mat)$p.value
}

for (n in 1:1000) {
	cat(n, pvalue(n), "\n")
}
