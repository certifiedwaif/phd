library(Matrix)

plot_matrix <- function(p, q, m, chol=FALSE)
{
	mat <- matrix(0, p + q * m, p + q * m)
	mat[1:p, ] <- 1/8
	mat[, 1:p] <- 1/8
	mat[1:p, 1:p] <- 1/4
	idx <- p + 1
	for (i in 1:q) {
		# inds <- (p + b * (m-1) + 1):((p + b * m))
		# cat("idx", idx, "\n")
		for (j in 1:m) {
			for (k in 1:m) {
				mat[idx+j-1, idx+k-1] <- 1/4
			}
		}
		idx <- idx + m
	}
	diag(mat) <- 1
	if (chol)
		mat <- t(chol(mat))
	image(Matrix(mat))
}

plot_matrix2 <- function(p, q, m, chol=FALSE)
{
	mat <- matrix(0, p + q * m, p + q * m)
	mat[q * m + 1:p, ] <- 1/8
	mat[, q * m + 1:p] <- 1/8
	mat[q * m + 1:p, q * m + 1:p] <- 1/4
	idx <- 1
	for (i in 1:q) {
		# inds <- (p + b * (m-1) + 1):((p + b * m))
		# cat("idx", idx, "\n")
		for (j in 1:m) {
			for (k in 1:m) {
				mat[idx+j-1, idx+k-1] <- 1/4
			}
		}
		idx <- idx + m
	}
	diag(mat) <- 1
	if (chol)
		mat <- t(chol(mat))
	image(Matrix(mat))
}

p <- 2
q <- 10
m <- 2
pdf("mX_mZ_mLambda.pdf")
plot_matrix(p, q, m)
dev.off()
pdf("mX_mZ_cholesky.pdf")
plot_matrix(p, q, m, chol=TRUE)
dev.off()
pdf("mZ_mX_mLambda.pdf")
plot_matrix2(p, q, m)
dev.off()
pdf("mZ_mX_cholesky.pdf")
plot_matrix2(p, q, m, chol=TRUE)
dev.off()

# PLot zero-inflated distribution
x <- c(0,7,3,4,5,3,2,6,5,0,0,1,0,0,5,0,2,3,6,4,0,5,4,0,
			 7,0,0,0,7,0,6,6,0,3,0,5,0,4,0,0,0,2,3,0,3,4,5,0,
			 8,0)
pdf("~/Dropbox/phd/poisson_does_not_fit.pdf", width=3, height=3)
hist(x, , prob=TRUE)
lambda_hat <- mean(x)
points(0:8 + .5, dpois(0:8, lambda_hat), col = "red", cex = 1.25, pch = 19)
dev.off()
