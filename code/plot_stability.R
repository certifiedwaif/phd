setwd("~/Dropbox/phd/code")
stability_intercept <- read.csv(file = "stability_intercept.csv")
stability_intercept <- stability_intercept[, 2:25]
pdf("stability_intercept.pdf")
par(mfrow = c(2, 2))
plot(stability_intercept[, 1], stability_intercept[, 2], type = "l", xlab = "Threshold", 
    ylab = "Accuracy", main = expression(paste("Accuracy of ", bold(beta[1]), " versus threshold")))
plot(stability_intercept[, 1], stability_intercept[, 3], type = "l", xlab = "Threshold", 
    ylab = "Accuracy", main = expression(paste("Accuracy of ", bold(u[1]), " versus threshold")))
plot(stability_intercept[, 1], stability_intercept[, 23], type = "l", xlab = "Threshold", 
    ylab = "Accuracy", main = expression(paste("Accuracy of ", bold(sigma[u[1]]^2), 
        " versus threshold")))
plot(stability_intercept[, 1], stability_intercept[, 24], type = "l", xlab = "Threshold", 
    ylab = "Accuracy", main = expression(paste("Accuracy of ", rho, " versus threshold")))
dev.off()
