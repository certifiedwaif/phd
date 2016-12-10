# plot_accuracy.R
accuracy <- read.csv(file="accuracy.csv")
pdf("results/accuracy_of_parameter_estimation.pdf")
plot(accuracy[,2], ylim=c(0, 1), type="l", ylab="Accuracy", xlab="Parameter", xaxt="n")
xaxis_labels <- c(expression(beta[0]),
									expression(beta[1]),
									expression(beta[2]),
									expression(u[1]),
									expression(u[2]),
									expression(u[3]),
									expression(u[4]),
									expression(u[5]),
									expression(u[6]),
									expression(u[7]),
									expression(u[8]),
									expression(u[9]),
									expression(u[10]),
									expression(u[11]),
									expression(u[12]),
									expression(sigma[u^2]),
									expression(rho))
axis(1, at=1:length(xaxis_labels), labels = xaxis_labels)
lines(accuracy[,3], ylim=c(0, 1), type="l", col="red")
lines(accuracy[,4], ylim=c(0, 1), type="l", col="blue")
# lines(accuracy[,5], ylim=c(0, 1), type="l", col="blue")
legend("bottomright",
				c("Laplace", "GVA covariance", "GVA precision"),
				fill=c("black", "red", "blue"))
dev.off()
