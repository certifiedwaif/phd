times <- c(45.36, 30.18, 25.57, 21.47, 21.38, 20.35, 20.19, 19.37)
cpu <- c(99, 173, 245, 313, 324, 348, 355, 370)

pdf("threads_versus_time.pdf")
plot(1:8, times, xlab = "Threads", ylab = "Time in seconds", type = "l")
dev.off()
pdf("threads_versus_utilisation.pdf")
plot(1:8, cpu, xlab = "Threads", ylab = "CPU utilisation", type = "l")
dev.off()
