plot_vw <- function(vw1) {
    vw1_mat <- as.matrix(vw1)
    n <- ncol(vw1_mat)
    r <- hist(vw1_mat, breaks = n, axes = FALSE, prob = TRUE, main = "", xlab = "")
    axis(1, r$mid, as.character(1:length(r$mid)))
}

vw1 <- read.csv("Hitters_vw1.csv", header = FALSE)
vw2 <- read.csv("Hitters_vw2.csv", header = FALSE)

pdf("Hitters_histogram.pdf")
par(mfrow = c(2, 1))
plot_vw(vw1)
plot_vw(vw2)
dev.off()

vw1 <- read.csv("bodyfat_vw1.csv", header = FALSE)
vw2 <- read.csv("bodyfat_vw2.csv", header = FALSE)

pdf("bodyfat_histogram.pdf")
par(mfrow = c(2, 1))
plot_vw(vw1)
plot_vw(vw2)
dev.off()

vw1 <- read.csv("GradRate_vw1.csv", header = FALSE)
vw2 <- read.csv("GradRate_vw2.csv", header = FALSE)

pdf("GradRate_histogram.pdf")
par(mfrow = c(2, 1))
plot_vw(vw1)
plot_vw(vw2)
dev.off()

vw1 <- read.csv("USCrime_vw1.csv", header = FALSE)
vw2 <- read.csv("USCrime_vw2.csv", header = FALSE)

pdf("USCrime_histogram.pdf")
par(mfrow = c(2, 1))
plot_vw(vw1)
plot_vw(vw2)
dev.off()

vw1 <- read.csv("Wage_vw1.csv", header = FALSE)
vw2 <- read.csv("Wage_vw2.csv", header = FALSE)

pdf("Wage_histogram.pdf")
par(mfrow = c(2, 1))
plot_vw(vw1)
plot_vw(vw2)
dev.off()
