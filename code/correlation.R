# correlation.R

# Anscombe's quartet
x1 <- c(10.0, 8.0, 13.0, 9.0, 11.0, 14.0, 6.0, 4.0, 12.0, 7.0, 5.0)
y1 <- c(8.04, 6.95, 7.58, 8.81, 8.33, 9.96, 7.24, 4.26, 10.84, 4.82, 5.68)
y2 <- c(9.14, 8.14, 8.74, 8.77, 9.26, 8.1, 6.13, 3.1, 9.13, 7.26, 4.74)
y3 <- c(7.46, 6.77, 12.74, 7.11, 7.81, 8.84, 6.08, 5.39, 8.15, 6.42, 5.73)
x2 <- c(rep(8.0, 7), 19, rep(8.0, 3))
y4 <- c(6.58, 5.76, 7.71, 8.84, 8.47, 7.04, 5.25, 12.50, 5.56, 7.91, 6.89)

anscombe_df <- data.frame(x1=x1, y1=y1, y2=y2, y3=y3, x2=x2, y4=y4)
mat <- scale(anscombe_df)
anscombe_df <- data.frame(x1=mat[, 1], y1=mat[, 2], y2=mat[, 3], y3=mat[, 4], x2=mat[, 5], y4=mat[, 6])

par(mfrow=c(2, 2))
plot(anscombe_df$x1, anscombe_df$y1)
fit <- lm(anscombe_df$y1~anscombe_df$x1)
abline(fit)

plot(anscombe_df$x1, anscombe_df$y2)
fit <- lm(anscombe_df$y2~anscombe_df$x1)
abline(fit)

plot(anscombe_df$x1, anscombe_df$y3)
fit <- lm(anscombe_df$y3~anscombe_df$x1)
abline(fit)

plot(x2, anscombe_df$y4)
fit <- lm(anscombe_df$y4~x2)
abline(fit)


fit2 <- lm(x1~y2+y3-1)
summary(fit2)
vy <- anscombe_df$x1
mX <- as.matrix(anscombe_df[, c("y2", "y3")])
vbeta <- coef(fit2)
cor(vy, mX %*% vbeta)
cor(vy, mX %*% vbeta)^2

write.table(file="anscombes_quartet.csv", anscombe_df, row.names=FALSE, col.names=FALSE, quote=FALSE, sep=",")
