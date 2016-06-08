# week2.R
x <- c(9.82, 9.85, 9.87, 9.90, 9.91, 9.94, 9.96, 10.02)
y <- c(9.79, 9.75, 9.63, 9.80, 9.94, 11.99, 9.88, 9.98)
race <- data.frame(x=x, y=y)
race

plot(x, y)
plot(x[-6], y[-6])
fit <- lm(y~x, data=race)
summary(fit)

79.27 / 8 * 3.571 + -25.287
[1] 10.09715
80.76 / 8
[1] 10.095
sum(resid(fit)^2)
[1] 3.821477
plot(fit)
