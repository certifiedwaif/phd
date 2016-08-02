# stat3012_week3.R

H = c(42.8,63.5,37.5,39.5,45.5,38.5,43,22.5,37,23.5,33,58)
W = c(40,93.5,35.5,30,52,17,38.5,8.5,33,9.5,21,79)
L = c(37,49.5,34.5,36,43,28,37,20,33.5,30.5,38.5,47)
dat = data.frame(H,W,L)
dat
lm1 = lm(L ~ . , data=dat)
summary(lm1)
sum(resid(lm1)^2)

par(mfrow=c(2,2))
boxplot(lm1$residuals)
plot(lm1,which=c(1,2,5),add.smooth=FALSE)
round(lm.influence(lm1)$h,3)
round(cooks.distance(lm1),4)
lm2 = lm(L~1+W, data=dat)
summary(lm2)
par(mfrow=c(2,2))
boxplot(lm2$residuals)
plot(lm2,which=c(1,2,5),add.smooth=FALSE)
predict(lm2,data.frame(W=90),se.fit=T)
