# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# TUTORIAL & LAB  --  WEEK 12                                       #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# set working directory to where the data is stored
setwd("/Users/mueller/Dropbox/usyd/teaching/S3012/data")  


# Question 02
dat = read.table(file="coil.txt",header=TRUE)
dat$mach = factor(dat$mach)
attach(dat)
# (a)
tapply(Y,mach,summary)
# (b)
summary(aov(Y~mach))
summary(aov(Y~Error(mach)))
# (c)
M1 = aov(Y~mach)
par(mfrow=c(2,2))
for (i in 1:4){
  qqnorm(M1$residuals[mach==i])
  print(shapiro.test(M1$residuals[mach==i]))
  }
# (d)
r=10
MSA = 200.83
df.ssa = 3
MSE = 7.15
df.rss = 36
f.l90 = qf(1-0.1/2,df.ssa,df.rss)
f.u90 = qf(0.1/2,df.ssa,df.rss)
f.l98 = qf(1-0.02/2,df.ssa,df.rss)
f.u98 = qf(0.02/2,df.ssa,df.rss)
L90 = 1/r * ( MSA/MSE * ( 1/f.l90) - 1 )
U90 = 1/r * ( MSA/MSE * ( 1/f.u90) - 1 )
CI90=c(L90/(1+L90), U90/(1+U90))
CI90
L98 = 1/r * ( MSA/MSE * ( 1/f.l98) - 1 )
U98 = 1/r * ( MSA/MSE * ( 1/f.u98) - 1 )
CI98=c(L98/(1+L98), U98/(1+U98))
CI98
(CI98[2]-CI98[1]) / (CI90[2]-CI90[1])

# (e)
# Y is a response vector of length n
# A is a balanced factor of length n
my1wayAOVII = function(Y,A,alpha=0.05){
  N = length(Y)
  t = dim(table(A))
  r = N/t
  Syy = (N-1)*var(Y)
  SSA = sum(tapply(Y,A,sum)**2)/r - (sum(Y))**2/N
  MSA = SSA/(t-1)
  RSS = Syy-SSA
  MSE = RSS/(N-t)
  f   = MSA/MSE
  p   = 1-pf(f,t-1,N-t)
  mu.hat = mean(Y)
  MSA = SSA/(t-1)
  CI.mu = c( mu.hat - qt(1-alpha/2,t-1)*sqrt(MSA/N),
             mu.hat + qt(1-alpha/2,t-1)*sqrt(MSA/N) )
  AOVtable = rbind(c(t-1,SSA,MSA,f,p),
                   c(N-t,RSS,MSE,NA,NA) )
  colnames(AOVtable) = c("df", "SS", "MeanSS", "F-test","p-value")
  rownames(AOVtable) = c("Factor", "Residuals")
  return(list(table = AOVtable,CI4mu = CI.mu))
} # end my1wayAOVII

# (f)
library(lme4)
lmer.out = lmer(Y~1+(1|mach),dat)
summary(lmer.out)

library(lattice)
densityplot(~Y, groups=mach, auto.key = list(columns=4))
densityplot(~Y|reorder(mach,Y),layout=c(1,4))
bwplot(reorder(mach,Y)~Y)