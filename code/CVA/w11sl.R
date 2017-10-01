# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #
# TUTORIAL & LAB  --  WEEK 11                                       #
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

# set working directory to where the data is stored
setwd("/Users/mueller/Desktop/usyd/teaching/S3012/data")  


# Question 01: 
# function that calculates Yates Algorithm iff y (response) is in standard order)
YATES = function(y){
  N = length(y)
  n = round(log(N)/log(2),0)
  Y = matrix(rep(0,(n+1)*N),ncol=(n+1))
  Y[,1] = y
  for (k in 2:(n+1)){
     for (j in 1:(N/2)){
     j1=2*j-1; j2=2*j;
     Y[j,k] = Y[j1,k-1] + Y[j2,k-1]
     Y[j+N/2,k] = Y[j2,k-1] - Y[j1,k-1]
     }
  }
  return(Y);
}

# Yates' Algorithm Example:
yi. = c(29,28,32,33,48,51,48,56)
YATES(yi.)

# b
sum(yi.**2)/4 - 325^2/32
y.j = c(55,82,98,90)
sum(y.j**2)/8 - 325^2/32

# c)
A = factor(c(0,1,0,1,0,1,0,1)); A = rep(A,4)
B = factor(c(0,0,1,1,0,0,1,1)); B = rep(B,4)
C = factor(c(0,0,0,0,1,1,1,1)); C = rep(C,4)
Block = factor(rep(1:4,c(8,8,8,8)))
Y = c(7, 3, 4, 5, 9, 8,11, 8,
      8, 7, 9, 8,13,14, 9,14,
      8, 7,11,12,12,15,16,17,
      6,11, 8, 8,14,14,12,17)
M1 = aov(Y~Block+(A*B*C))
summary(M1)

# d)
M2 = aov(Y~Block+C*A*B+Error(Block/C))
summary(M2)



# Question 03: 
#              V01 V01 ... V08
#              - - - - - - - - 
#  BLOCK I:    F01 F02 ... F08 
#  BLOCK II:   F09 F20 ... F16
#  BLOCK III:  F17 F18 ... F24

#(a)
dat = read.table(file="rubber.txt",header=TRUE) 
summary(dat)

names(dat)
attach(dat)
tapply(plants,variety,sum)
sum(tapply(plants,variety,sum)^2)
tapply(plants,treatment,sum)
sum(tapply(plants,treatment,sum)^2)
tapply(plants,flats,sum)
sum(tapply(plants,flats,sum)^2)
tapply(plants,list(variety,treatment),sum)
sum(tapply(plants,list(variety,treatment),sum)^2)
sum(plants^2)
sum(plants)

#(b)
M1 = lm(plants~treatment*variety)
anova(M1)
attach(dat)
interaction.plot(treatment,variety,plants)

#(c)
M2 = lm(plants~variety/flats)
anova(M2)

#(d)
M3 = lm(plants~flats+treatment)
anova(M3)

#(e)
# note that the order how variables are entered matters, e.g. lm(plants~flats + variety*treatment)
# produces a different output to the one below, namely 23 df are spent for the flats and there are
# no left for modelling main effect of variety
M4 = lm(plants~variety*treatment + flats)
anova(M4)
M4 = lm(plants~variety+variety:flats + treatment + variety:treatment)
anova(M4)

# sum of squares by hand:
#(f)i
sum(tapply(plants,variety,sum)^2)/12 - sum(plants)**2/96

#(f)ii
sum(tapply(plants,treatment,sum)^2)/24 - sum(plants)**2/96

#(f)iii
sum(tapply(plants,list(variety,treatment),sum)^2)/3 - sum(tapply(plants,treatment,sum)^2)/24 -sum(tapply(plants,variety,sum)^2)/12 + sum(plants)**2/96


