### BELIEFS DATA QUALITY
### Some experimentation
### Simon, 05/19/12

rm(list=ls(all=TRUE))

library(MASS)
library(MCMCpack)

### CREATE DATA SET
# model: y= b0 + b1 x1 + b2 x2 + b3 x3 + e
# x1 is binomial and of bad data quality

n <- 1000

# Design matrix X

#correlations between the variables
corr12 <- 0.3
corr13 <- 0.8
corr23 <- -0.2

mu <- c(1,-4,3)
Sigma <- matrix(c(1,corr12,corr13,corr12,1,corr23,corr13,corr23,1),3,3)

X <- mvrnorm(n, mu, Sigma)
X1old <- X[,1]
makebin <- ifelse(X[,1]>quantile(X[,1],0.4),1,0)
X[,1] <- makebin
X <- cbind(rnorm(n,1,1),X)

# Dependent Variable Y

# the regression coefficients
b1 <- 2
b2 <- -0.4
b3 <- 0.2

Y <- X[,1] + b1*X[,2] + b2*X[,3] + b3*X[,4] + rnorm(n,0,1)

data <- cbind(Y,X[,2],X[,3],X[,4])
colnames(data) <- c("dv","v1","v2","v3")
data <- as.data.frame(data)

m1 <- MCMCregress(dv ~ 1 + v1 + v2 + v3, data=data)
summary(m1)



### QUANTIFY BELIEFS ABOUT X1
# pi_i = probability x1=1 if t1=1
# model: logit(pi_1) = a0 + a1 w1

m <- 2

w1 <- X1old + rnorm(n,1,1)

a0 <- 0
a1 <- rnorm(m,-1,1)

pimat <- NULL
for(i in 1:m){
	pihelp <- a0 + a1[i]*w1
	pi <- exp(pihelp)/(1+exp(pihelp))
	pimat <- cbind(pimat,pi)
}


### CHANGING X1 ACCORDING TO BELIEFS

newpost <- NULL

for(i in 1:m){
	relevantpi <- pimat[,i]
	replace <- rbinom(n,1,relevantpi)
	v1change <- ifelse(data$v1==0,0,replace)
	data$v1change <- v1change
	m2 <- MCMCregress(dv ~ 1 + v1change + v2 + v3, data=data, burnin=1000, mcmc=5000, thin=50)
	newpost <- rbind(newpost,m2)
}

newpost <- as.data.frame(newpost)
mean(newpost$v1change)



