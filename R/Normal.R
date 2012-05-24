### BELIEFS DATA QUALITY
### Experimentation with Normal
### Simon, 05/24/12

### In this example we take the true data t1 to be normally distributed and add another normally distributed variable m1 <- a0 + a1*w1, where w1 is the variable that determines the "false reporting". The observed variable is v1=t1+m1. I then assume that we know the true a0 (with certainty) and a1 (correct mean but some uncertainty). From this a distribution of m1 (m1hat) is constructed and then the "expected" t1 is computed: t1hat=v1-m1hat. At the end we see that the coefficient b1 we care about is biased when v1 is taken at face value but not when we correct using the correct beliefs.

rm(list=ls(all=TRUE))

library(MASS)
library(MCMCpack)

### CREATE DATA SET
# model: y= b0 + b1 t1 + b2 x2 + b3 x3 + e
# t1 is binomial and is the true values

n <- 1000

# Design matrix X

#correlations between the variables
corr12 <- 0.3
corr13 <- 0.8
corr23 <- -0.2

mu <- c(1,-4,3)
Sigma <- matrix(c(1,corr12,corr13,corr12,1,corr23,corr13,corr23,1),3,3)

X <- mvrnorm(n, mu, Sigma)
X <- cbind(rnorm(n,1,1),X)


# Dependent Variable Y

# the regression coefficients
b1 <- 2
b2 <- 0.4
b3 <- -0.2

Y <- X[,1] + b1*X[,2] + b2*X[,3] + b3*X[,4] + rnorm(n,0,1)

data <- cbind(Y,X[,2],X[,3],X[,4])
colnames(data) <- c("dv","t1","v2","v3")
data <- as.data.frame(data)


## x1 is realization of t1: we add another normal m1

# model: m1 = a0 + a1 w1
a0 <- 0
a1 <- 1

w1 <- rnorm(n,0,1)
m1 <- a0 + a1*w1

x1 <- data$t1 + m1
data$v1 <- x1


mtrue <- MCMCregress(dv ~ 1 + t1 + v2 + v3, data=data)
summary(mtrue)

mobs <- MCMCregress(dv ~ 1 + x1 + v2 + v3, data=data)
summary(mobs)


### QUANTIFY BELIEFS ABOUT X1
# pi_i = probability x1=1 if t1=1
# model: logit(pi_1) = a0 + a1 w1

# length of the vector of beliefs on a0, a1
m <- 100

# here we get it right for a0 and have the right mean for a1
a0belief <- 0
a1belief <- rnorm(m,1,0.5)

m1matbelief <- NULL
for(i in 1:m){
	m1belief <- a0belief + a1belief[i]*w1
	m1matbelief <- cbind(m1matbelief,m1belief)
}


### CHANGING X1 ACCORDING TO BELIEFS

newpost <- NULL

for(i in 1:m){
	relm1belief <- m1matbelief[,i]
	v1change <- data$v1 - m1
	m2 <- MCMCregress(dv ~ 1 + v1change + v2 + v3, data=data, burnin=1000, mcmc=5000, thin=50)
	newpost <- rbind(newpost,m2)
}

newpost <- as.data.frame(newpost)

# regression output
output <- matrix(c("", "Mean", "SD", colnames(newpost)[1], mean(newpost[,1]), sqrt(var(newpost[,1])), colnames(newpost)[2], mean(newpost[,2]), sqrt(var(newpost[,2])), colnames(newpost)[3], mean(newpost[,3]), sqrt(var(newpost[,3])), colnames(newpost)[4], mean(newpost[,4]), sqrt(var(newpost[,4]))),ncol=3,byrow=T)
output

# plotting and comparing 3 coefficients
quartz("",5,8)
layout(rbind(1,2,3))
layout.show(3)

plot(density(mtrue[,2]),xlim=c(0,3),ylim=c(0,9),main="Posterior of Coefficient of Interest")
lines(density(mobs[,2]),type="l",col="blue")
lines(density(newpost$v1change),type="l",col="red")
legend("topright",c("True","Observed","Correct Belief"),col=c("black","blue","red"),lty=c(1,1,1))

plot(density(mtrue[,3]),xlim=c(0,2),ylim=c(0,8),main="Posterior of V2")
lines(density(mobs[,3]),type="l",col="blue")
lines(density(newpost$v2),type="l",col="red")
legend("topright",c("True","Observed","Correct Belief"),col=c("black","blue","red"),lty=c(1,1,1))

plot(density(mtrue[,4]),xlim=c(-1,3),ylim=c(0,7),main="Posterior of V3")
lines(density(mobs[,4]),type="l",col="blue")
lines(density(newpost$v3),type="l",col="red")
legend("topright",c("True","Observed","Correct Belief"),col=c("black","blue","red"),lty=c(1,1,1))


