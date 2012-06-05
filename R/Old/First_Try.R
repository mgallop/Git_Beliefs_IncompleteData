### BELIEFS DATA QUALITY
### Some experimentation
### Simon, 05/19/12

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
X1old <- X[,1]
makebin <- ifelse(X[,1]>quantile(X[,1],0.4),1,0)
X[,1] <- makebin
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


## x1 is realization of t1: x1=0 if t1=0, if t1=1 then x1=1 w/ prob pi_i

# model: logit(pi_1) = a0 + a1 w1
a0 <- 0
a1 <- 1

w1 <- X1old + rnorm(n,0,0)
pihelp <- a0 + a1*w1
pi <- exp(pihelp)/(1+exp(pihelp))


replace <- rbinom(n,1,pi)
x1 <- ifelse(data$t1==1,replace,0)
data$v1 <- x1


mtrue <- MCMCregress(dv ~ 1 + t1 + v2 + v3, data=data)
summary(mtrue)

mobs <- MCMCregress(dv ~ 1 + x1 + v2 + v3, data=data)
summary(mobs)


### QUANTIFY BELIEFS ABOUT X1
# pi_i = probability x1=1 if t1=1
# model: logit(pi_1) = a0 + a1 w1

# length of the vector of beliefs on a0, a1
m <- 50

# here we get it right for a0 and have the right mean for a1
a0belief <- 0
a1belief <- rnorm(m,1,0)

pimatbelief <- NULL
for(i in 1:m){
	pihelpbelief <- a0belief + a1belief[i]*w1
	pibelief <- exp(pihelpbelief)/(1+exp(pihelpbelief))
	pimatbelief <- cbind(pimatbelief,pibelief)
}


### CHANGING X1 ACCORDING TO BELIEFS

newpost <- NULL

for(i in 1:m){
	relpibelief <- pimatbelief[,i]
	# P(t=1 | x=0)
	gamma <- mean(data$v1)/mean(relpibelief)
	postbelief <- ((1-relpibelief)*gamma)/(1-relpibelief*gamma)
	replace <- rbinom(n,1,postbelief)
	v1change <- ifelse(data$v1==1,1,replace)
	#v1change <- ifelse(data$v1==1,1,postbelief)
	data$v1change <- v1change
	m2 <- MCMCregress(dv ~ 1 + v1change + v2 + v3, data=data, burnin=1000, mcmc=5000, thin=50)
	newpost <- rbind(newpost,m2)
}

newpost <- as.data.frame(newpost)
mean(newpost$v1change)

plot(density(mtrue[,2]),xlim=c(0,3))
lines(density(mobs[,2]),type="l",col="blue")
lines(density(newpost$v1change),type="l",col="red")
