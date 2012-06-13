rm(list=ls(all=TRUE))

library(MASS)
library(MCMCpack)

source("/Users/simonweschle/Political_Science/Research/Beliefs_Incomplete Data/Git_Beliefs_IncompleteData/R/Functions.R")

setwd("/Users/simonweschle/Political_Science/Research/Beliefs_Incomplete Data/Git_Beliefs_IncompleteData/Thoughts/graphs")

### CREATE DATA SET
# model: y= b0 + b1 t1 + b2 x2 + b3 x3 + e

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


## x1 is realization of t1: t1=x1*m1, so x1=t1/m1

# model: m1 = a0 + a1 w1
a0 <- 1.2
a1 <- 0.05

w1 <- rnorm(n,0,1)
m1 <- a0 + a1*w1

x1 <- data$t1/m1
data$v1 <- x1


mtrue <- MCMCregress(dv ~ 1 + t1 + v2 + v3, data=data)
summary(mtrue)

mobs <- MCMCregress(dv ~ 1 + v1 + v2 + v3, data=data)
summary(mobs)


### QUANTIFY BELIEFS ABOUT X1

# length of the vector of beliefs on a0, a1
m <- 100

# here we get it right for a0 and have the right mean for a1
	# the variance needs to be pretty low

a0belief <- rnorm(m,1.2,0.05)
a1belief <- rnorm(m,0.05,0.01)
beliefs <- list(a0belief, a1belief)

m1matbelief <- NULL
for(i in 1:m){
	m1belief <- a0belief + a1belief[i]*w1
	m1matbelief <- cbind(m1matbelief,m1belief)
}

bordermat <- NULL
for(i in 1:n){
	borders <- quantile(m1matbelief[i,],c(0.025,0.5,0.975))
	bordermat <- rbind(bordermat,borders)
}

#pdf(file="beliefplot.pdf",onefile=FALSE, paper="special", height=7, width=7*1.62)
plot(sort(w1),sort(bordermat[,2]),type="p",ylim=c(min(bordermat[,1]),max(bordermat[,3])))
points(sort(w1),sort(bordermat[,1]),type="l",lty=2)
points(sort(w1),sort(bordermat[,3]),type="l",lty=2)
#dev.off()



m1 <- simax.regress(dv ~ 1 + change + v2 + v3, data=data, modeltype="Gaussian", changevar=data$v1, polyorder=1, beliefslist=beliefs)
summary(m1)


a1range <- seq(-0.1,0.1,0.01)
a0range <- rep(1.2,length(a1range))
range <- list(a0range, a1range)
simax.sweep(dv ~ 1 + change + v2 + v3, data=data, modeltype="Gaussian", changevar=data$v1, polyorder=1, beliefsrange=range, coefinterest="change")



