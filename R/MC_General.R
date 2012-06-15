
####################################
###      SENSITIVITY BELIEFS     ###
###		      MC SIMON           ###
###          06/14/2012          ###
####################################


rm(list=ls(all=TRUE))

library(MASS)
library(MCMCpack)
library(FNN)

source("/Users/simonweschle/Political_Science/Research/Beliefs_Incomplete Data/Git_Beliefs_IncompleteData/R/Functions.R")

setwd("/Users/simonweschle/Political_Science/Research/Beliefs_Incomplete Data/Git_Beliefs_IncompleteData/Paper/graphs")


### THINGS THAT CAN CHANGE

# Number of observations
n <- 1000

#correlations between the variables
corr12 <- 0.3
corr13 <- -0.3
corr23 <- 0

# length of the vector of beliefs on a0, a1
m <- 100

# beliefs
a0belief <- rnorm(m,1,0.05)
a1belief <- rnorm(m,0.05,0.01)



### CREATE DATA SET
# model: y= b0 b1 x1 + b2 x2 + b3 x3 + e
# x3 is categorical

# Design matrix X

mu <- c(0,0,0)
Sigma <- matrix(c(1,corr12,corr13,corr12,1,corr23,corr13,corr23,1),3,3)

X <- mvrnorm(n, mu, Sigma)
cuts <- quantile(X[,3],c(1/3,2/3))
cat <- ifelse(X[,3]<cuts[[1]],1,2)
cat <- ifelse(X[,3]>cuts[[2]],3,cat)
X <- cbind(X[,1],X[,2],cat)

# Dependent Variable Y

# the regression coefficients
b1 <- 1
b2 <- 1
b3 <- 1

Y <- b1*X[,1] + b2*X[,2] + b3*X[,3] + rnorm(n,0,1)

data <- cbind(Y,X[,1],X[,2],X[,3])
colnames(data) <- c("dv","x1","x2","x3")
data <- as.data.frame(data)




### CASE 1: x1 is badly measured (t1 is true), w is continuous
## x1 is realization of t1: t1=x1*m1, so x1=t1/m1

# model: m1 = a0 + a1 w1
a0 <- 1
a1 <- 0.05

w1 <- rnorm(n,5,1)
m1 <- a0 + a1*w1
data$x1measured <- data$x1/m1

mtrue <- MCMCregress(dv ~ 1 + x1 + x2 + x3, data=data)
summary(mtrue)

mobs <- MCMCregress(dv ~ 1 + x1measured + x2 + x3, data=data)
summary(mobs)


# beliefs
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

#pdf(file="beliefplotS1.pdf",onefile=FALSE, paper="special", height=7, width=7*1.62)
plot(sort(w1),sort(bordermat[,2]),type="p",ylim=c(min(bordermat[,1]),max(bordermat[,3])))
points(sort(w1),sort(bordermat[,1]),type="l",lty=2)
points(sort(w1),sort(bordermat[,3]),type="l",lty=2)
#dev.off()


beliefslist <- beliefs
changevar <- data$x1measured
beliefsdata <- matrix(NA,nrow=dim(data)[1],ncol=length(beliefslist[[1]]))
for(i in 1:length(beliefslist[[1]])){
	beliefsdata[,i] <- changevar*(beliefslist[[1]][i] + w1*beliefslist[[2]][i])
}

posterior <- NULL
for(i in 1:dim(beliefsdata)[2]){
	data$change <- beliefsdata[,i]
	model <- MCMCregress(dv ~ 1 + change + x2 + x3, data=data)
	posterior <- rbind(posterior,model)
}

posterior.mcmc <- as.mcmc(posterior)
posterior.data <- as.data.frame(posterior)

summary(posterior.mcmc)

posterior.reduce <- posterior[seq(1,dim(posterior)[1],100),]


quartz("",5,8)
layout(rbind(1,2,3))
layout.show(3)

plot(density(mtrue[,2]),xlim=c(0.5,1.5),ylim=c(0,13),main="x1")
lines(density(mobs[,2]),type="l",col="blue")
lines(density(posterior.reduce[,2]),type="l",col="red")
legend("topright",c("True","Observed","Correct Belief"),col=c("black","blue","red"),lty=c(1,1,1))

plot(density(mtrue[,3]),xlim=c(0.5,1.5),ylim=c(0,13),main="x2")
lines(density(mobs[,3]),type="l",col="blue")
lines(density(posterior.reduce[,3]),type="l",col="red")
legend("topright",c("True","Observed","Correct Belief"),col=c("black","blue","red"),lty=c(1,1,1))

plot(density(mtrue[,4]),xlim=c(0.5,1.5),ylim=c(0,13),main="x3")
lines(density(mobs[,4]),type="l",col="blue")
lines(density(posterior.reduce[,4]),type="l",col="red")
legend("topright",c("True","Observed","Correct Belief"),col=c("black","blue","red"),lty=c(1,1,1))


div1 <- c(mean(KL.divergence(mtrue[,2], mobs[,2], k=5)), mean(KL.divergence(mtrue[,2], posterior.reduce[,2], k=5)))
div2 <- c(mean(KL.divergence(mtrue[,3], mobs[,3], k=5)), mean(KL.divergence(mtrue[,3], posterior.reduce[,3], k=5)))
div3 <- c(mean(KL.divergence(mtrue[,4], mobs[,4], k=5)), mean(KL.divergence(mtrue[,4], posterior.reduce[,4], k=5)))

div <- rbind(div1,div2,div3)
colnames(div) <- c("observed","belief-corrected")
div

