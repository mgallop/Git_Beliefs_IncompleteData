
####################################
###      SENSITIVITY BELIEFS     ###
###		  MC CASE 1          ###
###          06/15/2012          ###
####################################


rm(list=ls(all=TRUE))

library(MASS)
library(MCMCpack)
library(FNN)
library(Matrix)


setwd("C:/Users/sww7.WIN/Documents/Beliefs")


### THINGS THAT CAN CHANGE

# 1. number of observations
nvec <- c(100,500,1000)

# 2. corr12
corr12vec <- seq(-1,1,0.2)

# 3. corr13
corr13vec <- seq(-1,1,0.2)

# 4. corr23
corr23vec <- seq(-1,1,0.2)

# 5. beliefs distribution size
mvec <- c(10,25,50,100,250)

# 6. mean of a1 distribution
a1meanvec <- seq(0,1,0.01)

# 7. var of a1 distribution
a1varvec <- seq(0,0.05,0.01)


### THE LOOP

outmat <- NULL

for(l1 in 1: length(nvec)){
for(l2 in 1: length(corr12vec)){
for(l3 in 1: length(corr13vec)){
for(l4 in 1: length(corr23vec)){
for(l5 in 1: length(mvec)){
for(l6 in 1: length(a1meanvec)){
for(l7 in 1: length(a1varvec)){

### THE DRAWS

# 1. number of observations
n <- nvec[l1]

# 2. corr12
corr12 <- corr12vec[l2]

# 3. corr13
corr13 <- corr13vec[l3]

# 4. corr23
corr23 <- corr23vec[l4]

# 5. beliefs distribution size
m <- mvec[l5]

# 6. mean of a1 distribution
a1mean <- a1meanvec[l6]

# 7. var of a1 distribution
a1var <- a1varvec[l7]

# beliefs
a0belief <- rnorm(m,1,0.01)
a1belief <- rnorm(m,a1mean,a1var)


### CREATE DATA SET
# model: y= b0 b1 x1 + b2 x2 + b3 x3 + e
# x3 is categorical

# Design matrix X

mu <- c(0,0,0)
corrs <- matrix(c(1,corr12,corr13,corr12,1,corr23,corr13,corr23,1),3,3)
Sigma <- nearPD(corrs,cor=T)$mat

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


### CORRUPT THE VARIABLE

# model: m1 = a0 + a1 w1
a0 <- 1
a1 <- 0.05

w1 <- rnorm(n,3,1)
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
	m1belief <- a0belief[i] + a1belief[i]*w1
	m1matbelief <- cbind(m1matbelief,m1belief)
}

bordermat <- NULL
for(i in 1:n){
	borders <- quantile(m1matbelief[i,],c(0.025,0.5,0.975))
	bordermat <- rbind(bordermat,borders)
}

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

posterior.reduce <- posterior[seq(1,dim(posterior)[1],dim(posterior)[1]/10000),]

div1 <- c(mean(KL.divergence(mtrue[,2], mobs[,2], k=5),na.rm=T), mean(KL.divergence(mtrue[,2], posterior.reduce[,2], k=5),na.rm=T))
div2 <- c(mean(KL.divergence(mtrue[,3], mobs[,3], k=5),na.rm=T), mean(KL.divergence(mtrue[,3], posterior.reduce[,3], k=5),na.rm=T))
div3 <- c(mean(KL.divergence(mtrue[,4], mobs[,4], k=5),na.rm=T), mean(KL.divergence(mtrue[,4], posterior.reduce[,4], k=5),na.rm=T))

out <- c(n, Sigma[1,2], Sigma[1,3], Sigma[2,3], m, a1mean, a1var, div1[1], div1[2], div2[1], div2[2], div3[1], div3[2])
outmat <- rbind(outmat,out)

}}}}}}}

colnames(outmat) <- c("n","corr12","corr13","corr23","m","a1mean","a1var","KL.x1.obs","KL.x1.bel","KL.x2.obs","KL.x2.bel","KL.x2.obs","KL.x2.bel")

save(outmat, file="MC_Case1.rda")	

