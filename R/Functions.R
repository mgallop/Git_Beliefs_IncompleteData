simax.regress <- function(formula, data, modeltype, changevar, polyorder, beliefslist){
	
	# change the data according to the beliefs
	beliefsdata <- matrix(NA,nrow=dim(data)[1],ncol=length(beliefslist[[1]]))
	if(polyorder==1){
		for(i in 1:length(beliefslist[[1]])){beliefsdata[,i] <- changevar*(beliefslist[[1]][i] + changevar*beliefslist[[2]][i])}
	}
	if(polyorder==2){
		for(i in 1:length(beliefslist[[1]])){beliefsdata[,i] <- changevar*(beliefslist[[1]][i] + changevar*beliefslist[[2]][i] + changevar^2*beliefslist[[3]][i])}
	}
	if(polyorder==3){
		for(i in 1:length(beliefslist[[1]])){beliefsdata[,i] <- changevar*(beliefslist[[1]][i] + changevar*beliefslist[[2]][i] + changevar^2*beliefslist[[3]][i] + changevar^3*beliefslist[[4]][i])}
	}
	
	# run the model
	if(modeltype=="Gaussian"){
		posterior <- NULL
		for(i in 1:dim(beliefsdata)[2]){
			data$change <- beliefsdata[,i]
			model <- MCMCregress(dv ~ 1 + change + v2 + v3, data=data, burnin=1000, mcmc=5000, thin=50)
			posterior <- rbind(posterior,model)
		}
	}
	posterior <- as.mcmc(posterior)
	invisible(return(posterior))
}


simax.sweep <- function(formula, data, modeltype, changevar, polyorder, beliefsrange, coefinterest){
	
	# change the data according to the beliefs
	beliefsdata <- matrix(NA,nrow=dim(data)[1],ncol=length(beliefsrange[[1]]))
	if(polyorder==1){
		for(i in 1:length(beliefsrange[[1]])){beliefsdata[,i] <- changevar*(beliefsrange[[1]][i] + changevar*beliefsrange[[2]][i])}
	}
	if(polyorder==2){
		for(i in 1:length(beliefsrange[[1]])){beliefsdata[,i] <- changevar*(beliefsrange[[1]][i] + changevar*beliefsrange[[2]][i] + changevar^2*beliefsrange[[3]][i])}
	}
	if(polyorder==3){
		for(i in 1:length(beliefsrange[[1]])){beliefsdata[,i] <- changevar*(beliefsrange[[1]][i] + changevar*beliefsrange[[2]][i] + changevar^2*beliefsrange[[3]][i] + changevar^3*beliefsrange[[4]][i])}
	}
	
	# run the model
	if(modeltype=="Gaussian"){
		pointest <- lowerest <- upperest <- NULL
		for(i in 1:dim(beliefsdata)[2]){
			data$change <- beliefsdata[,i]
			model <- MCMCregress(dv ~ 1 + change + v2 + v3, data=data, burnin=1000, mcmc=10000, thin=10)
			addmed <- apply(model,2,median)
			inter <- HPDinterval(model,0.95)
			pointest <- rbind(pointest,addmed)
			lowerest <- rbind(lowerest,inter[,1])
			upperest <- rbind(upperest,inter[,2])
		}
	}
	NameVector <- beliefsrange[[2]]	
	howlong <- length(colnames(model))-1
	par(mfrow=c(howlong,1))
	for(i in 1:howlong){
		plot(NameVector,pointest[,i],ylim=c(min(lowerest[,i]),max(upperest[,i])),pch=16)
		for(j in 1:length(NameVector)){
			lines(x=c(NameVector[j],NameVector[j]),y=c(lowerest[j,i],upperest[j,i]))
		}
	}
}


