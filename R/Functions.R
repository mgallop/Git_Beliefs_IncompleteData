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
