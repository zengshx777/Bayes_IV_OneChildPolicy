#First Stage Plotting
#Produce Figure 3 in the paper
Subgroup.name<-c("Rural Females","Rural Males","Urban Females","Urban Males")
Response.name<-c("Confidence","Anxiety","Desperation")

newz.p=matrix(0,nrow=1000,ncol=4)
newz.p.low=matrix(0,nrow=1000,ncol=4)
newz.p.up=matrix(0,nrow=1000,ncol=4)
newz=seq(max(ifppr$ifppr),min(ifppr$ifppr),length=1000)
res=3
source("DataLoading.R")
for (k in 1:4){
#Extract MCMC Samples
MCMCSample<-get(paste(response_index[res],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;
  
#Mean of Threshold
mean_d<-MCMCSample[,"tinter"]+
  MCMCSample[,paste("theta[",1:4,"]",sep="")]%*%apply(data.jags[[k]][,covariate_index],2,mean)
#
p_treated<-sapply(newz,FUN=function(x){pnorm(-x-mean_d,mean=0,sd=sqrt(1/MCMCSample[,"taud.prec"]))})
newz.p[,k]<-apply(p_treated,2,mean)
newz.p.low[,k]<-apply(p_treated,2,FUN=function(x){quantile(x,0.05)})
newz.p.up[,k]<-apply(p_treated,2,FUN=function(x){quantile(x,0.95)})
}

