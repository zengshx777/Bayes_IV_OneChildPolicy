#First Stage Plotting
#Produce Figure 2 in the paper
load("qq604_result_monotonic_weight_model.RData")

Subgroup.name<-c("Rural Females","Rural Males","Urban Females","Urban Males")
Response.name<-c("Confidence","Nervous","Desperate")

newz.p=matrix(0,nrow=1000,ncol=4)
newz.p.low=matrix(0,nrow=1000,ncol=4)
newz.p.up=matrix(0,nrow=1000,ncol=4)
newz=seq(max(ifppr$ifppr),min(ifppr$ifppr),length=1000)
source("Data_Loading.R")
for (k in 1:4){
#Extract MCMC Samples
MCMCSample<-get(paste(response_index[res.id],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;
  
#Mean of Threshold
mean_d<-MCMCSample[,"tinter"]+
  MCMCSample[,paste("theta[",1:4,"]",sep="")]%*%apply(data.jags[[k]][,covariate_index],2,mean)
#
p_treated<-sapply(newz,FUN=function(x){pnorm(-x-mean_d,mean=0,sd=sqrt(1/MCMCSample[,"taud.prec"]))})
newz.p[,k]<-apply(p_treated,2,mean)
newz.p.low[,k]<-apply(p_treated,2,FUN=function(x){quantile(x,0.05)})
newz.p.up[,k]<-apply(p_treated,2,FUN=function(x){quantile(x,0.95)})
}


pdf("Firstage_Prob_Quantile.pdf",height=5,width=6)
index=seq(1,1000,by=100)
z.quantile=seq(0,1,length=1000)
plot(z.quantile,newz.p[,1],ylim=range(newz.p),type='l',lty=1,
     lwd=2,xlab="Quantile of IV values (IFPPR)",ylab="Probability of Being the Single Child",
     main="Probability of Being the Single Child\n 4 Subgroups Comparison")
points(z.quantile[index],newz.p[index,1],pch=1,cex=1)
lines(z.quantile,newz.p[,2],type='l',lty=1,lwd=2)
points(z.quantile[index],newz.p[index,2],pch=2,cex=1)
lines(z.quantile,newz.p[,3],type='l',lty=1,lwd=2)
points(z.quantile[index],newz.p[index,3],pch=5,cex=1)
lines(z.quantile,newz.p[,4],type='l',lty=1,lwd=2)
points(z.quantile[index],newz.p[index,4],pch=6,cex=1)
# 
# lines(z.quantile,newz.p.up[,1],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.low[,1],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.up[,2],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.low[,2],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.up[,3],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.low[,3],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.up[,4],type='l',lty=2,lwd=1)
# lines(z.quantile,newz.p.low[,4],type='l',lty=2,lwd=1)
legend("topleft",legend=Subgroup.name,lty=1,pch=c(1,2,5,6),cex=0.7)
dev.off()

pdf("Firstage_Prob.pdf",height=5,width=6)
z.grid=140-newz
index=seq(1,1000,by=100)
z.quantile=seq(0,1,length=1000)
plot(z.grid,newz.p[,1],ylim=range(newz.p),type='l',lty=1,
     lwd=2,xlab="IV values (IFPPR)",ylab="Probability of Being the Single Child",
     main="Probability of Being the Single Child\n 4 Subgroups Comparison")
points(z.grid[index],newz.p[index,1],pch=1,cex=1.3)
lines(z.grid,newz.p[,2],type='l',lty=1,lwd=2)
points(z.grid[index],newz.p[index,2],pch=2,cex=1.3)
lines(z.grid,newz.p[,3],type='l',lty=1,lwd=2)
points(z.grid[index],newz.p[index,3],pch=5,cex=1.3)
lines(z.grid,newz.p[,4],type='l',lty=1,lwd=2)
points(z.grid[index],newz.p[index,4],pch=6,cex=1.3)
legend("topleft",legend=Subgroup.name,lty=1,pch=c(1,2,5,6),cex=0.7)
dev.off()