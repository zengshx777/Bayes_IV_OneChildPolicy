library(msm)
load("qq603_result_monotonic_weight_model.RData")

logit<-function(x){1/(1+exp(-x))}
Response.name<-c("Confidence","Nervous","Nervous","Desperate","Tension",
                 "Difficulty","Meaningless","Confidence","Confidence")
Subgroup.name=c("Rural Female","Rural Male","Urban Female","Urban Male")
ncov=length(covariate_index)

zmin=-range(data1$ifppr)[2];zmax=-range(data1$ifppr)[1]
#Generating Data
source("Data_Loading.R")
source("Average_Utility_Func.R")
num.grid=1000
for(k in 1:4)
{
#Extract MCMC Samples
MCMCSample<-get(paste(response_index[res.id],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;

#Level of Response
res.level=length(unique(data.jags[[k]][,response_index[res.id]]))

num.clu<-length(unique(data.jags[[k]][,"province"]))
num.obs<-nrow(data.jags[[k]])
thres.hold<-MCMCSample[,paste("d[",1:num.obs,"]",sep="")]


#Support of IV
suppz<-seq(range(data1$ifppr)[1],range(data1$ifppr)[2],length=num.grid)

# #Mean Covariates Value
cov.mean<-apply(data.jags[[k]][,covariate_index],2,mean)
#Adjust to Han ethnicity
cov.mean[4]=1
#MTE Effect
#Grid of Threshold Based on Logistic Distribution
# {
thres_mean<-median(MCMCSample[,"tinter"]+MCMCSample[,paste("theta[",1:ncov,"]",sep="")]%*%
                     cov.mean)
thres_sd<-median(sqrt(1/MCMCSample[,"taud.prec"]))
thres.grid<-qnorm(seq(0.01,0.99,length=num.grid),thres_mean,thres_sd)

#thres.grid<-seq(zmin,zmax,length=num.grid)
#Effect on Expectation Y
#Latent Variable
Y_1_latent<-apply(MCMCSample[,"gamma1"]%*%t(thres.grid),2,FUN=function(x){x+
    MCMCSample[,paste("beta1[",1:ncov,"]",sep="")]%*%cov.mean})


Y_0_latent<-apply(MCMCSample[,"gamma0"]%*%t(thres.grid),2,FUN=function(x){x+
     MCMCSample[,paste("beta0[",1:ncov,"]",sep="")]%*%cov.mean})

#Calculate Expected Value
Y_1_expect<-apply(Y_1_latent,2,FUN=function(x){res.level-rowSums(
apply(MCMCSample[,paste("c1[",1:(res.level-1),"]",sep="")],2,FUN=function(y){logit(y+x)}))})

Y_0_expect<-apply(Y_0_latent,2,FUN=function(x){res.level-rowSums(
  apply(MCMCSample[,paste("c0[",1:(res.level-1),"]",sep="")],2,FUN=function(y){logit(y+x)}))})



#Expected Value MTE Plotting
mte_expect<-Y_1_expect-Y_0_expect
#pdf(file = paste(response_index[res.id],"Expected",k,"MTE.pdf"
#                 ,sep="_"),height=5,width=6)
plot(seq(0.01,0.99,length=num.grid),apply(mte_expect,2,mean),ylim=c(-4,4),type='l',ylab='MTE',lwd=2,
     xlab="Principal Strata (Quantile)",main=paste("Subgroup",Subgroup.name[k],sep=" "))
lines(seq(0.01,0.99,length=num.grid),apply(mte_expect,2,FUN=function(x){quantile(x,0.025)}),lty=2,lwd=2)
lines(seq(0.01,0.99,length=num.grid),apply(mte_expect,2,FUN=function(x){quantile(x,0.975)}),lty=2,lwd=2)
abline(h=0,lty=4,lwd=1)
legend("topright",legend=c("Posterior Mean","95% Credible Intervals"),lty=c(1,2))
}
#dev.off()
# }
#Calculating the posterior distribution PRTE and ATT


library(doParallel)
library(foreach)
#Register Multiple Cluster
registerDoParallel(4)
Aggre_result=vector("list",length=4)
#Extract Treated ID
for(k in 1:4){
#Extract MCMC Samples
MCMCSample<-get(paste(response_index[res.id],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;
  
#Level of Response
res.level=length(unique(data.jags[[k]][,response_index[res.id]]))
treated.id=which(data.jags[[k]]$onechild==1)


samplesize=nrow(MCMCSample)

print("==parallel start==")
#Parallel Computing
sample.index=sample(1:samplesize,10000,replace=T)
 Aggre_result[[k]]<-foreach(j=1:samplesize,.combine=rbind)%dopar%
   GetAverage(j,MCMCSample,data=data.jags[[k]],treated.id=treated.id,ncov=ncov,
              zmax=zmax,zmin=zmin,res.level = res.level)
# Aggre_result[[k]]<-t(sapply(1:samplesize,function(j) GetAverage(j,MCMCSample,data=data.jags[[k]],treated.id=treated.id,ncov=ncov,
#          zmax=zmax,zmin=zmin,res.level = res.level)))
  

assign(paste(response_index[res.id],k,"LATE",sep="_"),Aggre_result[[k]][,1])
assign(paste(response_index[res.id],k,"ATT",sep="_"),Aggre_result[[k]][,2])

print(paste("==",res.id,"==finish",k))
}
save(Aggre_result,file=paste(response_index[res.id],"_Average_ATTLATE.RData"))
stopImplicitCluster()


r_table<-NULL
for (k in 1:4)
{
  r_table<-rbind(r_table,c(mean(Aggre_result[[k]][,1]),
                           quantile(Aggre_result[[k]][,1],c(0.025,0.975))))
}
for (k in 1:4)
{
  r_table<-rbind(r_table,c(mean(Aggre_result[[k]][,2]),
                           quantile(Aggre_result[[k]][,2],c(0.025,0.975))))
}
library(xtable)
xtable(r_table,digits=3)
