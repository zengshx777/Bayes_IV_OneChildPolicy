#This Script Calculate TT PRTE from MCMCSample
#This Script Produce Table 3 and 4 in the Paper
library(msm)
#load("qq603_result_monotonic_weight_model.RData")

logit<-function(x){1/(1+exp(-x))}
Response.name<-c("Confidence","Nervous","Nervous","Desperate","Tension",
                 "Difficulty","Meaningless","Confidence","Confidence")
Subgroup.name=c("Rural Female","Rural Male","Urban Female","Urban Male")
ncov=length(covariate_index)

zmin=-range(data1$ifppr)[2];zmax=-range(data1$ifppr)[1]

source("Ordinal_Utility_Func.R")


source("Data_Loading.R")
#List to Store Results


library(doParallel)
library(foreach)
#Register Multiple Cluster
registerDoParallel(4)
Aggre_result=vector("list",length=4)
for(k in 1:4)
{
  #Extract MCMC Samples
  MCMCSample<-get(paste(response_index[res.id],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;
  
  #Level of Response
  res.level=length(unique(data.jags[[k]][,response_index[res.id]]))
  treated.id=which(data.jags[[k]]$onechild==1)
  
  
  samplesize=nrow(MCMCSample)
  

#ATT and PRTE Calculation
#Extract Treated ID
treated.id=which(data.jags[[k]]$onechild==1)
#Use the empirical distribution of X and Threshold Value
print("==parallel start==")
Aggre_result[[k]]<-foreach(j=1:samplesize,.combine=rbind)%dopar%
  GetBound(j,MCMCSample,data=data.jags[[k]],treated.id=treated.id,ncov=ncov,
           zmax=zmax,zmin=zmin,res.level = res.level)
print(paste("==",res.id,"==finish",k))}



save(Aggre_result,file=paste(response_index[res.id],"_ATTLATE.RData"))
# }
stopImplicitCluster()

