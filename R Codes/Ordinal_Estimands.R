#This Script Calculates TT PRTE (Ordinal) from MCMCSample
#Result for Figure 5
logit<-function(x){1/(1+exp(-x))}
zmin=-range(data1$ifppr)[2];zmax=-range(data1$ifppr)[1]

#Utility Functions for Ordinal Estimands
source("Utility_Functions.R")
#Store Results
RHO<-ETA<-NULL

library(doParallel)
library(foreach)
#Register Multiple Cluster
registerDoParallel(4)
for(res in 1:3){
source("DataLoading.R")
#List to Store Results
Aggre_result=vector("list",length=4)
for(k in 1:4)
{
  #Extract MCMC Samples
  MCMCSample<-get(paste(response_index[res],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;
  
  #Level of Response
  res.level=length(unique(data.jags[[k]][,response_index[res]]))
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


#Calculate ETA
ETA <- rbind(ETA, c(
  #TT Median of Lower and Upper Bound
  1 - apply(
    Aggre_result[[k]][, c(8, 7)],
    2,
    FUN = function(x) {
      median(x, na.rm = T)
    }
  ),
  #TT Lower and Upper CI
  1 - eval_credible(Aggre_result[[k]][, 8], Aggre_result[[k]][, 7]),
  #PRTE Median of Lower and Upper Bound
  1 - apply(
    Aggre_result[[k]][, c(2, 1)],
    2,
    FUN = function(x) {
      median(x, na.rm = T)
    }
  ),
  #PRTE Lower and Upper CI
  1 - eval_credible(Aggre_result[[k]][, 2], Aggre_result[[k]][, 1])
))

#RHO
RHO <- rbind(RHO, c(
  #TT Median of Lower and Upper Bound
  1 - apply(
    Aggre_result[[k]][, c(11, 10)],
    2,
    FUN = function(x) {
      median(x, na.rm = T)
    }
  ),
  #TT CI of Lower and Upper Bound
  1 - eval_credible(Aggre_result[[k]][, 11], Aggre_result[[k]][, 10]),
  #PRTE Median of Lower and Upper Bound
  1 - apply(
    Aggre_result[[k]][, c(5, 4)],
    2,
    FUN = function(x) {
      median(x, na.rm = T)
    }
  ),
  #PRTE CI of Lower and Upper Bound
  1 - eval_credible(Aggre_result[[k]][, 5], Aggre_result[[k]][, 4])
))

print(paste("==",res,"==finish",k))}
}

stopImplicitCluster()

