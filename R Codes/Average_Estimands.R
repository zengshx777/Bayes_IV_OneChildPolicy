#This Script Calculates TT PRTE (Average) from MCMCSample
#Result for Figure 4 and 5
logit<-function(x){1/(1+exp(-x))}
Response.name<-c("Confidence","Anxiety","Desperation")
Subgroup.name=c("Rural Female","Rural Male","Urban Female","Urban Male")
ncov=length(covariate_index)

zmin=-range(data1$ifppr)[2];zmax=-range(data1$ifppr)[1]

#Load the Utility function for calculating Average Estimands
source("Utility_Functions.R")
#Number of Grid to Evaluate the MTE
num.grid=1000
#Store the MTE Results
mte=vector("list",length=4)
TAU=NULL

#Calculating the posterior distribution PRTE and ATT
#Load Package for Parallel Computation
library(doParallel)
library(foreach)
#Register Multiple Cluster
registerDoParallel(4)

for(res in 1:3)
{
  #Store MTE Result
  mte[[res]]=vector("list",length=4)
  #Store PRTE,TT Result
  Aggre_result=vector("list",length=4)
  #Generating Data
  source("DataLoading.R")
for(k in 1:4)
{
#Extract MCMC Samples
MCMCSample<-get(paste(response_index[res],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;

#Level of Response(4 or 5)
res.level=length(unique(data.jags[[k]][,response_index[res]]))

#Sample size
nsample<-nrow(data.jags[[k]])

#Support of IV
suppz<-seq(range(data1$ifppr)[1],range(data1$ifppr)[2],length=num.grid)

# #Mean Covariates Value
cov.mean<-apply(data.jags[[k]][,covariate_index],2,mean)
#Adjust to Han ethnicity
cov.mean[4]=1
#MTE Effect
#Grid of Threshold Based on Logistic Distribution
thres_mean<-median(MCMCSample[,"tinter"]+MCMCSample[,paste("theta[",1:ncov,"]",sep="")]%*%
                     cov.mean)
thres_sd<-median(sqrt(1/MCMCSample[,"taud.prec"]))
#Quantile of Principal Strate Given mean covariates value
thres.grid<-qnorm(seq(0.01,0.99,length=num.grid),thres_mean,thres_sd)

#Effect on Expectation Y
#Latent Variable
Y_1_latent <-
  apply(
    MCMCSample[, "gamma1"] %*% t(thres.grid),
    2,
    FUN = function(x) {
      x +
        MCMCSample[, paste("beta1[", 1:ncov, "]", sep = "")] %*% cov.mean
    }
  )

Y_0_latent <-
  apply(
    MCMCSample[, "gamma0"] %*% t(thres.grid),
    2,
    FUN = function(x) {
      x +
        MCMCSample[, paste("beta0[", 1:ncov, "]", sep = "")] %*% cov.mean
    }
  )

#Calculate Expected Value
Y_1_expect <- apply(
  Y_1_latent,
  2,
  FUN = function(x) {
    res.level - rowSums(apply(
      MCMCSample[, paste("c1[", 1:(res.level - 1), "]", sep = "")] + 0.2,
      2,
      FUN = function(y) {
        logit(y + x)
      }
    ))
  }
)

Y_0_expect <- apply(
  Y_0_latent,
  2,
  FUN = function(x) {
    res.level - rowSums(apply(
      MCMCSample[, paste("c0[", 1:(res.level - 1), "]", sep = "")],
      2,
      FUN = function(y) {
        logit(y + x)
      }
    ))
  }
)

#Expected Value MTE 
mte_expect<-Y_1_expect-Y_0_expect
#Store the Result for Subgroup and Measures
mte[[res]][[k]]=rbind(apply(mte_expect,2,mean),
                      apply(mte_expect,2,FUN=function(x){quantile(x,0.975)}),
                      apply(mte_expect,2,FUN=function(x){quantile(x,0.025)}))

#Calculate TT and PRTE 
treated.id=which(data.jags[[k]]$onechild==1)
#MCMC Sample Size
samplesize=nrow(MCMCSample)

print("==parallel start==")
#Parallel Computing
Aggre_result[[k]]<-foreach(j=1:samplesize,.combine=rbind)%dopar%
  GetAverage(j,MCMCSample,data=data.jags[[k]],treated.id=treated.id,ncov=ncov,
             zmax=zmax,zmin=zmin,res.level = res.level)

#Construct TAU
PRTE=Aggre_result[[k]][,1]
TT=Aggre_result[[k]][,2]
TAU <- rbind(TAU, c(
  median(PRTE, na.rm = T),
  quantile(PRTE, c(0.025, 0.975), na.rm = T),
  c(median(TT, na.rm = T),
    quantile(TT, c(0.025, 0.975), na.rm = T))
))

print(paste("==",res,"==finish",k))
}
}
#Quit Parallel Computation
stopImplicitCluster()


