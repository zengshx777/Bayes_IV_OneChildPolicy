#Model Checking Template Sample Predictive P-value 
#Produce Table 5 in the paper
#This is paralleled version
#Extract MCMC Samples
#load("qq604_result_monotonic_weight_model.RData")
ncov=length(covariate_index)
sppv_comparison<-function(j,rep.time,Checking_Sample,Original_Data)
{
  tryCatch({
  #Sample Size
  nsize=nrow(Original_Data)
  #Covariates
  X.obs<-as.matrix(Original_Data[,covariate_index])
  #Response
  yobs<-Original_Data[,response_index[res.id]];res.level=length(unique(yobs))
  
  #IV Values
  zvalues=-Original_Data$ifppr
  #Dichtomize into two categories
  low.iv.ind<-which(zvalues<=median(zvalues))
  high.iv.ind<-which(zvalues>median(zvalues))
  
  #Temp Discrepancy Restore for replicated data set
  temp.signal.rep.nt<-temp.signal.rep.at<-temp.signal.rep.cp<-numeric(rep.time)
  temp.noise.rep.nt<-temp.noise.rep.at<-temp.noise.rep.cp<-numeric(rep.time)
  temp.snr.rep.nt<-temp.snr.rep.at<-temp.snr.rep.cp<-numeric(rep.time)
  
  
  #IV Value
  zsupport<-range(zvalues)
  
  #Threshold Value
  c1=Checking_Sample[j,paste("c1[", 1:(res.level-1), "]",sep="")]
  c0=Checking_Sample[j,paste("c0[", 1:(res.level-1), "]",sep="")]
  
  
  #Observed Data Measure Calculation
  #Imputated Principal Strata Based on Observed Data
  #Here we use the strata directly to avoid potential 
  #conflications with observed Treatment status
  obs.d<-Checking_Sample[j,paste("d[",1:nsize,"]",sep="")]
  nt.obs.id<-which(obs.d>zsupport[2])
  at.obs.id<-which(obs.d<zsupport[1])
  cp.obs.id<-which(obs.d<=zsupport[2]&obs.d>=zsupport[1])
  
  obs.setA<-intersect(nt.obs.id,low.iv.ind);sizeA=length(obs.setA)
  obs.setB<-intersect(nt.obs.id,high.iv.ind);sizeB=length(obs.setB)
  
  signal.obs.nt<-mean(yobs[obs.setA])-mean(yobs[obs.setB])
  noise.obs.nt<-sqrt(var(yobs[obs.setA])/sizeA+var(yobs[obs.setB])/sizeB)
  snr.obs.nt<-signal.obs.nt/noise.obs.nt
  
  obs.setA<-intersect(at.obs.id,low.iv.ind);sizeA=length(obs.setA)
  obs.setB<-intersect(at.obs.id,high.iv.ind);sizeB=length(obs.setB)
  
  signal.obs.at<-mean(yobs[obs.setA])-mean(yobs[obs.setB])
  noise.obs.at<-sqrt(var(yobs[obs.setA])/sizeA+var(yobs[obs.setB])/sizeB)
  snr.obs.at<-signal.obs.at/noise.obs.at
  
  obs.setA<-intersect(cp.obs.id,low.iv.ind);sizeA=length(obs.setA)
  obs.setB<-intersect(cp.obs.id,high.iv.ind);sizeB=length(obs.setB)
  
  signal.obs.cp<-mean(yobs[obs.setA])-mean(yobs[obs.setB])
  noise.obs.cp<-sqrt(var(yobs[obs.setA])/sizeA+var(yobs[obs.setB])/sizeB)
  snr.obs.cp<-signal.obs.cp/noise.obs.cp
  
  #Mean Value for Principal Strata
  #Mean Value for Principal Strata
  mud<-Checking_Sample[j,"tinter"]+
    X.obs%*%Checking_Sample[j,paste("theta[", 1:ncov, "]",sep="")]+
    Checking_Sample[j,paste("t.tau[",Original_Data$rd.ind,"]",sep="")]     
  std<-1/sqrt(Checking_Sample[j,"taud.prec"]) 
  
  
  for (k in 1:rep.time){  
    #Generate Principal Strata
    simu.d<-rnorm(nsize,mean=mud,sd=std)
    
    #Determine Treatment Indicator
    treat.status<-(zvalues>simu.d)
    
    #Determine Potential Outcome
    #Latent Values
    y1<-X.obs%*%Checking_Sample[j,paste("beta1[", 1:ncov, "]",sep="")]+
      Checking_Sample[j,"gamma1"]*simu.d+Checking_Sample[j,paste("y1.tau[",Original_Data$rd.ind,"]",sep="")] 
    
    y0<-X.obs%*%Checking_Sample[j,paste("beta0[", 1:ncov, "]",sep="")]+
      Checking_Sample[j,"gamma0"]*simu.d+Checking_Sample[j,paste("y0.tau[",Original_Data$rd.ind,"]",sep="")] 
    
    ylatent=y1*treat.status+y0*(1-treat.status)
    
    #Determine Category
    cate_p1<-sapply(ylatent,FUN=function(x){logit(x+c1)})
    cate_p0<-sapply(ylatent,FUN=function(x){logit(x+c0)})
    cate_p=cate_p0
    deter_p<-matrix(0,nrow=nsize,ncol=res.level);
    deter_p[,res.level]=1-cate_p[res.level-1,]
    for (m in 2:(res.level-1))
    {
      deter_p[,m]=cate_p[m,]-cate_p[m-1,]
    }
    deter_p[,1]=cate_p[1,]
    
    #Sample Observe Outcome
    yobs_simu<-apply(deter_p,1,FUN=function(x){sample((6-res.level):5,1,prob=x)})
    
    #Replicate Data Measure Calculation
    nt.rep.id<-which(simu.d>zsupport[2])
    at.rep.id<-which(simu.d<zsupport[1])
    cp.rep.id<-which(simu.d<=zsupport[2]&simu.d>=zsupport[1])
    
    rep.setA<-intersect(nt.rep.id,low.iv.ind);sizeA=length(rep.setA)
    rep.setB<-intersect(nt.rep.id,high.iv.ind);sizeB=length(rep.setB)
    
    temp.signal.rep.nt[k]<-mean(yobs_simu[rep.setA])-mean(yobs_simu[rep.setB])
    temp.noise.rep.nt[k]<-sqrt(var(yobs_simu[rep.setA])/sizeA+var(yobs_simu[rep.setB])/sizeB)
    temp.snr.rep.nt[k]<- temp.signal.rep.nt[k]/ temp.noise.rep.nt[k]
    
    rep.setA<-intersect(at.rep.id,low.iv.ind);sizeA=length(rep.setA)
    rep.setB<-intersect(at.rep.id,high.iv.ind);sizeB=length(rep.setB)
    
    temp.signal.rep.at[k]<-mean(yobs_simu[rep.setA])-mean(yobs_simu[rep.setB])
    temp.noise.rep.at[k]<-sqrt(var(yobs_simu[rep.setA])/sizeA+var(yobs_simu[rep.setB])/sizeB)
    temp.snr.rep.at[k]<-temp.signal.rep.at[k]/temp.noise.rep.at[k]
    
    rep.setA<-intersect(cp.rep.id,low.iv.ind);sizeA=length(rep.setA)
    rep.setB<-intersect(cp.rep.id,high.iv.ind);sizeB=length(rep.setB)
    
    temp.signal.rep.cp[k]<-mean(yobs_simu[rep.setA])-mean(yobs_simu[rep.setB])
    temp.noise.rep.cp[k]<-sqrt(var(yobs_simu[rep.setA])/sizeA+var(yobs_simu[rep.setB])/sizeB)
    temp.snr.rep.cp[k]<-temp.signal.rep.cp[k]/temp.noise.rep.cp[k]
  }
  vare<-runif(1)
  #Parameter for Beta Distribution
  #Nevertakers
  a=sum(temp.signal.rep.nt>signal.obs.nt,na.rm=T)+
    vare*sum(temp.signal.rep.nt==signal.obs.nt,na.rm=T)
  b=sum(temp.signal.rep.nt<signal.obs.nt,na.rm=T)+
    (1-vare)*sum(temp.signal.rep.nt==signal.obs.nt,na.rm=T)
  signal.nt.sp<-rbeta(1,a+1,b+1)
  
  vare<-runif(1)
  a=sum(temp.noise.rep.nt>noise.obs.nt,na.rm=T)+
    vare*sum(temp.noise.rep.nt==noise.obs.nt,na.rm=T)
  b=sum(temp.noise.rep.nt<noise.obs.nt,na.rm=T)+
    (1-vare)*sum(temp.noise.rep.nt==noise.obs.nt,na.rm=T)
  noise.nt.sp<-rbeta(1,a+1,b+1)
  
  vare<-runif(1)
  a=sum(temp.snr.rep.nt>snr.obs.nt,na.rm=T)+
    vare*sum(temp.snr.rep.nt==snr.obs.nt,na.rm=T)
  b=sum(temp.snr.rep.nt<snr.obs.nt,na.rm=T)+
    (1-vare)*sum(temp.snr.rep.nt==snr.obs.nt,na.rm=T)
  snr.nt.sp<-rbeta(1,a+1,b+1)
  
  #AlwaysTakers
  vare<-runif(1)
  a=sum(temp.signal.rep.at>signal.obs.at,na.rm=T)+
    vare*sum(temp.signal.rep.at==signal.obs.at,na.rm=T)
  b=sum(temp.signal.rep.at<signal.obs.at,na.rm=T)+
    (1-vare)*sum(temp.signal.rep.at==signal.obs.at,na.rm=T)
  signal.at.sp<-rbeta(1,a+1,b+1)
  
  vare<-runif(1)
  a=sum(temp.noise.rep.at>noise.obs.at,na.rm=T)+
    vare*sum(temp.noise.rep.at==noise.obs.at,na.rm=T)
  b=sum(temp.noise.rep.at<noise.obs.at,na.rm=T)+
    (1-vare)*sum(temp.noise.rep.at==noise.obs.at,na.rm=T)
  noise.at.sp<-rbeta(1,a+1,b+1)
  
  vare<-runif(1)
  a=sum(temp.snr.rep.at>snr.obs.at,na.rm=T)+
    vare*sum(temp.snr.rep.at==snr.obs.at,na.rm=T)
  b=sum(temp.snr.rep.at<snr.obs.at,na.rm=T)+
    (1-vare)*sum(temp.snr.rep.at==snr.obs.at,na.rm=T)
  snr.at.sp<-rbeta(1,a+1,b+1)
  
  #Compliers
  vare<-runif(1)
  a=sum(temp.signal.rep.cp>signal.obs.cp,na.rm=T)+
    vare*sum(temp.signal.rep.cp==signal.obs.cp,na.rm=T)
  b=sum(temp.signal.rep.cp<signal.obs.cp,na.rm=T)+
    (1-vare)*sum(temp.signal.rep.cp==signal.obs.cp,na.rm=T)
  signal.cp.sp<-rbeta(1,a+1,b+1)
  
  vare<-runif(1)
  a=sum(temp.noise.rep.cp>noise.obs.cp,na.rm=T)+
    vare*sum(temp.noise.rep.cp==noise.obs.cp,na.rm=T)
  b=sum(temp.noise.rep.cp<noise.obs.cp,na.rm=T)+
    (1-vare)*sum(temp.noise.rep.cp==noise.obs.cp,na.rm=T)
  noise.cp.sp<-rbeta(1,a+1,b+1)
  
  vare<-runif(1)
  a=sum(temp.snr.rep.cp>snr.obs.cp,na.rm=T)+
    vare*sum(temp.snr.rep.cp==snr.obs.cp,na.rm=T)
  b=sum(temp.snr.rep.cp<snr.obs.cp,na.rm=T)+
    (1-vare)*sum(temp.snr.rep.cp==snr.obs.cp,na.rm=T)
  snr.cp.sp<-rbeta(1,a+1,b+1)

  return(c(signal.nt.sp,noise.nt.sp,snr.nt.sp,signal.cp.sp,noise.cp.sp,snr.cp.sp,
           signal.at.sp,noise.at.sp,snr.at.sp))
  },error=function(e){
    return(rep(NA,9))
  })
}

library(doParallel)
library(foreach)
num_Cores<-detectCores()
#Register Multiple Cluster
registerDoParallel(num_Cores)

#for (res.id in c(3,4,7)){
source("DataLoading.R")
#Store Sample
p_result<-vector("list",length=4)
for (k in 1:4){
#Extract MCMC Samples
Checking_Sample<-get(paste(response_index[res.id],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;
  
#Extract Original Data
Original_Data<-data.jags[[k]]
samplesize=nrow(Checking_Sample)
#Replicate Time Within Each Posterior Draw
rep.time=50

#Parallel Implementation
p_result[[k]]<-foreach(j=1:samplesize,.combine=rbind)%dopar%
sppv_comparison(j,rep.time,Checking_Sample,Original_Data)
assign(paste(response_index[res.id],k,"SPPV",sep="_"),p_result)
print(paste("==",res.id,"==finish",k))
save(p_result,file=paste(response_index[res.id],"_SPPV.RData"))
}
#}
stopImplicitCluster()



#Calculate Empirical P-value
response_index<-c("qq601","qq602","qq603","qq604","qq605","qq606",
                  "qm404","qk801","qk802","qk803","qk804")
pppv_result<-NULL
for(k in 1:4){
U<-runif(9,0,1)
assign(paste(response_index[res.id],k,"SPPV_NT",sep="_"),
       c(mean(p_result[[k]][,1]>U[1],na.rm=T),
         mean(p_result[[k]][,2]>U[2],na.rm=T),
         mean(p_result[[k]][,3]>U[3],na.rm=T)))

assign(paste(response_index[res.id],k,"SPPV_CP",sep="_"),
       c(mean(p_result[[k]][,4]>U[4],na.rm=T),
         mean(p_result[[k]][,5]>U[5],na.rm=T),
         mean(p_result[[k]][,6]>U[6],na.rm=T)))

assign(paste(response_index[res.id],k,"SPPV_AT",sep="_"),
       c(mean(p_result[[k]][,7]>U[7],na.rm=T),
         mean(p_result[[k]][,8]>U[8],na.rm=T),
         mean(p_result[[k]][,9]>U[9],na.rm=T)))

#Restore Discrepancy Measure PPPV
assign(paste(response_index[res.id],k,"SPPV",sep="_"),
       rbind(get(paste(response_index[res.id],k,"SPPV_NT",sep="_")),
             get(paste(response_index[res.id],k,"SPPV_CP",sep="_")),
             get(paste(response_index[res.id],k,"SPPV_AT",sep="_"))))
pppv_result<-rbind(pppv_result,
                   as.vector(get(paste(response_index[res.id],k,"SPPV",sep="_"))))
}
require(xtable)
xtable(pppv_result,digits=3)

