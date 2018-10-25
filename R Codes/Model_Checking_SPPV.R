#Model Checking Template Sample Predictive P-value \
#Produce Table 5 in the paper
#This is non-paralleled version
#Extract MCMC Samples
#Example
load("qq604_result_monotonic_weight_model.RData")
ncov=length(covariate_index)

for (k in 1:4){
Checking_Sample<-get(paste(response_index[res.id],k,"JAGSObject",sep="_"))$BUGSoutput$sims.matrix;

#Extract Original Data
Original_Data<-data.jags[[k]]

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

samplesize=nrow(Checking_Sample)

#Adjusted p value
signal.nt.sp<-signal.at.sp<-signal.cp.sp<-numeric(samplesize)
noise.nt.sp<-noise.at.sp<-noise.cp.sp<-numeric(samplesize)
snr.nt.sp<-snr.at.sp<-snr.cp.sp<-numeric(samplesize)


signal.obs.nt<-signal.obs.at<-signal.obs.cp<-numeric(samplesize)
noise.obs.nt<-noise.obs.at<-noise.obs.cp<-numeric(samplesize)
snr.obs.nt<-snr.obs.at<-snr.obs.cp<-numeric(samplesize)

#Replicate Time Within Each Posterior Draw
rep.time=50
#Replicate Once for each Posterior Draw
for (j in 1:samplesize){
  
  #Temp Discrepancy Restore
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
  
  signal.obs.nt[j]<-mean(yobs[obs.setA])-mean(yobs[obs.setB])
  noise.obs.nt[j]<-sqrt(var(yobs[obs.setA])/sizeA+var(yobs[obs.setB])/sizeB)
  snr.obs.nt[j]<-signal.obs.nt[j]/noise.obs.nt[j]
  
  obs.setA<-intersect(at.obs.id,low.iv.ind);sizeA=length(obs.setA)
  obs.setB<-intersect(at.obs.id,high.iv.ind);sizeB=length(obs.setB)
  
  signal.obs.at[j]<-mean(yobs[obs.setA])-mean(yobs[obs.setB])
  noise.obs.at[j]<-sqrt(var(yobs[obs.setA])/sizeA+var(yobs[obs.setB])/sizeB)
  snr.obs.at[j]<-signal.obs.at[j]/noise.obs.at[j]
  
  obs.setA<-intersect(cp.obs.id,low.iv.ind);sizeA=length(obs.setA)
  obs.setB<-intersect(cp.obs.id,high.iv.ind);sizeB=length(obs.setB)
  
  signal.obs.cp[j]<-mean(yobs[obs.setA])-mean(yobs[obs.setB])
  noise.obs.cp[j]<-sqrt(var(yobs[obs.setA])/sizeA+var(yobs[obs.setB])/sizeB)
  snr.obs.cp[j]<-signal.obs.cp[j]/noise.obs.cp[j]
  
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
  cate_p[,which(treat.status==1)]=cate_p1[,which(treat.status==1)]
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
  a=sum(temp.signal.rep.nt>signal.obs.nt[j],na.rm=T)+
    vare*sum(temp.signal.rep.nt==signal.obs.nt[j],na.rm=T)
  b=sum(temp.signal.rep.nt<signal.obs.nt[j],na.rm=T)+
    (1-vare)*sum(temp.signal.rep.nt==signal.obs.nt[j],na.rm=T)
  signal.nt.sp[j]<-rbeta(1,a+1,b+1)
  
  a=sum(temp.noise.rep.nt>noise.obs.nt[j],na.rm=T)+
    vare*sum(temp.noise.rep.nt==noise.obs.nt[j],na.rm=T)
  b=sum(temp.noise.rep.nt<noise.obs.nt[j],na.rm=T)+
    (1-vare)*sum(temp.noise.rep.nt==noise.obs.nt[j],na.rm=T)
  noise.nt.sp[j]<-rbeta(1,a+1,b+1)
  
  a=sum(temp.snr.rep.nt>snr.obs.nt[j],na.rm=T)+
    vare*sum(temp.snr.rep.nt==snr.obs.nt[j],na.rm=T)
  b=sum(temp.snr.rep.nt<snr.obs.nt[j],na.rm=T)+
    (1-vare)*sum(temp.snr.rep.nt==snr.obs.nt[j],na.rm=T)
  snr.nt.sp[j]<-rbeta(1,a+1,b+1)

  #AlwaysTakers
  a=sum(temp.signal.rep.at>signal.obs.at[j],na.rm=T)+
    vare*sum(temp.signal.rep.at==signal.obs.at[j],na.rm=T)
  b=sum(temp.signal.rep.at<signal.obs.at[j],na.rm=T)+
    (1-vare)*sum(temp.signal.rep.at==signal.obs.at[j],na.rm=T)
  signal.at.sp[j]<-rbeta(1,a+1,b+1)
  
  a=sum(temp.noise.rep.at>noise.obs.at[j],na.rm=T)+
    vare*sum(temp.noise.rep.at==noise.obs.at[j],na.rm=T)
  b=sum(temp.noise.rep.at<noise.obs.at[j],na.rm=T)+
    (1-vare)*sum(temp.noise.rep.at==noise.obs.at[j],na.rm=T)
  noise.at.sp[j]<-rbeta(1,a+1,b+1)
  
  a=sum(temp.snr.rep.at>snr.obs.at[j],na.rm=T)+
    vare*sum(temp.snr.rep.at==snr.obs.at[j],na.rm=T)
  b=sum(temp.snr.rep.at<snr.obs.at[j],na.rm=T)+
    (1-vare)*sum(temp.snr.rep.at==snr.obs.at[j],na.rm=T)
  snr.at.sp[j]<-rbeta(1,a+1,b+1)
  
  #Compliers
  a=sum(temp.signal.rep.cp>signal.obs.cp[j],na.rm=T)+
    vare*sum(temp.signal.rep.cp==signal.obs.cp[j],na.rm=T)
  b=sum(temp.signal.rep.cp<signal.obs.cp[j],na.rm=T)+
    (1-vare)*sum(temp.signal.rep.cp==signal.obs.cp[j],na.rm=T)
  signal.cp.sp[j]<-rbeta(1,a+1,b+1)
  
  a=sum(temp.noise.rep.cp>noise.obs.cp[j],na.rm=T)+
    vare*sum(temp.noise.rep.cp==noise.obs.cp[j],na.rm=T)
  b=sum(temp.noise.rep.cp<noise.obs.cp[j],na.rm=T)+
    (1-vare)*sum(temp.noise.rep.cp==noise.obs.cp[j],na.rm=T)
  noise.cp.sp[j]<-rbeta(1,a+1,b+1)
  
  a=sum(temp.snr.rep.cp>snr.obs.cp[j],na.rm=T)+
    vare*sum(temp.snr.rep.cp==snr.obs.cp[j],na.rm=T)
  b=sum(temp.snr.rep.cp<snr.obs.cp[j],na.rm=T)+
    (1-vare)*sum(temp.snr.rep.cp==snr.obs.cp[j],na.rm=T)
  snr.cp.sp[j]<-rbeta(1,a+1,b+1)
  print(paste("==",j,"=="))

#Calculate Empirical P-value

assign(paste(response_index[res.id],k,"SPPV_NT",sep="_"),
       c(mean(signal.nt.sp>runif(1)),
         mean(noise.nt.sp>runif(1)),
         mean(snr.nt.sp>runif(1))))

assign(paste(response_index[res.id],k,"SPPV_AT",sep="_"),
       c(mean(signal.at.sp>runif(1)),
         mean(noise.at.sp>runif(1)),
         mean(snr.at.sp>runif(1))))

assign(paste(response_index[res.id],k,"SPPV_CP",sep="_"),
       c(mean(signal.cp.sp>runif(1)),
         mean(noise.cp.sp>runif(1)),
         mean(snr.cp.sp>runif(1))))

#Restore Discrepancy Measure PPPV
assign(paste(response_index[res.id],k,"SPPV",sep="_"),
       rbind(get(paste(response_index[res.id],k,"SPPV_NT",sep="_")),
             get(paste(response_index[res.id],k,"SPPV_AT",sep="_")),
             get(paste(response_index[res.id],k,"SPPV_CP",sep="_"))))
}}