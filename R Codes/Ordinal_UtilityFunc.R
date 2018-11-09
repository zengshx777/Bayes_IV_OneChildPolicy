#This Scripts Include the Functions to Calculate Ordinal Estimands.
library(msm)

Largerthan<-function(a,b,strict)
{
  #Function Calculate the Probability of X>Y or X>=Y
  #a is the probability cell for X
  #b is the probability cell for Y
  #strict controls whether the inequality is strict  
  sum=numeric(nrow(a))
  if(strict==0)
  {
    for (k in 1:nrow(a))
    {
      for (i in 1:ncol(a))
      {
        sum[k]=sum[k]+sum(a[k,i]*b[k,i:ncol(a)])
      }
    }
  }
  else{
    for (k in 1:nrow(a))
    {
      for (i in 1:(ncol(a)-1))
      {
        sum[k]=sum[k]+sum(a[k,i]*b[k,(i+1):ncol(a)])
      }
    }
  }
  return(sum)
}

evalbin<-function(A,B,alpha)
{
  #This function returns the interval based on the posterior sample 
  #of sharp bounds. A is the upper bound, B is the lower bound
  #alpha is the target covarage probability
  t=numeric(ncol(A))
  for (k in 1:ncol(A))
  {
    for (z in seq(0,1,length=1000))
    {
      if (mean(B[,k]>(mean(B[,k])-z)&(A[,k]<mean(A[,k])+z))>=alpha)
      {
        t[k]=z
        break
      }
    }
  }
  return(t)
}


GetBound<-function(t,MCMCSample,data,treated.id,ncov,zmax,zmin,res.level,treat.ind=1,control.ind=0)
{
  #This function derive the PRTE and TT for a given posterior sample
  #It basically follows the same procedure in "Ordinal MTE" but perform it
  #for each posterior draw and take a weighted average
  
  tryCatch({
  #Use the empirical distribution of X and Threshold Value
  #Extract Threshold Values
  ps<-MCMCSample[t,grep("d",colnames(MCMCSample))]
  y1_mean<-MCMCSample[t,paste("beta1[",1:ncov,"]",sep="")]%*%t(data[,covariate_index])+
    MCMCSample[t,paste("y1.tau[",data$rd.ind,"]",sep="")]
  y0_mean<-MCMCSample[t,paste("beta0[",1:ncov,"]",sep="")]%*%t(data[,covariate_index])+
    MCMCSample[t,paste("y0.tau[",data$rd.ind,"]",sep="")]
  
  #Obtain the individual with principal strata within the range of IV
  #PRTE takes the average over those samples
  c.id<-which(zmin<ps&ps<zmax)
  Y_1_latent<-y1_mean[c.id]+MCMCSample[t,"gamma1"]*ps[c.id]
  Y_0_latent<-y0_mean[c.id]+MCMCSample[t,"gamma0"]*ps[c.id]
  
  
  
  
  
  p1_cell <- t(sapply(
    Y_1_latent,
    FUN = function(y) {
      logit(y+MCMCSample[t, paste("c1[", 1:(res.level - 1), "]", sep ="")])
    }
  ))
  
  p0_cell <- t(sapply(
    Y_0_latent,
    FUN = function(y) {
      logit(y+MCMCSample[t, paste("c0[", 1:(res.level - 1), "]", sep = "")])
    }
  ))
  
  mass_p1<-mass_p0<-matrix(0,nrow=nrow(p1_cell),ncol=res.level)
  mass_p1[,1]=p1_cell[,1];mass_p0[,1]=p0_cell[,1]
  for (k in 2:(res.level-1))
  {
    mass_p1[,k]=p1_cell[,k]-p1_cell[,k-1]
    mass_p0[,k]=p0_cell[,k]-p0_cell[,k-1]
  }
  mass_p1[,res.level]=1-p1_cell[,res.level-1]
  mass_p0[,res.level]=1-p0_cell[,res.level-1]
  aggr.p1=mass_p1;aggr.p0=mass_p0
  for (k in (res.level-1):1)
  {
    aggr.p1[,k]=aggr.p1[,k+1]+aggr.p1[,k]
    aggr.p0[,k]=aggr.p0[,k+1]+aggr.p0[,k]
  }
  #For Sharp Bounds Construction, we refer to the paper on 
  #treatment effect on ordinal outcome by Lu, Ding and Dasgupta
  #P(Y1>=Y0)
  #delta=p(y1>=j)-p(y0>=j)
  delta=aggr.p1-aggr.p0;delta[,1]=0
  late_tau_low<-mean(apply(delta+mass_p0,1,max))
  late_tau_up<-mean(apply(1+delta,1,min))
  late_tau_indi<-mean(Largerthan(mass_p0,mass_p1,0))
  late_eta_low<-mean(apply(delta,1,max))
  late_eta_up<-mean(1+apply(delta-mass_p1,1,min))
  late_eta_indi<-mean(Largerthan(mass_p0,mass_p1,1))
  
  #ATT Calculation
  #Generate Threshold from Truncate Normal
  #to approximate the distribution of principal strata for those being treated
  
  #The calculation procedure below is similar to the one for PRTE
  #ATT
  Y_1_latent<-y1_mean[treated.id]+MCMCSample[t,"gamma1"]*ps[treated.id]
  Y_0_latent<-y0_mean[treated.id]+MCMCSample[t,"gamma0"]*ps[treated.id]
  
  
  p1_cell <- t(sapply(
    Y_1_latent,
    FUN = function(y) {
      logit(y+MCMCSample[t, paste("c1[", 1:(res.level - 1), "]", sep ="")]+rep(treat.ind,res.level-1))
    }
  ))
  
  p0_cell <- t(sapply(
    Y_0_latent,
    FUN = function(y) {
      logit(y+MCMCSample[t, paste("c0[", 1:(res.level - 1), "]", sep = "")]+rep(control.ind,res.level-1))
    }
  ))
  
  mass_p1<-mass_p0<-matrix(0,nrow=nrow(p1_cell),ncol=res.level)
  mass_p1[,1]=p1_cell[,1];mass_p0[,1]=p0_cell[,1]
  for (k in 2:(res.level-1))
  {
    mass_p1[,k]=p1_cell[,k]-p1_cell[,k-1]
    mass_p0[,k]=p0_cell[,k]-p0_cell[,k-1]
  }
  mass_p1[,res.level]=1-p1_cell[,res.level-1]
  mass_p0[,res.level]=1-p0_cell[,res.level-1]
  aggr.p1=mass_p1;aggr.p0=mass_p0
  for (k in (res.level-1):1)
  {
    aggr.p1[,k]=aggr.p1[,k+1]+aggr.p1[,k]
    aggr.p0[,k]=aggr.p0[,k+1]+aggr.p0[,k]
  }  
  #For Sharp Bounds Construction, we refer to the paper on 
  #treatment effect on ordinal outcome by Lu, Ding and Dasgupta
  #P(Y1>=Y0)
  #delta=p(y1>=j)-p(y0>=j)
  delta=aggr.p1-aggr.p0;delta[,1]=0
  att_tau_low<-mean(apply(delta+mass_p0,1,max))
  att_tau_up<-mean(apply(1+delta,1,min))
  att_tau_indi<-mean(Largerthan(mass_p0,mass_p1,0))
  att_eta_low<-mean(apply(delta,1,max))
  att_eta_up<-mean(1+apply(delta-mass_p1,1,min))
  att_eta_indi<-mean(Largerthan(mass_p0,mass_p1,1))
  
  
  
  #Return the mean of values of interests.
  return(
    c(late_tau_low,
      late_tau_up,
      late_tau_indi,
      late_eta_low,
      late_eta_up,
      late_eta_indi,
      att_tau_low,
      att_tau_up,
      att_tau_indi,
      att_eta_low,
      att_eta_up,
      att_eta_indi)
  )} 
, error = function(e) {
  
  print(paste("==",t,"==","Find Errors"))
  return(rep(NA,12))
})
}
ncov=length(covariate_index)

#Evaluate the CI of Upper and Lower Bound
eval_credible<-function(A,B,alpha=0.95)
{
  return(c(quantile(A,1-(1-alpha)/2,na.rm=T),quantile(B,(1-alpha)/2,na.rm=T)))
}