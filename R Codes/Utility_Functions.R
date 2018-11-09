library(msm)
#Utility Functions to Calculate Estimands on Average 
GetAverage<-function(t,MCMCSample,data,treated.id,ncov,zmax,zmin,res.level,treat.ind=1,control.ind=0)
{
  tryCatch({
    #Use the empirical distribution of X and Threshold Value
    #Extract Threshold Values
    ps<-MCMCSample[t,grep("d",colnames(MCMCSample))]
    y1_mean<-MCMCSample[t,paste("beta1[",1:ncov,"]",sep="")]%*%t(data[,covariate_index])+
      MCMCSample[t,paste("y1.tau[",data$rd.ind,"]",sep="")]
    y0_mean<-MCMCSample[t,paste("beta0[",1:ncov,"]",sep="")]%*%t(data[,covariate_index])+
      MCMCSample[t,paste("y0.tau[",data$rd.ind,"]",sep="")]
    
    #Obtain the individual with principal strata within the range of IV
    c.id<-which(zmin<ps&ps<zmax)
    Y_1_latent<-y1_mean[c.id]+MCMCSample[t,"gamma1"]*ps[c.id]
    Y_0_latent<-y0_mean[c.id]+MCMCSample[t,"gamma0"]*ps[c.id]
    
    Y_1_expect <-
      sapply(
        Y_1_latent,
        FUN = function(x) {
          - sum(sapply(MCMCSample[t, paste("c1[", 1:(res.level - 1), "]", sep = "")]+rep(treat.ind,res.level-1),
                       FUN = function(y) { logit(y + x)}
          ))
        }
      )
    
    Y_0_expect <-
      sapply(
        Y_0_latent,
        FUN = function(x) {
          - sum(sapply(
            MCMCSample[t, paste("c0[", 1:(res.level - 1), "]", sep = "")]+rep(control.ind,res.level-1),
            FUN = function(y) {
              logit(y + x)
            }
          ))
        }
      )
    prte=mean(Y_1_expect-Y_0_expect)
    #ATT
    Y_1_latent<-y1_mean[treated.id]+MCMCSample[t,"gamma1"]*ps[treated.id]
    Y_0_latent<-y0_mean[treated.id]+MCMCSample[t,"gamma0"]*ps[treated.id]
    
    Y_1_expect <-
      sapply(
        Y_1_latent,
        FUN = function(x) {
          -sum(sapply(
            MCMCSample[t, paste("c1[", 1:(res.level - 1), "]", sep = "")]+rep(treat.ind,res.level-1),
            FUN = function(y) {
              logit(y + x)
            }
          ))
        }
      )
    
    Y_0_expect <-
      sapply(
        Y_0_latent,
        FUN = function(x) {
          -sum(sapply(
            MCMCSample[t, paste("c0[", 1:(res.level - 1), "]", sep = "")]+rep(control.ind,res.level-1),
            FUN = function(y) {
              logit(y + x)
            }
          ))
        }
      )
    att=mean(Y_1_expect-Y_0_expect)
    
    return(c(prte,att))
  },error = function(e) {print(paste("==",t,"==","With Error")) 
    return(c(NA,NA))}  )
}


#Utility Functions to Calculate Ordinal Estimands.


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





#Utility Function for Other Methods Result Comparison

#Get Fitted Value
predict.expout<-function(model,data,ylevel)
{
  p_fitted<-matrix(0,nrow=nrow(data),ncol=ylevel)
  fitdata<-data
  for (k in (6-ylevel):5)
  {
    fitdata$y=factor(k,levels=1:ylevel)
    p_fitted[,k]=predict(model,newdata=fitdata)
  }
  return(apply(p_fitted,1,FUN=function(x){sum(x*((6-ylevel):5))}))
}

#ATT
#Obtain ATT from Outcome Modeling
#Temp.data is the input population Data
#ylevel is the number of unique values for the outcome
Dir_ATT<-function(model,temp.data,ylevel){
  treated.temp.data<-subset(temp.data,temp.data$onechild==1)
  treated.outcome=as.numeric(treated.temp.data$y)
  treated.temp.data$onechild=0
  impute<-predict.expout(model,treated.temp.data,ylevel)
  ATT_dir<-mean(treated.outcome-impute)
  return(ATT_dir)
}

#Get Overlapped Data, with Propensity Score Estimated
overlap_data<-function(data)
{
  low_p<-max(min(data[data$onechild==1,]$p),min(data[data$onechild==0,]$p))
  up_p<-min(max(data[data$onechild==1,]$p),max(data[data$onechild==0,]$p))
  subset(data,p>=low_p&p<=up_p)
}

#Double Robust Calculation
#Augmented model with the fitted propensity score
DR_ATT<-function(model,data,ylevel){
  #temp.data should contain fitted propensity score
  treated.temp.data<-subset(data,onechild==1)
  treated.outcome=as.numeric(treated.temp.data$y)
  control.data<-data;control.data$onechild=0
  impute<-predict.expout(model,control.data,ylevel)
  ATT_dr<-mean(treated.outcome)-sum((data[,response_index[res]]*data$p*(1-data$onechild)+
                                       impute*(data$onechild-data$p))/(1-data$p))/length(treated.outcome)   
  return(ATT_dr)
}

#nonparametric Bootstrap
#Bootstrap within each cluster, replace sampling within each cluster
#Keep the number of clusters, cluster size fixed
fix_cluster_boot<-function(data){
  group_size<-table(data$group)
  bdata<-NULL
  for (k in 1:length(group_size))
  {
    bdata<-rbind(bdata,subset(data,group==k)[
      sample(1:group_size[k],group_size[k],replace=T),])
  }
  return(bdata)
}
#Two Stages Cluster Sample
two_stage_cluster_boot<-function(data)
{
  group_size<-table(data$group)
  #Sample Cluster according to cluster size
  cluster.ind=sample(1:length(group_size),length(group_size),
                     prob=group_size,replace=T)
  bdata<-NULL
  #Within Each Cluster, Sample with Replacement
  for(k in cluster.ind)
  {
    bdata<-rbind(bdata,subset(data,group==k)[
      sample(1:group_size[k],group_size[k],replace=T),])
  }
  return(bdata)
}
#A compounding function to fit OLS, IPW,DR
#Applying Bootstrap to get the variance
Fit_Boot_method<-function(data,res,B=250,model,choice_model,ylevel)
{
  
  
  
  #Get Overlapping Data
  olp<-overlap_data(data)
  
  #Regression ATT
  att_reg_est<-Dir_ATT(model,olp,ylevel)
  
  #Matching ATT
  m<-Match(Y=data[,response_index[res]],M=5,
           Tr=data$onechild,X=data$p,
           estimand = "ATT")
  att_matching_est=m$est
  
  #IPW ATT
  att_ipw_est<-sum(olp[,response_index[res]]*olp$onechild)/sum(olp$onechild)-
    sum(olp[,response_index[res]]*olp$p*(1-olp$onechild)/(1-olp$p))/sum(olp$p*(1-olp$onechild)/(1-olp$p))
  
  #DR ATT
  att_dr_est<-DR_ATT(model,olp,ylevel)
  
  #IV
  data$num.y=as.numeric(data$y)
  IVfit<-ivreg(num.y~meduy+feduy+age+onechild|age+meduy+feduy+han+ifppr,data=data)
  iv_est=IVfit$coefficients["onechild"]
  
  f1<-formula(choice_model)
  
  att_reg<-rep(NA,B)
  att_matching<-rep(NA,B)
  att_ipw<-rep(NA,B)
  att_dr<-rep(NA,B)
  att_iv<-rep(NA,B)
  print("==Bootstrap Start==")
  for (i in 1:B)
  {
    #Bootstrap Data
    #bdata<-fix_cluster_boot(data)
    bdata<-two_stage_cluster_boot(data)
    tryCatch({
      suppressMessages(
        #Refit the Propensity Score Model
        choice<-glm(f1,data=bdata,family=binomial(link="logit")))
      bdata$p<-as.vector(fitted(choice,data=bdata))
      
      #Matching
      att_matching[i]<-Match(Y=bdata[,response_index[res]],M=5,
                             Tr=bdata$onechild,X=bdata$p,
                             estimand = "ATT")$est
      
      #IPW Calculation
      bolp<-overlap_data(bdata)
      att_ipw[i]<-sum(bolp[,response_index[res]]*bolp$onechild)/sum(bolp$onechild)-
        sum(bolp[,response_index[res]]*bolp$p*(1-bolp$onechild)/(1-bolp$p))/
        sum(bolp$p*(1-bolp$onechild)/(1-bolp$p))
      
      #Direct Modeling
      bmodel_dir<-clmm2(y~meduy+feduy+age+onechild+han+onechild*meduy+
                          onechild*feduy+onechild*age+onechild*han,random=group,
                        data=bdata,link="logistic",control=clmm2.control(method="nlminb"))
      #Fitted ATT
      bylevel=length(unique(bdata$y))
      att_reg[i]<-Dir_ATT(bmodel_dir,bdata,bylevel)
      
      #ATT for Matching
      
      #ATT for Double Robust
      att_dr[i]<-DR_ATT(bmodel_dir,bolp,bylevel)
      
      #Est for IV
      b_IVfit<-ivreg(num.y~meduy+feduy+age+onechild|age+meduy+feduy+han+ifppr,data=bdata)
      att_iv[i]=b_IVfit$coefficients["onechild"]
      
      #      if(abs(att_ipw[i])>1){att_ipw[i]=NA}
      #      if(abs(att_reg[i])>1){att_reg[i]=NA}
      #      if(abs(att_dr[i])>1){att_dr[i]=NA}
      #      if(abs(att_matching[i])>1){att_dr[i]=NA}
      if(abs(att_iv[i])>1){att_dr[i]=NA}
    },error=function(e){att_ipw[i]<-NA;att_reg[i]<-NA;att_dr[i]<-NA;att_matching[i]=NA})
    print(paste("==",i,"=="))
  }
  return(list(reg.est=att_reg_est,
              match.est=att_matching_est,
              ipw.est=att_ipw_est,
              dr.est=att_dr_est,
              iv.est=iv_est,
              b_reg=att_reg,b_ipw=att_ipw,b_dr=att_dr,
              b_matching=att_matching,
              b_iv=att_iv))
}

