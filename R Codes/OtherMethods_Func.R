#Function for Result Comparison

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
  ATT_dr<-mean(treated.outcome)-sum((data[,response_index[res.id]]*data$p*(1-data$onechild)+
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


#A compounding function to fit OLS, IPW,DR
#Applying Bootstrap to get the variance
Fit_Boot_method<-function(data,res.id,B=250,model,choice_model,ylevel)
{
  #Direct Modeling
  att_reg_est<-Dir_ATT(model,data,ylevel)
  
  #Get Overlapping Data
  olp<-overlap_data(data)
  #IPW ATT
  att_ipw_est<-sum(olp[,response_index[res.id]]*olp$onechild)/sum(olp$onechild)-
    sum(olp[,response_index[res.id]]*olp$p*(1-olp$onechild)/(1-olp$p))/sum(olp$p*(1-olp$onechild)/(1-olp$p))
  #DR ATT
  att_dr_est<-DR_ATT(model,olp,ylevel)
  
  f1<-formula(choice_model)

  att_reg<-rep(NA,B)
  att_ipw<-rep(NA,B)
  att_dr<-rep(NA,B)
  for (i in 1:B)
  {
    #Bootstrap Data
    bdata<-fix_cluster_boot(data)
    tryCatch({
      suppressMessages(
        choice<-glm(f1,data=bdata,family=binomial(link="logit")))
      bdata$p<-as.vector(fitted(choice,data=bdata))
      #IPW Calculation
      bolp<-overlap_data(bdata)
      att_ipw[i]<-sum(bolp[,response_index[res.id]]*bolp$onechild)/sum(bolp$onechild)-
        sum(bolp[,response_index[res.id]]*bolp$p*(1-bolp$onechild)/(1-bolp$p))/
        sum(bolp$p*(1-bolp$onechild)/(1-bolp$p))
      
      #Direct Modeling
      bmodel_dir<-clmm2(y~meduy+feduy+age+onechild+han+onechild*meduy+
                          onechild*feduy+onechild*age+onechild*han,random=group,
                        data=bdata,link="logistic",control=clmm2.control(method="nlminb"))
      #Fitted ATT
      bylevel=length(unique(bdata$y))
      att_reg[i]<-Dir_ATT(bmodel_dir,bdata,bylevel)

      #ATT for Double Robust
      att_dr[i]<-DR_ATT(bmodel_dir,bolp,bylevel)
      if(abs(att_ipw[i])>1){att_ipw[i]=NA}
      if(abs(att_reg[i])>1){att_reg[i]=NA}
      if(abs(att_dr[i])>1){att_dr[i]=NA}
    },
    error=function(e){att_ipw[i]<-NA;att_reg[i]<-NA;att_dr[i]<-NA;
  })
    print(paste("==",i,"=="))
  }
  return(list(reg.est=att_reg_est,
              ipw.est=att_ipw_est,
              dr.est=att_dr_est,
              b_reg=att_reg,b_ipw=att_ipw,b_dr=att_dr))
}

