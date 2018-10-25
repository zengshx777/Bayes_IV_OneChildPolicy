source("Data_Reading.R")
library(Matching)
library(ordinal)
library(AER)

#Other Methods Comparison
source("OtherMethods_Func.R")


for (res.id in c(3,4,7))
{
  source("Data_Loading.R") 
  for (i in 1:4)
{
temp.data<-data.jags[[i]]
temp.data$y<-as.factor(temp.data[,response_index[res.id]])
temp.data$group<-as.factor(temp.data$rd.ind)
ylevel=length(unique(temp.data$y))

#ATT Estimation
tryCatch({
#Direct Modeling (Ordinal Probit)
#Outcome Formula
out.formula=as.formula("y~meduy+feduy+age+onechild+han+onechild*meduy+
                   onechild*feduy+onechild*age+onechild*han")

model_dir<-clmm2(out.formula,random=group,
                 data=temp.data,link="logistic",
                 control=clmm2.control(method="nlminb"))


#IPW
#Fitted First Stage Modeling
first_stage<-glm(onechild~age+meduy+feduy+han,family = binomial(link = "probit"),data=temp.data)
temp.data$p<-fitted(first_stage,type="response")

#Get the Point Estimation of OLS, IPW and Double Robust
REG_IPW_Result<-Fit_Boot_method(temp.data,res.id,B=250,model=model_dir,choice_model=first_stage,ylevel)

#Matching
m<-Match(Y=temp.data[,response_index[res.id]],M=5,
         Tr=temp.data$onechild,X=cbind(temp.data[,covariate_index],temp.data$p),
         estimand = "ATT")
Matching_est=m$est
Matching_se=m$se


#Using 2SLS to get LATE
temp.data$num.y=as.numeric(temp.data$y)
IVfit<-ivreg(num.y~meduy+feduy+age+onechild|age+meduy+feduy+han+ifppr,data=temp.data)
s<-summary(IVfit,vcov=sandwich)
LATE_EST<-s$coefficients["onechild",1]
LATE_SE<-s$coefficients["onechild",2]

#Point Estimation
Est_result<-c(REG_IPW_Result$reg.est,
              Matching_est,
              REG_IPW_Result$ipw.est,REG_IPW_Result$dr.est,
              LATE_EST)

#Standard Deviation
Se_result<-c(sd(REG_IPW_Result$b_reg,na.rm=T),
             Matching_se,
             sd(REG_IPW_Result$b_ipw,na.rm=T),
             sd(REG_IPW_Result$b_dr,na.rm=T),
              LATE_SE)

#Upper Credible Interval
UCI_result<-c(quantile(REG_IPW_Result$b_reg,0.975,na.rm=T),
              Matching_est+qnorm(0.975)*Matching_se,
              quantile(REG_IPW_Result$b_ipw,0.975,na.rm=T),
              quantile(REG_IPW_Result$b_dr,0.975,na.rm=T),
              LATE_EST+qnorm(0.975)*LATE_SE) 

#Lower Credible Interval
LCI_result<-c(quantile(REG_IPW_Result$b_reg,0.025,na.rm=T),
              Matching_est-qnorm(0.975)*Matching_se,
              quantile(REG_IPW_Result$b_ipw,0.025,na.rm=T),
              quantile(REG_IPW_Result$b_dr,0.025,na.rm=T),
              LATE_EST-qnorm(0.975)*LATE_SE) 
},
error=function(e){
  print("==ERROR OCCUR==")
  Est_result<-rep(NA,5);Se_result<-rep(NA,5)
  LCI_result<-rep(NA,5);UCI_result<-rep(NA,5)
})
#Store Values
print(paste("==",res.id,"==",i,"==",seq=""))
assign(paste(response_index[res.id],i,"Estimation",seq="_"),Est_result)
assign(paste(response_index[res.id],i,"SE",seq="_"),Se_result)
assign(paste(response_index[res.id],i,"UCI",seq="_"),UCI_result)
assign(paste(response_index[res.id],i,"LCI",seq="_"),LCI_result)
assign(paste(response_index[res.id],i,"Bootsample",seq="_"),REG_IPW_Result)
}
}


save.image("Other Methods Results.RData")

#Result Extraction to Produce Table 6
T_MATRIX<-NULL
for (i in 1:4)
{
  R_MATRIX<-NULL
  for(res.id in c(7,3,4))
  {
    R_MATRIX<-cbind(R_MATRIX,
    t(rbind(get(paste(response_index[res.id],i,"Estimation",seq="_")),
            get(paste(response_index[res.id],i,"LCI",seq="_")),
            get(paste(response_index[res.id],i,"UCI",seq="_")))))
  }
  row.names(R_MATRIX)=c("OLS","Matching","IPW","Double Robust","2SLS")
  T_MATRIX<-rbind(T_MATRIX,R_MATRIX)
}
library(xtable)
xtable(T_MATRIX,digits=3)