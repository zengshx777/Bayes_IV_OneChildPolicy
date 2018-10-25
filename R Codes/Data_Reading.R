#This Script preprocess the data for subsequent analysis
load("CFPS_2010_Adult_Demo.RData")
#read the data of instrument variable
ifppr<-read.table("IFPPR.txt",header=TRUE)
#Unify the name of province variable
ifppr$provcd_born=ifppr$provcd

#Extract Range
start=1979;end=2010;iend=140

#combine the data set
merge<-merge(adult2010,ifppr,by=c("provcd_born"))
full_data<-subset(merge,ifppr>0&meduy>=0&feduy>=0&
                    age>=0&gender>=0&urban_3>=0)

#Extract the data within the range specified
data1<-subset(full_data,age>=(2010-end)&age<=(2010-start)&ifppr<iend)

#Divide Data into 4 Subgroups by Hukou Status and Gender
data_00=subset(data1,urban_3==0&gender==0)
data_10=subset(data1,urban_3==1&gender==0)
data_01=subset(data1,urban_3==0&gender==1)
data_11=subset(data1,urban_3==1&gender==1)

#Codebook for Variable Name
#Notice that sampling weight (national combined weight) is included
#covariate_index<-c("meduy","feduy","age","han","rswt_nat")
covariate_index<-c("meduy","feduy","age","han")
response_index<-c("qq601","qq602","qq603","qq604","qq605","qq606",
                  "qm404","qk801","qk802","qk803","qk804")
treatment_index<-"onechild"
instrument_index<-"ifppr"
cluster_index<-"provcd_born"

#Create a List to Combine Differ Data Set
data.combine<-list(data_00,data_01,data_10,data_11)
data.jags=data.combine

