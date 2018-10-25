#This Script reformats the data frame to the one for JAGS Sampling
#Need to supply res.id:which response to examine
for (m in 1:4)
{
#Create dataframe for each subgroup
data.try<-data.combine[[m]][,c(response_index[res.id],treatment_index,
                     covariate_index,instrument_index,cluster_index)]
#Extract Complete Data
data.try<-data.try[complete.cases(data.try),]
data.try<-subset(data.try,data.try[,response_index[res.id]]>0)
#extract the cluster with size greater than 10
data.try<-subset(data.try,provcd_born%in%names(table(data.try$provcd_born))[
  table(data.try$provcd_born)>10])
data.try$province<-factor(data.try$provcd_born)
#Sort by Province
data.try<-data.try[order(data.try$province),]

#Create Cluster Indicator Variale
rd.ind<-NULL
rd.info<-table(data.try$province)
for (k in 1:length(rd.info))
{
  rd.ind<-c(rd.ind,rep(k,rd.info[k]))
}
data.try$rd.ind<-rd.ind

#Normalize the Weight
#data.try$rswt_nat<-data.try$rswt_nat/mean(data.try$rswt_nat)

data.jags[[m]]=data.try
print(paste("==",m,"=="))
}
