#Produce Figure 1 in the paper, Heuristic Comparison
library(ggplot2)
load("Figure.RData")
#load("Figure1.RData")
covariate_index<-c("meduy","feduy","age","han")
response_index<-c("qq601","qq602","qq603","qq604","qq605","qq606",
                  "qm404","qk801","qk802","qk803","qk804")

#Summary Stats for Rural Females
data=data_00
data.t<-subset(data,onechild==1)
data.c<-subset(data,onechild==0)
dif1<-apply(data.t[,c(covariate_index,response_index[c(7,3,4)])],
      2,FUN=function(x){mean(x[x>=0],na.rm=T)})-
apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
      2,FUN=function(x){mean(x[x>=0],na.rm=T)})
sd1<-sqrt(apply(data.t[,c(covariate_index,response_index[c(7,3,4)])],
           2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.t)+
  apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
        2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.c))
ASD_1<-dif1/sd1
names=c("MEDUY","FEDUY","AGE","ETHINICITY","CONFIDENCE","ANXIETY","DESPERATION")
names<-factor(names,levels=names)
ASD_1<-data.frame(asd=ASD_1,xnames=names,var=c(rep(1,4),rep(2,3)))

#Summary Stats for Rural Males
data=data_01
data.t<-subset(data,onechild==1)
data.c<-subset(data,onechild==0)
dif1<-apply(data.t[,c(covariate_index,response_index[c(7,3,4)])],
            2,FUN=function(x){mean(x[x>=0],na.rm=T)})-
  apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
        2,FUN=function(x){mean(x[x>=0],na.rm=T)})
sd1<-sqrt(apply(data.t[,c(covariate_index,response_index[c(7,3,4)])],
                2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.t)+
            apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
                  2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.c))
ASD_2<-dif1/sd1
names=c("MEDUY","FEDUY","AGE","ETHINICITY","CONFIDENCE","ANXIETY","DESPERATION")
names<-factor(names,levels=names)
ASD_2<-data.frame(asd=ASD_2,xnames=names,var=c(rep(1,4),rep(2,3)))

#Summary Stats for Urban Females
data=data_10
data.t<-subset(data,onechild==1)
data.c<-subset(data,onechild==0)
dif1<-apply(data.t[,c(covariate_index,response_index[c(7,3,3)])],
            2,FUN=function(x){mean(x[x>=0],na.rm=T)})-
  apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
        2,FUN=function(x){mean(x[x>=0],na.rm=T)})
sd1<-sqrt(apply(data.t[,c(covariate_index,response_index[c(7,3,3)])],
                2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.t)+
            apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
                  2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.c))
ASD_3<-dif1/sd1
names=c("MEDUY","FEDUY","AGE","ETHINICITY","CONFIDENCE","ANXIETY","DESPERATION")
names<-factor(names,levels=names)
ASD_3<-data.frame(asd=ASD_3,xnames=names,var=c(rep(1,4),rep(2,3)))

#Summary Stats for Urban Males
data=data_11
data.t<-subset(data,onechild==1)
data.c<-subset(data,onechild==0)
dif1<-apply(data.t[,c(covariate_index,response_index[c(7,3,4)])],
            2,FUN=function(x){mean(x[x>=0],na.rm=T)})-
  apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
        2,FUN=function(x){mean(x[x>=0],na.rm=T)})
sd1<-sqrt(apply(data.t[,c(covariate_index,response_index[c(7,3,4)])],
                2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.t)+
            apply(data.c[,c(covariate_index,response_index[c(7,3,4)])],
                  2,FUN=function(x){var(x[x>=0],na.rm=T)})/nrow(data.c))
ASD_4<-dif1/sd1
names=c("MEDUY","FEDUY","AGE","ETHINICITY","CONFIDENCE","ANXIETY","DESPERATION")
names<-factor(names,levels=names)
ASD_4<-data.frame(asd=ASD_4,xnames=names,var=c(rep(1,4),rep(2,3)))

#Combine the Data Set
ASD_Combine<-rbind(ASD_1,ASD_2,ASD_3,ASD_4)
ASD_Combine$ruralgroup=rep(c("Rural Female","Rural Male",
                             "Urban Female","Urban Male"),each=7)
#assign(paste("p",k),
p1<-ggplot(data = ASD_Combine,
           aes(
             x = xnames,
             y = asd,
             shape = factor(var),
             ymin = asd,
             ymax = asd
           )) + geom_pointrange(fatten = 8) + ylim(c(-10, 10)) +
  facet_grid(.~ruralgroup)+
  geom_hline(yintercept = c(-1.96, 0, 1.96), lty = 2) + coord_flip()  + xlab("")+
  ylab("Standard Difference") +scale_shape_discrete(solid=F)+
   theme_bw()+theme(plot.title = element_text(hjust=0.5),
                                            legend.position = "none")


#pdf("DescriptiveASD.pdf",height=5.5,width=16)
p1
#dev.off()
