#Get Average Function
library(msm)

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
