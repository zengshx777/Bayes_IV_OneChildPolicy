#Define Model
rr.model<-function()
{
  for (j in 1:nsample)
  {
    #Modeling for Principal strata:combination of covariates 
    #and random effect
    mud[j]=tinter+inprod(theta,X[j,])+t.tau[rd.ind[j]]
    #Draw Principal strata from normal distribution
    d[j]~dnorm(mud[j],1)
    #Compare the value of strata and IV to determine the treatment 
    #status
    determine[j]=-Z[j]-d[j]
    t[j]~dinterval(determine[j],0)
    
    #Process for Outcome Process
    #Latent Variable for cumulative logistic
    y1[j]<-inprod(beta1,X[j,])+gamma1*d[j]+y1.tau[rd.ind[j]]
    y0[j]<-inprod(beta0,X[j,])+gamma0*d[j]+y0.tau[rd.ind[j]]
    
    #Categorical Model
    #Caculate the prob in each cell with inverse logit function
    g1[j,1]=ilogit(c1[1]+y1[j])
    g1[j,2]=ilogit(c1[2]+y1[j])
    g1[j,3]=ilogit(c1[3]+y1[j])
    g1[j,4]=ilogit(c1[4]+y1[j])

    py1[j,1]=g1[j,1]
    py1[j,2]=g1[j,2]-g1[j,1]
    py1[j,3]=g1[j,3]-g1[j,2]
    py1[j,4]=g1[j,4]-g1[j,3]
    py1[j,5]=1-g1[j,4]
    #Draw potential outcome given the probability cell
    y1_obs[j]~dcat(py1[j,])
    
    #Similar Procedure for the control
    g0[j,1]=ilogit(c0[1]+y0[j])
    g0[j,2]=ilogit(c0[2]+y0[j])
    g0[j,3]=ilogit(c0[3]+y0[j])
    g0[j,4]=ilogit(c0[4]+y0[j])
    
    py0[j,1]=g0[j,1]
    py0[j,2]=g0[j,2]-g0[j,1]
    py0[j,3]=g0[j,3]-g0[j,2]
    py0[j,4]=g0[j,4]-g0[j,3]
    py0[j,5]=1-g0[j,4]
    y0_obs[j]~dcat(py0[j,])
    #Observed outcome is one of the realization.
    y[j]~dsum(y0_obs[j]*(1-t[j]),y1_obs[j]*t[j])
    
  }
  #Draw province random effect
  #Random effect in y1,y0 and principal strata
  for (i in 1:N)
  {
    y1.tau[i]~dnorm(0,y1.tauprec)
    y0.tau[i]~dnorm(0,y0.tauprec)
    t.tau[i]~dnorm(0,t.tauprec)
  }
  
  #Random Effect Precision Parameter Prior
  y1.tauprec~dgamma(1.0E-6,1.0E-6)
  y0.tauprec~dgamma(1.0E-6,1.0E-6)
  t.tauprec~dgamma(1.0E-6,1.0E-6)
  taud.prec~dgamma(1.0E-6,1.0E-6)
  
  #Coefficents for Covariates Prior
  theta~dmnorm(rep(0,ncov),RX)
  beta1~dmnorm(rep(0,ncov),RX)
  beta0~dmnorm(rep(0,ncov),RX)
  
  #Other Parameters Prior
  #Intercept for principal strata and outcome
  tinter~dnorm(0,1/100)
  gamma1~dnorm(0,1/100)
  gamma0~dnorm(0,1/100)
  
  #Threshold Prior
  c1~dmnorm(c(0,0,0,0),RT)
  c0~dmnorm(c(0,0,0,0),RT)
}

#Exact Same JAGS model except with only 4 categories.
#Maybe used if we only observe 4 categories in real data.
rr.model.den<-function()
{
  for (j in 1:nsample)
  {
    #Process for Treatment Assignment
    mud[j]=tinter+inprod(theta,X[j,])+t.tau[rd.ind[j]]
    d[j]~dnorm(mud[j], 1)
    determine[j]=-Z[j]-d[j]
    t[j]~dinterval(determine[j],0)
    
    #Process for Outcome Process
    y1[j]<-inprod(beta1,X[j,])+gamma1*d[j]+y1.tau[rd.ind[j]]
    y0[j]<-inprod(beta0,X[j,])+gamma0*d[j]+y0.tau[rd.ind[j]]
#    yobs_latent[j]=y1[j]*t[j]+y0[j]*(1-t[j])
    
    # #Categorical Model
    # yobs_cont[j]~dlogis(yobs_latent[j],1)
    # y[j]~dinterval(yobs_cont[j],c)
    g1[j,1]=ilogit(c1[1]+y1[j])
    g1[j,2]=ilogit(c1[2]+y1[j])
    g1[j,3]=ilogit(c1[3]+y1[j])
    
    py1[j,1]=g1[j,1]
    py1[j,2]=g1[j,2]-g1[j,1]
    py1[j,3]=g1[j,3]-g1[j,2]
    py1[j,4]=1-g1[j,3]
    y1_obs[j]~dcat(py1[j,])
    
    g0[j,1]=ilogit(c0[1]+y0[j])
    g0[j,2]=ilogit(c0[2]+y0[j])
    g0[j,3]=ilogit(c0[3]+y0[j])

    
    py0[j,1]=g0[j,1]
    py0[j,2]=g0[j,2]-g0[j,1]
    py0[j,3]=g0[j,3]-g0[j,2]
    py0[j,4]=1-g0[j,3]
    y0_obs[j]~dcat(py0[j,])
    y[j]~dsum(y0_obs[j]*(1-t[j]),y1_obs[j]*t[j])
  }
  
  for (i in 1:N)
  {
    y1.tau[i]~dnorm(0,y1.tauprec)
    y0.tau[i]~dnorm(0,y0.tauprec)
    t.tau[i]~dnorm(0,t.tauprec)
  }
  
  #Random Effect Precision
  y1.tauprec~dgamma(1.0E-6,1.0E-6)
  y0.tauprec~dgamma(1.0E-6,1.0E-6)
  t.tauprec~dgamma(1.0E-6,1.0E-6)
  taud.prec~dgamma(1.0E-6,1.0E-6)
  
  #Coefficents for Covariates
  theta~dmnorm(rep(0,ncov),RX)
  beta1~dmnorm(rep(0,ncov),RX)
  beta0~dmnorm(rep(0,ncov),RX)
  
  #Other Parameters
  tinter~dnorm(0,1/100)
  gamma1~dnorm(0,1/100)
  gamma0~dnorm(0,1/100)
  
  #Threshold 
  c1~dmnorm(c(0,0,0),RT)
  c0~dmnorm(c(0,0,0),RT)
}

#Parameters to Store After Sampling
parameters=c("gamma0","gamma1","theta","beta1","tinter","y1","y0",
             "beta0","d","t.tau","y1.tau","y0.tau","c1","c0","taud.prec")
