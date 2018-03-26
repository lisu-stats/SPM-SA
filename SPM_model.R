###############################################################
#### WinBUGS code to fit the SPM to the HERS data
#### author: Li Su; Qiuju Li


model
{
  
  for(i in 1:N)  ## N is the total number of longitudinal measurements;
   {
    lonyt[i]~dnorm(mu[i],tauepsilon)
    
    mu[i]<-lonalpha[1]+lonalpha[2]*lonx1[i]+lonalpha[3]*lonx2[i]+lonalpha[4]*lonx3[i]+lonalpha[5]*lonx4[i]+lonalpha[6]*lonx5[i]+lonalpha[7]*lonx6[i]+lonalpha[8]*lonx7[i]+lonalpha[9]*lonx8[i]+lonalpha[10]*lonx9[i]+lonalpha[11]*lonx10[i]+lonalpha[12]*lonx11[i]+randombi[numsubjid[i],1]+randombi[numsubjid[i],2]*lonx1[i]
    }
  
  
   ### drop-out  process;
   for(j in 1:M) ### M is the number of subjects;
   {
     
     ### drop-out process;
      onesD[j]<-1                ### ones trick;
      for(kd in 1:numK)
      {
        quantD[j,kd]<-phi(lambdaD2[j,kd])
       ##### NB: 1-step() to avoid  numerical problems in probit function
       lambdaD2[j,kd]<-lambdaD[j,kd]*step(lambdaD[j,kd]+7)*step(7-lambdaD[j,kd])+(-7)*(1-step(lambdaD[j,kd]+7))+7*(1-step(7-lambdaD[j,kd])) 
       lambdaD[j,kd]<-(survalphaD[1]+survalphaD[2]*survDx1[j]+survalphaD[3]*survDx2[j]+survalphaD[4]*survDx3[j]+survalphaD[5]*survDx4[j]+survalphaD[6]*survDx5[j]+survalphaD[7]*survDvisit[kd]+survalphaD[8]*pow(survDvisit[kd],2)+survgammaD[1]*randombi[j,1]+survgammaD[2]*randombi[j,2])
       }
       probD[j]<-prod(quantD[j,1:survD[j]])*pow((1-quantD[j,survD[j]]),survdeltaD[j])/pow(quantD[j,(survD[j])],survdeltaD[j])*survdeltaD[j]+(1-survdeltaD[j])*prod(quantD[j,1:11])
  
      onesD[j]~dbern(probD[j])
  
     
  
  #### random effects;
  randombi[j,1]~dnorm(0,taub1)  ## random intercept;
  randombi[j,2]<-theta12*randombi[j,1]+e2[j]  ## random slope: var(b2)=theta21^2*var(b1)+var(e2)
  e2[j]~dnorm(0,taub2)
}


#####################################
### prior;
tauepsilon~dgamma(0.001,0.001)
### lonalpha;
for(p1 in 1:lonalphaP)
{lonalpha[p1]~dnorm(0,1.0E-2)}
### survD-alpha;
for(p2 in 1:survalphaDP)
{survalphaD[p2]~dnorm(0,0.25)}
### survD-gamma;
for(p3 in 1:survgammaDP)
{survgammaD[p3]~dnorm(0,0.25)}



taub1<-1/(sigmab11*sigmab11)
taub2<-1/(sigmab22*sigmab22)
sigmab11~dunif(0,5)
sigmab22~dunif(0,5)
theta12~dnorm(0,1.0E-2)

###### nuisance parameters;
sigma2epsilon<-1/tauepsilon
sigma2b11<-1/taub1
sigma2b22<-sigma2b11*pow(theta12,2)+1/taub2
#### covariance and correlation; 
sigma2b12<-theta12*sigma2b11
rho12<-theta12*sqrt(sigma2b11)/sqrt(sigma2b22)
} 


