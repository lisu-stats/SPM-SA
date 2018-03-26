
################## function to generate samples to perform G-computation################################################# author: Li Su


Gcomp<-function(j, Nall=82700, workdir='C://')
{
 
  ### j = which of the posterior samples
  ### Nall= number of samples needed
  ### workdir= working directory, all code and posterior samples files are saved here.
  
  setwd(workdir)
  getwd()
  
  ## load posterior samples
  load(file='postsample100.Rdata') 
  library(nlme)
  library(mvtnorm) 
  library(mc2d)
  source(file="EV_skewednormaldist.R")

  dir<-workdir  ## directory to save G-computation samples
 
  ni=12
  ri<-(1:12-1)/11
  means=c(0,0)
  ti=ri

####################### load posterior sample values               
gammaD<-gammaD[j,]
alphaD<-alphaD[j,]

alpha<-alpha[j,]
sigma=matrix(sigma[j,], ncol=2)
sigmae=sigmae[j]

#a<--runif(1)  # sensitivity parameter
a<-rtriang(1, -2, -1, 0)  # triangular distribution

### covariate values for replicated data
survmatrixall=matrix(
            c(1, 0, 0, 0, 0, 0,
             1, 0, 0, 0, 0, 1, 
             1, 1, 0, 0, 0, 0,
             1, 1, 0, 0, 0, 1,
             1, 0, 1, 0, 0, 0,
             1, 0, 1, 0, 0, 1,
             1, 0, 0, 1, 0, 0,
             1, 0, 0, 1, 0, 1
             ), byrow=T, ncol=6)

for(i in 1:8)   ### i=which combination of the baseline covariates; there are 8 combinations
{
  ### which combination of the baseline covariates
    longmatrix=cbind(1,ti,t(t(rep(1,ni)))%*%(t(survmatrixall[i,2:6])), (t(t(ti))%*%(t(survmatrixall[i,2:6]))))
    survmatrix<-survmatrixall[i,]


   alldata<-NULL

for (k in 1:Nall)
{
  re=rmvnorm(1,means,sigma)##random effects
  e=rnorm(ni,0,sigmae)
  
  XX<-longmatrix
  ZZ<-longmatrix[,1:2]
  mu<-XX%*%alpha
  
  yi=mu+re[1]+re[2]*ti+e
    
  etaDi<-alphaD[1]+alphaD[2]*survmatrix[2]+alphaD[3]*survmatrix[3]+alphaD[4]*survmatrix[4] +alphaD[5]*survmatrix[5]+alphaD[6]*survmatrix[6]+alphaD[7]*ri+alphaD[8]*ri^2+gammaD[1]*re[1]+gammaD[2]*re[2]
  
  lambdaDi=pnorm(etaDi)
 
  statusDi=as.numeric(runif(ni)<lambdaDi)
  
  Di=try(min(which(statusDi==0)), silent=T)
  if(Di==Inf){Di=13}
 
  whi<-which(ti*11+1<=Di)
  yiobs<-yi[whi]
  mi=length(whi)
 
  #############################################
  ### calculate posterior mean and var of random effects
  if (mi<12)
  {
    if (mi>1)
     {tmp.H <- ZZ[whi,]
     }else{
    tmp.H <- t(as.matrix(ZZ[whi,]))
     }
  
    Hsig.y2 <- t(tmp.H) %*% tmp.H 
    H <- Hsig.y2/sigmae^2
    Sig.u.inv <- solve(sigma)
    V <- chol2inv(chol(H+Sig.u.inv))
  
    yran <- yiobs - mu[whi]
    h <- sweep(tmp.H,1,yran,"*") 
    h <- apply(h,2,sum) / sigmae^2
  
    newmeans <- V%*%h 
  
    ##### calculate nu
    L=matrix(rep(gammaD, mi),byrow=T, ncol=length(gammaD))
    sXi=matrix(rep(c(survmatrix[1:6]), mi),byrow=T, ncol=length(alphaD)-2)
    sXi<-cbind(sXi, ri[whi], ri[whi]^2)


    L[Di,]<- -L[Di,]
    sXi[Di,]<--sXi[Di,]
  
    nu <- -sXi%*%alphaD - (L%*%newmeans)


    ###sample random effects
    Sig11 <- V
    Sig21 <- -L%*%V
    Sig22 <- diag(Di) + L%*%V%*%t(L)
    
    bhat<-re
    sam <- 1
    omega <- rep(1,Di)
    breakflag <- F
    while(any(omega>0))
    {
      omega <- rmvnorm(1,nu, Sig22) 
      sam <- sam+1
      if((sam>10000) & any(omega>0))
      {
        breakflag <- T
        break
      }
    }
    if(!breakflag)
    {
      mucond <- newmeans + t(Sig21) %*% solve(Sig22) %*% (t(omega) - nu) 
      Sigcond <- Sig11 - t(Sig21) %*% solve(Sig22) %*% Sig21
      bhat<- rmvnorm(1, mucond, Sigcond)
    }
   
    
   
    
    out<-EV_skewdnormal(mu=newmeans,Sigma=V,D=L,nu=nu,Delta=diag(mi))
    newsigma=out$variance
    
    
    #######default extrapolation#####################
    newei=rnorm(ni,0,sigmae)
    newyi=newei+mu+ZZ%*%t(bhat)
    
    ## sensitivity analysis extrapolation
    sensiyi=newei+mu+ZZ%*%t(bhat)+sqrt(newsigma[4])*a*((12-mi)/11)*((1:12-mi)/11)*as.numeric(1:12>mi)
    
  
    alldata<-rbind(alldata, cbind(1:12, rep(mi, 12), yi, c(yiobs, newyi[-whi]), c(yiobs, sensiyi[-whi])))
  }else
  {alldata<-rbind(alldata, cbind(1:12, rep(mi, 12), yi, yi, yi))}
  }

   
   #### save G-computation samples
   save(alldata, file=sprintf(paste(dir,'gcomp_%d_%d.Rdata',sep=''), i,j) )
}
  
}

################# perform G-compuation for seperate posterior samples of SPM parameters in parallel 
library(parallel)
no_core<-32
system.time(mclapply(1:32, Gcomp, mc.cores=no_core)) 