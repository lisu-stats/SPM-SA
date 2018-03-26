#############################################################################;
## Mean and Variance for a multivariate skewed normal disrtribution;
## Author: Qiuju Li
## document name: EV_skewednormaldist/EV stands for Expectation and Variance;
#############################################################################;


#install.packages("mvtnorm")
library(mvtnorm)  ### pmvnorm() uo to 10000 dims; much slower;
#install.packages("mnormt")
library(mnormt)   #### pmnorm() up to 20 dims; in our case: 1<=dim<=24;

### nu: vector of length p;

dev_normalcdf<-function(G, nu, Psi){
nu<-as.numeric(nu)  # thus nu is set to be a vector;
p<-length(nu)  
if(p==1)  ### G, nu, Psi are scalar;
   {G<-as.numeric(G)  ## make sure G is a vector instaed of a maxtrix;
    Psi<-as.numeric(Psi) ## make sure Psi is a number;
    dev1<-dnorm(0,mean=nu,sd=sqrt(Psi))*G
    dev2<-dnorm(0,mean=nu,sd=sqrt(Psi))*nu*G%*%t(G)/Psi
    }
if(p==2) 
   {G<-as.matrix(G)  ## G has to be a matrix of 2 rows;
    if(dim(G)[1]!=2) stop("G has to be a matrix of 2 rows")
    ellmean1<-nu[2]-Psi[2,1]*nu[1]/Psi[1,1]
    ellcov1<-Psi[2,2]-Psi[2,1]*Psi[1,2]/Psi[1,1]
    ellmean2<-nu[1]-Psi[1,2]*nu[2]/Psi[2,2]
    ellcov2<-Psi[1,1]-Psi[1,2]*Psi[2,1]/Psi[2,2]
  
    ellprobvec1<-pnorm(0,mean=ellmean1,sd=sqrt(ellcov1))*dnorm(0,mean=nu[1],sd=sqrt(Psi[1,1]))*G[1,]
    ellprobvec2<-pnorm(0,mean=ellmean2,sd=sqrt(ellcov2))*dnorm(0,mean=nu[2],sd=sqrt(Psi[2,2]))*G[2,]
    
    ellG1<-G[2,]-Psi[2,1]*G[1,]/Psi[1,1] ## vector;
    ellG2<-G[1,]-Psi[1,2]*G[2,]/Psi[2,2]
    ### first derivative;
    dev1<-ellprobvec1+ellprobvec2
    ### second derivative; please note NEVER write one formula in two lines, error may occur!!!
    dev2<-dnorm(0,mean=ellmean1,sd=sqrt(ellcov1))*dnorm(0,mean=nu[1],sd=sqrt(Psi[1,1]))*G[1,]%*%t(ellG1)+dnorm(0,mean=ellmean2,sd=sqrt(ellcov2))*dnorm(0,mean=nu[2],sd=sqrt(Psi[2,2]))*G[2,]%*%t(ellG2)+ellprobvec1%*%t(G[1,])*nu[1]/Psi[1,1]+ellprobvec2%*%t(G[2,])*nu[2]/Psi[2,2]
    }
if(p>2)
   {dev1<-0
    dev2<-0
    G<-as.matrix(G)
    if(dim(G)[1]!=p) stop("please check the dimension of matrix G")
    for (ell in 1:p)
        {ellmean<-nu[-ell]-Psi[-ell,ell]*nu[ell]/Psi[ell,ell]
         ellcov<-Psi[-ell,-ell]-Psi[-ell,ell]%*%t(Psi[ell,-ell])/Psi[ell,ell]
         if(p<=21) ellprobvec<-pmnorm(rep(0,p-1),mean=ellmean,varcov=ellcov)*dnorm(0,mean=nu[ell],sd=sqrt(Psi[ell,ell]))*G[ell,] else
                  {ellprobvec<-pmvnorm(upper=rep(0,p-1),mean=ellmean,sigma=ellcov)*dnorm(0,mean=nu[ell],sd=sqrt(Psi[ell,ell]))*G[ell,]}         
         ### fisrt derivative;
         dev1<-dev1+ellprobvec
         ### second derivative;
         ellG<-G[-ell,]-Psi[-ell,ell]%*%t(G[ell,])/Psi[ell,ell]
         kterm<-0
         for(k in 1:(p-1))
           {kmean<-ellmean[-k]-ellcov[-k,k]*ellmean[k]/ellcov[k,k]
            kcov<-ellcov[-k,-k]-ellcov[-k,k]%*%t(ellcov[k,-k])/ellcov[k,k]
            kpdfvec<-dnorm(0,mean=ellmean[k],sd=sqrt(ellcov[k,k]))*ellG[k,]
            if(p==3) kterm<-kterm+pnorm(0,mean=kmean,sd=sqrt(as.numeric(kcov)))*kpdfvec else
                     {if(p<=22) kterm<-kterm+pmnorm(rep(0,p-2),mean=kmean,varcov=kcov)*kpdfvec else
                                {kterm<-kterm+pmvnorm(upper=rep(0,p-2),mean=kmean,sigma=kcov)*kpdfvec} 
                     }
            }
         dev2<-dev2+dnorm(0,mean=nu[ell],sd=sqrt(Psi[ell,ell]))*G[ell,]%*%t(kterm)+ellprobvec%*%t(G[ell,])*nu[ell]/Psi[ell,ell]
         }
    }
output<-list(dev1=dev1,dev2=dev2)
return(output)
}

### mu is a vector;
### Q: D has to be a matrix????
EV_skewdnormal<-function(mu,Sigma,D,nu,Delta){
nu<-as.numeric(nu) 
mu<-as.numeric(mu)
q<-length(mu)
p<-length(nu)
dencov<-as.matrix(Delta)+as.matrix(D)%*%as.matrix(Sigma)%*%as.matrix(t(D))  ### covariance matrix for the denominator part;
if(p==1) denprob<-pnorm(0,mean=nu,sd=sqrt(as.numeric(dencov))) else
  {if(p<=20) denprob<-pmnorm(x=rep(0,p),mean=nu,varcov=dencov) else denprob<-pmvnorm(upper=rep(0,p),mean=nu,sigma=dencov)}
Gmatrix<-as.matrix(D)%*%as.matrix(Sigma)
### derivatives output;
devout<-dev_normalcdf(G=Gmatrix, nu=nu, Psi=dencov)
### mean
outmean<-mu+devout$dev1/denprob
### variance
outvar<-Sigma+devout$dev2/denprob-(devout$dev1/denprob)%*%t(devout$dev1/denprob)
if(q==1) outvar<-as.numeric(outvar) ## if q==1, then variance is return in a number instead of matirx;
output<-list(mean=outmean,variance=outvar)  
return(output)
}


