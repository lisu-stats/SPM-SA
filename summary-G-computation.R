
##### summarize G-computation results
########### author: Li Su


dir<-'C://'
allmeans<-NULL
allmedian<-NULL
allq25<-allq975<-NULL

  for (i in 1:8)
  {   
    alldatasim=NULL
    for (j in 1:200)
    {
    load(file=sprintf(paste(dir,'gcomp_%d_%d.Rdata',sep=''), i,j) )
    alldatasim=rbind(alldatasim, alldata)
    }
    
    for (v in 1:12)
    {
      subdata<-alldatasim[alldatasim[,1]==v,-c(1:3)]
      
      allmeans<-rbind(allmeans, c(as.numeric(i%%2==0), (i+1)%/%2, v, apply(subdata, 2, mean)))
      allmedian<-rbind(allmedian, c( as.numeric(i%%2==0), (i+1)%/%2, v, apply(subdata, 2, median)))
      allq25<-rbind(allq25, c( as.numeric(i%%2==0), (i+1)%/%2, v, apply(subdata, 2, function(x) quantile(x, prob=0.025))))
      allq975<-rbind(allq975, c(as.numeric(i%%2==0), (i+1)%/%2, v, apply(subdata, 2, function(x) quantile(x,prob=0.975))))
    }

  }

### transform to the original square root CD4 count scale
allmeans[,4:5]<-allmeans[,4:5]*7+18
allmedian[,4:5]<-allmedian[,4:5]*7+18
allq25[,4:5]<-allq25[,4:5]*7+18
allq975[,4:5]<-allq975[,4:5]*7+18

### save results
save(allmeans, allmedian, allq25, allq975, file="C://Gsummary.Rdata")


load(file="C://Gsummary.Rdata")

###### plot predicted longitudinal profiles, mean only####
allmeans<-data.frame(allmeans)
cexsize2=0.8
x11()
par(mfrow=c(1,2))
plot(1:12, allmeans[allmeans$X1==0 & allmeans$X2==1,4], ylim=c(0,26), pch=19, col='red', type='n', ylab='square root of CD4 count',xlab='visit', main='No ART at baseline', xlim=c(-1,12))

lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==1,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==1,5], col='gray', pch=1,cex=cexsize2, lwd=2)

lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==2,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==2,5], col='gray', pch=1,cex=cexsize2, lwd=2)

lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==3,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==3,5], col='gray', pch=1,cex=cexsize2, lwd=2)

lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==4,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==0 & allmeans$X2==4,5], col='gray', pch=1,cex=cexsize2, lwd=2)

##points
points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==1,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==1,5], col='gray', pch=1,cex=cexsize2)

points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==2,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==2,5], col='gray', pch=1,cex=cexsize2)

points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==3,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==3,5], col='gray', pch=1,cex=cexsize2)

points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==4,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==0 & allmeans$X2==4,5], col='gray', pch=1,cex=cexsize2)

text(1.0,26,"viral load (0-500)",cex=cexsize2)
text(1.0,22.0,"viral load (500-5k)",cex=cexsize2)
text(1.0,18,"viral load (5k-30k)",cex=cexsize2)
text(1.0,15.0,"viral load (>30k)",cex=cexsize2)


plot(1:12, allmeans[allmeans$X1==1 & allmeans$X2==1,4], ylim=c(0,26), pch=19, col='red', type='n', ylab='square root of CD4 count',xlab='visit', main='with ART at baseline', xlim=c(-1,12))
lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==1,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==1,5], col='gray', pch=1,cex=cexsize2, lwd=2)

lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==2,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==2,5], col='gray', pch=1,cex=cexsize2, lwd=2)

lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==3,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==3,5], col='gray', pch=1,cex=cexsize2, lwd=2)

lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==4,4],pch=1,cex=cexsize2, lwd=2)
lines(1:12, allmeans[allmeans$X1==1 & allmeans$X2==4,5], col='gray', pch=1,cex=cexsize2, lwd=2)

points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==1,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==1,5], col='gray', pch=1,cex=cexsize2)

points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==2,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==2,5], col='gray', pch=1,cex=cexsize2)

points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==3,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==3,5], col='gray', pch=1,cex=cexsize2)

points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==4,4],pch=1,cex=cexsize2)
points(1:12, allmeans[allmeans$X1==1 & allmeans$X2==4,5], col='gray', pch=1,cex=cexsize2)


text(1.0,21,"viral load (0-500)",cex=cexsize2)
text(1.0,18,"viral load (500-5k)",cex=cexsize2)
text(1.0,14,"viral load (5k-30k)",cex=cexsize2)
text(1.0,11.0,"viral load (>30k)",cex=cexsize2)
dev.copy2pdf(device=x11,width=12,height=8,file="C://sensi_cd4.pdf")



