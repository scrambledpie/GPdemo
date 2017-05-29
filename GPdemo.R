#par(mfrow=1)
#####################################################
#a single normal poiont
while (TRUE){
  y=rnorm(1)
  plot(0:1,c(0,0),ylim=c(-2.5,2.5),type="l")
  lines(c(0.5,0.5),c(0,y))
  points(0.5,y,pch=19,col="red")
  Sys.sleep(1)
}


#####################################################
#2 covarying points
library(MASS)
mean=c(0,0); sigma=matrix(c(1,0.80,0.80,1),2,2)
while (TRUE){
  y=mvrnorm(1,mean,sigma)
  plot(0:1,c(0,0),ylim=c(-2.5,2.5),type="l")
  lines(c(0.25,0.25),c(0,y[1]))
  lines(c(0.75,0.75),c(0,y[2]))
  
  points(c(0.25,0.75),y,pch=19,col="red")
  Sys.sleep(1)
}

#####################################################
#points in 1D Instance Space = {1,...,50}
I=1:50
CovMatrix=matrix(0,50,50)
c("red","darkgreen","blue") ->C
for (i in I){for (j in I){CovMatrix[i,j]=exp( -0.5* (i-j)^2 /10^2 )}}
#for (i in I){CovMatrix[i,i]=1.1}
while(TRUE){
  plot(I,rep(0,length(I)),type="l",ylim=c(-3,3) )
  points(I,mvrnorm(1,rep(0,50),CovMatrix),col=C[1],pch=19)
  Sys.sleep(1)
}

#####################################################
#lines in 1D Instance Space = {1,...,50}
I=1:50
CovMatrix=matrix(0,50,50)
for (i in I){for (j in I){CovMatrix[i,j]=exp( -0.5* (i-j)^2 /5^2 )}}
while(TRUE){
  plot(I,rep(0,length(I)),type="l",ylim=c(-3,3) )
  for (i in 1:3){
  lines(I,mvrnorm(1,rep(0,50),CovMatrix),col=C[i],pch=19,lwd=4)
  }
  Sys.sleep(1)
}


####################################################

#interpolation
DI=c(10,25,40); FI1=c(-1,2,0.5)
c("red","darkgreen","blue") ->C

DDCov=matrix(0,3,3)
for (i in 1:3){for (j in 1:3){DDCov[i,j]=exp( -0.5* (I[DI[i]]-I[DI[j]])^2 /5^2 )}} 

DCov=matrix(0,3,50)
for (i in 1:3){for (j in 1:50){DCov[i,j]=exp( -0.5* (I[DI[i]]-I[j])^2 /5^2 )}}

fbar=t(DCov)%*%solve(DDCov)%*%FI1
FullCov=CovMatrix - t(DCov)%*%solve(DDCov)%*%(DCov)

library(MASS)
while(TRUE){
plot(I,rep(0,50),ylim=c(-3,3),type="l")
Y1=mvrnorm(1,fbar,FullCov); lines(I,Y1, col="blue",lwd=2)
Y2=mvrnorm(1,fbar,FullCov); lines(I,Y2, col="red",lwd=2)
Y3=mvrnorm(1,fbar,FullCov); lines(I,Y3, col="darkgreen",lwd=2)
points(DI,FI1,col="black",pch=19)
Sys.sleep(1)
}

####################################################

SD=2*sqrt(diag(FullCov))
# for(i in I){SD[i]=2*sqrt(FullCov[i,i])}
plot(I,fbar,lwd=3,col="black",type="l", ylim=c(-3,3))
lines(I,fbar+SD,col="red")
lines(I,fbar-SD,col="red")
points(DI,FI1,col="black",pch=19)
SD=0;


####################################################

#c(-1.5,2,0.5)
DDCov=matrix(0,3,3);for (i in 1:3){for (j in 1:3){DDCov[i,j]=exp( -0.5* (I[DI[i]]-I[DI[j]])^2 /5^2 )}} 
DCov=matrix(0,3,50);for (i in 1:3){for (j in 1:50){DCov[i,j]=exp( -0.5* (I[DI[i]]-I[j])^2 /5^2 )}}
FullCov=CovMatrix - t(DCov)%*%solve(DDCov)%*%(DCov)
c("red","darkgreen","blue") ->C
t=0
while(TRUE){
  t=t+0.1
  F1=1.5*sin(t+2*(1:3))
  fbar=t(DCov)%*%solve(DDCov)%*%F1
  
  plot(I,fbar,lwd=3,col="black",type="l", ylim=c(-3,3))
  points(DI,F1,col="black",pch=19)
  SD=0;
  for(i in I){SD[i]=2*sqrt(FullCov[i,i])}
  lines(I,fbar+SD,col="red")
  lines(I,fbar-SD,col="red")
  Sys.sleep(0.23)
}
####################################################

#2D animted MV samples
library(rgl)
N=10
l=2^2
testx=cbind(rep(1:N,each=N),rep(1:N,N))
Ntrain=3
trainx=cbind(rep(N*(-0.5+1:Ntrain)/Ntrain,Ntrain),rep(N*(-0.5+1:Ntrain)/Ntrain,each=Ntrain))
K=matrix(0,nrow(trainx),nrow(trainx))
for(i in 1:nrow(trainx)){
  for(j in 1:nrow(trainx)){
       K[i,j]=exp(-0.5*sum((trainx[i,]-trainx[j,])^2)/l)
  }
}
iK=ginv(K)
Ks=matrix(0,nrow(trainx),nrow(testx))
for(i in 1:nrow(trainx)){for(j in 1:nrow(testx)){
  Ks[i,j]=exp(-0.5*sum((trainx[i,]-testx[j,])^2)/l)
}}
KsiK=t(Ks)%*%(iK)
t=0
phase=2*pi*runif(nrow(trainx))
while(T){
  clear3d()
  trainy=1.5*sin(t+phase)
  meany=KsiK%*%trainy
  surface3d(1:N,1:N,meany,lit=T,col="red")
  points3d(trainx[,2],trainx[,1],trainy,col="blue",size=15)
  surface3d(1:N,1:N,rep(-1.5,N*N))
  box3d()
  t=t+0.2
  print(t)
  Sys.sleep(0.25)
}

