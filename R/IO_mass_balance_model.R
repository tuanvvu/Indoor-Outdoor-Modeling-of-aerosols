setwd("C:/MBmodel/")
workingDirectory<<-"C:/MBmodel/"
library (openair)
library(lubridate)
mydata<- read.csv("T1OI20.csv",header=TRUE)

###01. Run mass balance modelzz

dp<-seq(0.2,1,0.01) ### P=1 as the design of this experiments
dk<-seq(0,4,0.01)  ### k: deposition rate
ds<-seq(2,27,1)   ### s: size bin
sd<-matrix(nrow=26,ncol=401) ### residual between modelled and measured data
r2<-matrix(nrow=26,ncol=401 )
rsd<-seq(1,2,1)
### Note- O: Outdoor; I: Indoor; Im: Indoor modeled; Ip: Indoor modeled - Indoor measured

ks<-function(k,s) {
  mydata[1,s+52]<- mydata[1,s + 26]
  mydata[1, s + 78] <- 0
  for (i in 2:1061){
    mydata[i, s + 52]<- (0.5*1/(k+0.5))*(1-(exp(-(0.5+k)/3)))*mydata[i, s] + (exp(-(k+0.5)/3))*mydata[i-1, s + 52] 
    if (mydata[i,s] >0) {
    mydata[i, s + 78] <- mydata[i, s+ 52]/mydata [i,s]- mydata[i, s + 26]/mydata[i,s]
    } else { mydata[i,s+78] <- 0 } }
  rsd[1]<-sqrt(mean(mydata[,s+78]^2, na.rm=TRUE)) 
  rsd[2]<-cor(mydata[,s + 52], mydata[,s + 26],use="na.or.complete")^2
  rsd
  }

for (i in 1:26) {       ### Run from 9:46 and end at -- Run time 25 minutes # 
  for (j in 1:401) {
    s<-ds[i]
    k<-dk[j]
    sd[i,j]<-ks(k,s)[1]
    r2[i,j]<-ks(k,s)[2]}}
write.csv(sd,paste(workingDirectory,"sd.csv",sep=""))  ###sd: Residual between modelled and measured data.

KR<-ds ## KR: deposition rate optimize for each size
R2Max<-KR
Rs<-KR
Sdmin<-Rs
SDmin<-KR       #### SDmin: minimum of residual between modelled and measured data set.
for (i in 1:26) {
  KR[i] <-dk[which(sd[i,] == min(sd[i,]), arr.ind = TRUE)]
  s<- ds[i]
  k<-KR[i]
  Rs[i]<-ks(k,s)[2]
  Sdmin[i]<-ks(k,s)[1]
  R2Max[i]<-max(r2[i,])
}
KR
Rs
Sdmin

###03. Run mass balance model with added coagulation 

ks<-function(k,s) {
  mydata[1,s+52]<- mydata[1,s + 26]
  mydata[1, s + 78] <- 0
  for (i in 2:1061){
    mydata[i, s + 52]<- (0.5*1/(k+0.5))*(1-(exp(-(0.5+k)/3)))*mydata[i, s] + (exp(-(k+0.5)/3))*mydata[i-1, s + 52] + mydata[i,s+104]
    if (mydata[i,s] >0) {
      mydata[i, s + 78] <- mydata[i, s+ 52]/mydata [i,s]- mydata[i, s + 26]/mydata[i,s]
    } else { mydata[i,s+78] <- 0 } }
  rsd[1]<-sqrt(mean(mydata[,s+78]^2, na.rm=TRUE)) 
  rsd[2]<-cor(mydata[,s + 52], mydata[,s + 26],use="na.or.complete")^2
  rsd}

for (i in 1:26) {       ### Run from 9:46 and end at -- Run tim 25 minutes # 
  for (j in 1:401) {
    s<-ds[i]
    k<-dk[j]
    sd[i,j]<-ks(k,s)[1]
    r2[i,j]<-ks(k,s)[2]}}
write.csv(sd,paste(workingDirectory,"sd.csv",sep=""))  ###sd: Residual between modelled and measured data.

KR<-ds ## KR: deposition rate optimize for each size
R2Max<-KR
Rs<-KR
Sdmin<-Rs
SDmin<-KR       #### SDmin: minimum of residual between modelled and measured data set.
for (i in 1:26) {
  KR[i] <-dk[which(sd[i,] == min(sd[i,]), arr.ind = TRUE)]
  s<- ds[i]
  k<-KR[i]
  Rs[i]<-ks(k,s)[2]
  Sdmin[i]<-ks(k,s)[1]
  R2Max[i]<-max(r2[i,])
}
KR
Rs
Sdmin
