setwd("C:/MBmodel/")
workingDirectory<<-"C:/MBmodel/"
library (openair)
library(lubridate)
mydata<- read.csv("T1OI20.csv",header=TRUE)

### Run Coagulation model

nt<-matrix(nrow=1061, ncol=26)
vt<-matrix(nrow=1061, ncol=26)
dp<-c(1.43E-06, 1.65E-06, 1.91E-06, 2.21E-06, 2.55E-06, 2.94E-06,3.40E-06,3.92E-06, 4.53E-06,5.23E-06, 6.04E-06, 6.98E-06, 8.06E-06, 9.31E-06, 1.08E-05, 1.24E-05, 1.43E-05, 1.66E-05, 1.91E-05, 2.21E-05, 2.55E-05, 2.95E-05, 3.4E-05, 3.92E-05, 4.53E-05, 5.23E-05)
for (i in 2:1061){
for (s in 1:26) {
  nt[1,s]<- 0.5*mydata[1,s+1]
  nt[i,s] <-0.5*mydata[i,s + 1] + mydata[i-1,s+27]
  vt[,s] <- nt[,s]*(dp[s]^3)*pi/6
}}
write.csv(nt,paste(workingDirectory,"nt.csv",sep=""), row.names=FALSE)
write.csv(vt,paste(workingDirectory,"vt.csv",sep=""), row.names=FALSE)
nt<-read.csv("C:/MBmodel/nt.csv", header=TRUE) 
vt<-read.csv("C:/MBmodel/vt.csv", header=TRUE)

sizevolume<-read.csv("sizevolume.csv", header=TRUE)  ### Input file information of particle number, diameter and volume
v<-sizevolume$v               ### A Volume of a particle in size bin 1 to 26 of FMPS measurement
dn<-matrix(nrow=1061, ncol=26)               ### A array for partilce loss due to coagulation from size bin 1 to size bin 26 with 1061 measurements
hour<-rep(60,length.out=360)  ### Time step =60s with a simulation time is 6 hours
## nt      ### Number particles in each size bin
## vt      ### Total volume of particle concentration vt=nt*v
m1<-seq(1,26,1)               ### Set up array for denominator of coagualation fraction;
dk<-seq(1,26,1)               ### The particle remain after coagulation, dk <- nt + dn for 1 time step  
dk1<-matrix(nrow=26,ncol=360)  ### An matrix for dk at 360 time step

KBV<- read.csv("KBVE.csv", header=FALSE)   ### Input the coagulation kernel, KBVE: Brownian kernel with Convection enhance and Van der Waal/vicous force
volume<-matrix(nrow=26, ncol=26)   ### Volume fraction for new particle by sum of size bin i and j
for (i in 1:26) {
  for (j  in 1:26) {
    volume[i,j] <-v [i] + v[j]
  }}
write.csv(volume,paste(workingDirectory,"volume.csv",sep=""))

h<-60 ### Time step for coagulation 
y<-function (h,r) {
  ###01. q =1 at dp =16.5 nm. The first size bin
  for (j in 1:26) {
    m1[j]<- - h*KBV[1,j]*nt[r,1]*nt[r,j]
  }
  dn[r,1]<-sum(m1)
  
  ###02. q =2 at dp =16.5 nm. The 2nd size bin 
  f2<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 1
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[3] & volume [i,j] > v[2]) {
        f2[i,j] <- (v[3]- volume[i,j])*v[2]/((v[3]-v[2])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[2] & volume [i,j] > v[1]) {
        f2[i,j] <- 1- (v[2]- volume[i,j])*v[1]/((v[2]-v[1])*volume[i,j])}}}
  
  f2[is.na(f2)]<-0
  write.csv(f2,paste(workingDirectory,"f2.csv",sep="")) ### fi,j,q
  for (j in 1:26) {
    m1[j]<-(1-f2[2,j])*KBV[2,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=1, ncol=2)
  for (i in 1:1) {
    for (j in 1:2) {
      t1[i,j]<-f2[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan16.5<- (t- vt[r,2]*(m-1))/(m*v[2])
  dn[r,2]<-deltan16.5 
  
  ###03. q =3 at dp =19.1 nm. The 3rd size bin 
  f3<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 1
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[4] & volume [i,j] > v[3]) {
        f3[i,j] <- (v[4]- volume[i,j])*v[3]/((v[4]-v[3])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[3] & volume [i,j] > v[2]) {
        f3[i,j] <- 1- (v[3]- volume[i,j])*v[2]/((v[3]-v[2])*volume[i,j])}}}
  
  f3[is.na(f3)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f3[3,j])*KBV[3,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=2, ncol=3)
  for (i in 1:2) {
    for (j in 1:3) {
      t1[i,j]<-f3[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan19.1<- (t- vt[r,3]*(m-1))/(m*v[3])
  dn[r,3]<-deltan19.1
  
  ###04. q =4 at dp =22.1 nm. The 4th size bin 
  f4<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 4
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[5] & volume [i,j] > v[4]) {
        f4[i,j] <- (v[5]- volume[i,j])*v[4]/((v[5]-v[4])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[4] & volume [i,j] > v[3]) {
        f4[i,j] <- 1- (v[4]- volume[i,j])*v[3]/((v[4]-v[3])*volume[i,j])}}}
  
  f4[is.na(f4)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f4[4,j])*KBV[4,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=3, ncol=4)
  for (i in 1:3) {
    for (j in 1:4) {
      t1[i,j]<-f4[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan22.1<- (t- vt[r,4]*(m-1))/(m*v[4])
  dn[r,4]<-deltan22.1
  
  ###05. q =5 at dp =25.5 nm. The 5th size bin 
  f5<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 5
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[6] & volume [i,j] > v[5]) {
        f5[i,j] <- (v[6]- volume[i,j])*v[5]/((v[6]-v[5])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[5] & volume [i,j] > v[4]) {
        f5[i,j] <- 1- (v[5]- volume[i,j])*v[4]/((v[5]-v[4])*volume[i,j])}}}
  
  f5[is.na(f5)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f5[5,j])*KBV[5,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=4, ncol=5)
  for (i in 1:4) {
    for (j in 1:5) {
      t1[i,j]<-f5[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan25.5<- (t- vt[r,5]*(m-1))/(m*v[5])
  dn[r,5]<-deltan25.5
  
  ###06. q =6 at dp =29.4 nm. The 6th size bin 
  f6<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 6
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[7] & volume [i,j] > v[6]) {
        f6[i,j] <- (v[7]- volume[i,j])*v[6]/((v[7]-v[6])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[6] & volume [i,j] > v[5]) {
        f6[i,j] <- 1- (v[6]- volume[i,j])*v[5]/((v[6]-v[5])*volume[i,j])}}}
  
  f6[is.na(f6)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f6[6,j])*KBV[6,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=5, ncol=6)
  for (i in 1:5) {
    for (j in 1:6) {
      t1[i,j]<-f6[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan29.4<- (t- vt[r,6]*(m-1))/(m*v[6])
  dn[r,6]<-deltan29.4
  
  ###07. q =7 at dp =34 nm. The 7th size bin 
  f7<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 7
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[8] & volume [i,j] > v[7]) {
        f7[i,j] <- (v[8]- volume[i,j])*v[7]/((v[8]-v[7])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[7] & volume [i,j] > v[6]) {
        f7[i,j] <- 1- (v[7]- volume[i,j])*v[6]/((v[7]-v[6])*volume[i,j])}}}
  
  f7[is.na(f7)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f7[7,j])*KBV[7,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=6, ncol=7)
  for (i in 1:6) {
    for (j in 1:7) {
      t1[i,j]<-f7[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan34<- (t- vt[r,7]*(m-1))/(m*v[7])
  dn[r,7]<-deltan34
  
  ###08. q =8 at dp =39.2 nm. The 8th size bin 
  f8<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 8
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[9] & volume [i,j] > v[8]) {
        f8[i,j] <- (v[9]- volume[i,j])*v[8]/((v[9]-v[8])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[8] & volume [i,j] > v[7]) {
        f8[i,j] <- 1- (v[8]- volume[i,j])*v[7]/((v[8]-v[7])*volume[i,j])}}}
  
  f8[is.na(f8)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f8[8,j])*KBV[8,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=7, ncol=8)
  for (i in 1:7) {
    for (j in 1:8) {
      t1[i,j]<-f8[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan39.2<- (t- vt[r,8]*(m-1))/(m*v[8])
  dn[r,8]<-deltan39.2
  
  ###09. q =9 at dp =45.3 nm. The 9th size bin 
  f9<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 8
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[10] & volume [i,j] > v[9]) {
        f9[i,j] <- (v[10]- volume[i,j])*v[9]/((v[10]-v[9])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[9] & volume [i,j] > v[8]) {
        f9[i,j] <- 1- (v[9]- volume[i,j])*v[8]/((v[9]-v[8])*volume[i,j])}}}
  
  f9[is.na(f9)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f9[9,j])*KBV[9,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=8, ncol=9)
  for (i in 1:8) {
    for (j in 1:9) {
      t1[i,j]<-f9[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan45.3<- (t- vt[r,9]*(m-1))/(m*v[9])
  dn[r,9]<-deltan45.3
  
  ###10. q =10 at dp =52.3 nm. The 10th size bin 
  f10<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[11] & volume [i,j] > v[10]) {
        f10[i,j] <- (v[11]- volume[i,j])*v[10]/((v[11]-v[10])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[10] & volume [i,j] > v[9]) {
        f10[i,j] <- 1- (v[10]- volume[i,j])*v[9]/((v[10]-v[9])*volume[i,j])}}}
  
  f10[is.na(f10)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f10[10,j])*KBV[10,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=9, ncol=10)
  for (i in 1:9) {
    for (j in 1:10) {
      t1[i,j]<-f10[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan52.3<- (t- vt[r,10]*(m-1))/(m*v[10])
  dn[r,10]<-deltan52.3
  
  ###11. q =11 at dp =60.4 nm. The 11th size bin 
  f11<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 11
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[12] & volume [i,j] > v[11]) {
        f11[i,j] <- (v[12]- volume[i,j])*v[11]/((v[12]-v[11])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[11] & volume [i,j] > v[10]) {
        f11[i,j] <- 1- (v[11]- volume[i,j])*v[10]/((v[11]-v[10])*volume[i,j])}}}
  
  f11[is.na(f11)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f11[11,j])*KBV[11,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=10, ncol=11)
  for (i in 1:10) {
    for (j in 1:11) {
      t1[i,j]<-f11[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan60.4<- (t- vt[r,11]*(m-1))/(m*v[11])
  dn[r,11]<-deltan60.4
  
  ###12. q =12 at dp =69.8 nm. The 12th size bin 
  f12<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 12
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[13] & volume [i,j] > v[12]) {
        f12[i,j] <- (v[13]- volume[i,j])*v[12]/((v[13]-v[12])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[12] & volume [i,j] > v[11]) {
        f12[i,j] <- 1- (v[12]- volume[i,j])*v[11]/((v[12]-v[11])*volume[i,j])}}}
  
  f12[is.na(f12)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f12[12,j])*KBV[12,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=11, ncol=12)
  for (i in 1:11) {
    for (j in 1:12) {
      t1[i,j]<-f12[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan69.8<- (t- vt[r,12]*(m-1))/(m*v[12])
  dn[r,12]<-deltan69.8
  
  ###13. q =13 at dp =80.6 nm. The 13th size bin 
  f13<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[14] & volume [i,j] > v[13]) {
        f13[i,j] <- (v[14]- volume[i,j])*v[13]/((v[14]-v[13])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[13] & volume [i,j] > v[12]) {
        f13[i,j] <- 1- (v[13]- volume[i,j])*v[12]/((v[13]-v[12])*volume[i,j])}}}
  
  f13[is.na(f13)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f13[13,j])*KBV[13,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=12, ncol=13)
  for (i in 1:12) {
    for (j in 1:13) {
      t1[i,j]<-f13[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan80.6<- (t- vt[r,13]*(m-1))/(m*v[13])
  dn[r,13]<-deltan80.6
  
  ###14. q =14 at dp =93.1 nm. The 14th size bin 
  f14<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[15] & volume [i,j] > v[14]) {
        f14[i,j] <- (v[15]- volume[i,j])*v[14]/((v[15]-v[14])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[14] & volume [i,j] > v[13]) {
        f14[i,j] <- 1- (v[14]- volume[i,j])*v[13]/((v[14]-v[13])*volume[i,j])}}}
  
  f14[is.na(f14)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f14[14,j])*KBV[14,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=13, ncol=14)
  for (i in 1:13) {
    for (j in 1:14) {
      t1[i,j]<-f14[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan93.1<- (t- vt[r,14]*(m-1))/(m*v[14])
  dn[r,14]<-deltan93.1
  
  ###15. q =15 at dp =107.5 nm. The 15th size bin 
  f15<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[16] & volume [i,j] > v[15]) {
        f15[i,j] <- (v[16]- volume[i,j])*v[15]/((v[16]-v[15])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[15] & volume [i,j] > v[14]) {
        f15[i,j] <- 1- (v[15]- volume[i,j])*v[14]/((v[15]-v[14])*volume[i,j])}}}
  
  f15[is.na(f15)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f15[15,j])*KBV[15,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=14, ncol=15)
  for (i in 1:14) {
    for (j in 1:15) {
      t1[i,j]<-f15[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan107.5<- (t- vt[r,15]*(m-1))/(m*v[15])
  dn[r,15]<-deltan107.5
  
  ###16. q =16 at dp =124.1 nm. The 16th size bin 
  f16<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[17] & volume [i,j] > v[16]) {
        f16[i,j] <- (v[17]- volume[i,j])*v[16]/((v[17]-v[16])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[16] & volume [i,j] > v[15]) {
        f16[i,j] <- 1- (v[16]- volume[i,j])*v[15]/((v[16]-v[15])*volume[i,j])}}}
  
  f16[is.na(f16)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f16[16,j])*KBV[16,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=15, ncol=16)
  for (i in 1:15) {
    for (j in 1:16) {
      t1[i,j]<-f16[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan124.1<- (t- vt[r,16]*(m-1))/(m*v[16])
  dn[r,16]<-deltan124.1
  
  ###17. q =17 at dp =143.3 nm. The 17th size bin 
  f17<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[18] & volume [i,j] > v[17]) {
        f17[i,j] <- (v[18]- volume[i,j])*v[17]/((v[18]-v[17])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[17] & volume [i,j] > v[16]) {
        f17[i,j] <- 1- (v[17]- volume[i,j])*v[16]/((v[17]-v[16])*volume[i,j])}}}
  
  f17[is.na(f17)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f17[17,j])*KBV[17,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=16, ncol=17)
  for (i in 1:16) {
    for (j in 1:17) {
      t1[i,j]<-f17[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan143.3<- (t- vt[r,17]*(m-1))/(m*v[17])
  dn[r,17]<-deltan143.3
  
  ###18. q =18 at dp =165.3 nm. The 18th size bin 
  f18<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[19] & volume [i,j] > v[18]) {
        f18[i,j] <- (v[19]- volume[i,j])*v[18]/((v[19]-v[18])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[18] & volume [i,j] > v[17]) {
        f18[i,j] <- 1- (v[18]- volume[i,j])*v[17]/((v[18]-v[17])*volume[i,j])}}}
  
  f18[is.na(f18)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f18[18,j])*KBV[18,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=17, ncol=18)
  for (i in 1:17) {
    for (j in 1:18) {
      t1[i,j]<-f18[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan165.5<- (t- vt[r,18]*(m-1))/(m*v[18])
  dn[r,18]<-deltan165.5
  
  ###19. q =19 at dp =191.1 nm. The 19th size bin 
  f19<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[20] & volume [i,j] > v[19]) {
        f19[i,j] <- (v[20]- volume[i,j])*v[19]/((v[20]-v[19])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[19] & volume [i,j] > v[18]) {
        f19[i,j] <- 1- (v[19]- volume[i,j])*v[18]/((v[19]-v[18])*volume[i,j])}}}
  
  f19[is.na(f19)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f19[19,j])*KBV[19,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=18, ncol=19)
  for (i in 1:18) {
    for (j in 1:19) {
      t1[i,j]<-f19[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan191.1<- (t- vt[r,19]*(m-1))/(m*v[19])
  dn[r,19]<-deltan191.1
  
  ###20. q =20 at dp =220.7 nm. The 20th size bin 
  f20<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[21] & volume [i,j] > v[20]) {
        f20[i,j] <- (v[21]- volume[i,j])*v[20]/((v[21]-v[20])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[20] & volume [i,j] > v[19]) {
        f20[i,j] <- 1- (v[20]- volume[i,j])*v[19]/((v[20]-v[19])*volume[i,j])}}}
  
  f20[is.na(f20)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f20[20,j])*KBV[20,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=19, ncol=20)
  for (i in 1:19) {
    for (j in 1:20) {
      t1[i,j]<-f20[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan220.7<- (t- vt[r,20]*(m-1))/(m*v[20])
  dn[r,20]<-deltan220.7
  
  ###21. q =21 at dp =254.8 nm. The 21th size bin 
  f21<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[22] & volume [i,j] > v[21]) {
        f21[i,j] <- (v[22]- volume[i,j])*v[21]/((v[22]-v[21])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[21] & volume [i,j] > v[20]) {
        f21[i,j] <- 1- (v[21]- volume[i,j])*v[20]/((v[21]-v[20])*volume[i,j])}}}
  
  f21[is.na(f21)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f21[21,j])*KBV[21,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=20, ncol=21)
  for (i in 1:20) {
    for (j in 1:21) {
      t1[i,j]<-f21[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan254.8<- (t- vt[r,21]*(m-1))/(m*v[21])
  dn[r,21]<-deltan254.8
  
  ###22. q =22 at dp =294.4 nm. The 22th size bin 
  f22<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[23] & volume [i,j] > v[22]) {
        f22[i,j] <- (v[23]- volume[i,j])*v[22]/((v[23]-v[22])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[22] & volume [i,j] > v[21]) {
        f22[i,j] <- 1- (v[22]- volume[i,j])*v[21]/((v[22]-v[21])*volume[i,j])}}}
  
  f22[is.na(f22)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f22[22,j])*KBV[22,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=21, ncol=22)
  for (i in 1:21) {
    for (j in 1:22) {
      t1[i,j]<-f22[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan294.5<- (t- vt[r,22]*(m-1))/(m*v[22])
  dn[r,22]<-deltan294.5
  
  ###23. q =23 at dp =339.8 nm. The 23th size bin 
  f23<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[24] & volume [i,j] > v[23]) {
        f23[i,j] <- (v[24]- volume[i,j])*v[23]/((v[24]-v[23])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[23] & volume [i,j] > v[22]) {
        f23[i,j] <- 1- (v[23]- volume[i,j])*v[22]/((v[23]-v[22])*volume[i,j])}}}
  
  f23[is.na(f23)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f23[23,j])*KBV[23,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=22, ncol=23)
  for (i in 1:22) {
    for (j in 1:23) {
      t1[i,j]<-f23[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan339.8<- (t- vt[r,23]*(m-1))/(m*v[23])
  dn[r,23]<-deltan339.8
  
  ###24. q =24 at dp =392.4 nm. The 24th size bin 
  f24<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[25] & volume [i,j] > v[24]) {
        f24[i,j] <- (v[25]- volume[i,j])*v[24]/((v[25]-v[24])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[24] & volume [i,j] > v[23]) {
        f24[i,j] <- 1- (v[24]- volume[i,j])*v[23]/((v[24]-v[23])*volume[i,j])}}}
  
  f24[is.na(f24)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f24[24,j])*KBV[24,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=23, ncol=24)
  for (i in 1:23) {
    for (j in 1:24) {
      t1[i,j]<-f24[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan392.4<- (t- vt[r,24]*(m-1))/(m*v[24])
  dn[r,24]<-deltan392.4
  
  ###25. q =25 at dp =453.2 nm. The 24th size bin 
  f25<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] <= v[26] & volume [i,j] > v[25]) {
        f25[i,j] <- (v[26]- volume[i,j])*v[25]/((v[26]-v[25])*volume[i,j])}}}
  for (i in 1:26) {
    for (j  in 1:26) {
      if (volume[i,j] < v[25] & volume [i,j] > v[24]) {
        f25[i,j] <- 1- (v[25]- volume[i,j])*v[24]/((v[25]-v[24])*volume[i,j])}}}
  
  f25[is.na(f25)]<-0
  for (j in 1:26) {
    m1[j]<-(1-f25[25,j])*KBV[25,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=24, ncol=25)
  for (i in 1:24) {
    for (j in 1:25) {
      t1[i,j]<-f25[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan453.2<- (t- vt[r,25]*(m-1))/(m*v[25])
  dn[r,25]<-deltan453.2
  
  ###26. q =26 at dp =523.3 nm. The 26th size bin 
  f26<-matrix (nrow=26, ncol=26) ### Correction coagulation rate for size bin 10
  for (i in 1:26) {
    for (j  in 1:26) {  
      if (volume[i,j] >= v[26]) {
        f26[i,j] <- 1 }}}
  
  f26[is.na(f26)]<-0 
  
  for (j in 1:26) {
    m1[j]<-(1-f26[26,j])*KBV[26,j]*nt[r,j]
    m <- 1 + h* sum(m1)
  }
  t1<-matrix(nrow=25, ncol=26)
  for (i in 1:25) {
    for (j in 1:26) {
      t1[i,j]<-f26[i,j]*KBV[i,j]*vt[r,i]*nt[r,j]
    }
  }
  t<- h*sum(t1)
  deltan523.3<- (t- vt[r,26]*(m-1))/(m*v[26])
  dn[r,26]<-deltan523.3   
  dn
}

for (k in 1:360) {
  h<-hour2[k]
  for(i in 1:26) {
    dk1[i,k] <- nt[i] + y (h)[i]
    nt[i]<-dk1[i,k]
    vt[i]<-nt[i]*v[i]
    dk1[i,k]
  }}

dnn<-dn #### Particle lose
hour<-rep(1,length.out=120)  ### Time step =60s with a simulation time is 6 hours
for (k in 1:120)  {
     h<-hour[k]
for (r in 1:1061) {     
  dnn[r,]<-y(h,r)[r,]
  for (i in 1:26) {
  nt[r,i]<-nt[r,i] + dnn[r,i]
  vt[r,i] <- nt[r,i]*(dp[i]^3)*pi/6
  dnn
}}}
write.csv(dnn,paste(workingDirectory,"dnn10minutes2.csv",sep="")) 


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
