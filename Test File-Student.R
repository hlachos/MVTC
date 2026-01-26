################################################################################
#### Run this file to estimate the parameters in "simulated data"
#### The data can be generated under MVNC or MVTC models
##### Slash data is also allowed
################################################################################

rm(list=ls(all=TRUE))

library(matrixNormal) #must be installed 
library(MixMatrix) #must be installed 
library(mvtnorm)
library(MomTrunc) #must be installed 
library(mnormt) #must be installed 
library(Matrix) #must be installed 
library(LaplacesDemon)

source("Functions File -Student.R")


p <- 3   # dimension of Sigma
q <- 4   # dimension of Psi
n <- 100 # Sample Size

M <- matrix(data = c(1,2,1,2,
                     2,2,1,1,
                     3,3,2,3), p, q, byrow = T)*3
                     
Sigma <- ar1_cor(p, rho = 0.4)*1.5
Psi <- ar1_cor (q, rho = 0.2)
Psi<-Psi/det(Psi)^(1/q)  # constraint: determinant = 1

cens <- 0.20 # percentage of observations to be censored (number in [0;1])

#dist<-"Normal"  ## Normal data
#dist<-"SL"   ## Slash data
dist<-"t"   ## Student-t data
nu<-4
Ind<-3 # all censored data: use 1,2, or 3

Matrixdata<-GeraMatrix(n=n, cens=cens, Ind=Ind, M=M, U=Sigma, V=Psi, nu=nu, dist=dist)

###### Data entering in the EM functions
X.cens<-Matrixdata$X.cens
cc<- Matrixdata$cc
LS<- Matrixdata$LS

### Student-T fit: estimating nu. For fixed nu at 4, use nu=4
ResCt<-EM.MatrixtC.cens3(X.cens, cc, LS, nuI=TRUE, precision=1e-5, MaxIter=40)
ResCt$loglik
ResCt$mu
M
ResCt$Psi
Psi
ResCt$Sigma
Sigma

### Normal fit: 
ResCN<- EM.MatrixNC.cens3(X.cens, cc, LS, precision=1e-5, MaxIter=40)
ResCN$loglik
ResCN$mu
M




