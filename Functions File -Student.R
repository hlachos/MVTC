################################################################################
##  EM algorithm for estimating the parameters of MVTC and MVNC Models
##                             01/27/2026
##    This file contains all the necessary functions for running the EM
##          Authors: Victor Lachos, Carlos Diniz and Jongwoo Choi
################################################################################

################################################################################
#### log-likelihood: Normal Case
################################################################################

loglikel<-function(dados, muM, SigmaM, PsiM){

  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]
  suma1<-0

  for(j in 1:n){

    suma1<-suma1 + matrixNormal::dmatnorm(dados[,,j], muM, SigmaM, PsiM, tol = .Machine$double.eps^0.5, log = TRUE)

  }

  return(suma1)
}

################################################################################
#### log-likelihood: Student-t Case
################################################################################
loglikelt<-function(dados, muM, SigmaM, PsiM, nu){

  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]
  suma1<-0

  for(j in 1:n){
    suma1<-suma1 + dmvt(vec(dados[,,j]), vec(muM), kronecker(PsiM,SigmaM), df=nu, log = TRUE)
  }

  return(suma1)
}

################################################################################
#### Log-Likelihood censoring data :  MVNC case
################################################################################

logLikCens <- function(cc, LS, dados, muM, SigmaM, PsiM){

  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]
  ver <- matrix(0,n,1)

  Vari<- kronecker(PsiM,SigmaM)

  for (j in 1:n){

    cc1<-as.vector(matrixNormal::vec(cc[,,j]))
    LS1<-as.vector(matrixNormal::vec(LS[,,j]))
    y1<-as.vector(matrixNormal::vec(dados[,,j]))
    mu1<- as.vector(matrixNormal::vec(muM))

    if(sum(cc1)==0){

      ver[j]<-matrixNormal::dmatnorm(dados[,,j], muM, SigmaM, PsiM, tol = .Machine$double.eps^0.5, log = TRUE)
    }

    if(sum(cc1)>=1){
      if(sum(cc1)==p*q){
        ver[j]<- matrixNormal::pmatnorm(Lower = dados[,,j], Upper = LS[,,j], muM, SigmaM, PsiM, tol = .Machine$double.eps^0.5, log = TRUE)
      }
      else{
        muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0],tol = .Machine$double.eps^0.5)%*%(y1[cc1==0]-mu1[cc1==0])
        Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0],tol = .Machine$double.eps^0.5)%*%Vari[cc1==0,cc1==1]
        Sc  <- (Sc+t(Sc))/2

        auxM<- Vari[cc1==0,cc1==0]
        auxM<-(auxM+t(auxM))/2

        auxcdf<-mnormt::sadmvn(lower=as.vector(y1[cc1==1]), upper=as.vector(LS1[cc1==1]), as.vector(muc), Sc, abseps=.Machine$double.eps^0.5)
        if(auxcdf==0)  auxcdf <- .Machine$double.xmin

        ver[j]<- mnormt::dmnorm(as.vector(y1[cc1==0]),as.vector(mu1[cc1==0]),auxM, log = TRUE)+log(auxcdf)
       }
    }
  }


  if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin

  return(sum(ver))

}

################################################################################
#### Log-Likelihood censoring data :  MVTC case
################################################################################

logLikCenst <- function(nu, cc, LS, dados, muM, SigmaM, PsiM){
# GB = GenzBretz(maxpts = 5e4, abseps = 1e-9, releps = 0)
  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]
  ver <- matrix(0,n,1)
  Vari<- kronecker(PsiM,SigmaM)
  for (j in 1:n){
    cc1<-as.vector(matrixNormal::vec(cc[,,j]))
    LS1<-as.vector(matrixNormal::vec(LS[,,j]))
    y1<-as.vector(matrixNormal::vec(dados[,,j]))
    mu1<- as.vector(matrixNormal::vec(muM))
    if(sum(cc1)==0){
      ver[j]<- mvtnorm::dmvt(y1, delta=mu1, sigma=Vari, df=nu, log = TRUE)
    }
    if(sum(cc1)>=1){
      if(sum(cc1)==p*q){
        # arss<-mvtnorm::pmvt(lower = y1, upper = LS1, delta=mu1, sigma=Vari, df=nu)
         arss<-MomTrunc::pmvnormt(lower = y1,upper = LS1,mean = mu1,sigma = Vari, nu = nu)
        if(arss==0)  arss <- .Machine$double.xmin
        ver[j]<-log(arss)
      }
      else{
        nu1<-as.numeric(nu+length(cc1[cc1==0]))
        auxM<- Vari[cc1==0,cc1==0]
        auxM<-as.matrix((auxM+t(auxM))/2)
        muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(auxM)%*%(y1[cc1==0]-mu1[cc1==0])
        Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(auxM)%*%Vari[cc1==0,cc1==1]
        Sc  <- (Sc+t(Sc))/2
        Qy1<-t(y1[cc1==0]-mu1[cc1==0])%*%solve(auxM)%*%(y1[cc1==0]-mu1[cc1==0])
        auxcte<-as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
        Sc22<-auxcte*Sc
        muUi<-muc
        SigmaUi<-Sc22
        SigmaUi<-(SigmaUi+t(SigmaUi))/2
        #auxcdf<- mvtnorm::pmvt(lower = as.vector(y1[cc1==1]), upper = as.vector(LS1[cc1==1]), delta=as.vector(muUi), sigma=SigmaUi, df=nu1)
        auxcdf<- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi), sigma = SigmaUi, nu = nu1)
        if(auxcdf==0)  auxcdf <- .Machine$double.xmin
        aux123<- mvtnorm::dmvt(as.vector(y1[cc1==0]), delta=as.vector(mu1[cc1==0]), sigma=auxM, df=nu, log = TRUE)
        ver[j]<- aux123+log(auxcdf)
      }
    }
  }
  if(length(which(ver == 0)) > 0) ver[which(ver == 0)] <- .Machine$double.xmin
  return(sum(ver))
}

################################################################################
#### Auxiliary function for the the EM: Used in the Cholesky decomposition #####
################################################################################

somaL3<-function(L1,Sigma,Psi){
  n<-ncol(L1)
  p<-ncol(Sigma)
  q<-ncol(Psi)
  suma2<-matrix(0,q,q)
  suma3<-matrix(0,p,p)
  for(j in 1:n){
    aux<- matrix(L1[,j],p,q)
    suma2<-suma2+t(aux)%*%solve(Sigma)%*%aux
    suma3<-suma3+aux%*%solve(Psi)%*%t(aux)
  }
  return(list(sPsi=suma2,sSigma=suma3))
}

################################################################################
############ ECM algorithm for interval censored and missing data ##############
####################### NORMAL: MVNC ###########################################
################################################################################

EM.MatrixNC.cens3<-function(dados, cc, LS, precision=0.0000001, MaxIter=50){
################################################################################
## dados: if censoring it should contain the lower limit, if missing use 0 (array)
## cc: indicator matrix 1: censoring 0: observed    (array)
## LS: superior Limit of the interval censoring, if missing then use +Inf   (array)
## precision:  is the precision used in the convergence
## MaxIter: Maximum number of iterations
################################################################################
  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]

  loglik<-numeric()
  criterio<-1
  count<-0

  # initial values
  dados1<-dados
  dados2<-dados
  dados2[cc==1]<-mean(dados[cc==0])

  mu<-apply(dados2, c(1,2), mean)

  Psi<-diag(q)
  Sigma<-diag(p)

  Vari<- kronecker(Psi,Sigma)

  while(criterio > precision){

    count <- count + 1

    suma0<-matrix(0,p,q)
    suma2<-matrix(0,q,q)
    suma3<-matrix(0,p,p)

    for (j in 1:n){
      cc1<-as.vector(matrixNormal::vec(cc[,,j]))
      LS1<-as.vector(matrixNormal::vec(LS[,,j]))
      y1<-as.vector(matrixNormal::vec(dados[,,j]))
      mu1<- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){
        omega<- matrix(y1,p,q)
        suma0<- suma0 + omega
      }

      if(sum(cc1) >= 1){

        if(sum(cc1) == p*q) {
          muc  <- mu1
          Sc   <- Vari
          aux  <- MomTrunc::meanvarTMD(y1,LS1,muc,Sc,dist="normal")
          tuy  <- aux$mean
          dados1[,,j]<- matrix(tuy,p,q)

        } else {
          muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
          Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%Vari[cc1==0,cc1==1]
          Sc  <- (Sc+t(Sc))/2
          aux <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muc,Sc,dist="normal")
          tuy <- matrix(y1,p*q,1)
          tuy[cc1==1] <- aux$mean
          dados1[,,j]<- matrix(tuy,p,q)
        }

        suma0<-suma0+matrix(tuy,p,q)

      }

    }

    mu<- suma0/n

    for (j in 1:n){
      cc1<-as.vector(matrixNormal::vec(cc[,,j]))
      LS1<-as.vector(matrixNormal::vec(LS[,,j]))
      y1<-as.vector(matrixNormal::vec(dados[,,j]))
      mu1<- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){
        omega<- matrix(y1,p,q)
        suma2<-suma2+t(omega-mu)%*%solve(Sigma)%*%(omega - mu)

      }

      if(sum(cc1) >= 1){

        if(sum(cc1) == p*q) {
          muc  <- mu1
          Sc   <- Vari
          aux  <- MomTrunc::meanvarTMD(y1,LS1,muc,Sc,dist="normal")
          tuy  <- aux$mean
          tuyy <- aux$EYY
          dados1[,,j]<- matrix(tuy,p,q)

        } else {
          muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
          Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%Vari[cc1==0,cc1==1]
          Sc  <- (Sc+t(Sc))/2
          aux <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muc,Sc,dist="normal")
          tuy <- matrix(y1,p*q,1)
          tuy[cc1==1] <- aux$mean
          tuyy <- matrix(0,p*q,p*q)
          tuyy[cc1==1,cc1==1] <- aux$varcov
          tuyy <- tuyy+tuy%*%t(tuy)
          dados1[,,j]<- matrix(tuy,p,q)
        }

        omega1<- tuyy-(tuy)%*%t(mu1)-(mu1)%*%t(tuy)+mu1%*%t(mu1)
        omega<- as.matrix(Matrix::nearPD(omega1)$mat)
        omega<- (omega+t(omega))/2
        L1<- t(Matrix::chol(omega))
        aux<- somaL3(L1,Sigma,Psi)
        suma2<-suma2+aux$sPsi

      }

    }

    Psi <- suma2/det(suma2)^(1/q)

    Vari<- kronecker(Psi,Sigma)

    for (j in 1:n){
      cc1<-as.vector(matrixNormal::vec(cc[,,j]))
      LS1<-as.vector(matrixNormal::vec(LS[,,j]))
      y1<-as.vector(matrixNormal::vec(dados[,,j]))
      mu1<- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){
        omega<- matrix(y1,p,q)
        suma3<-suma3+(omega- mu)%*%solve(Psi)%*%t(omega- mu)

      }

      if(sum(cc1)>= 1){

        if(sum(cc1) == p*q) {
          muc  <- mu1
          Sc   <- Vari
          aux  <- MomTrunc::meanvarTMD(y1,LS1,muc,Sc,dist="normal")
          tuy  <- aux$mean
          tuyy <- aux$EYY
          dados1[,,j]<- matrix(tuy,p,q)

        } else {
          muc <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
          Sc  <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%Vari[cc1==0,cc1==1]
          Sc  <- (Sc+t(Sc))/2
          aux <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muc,Sc,dist="normal")
          tuy <- matrix(y1,p*q,1)
          tuy[cc1==1] <- aux$mean
          tuyy <- matrix(0,p*q,p*q)
          tuyy[cc1==1,cc1==1] <- aux$varcov
          tuyy <- tuyy+tuy%*%t(tuy)
          dados1[,,j]<- matrix(tuy,p,q)
        }

        omega1<- tuyy-(tuy)%*%t(mu1)-(mu1)%*%t(tuy)+mu1%*%t(mu1)
        omega<- as.matrix(Matrix::nearPD(omega1)$mat)
        omega<- (omega+t(omega))/2
        L1<- t(Matrix::chol(omega))
        aux<- somaL3(L1,Sigma,Psi)
        suma3<-suma3+aux$sSigma

      }

    }

    Sigma <- suma3/(q*n)

    Vari<- kronecker(Psi,Sigma)

    loglik[count]<- logLikCens(cc, LS, dados, mu, Sigma, Psi)

    if (count>2){
      at<- (loglik[count]-loglik[count-1])/(loglik[count-1]-loglik[count-2])
      criterio<-(loglik[count]-loglik[count-1])/(1-at)
    }

    if (count==MaxIter){
      criterio <- 0.000000000000001
    }


  }

  npar <- (p*q)+(p*(p+1)/2)+(q*(q+1)/2)-1
  BIC   <- -2*loglik[count] + npar*log(n)    # to be minimized

  obj.out <- list(mu=mu,Sigma=Sigma, Psi=Psi, dadosPred=dados1, loglik=loglik[count], BIC = BIC, iter = count)

  class(obj.out) <- "EM.MatrixNC.cens3"
  return(obj.out)

}


################################################################################
############ ECM algorithm for interval censored and missing data ##############
####################### Student-t: MVTC ###########################################
################################################################################


EM.MatrixtC.cens3<-function(dados, cc, LS, nuI=TRUE, precision=0.0000001, MaxIter=50){
################################################################################
## dados: if censoring it should contain the lower limit, if missing use 0 (array)
## cc: indicator matrix 1: censoring 0: observed    (array)
## LS: superior Limit of the interval censoring, if missing then use +Inf   (array)
## precision:  is the precision used in the convergence
## MaxIter: Maximum number of iterations
## if nuI=TRUE, then nu is estimated. For fixed, use any positive real value>2.1 otherwise
################################################################################  
  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]

  loglik<-numeric()
  criterio<-1
  count<-0

  # initial values
  dados1<-dados
  dados2<-dados
  dados2[cc==1]<-mean(dados[cc==0])
  IV.cens<-EM.MatrixNC2(dados2, precision=0.000001, MaxIter=30)
  #mu<-apply(dados2, c(1,2), mean)
  #Psi<-diag(q) 
  #Sigma<-diag(p)
  
  mu<- IV.cens$mu
  Psi<-IV.cens$Psi
  Sigma<-IV.cens$Sigma
  Vari<- kronecker(Psi,Sigma)
  
  if(nuI==TRUE){nu=3}
  else {nu=nuI}
  while(criterio > precision){

    count <- count + 1
    print(count)
    suma<-0
    suma0<-matrix(0,p,q)
    suma2<-matrix(0,q,q)
    suma3<-matrix(0,p,p)

    for (j in 1:n){
        
      cc1<-as.vector(matrixNormal::vec(cc[,,j]))
      LS1<-as.vector(matrixNormal::vec(LS[,,j]))
      y1<-as.vector(matrixNormal::vec(dados[,,j]))
      mu1<- as.vector(matrixNormal::vec(mu))


      if(sum(cc1)==0){
           dm  <- t(y1-mu1)%*%solve(Vari)%*%(y1-mu1)
            cdm <- as.numeric((nu+p*q)/(nu+dm))
           omega<- matrix(y1,p,q)
         suma0<- suma0 + omega*cdm
         suma<-suma+cdm
      }

      if(sum(cc1) >= 1){

        if(sum(cc1) == p*q) {

              muUi     <- mu1
              SigmaUi  <- Vari
              SigmaUiA <- SigmaUi*nu/(nu+2)
              auxU1    <- MomTrunc::pmvnormt(lower = as.vector(y1),upper = as.vector(LS1), mean = as.vector(muUi),sigma = SigmaUiA, nu = nu+2)
              auxU2    <- MomTrunc::pmvnormt(lower = as.vector(y1),upper = as.vector(LS1), mean = as.vector(muUi),sigma = SigmaUi, nu = nu)
              if(auxU2==0) {auxU2 <- .Machine$double.xmin}
              MoMT     <- MomTrunc::meanvarTMD(y1,LS1,muUi,SigmaUiA,nu=nu+2, dist="t")
              U0       <- as.numeric(auxU1/auxU2)
              U1       <- auxU1/auxU2*MoMT$mean
              U2       <- auxU1/auxU2*MoMT$EYY
              ty       <- MoMT$mean
              tuy      <- U1
              tuyy     <- U2
              tu       <- U0
              dados1[,,j]<- matrix(ty,p,q)

        } else {

              PsiA <- Vari*nu/(nu+2)
              nu1  <- (nu+length(cc1[cc1==0]))

              muc  <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
              Sc   <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%Vari[cc1==0,cc1==1]
              Sc   <- (Sc+t(Sc))/2
              ScA  <- nu/(nu+2)*Sc

              Qy1  <- t(y1[cc1==0]-mu1[cc1==0])%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
              Qy2  <- t(y1[cc1==0]-mu1[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])

              auxcte  <- as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
              auxcte1 <- as.numeric((nu+2+Qy2)/(nu+2+length(cc1[cc1==0])))

              Sc22 <- auxcte*Sc

              muUi    <- muc
              SigmaUi <- Sc22

              SigmaUiA <- auxcte1*ScA
              auxU1 <- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi),sigma = SigmaUiA, nu = nu1+2)
              auxU2<- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi),sigma = SigmaUi, nu = nu1)        
               if(auxU2==0) { auxU2 <- .Machine$double.xmin }
              MoMT <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muUi,SigmaUiA,nu=nu1+2, dist="t")
              U0 <- as.numeric(auxU1/auxU2)/auxcte
              U1 <- (U0)*( MoMT$mean)
              U2 <- (U0)*(MoMT$EYY)

              Auxtuy <- (matrix(y1,p*q,1))
              ty<- Auxtuy
              ty[cc1==1]<-MoMT$mean
              tuy <- Auxtuy*U0
              tuy[cc1==1]<- U1

              tuyy <- (Auxtuy%*%t(Auxtuy))

              AAx <- tuyy[cc1==0,cc1==0]*U0
              ABx <- Auxtuy[cc1==0]%*%t(U1)
              BAx <- t(ABx)
              BBx <- U2

              tuyy[cc1==0,cc1==0] <- AAx
              tuyy[cc1==0,cc1==1] <- ABx
              tuyy[cc1==1,cc1==0] <- BAx
              tuyy[cc1==1,cc1==1] <- BBx


              tu <- U0


          dados1[,,j]<- matrix(ty,p,q)
        }

        suma0<-suma0+matrix(tuy,p,q)
        suma<-suma+tu
      }

    }

    mu<- suma0/suma

    for (j in 1:n){
      cc1<-as.vector(matrixNormal::vec(cc[,,j]))
      LS1<-as.vector(matrixNormal::vec(LS[,,j]))
      y1<-as.vector(matrixNormal::vec(dados[,,j]))
      mu1<- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){

            dm  <- t(y1-mu1)%*%solve(Vari)%*%(y1-mu1)
            cdm <- as.numeric((nu+p*q)/(nu+dm))
           omega<- matrix(y1,p,q)
          suma2  <-  suma2+(t(omega-mu)%*%solve(Sigma)%*%(omega - mu))*cdm

      }

      if(sum(cc1) >= 1){

        if(sum(cc1) == p*q) {

              muUi     <- mu1
              SigmaUi  <- Vari
              SigmaUiA <- SigmaUi*nu/(nu+2)

              auxU1    <- MomTrunc::pmvnormt(lower = as.vector(y1),upper = as.vector(LS1), mean = as.vector(muUi),sigma = SigmaUiA, nu = nu+2)
              auxU2    <- MomTrunc::pmvnormt(lower = as.vector(y1),upper = as.vector(LS1), mean = as.vector(muUi),sigma = SigmaUi, nu = nu)
              if(auxU2==0){  auxU2 <- .Machine$double.xmin}
              MoMT     <- MomTrunc::meanvarTMD(y1,LS1,muUi,SigmaUiA,nu=nu+2, dist="t")
              U0       <- as.numeric(auxU1/auxU2)
              U1       <- auxU1/auxU2*MoMT$mean
              U2       <- auxU1/auxU2*MoMT$EYY
              ty       <- MoMT$mean
              tuy      <- U1
              tuyy     <- U2
              tu       <- U0
              dados1[,,j]<- matrix(ty,p,q)

        } else {
            PsiA <- Vari*nu/(nu+2)
              nu1  <- (nu+length(cc1[cc1==0]))

              muc  <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
              Sc   <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%Vari[cc1==0,cc1==1]
              Sc   <- (Sc+t(Sc))/2
              ScA  <- nu/(nu+2)*Sc

              Qy1  <- t(y1[cc1==0]-mu1[cc1==0])%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
              Qy2  <- t(y1[cc1==0]-mu1[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])

              auxcte  <- as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
              auxcte1 <- as.numeric((nu+2+Qy2)/(nu+2+length(cc1[cc1==0])))

              Sc22 <- auxcte*Sc

              muUi    <- muc
              SigmaUi <- Sc22

              SigmaUiA <- auxcte1*ScA

              auxU1 <- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi),sigma = SigmaUiA, nu = nu1+2)
              auxU2<- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi),sigma = SigmaUi, nu = nu1)          
              if(auxU2==0){  auxU2 <- .Machine$double.xmin}
              MoMT <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1],muUi,SigmaUiA,nu=nu1+2, dist="t")
              U0 <- as.numeric(auxU1/auxU2)/auxcte
              U1 <- (U0)*( MoMT$mean)
              U2 <- (U0)*(MoMT$EYY)

              Auxtuy <- (matrix(y1,p*q,1))
              ty<- Auxtuy
              ty[cc1==1]<-MoMT$mean
              tuy <- Auxtuy*U0
              tuy[cc1==1]<- U1

              tuyy <- (Auxtuy%*%t(Auxtuy))

              AAx <- tuyy[cc1==0,cc1==0]*U0
              ABx <- Auxtuy[cc1==0]%*%t(U1)
              BAx <- t(ABx)
              BBx <- U2

              tuyy[cc1==0,cc1==0] <- AAx
              tuyy[cc1==0,cc1==1] <- ABx
              tuyy[cc1==1,cc1==0] <- BAx
              tuyy[cc1==1,cc1==1] <- BBx
              tu <- as.numeric(U0)
          dados1[,,j]<- matrix(ty,p,q)
        }

        omega1<- tuyy-(tuy)%*%t(mu1)-(mu1)%*%t(tuy)+(mu1%*%t(mu1))*tu
        omega<- as.matrix(Matrix::nearPD(omega1)$mat)
        omega<- (omega+t(omega))/2
        L1<- t(Matrix::chol(omega))
        aux<- somaL3(L1,Sigma,Psi)
        suma2<-suma2+aux$sPsi

      }

    }

    Psi <- suma2/det(suma2)^(1/q)
     Psi<- (Psi+t(Psi))/2
    Vari<- kronecker(Psi,Sigma)

    for (j in 1:n){
      cc1<-as.vector(matrixNormal::vec(cc[,,j]))
      LS1<-as.vector(matrixNormal::vec(LS[,,j]))
      y1<-as.vector(matrixNormal::vec(dados[,,j]))
      mu1<- as.vector(matrixNormal::vec(mu))

      if(sum(cc1)==0){

            dm  <- t(y1-mu1)%*%solve(Vari)%*%(y1-mu1)
            cdm <- as.numeric((nu+p*q)/(nu+dm))
           omega<- matrix(y1,p,q)
         suma3<-suma3+((omega- mu)%*%solve(Psi)%*%t(omega- mu))*cdm

      }

      if(sum(cc1)>= 1){

        if(sum(cc1) == p*q) {
              muUi     <- mu1
              SigmaUi  <- Vari
              SigmaUiA <- SigmaUi*nu/(nu+2)

              auxU1    <- MomTrunc::pmvnormt(lower = as.vector(y1),upper = as.vector(LS1), mean = as.vector(muUi),sigma = SigmaUiA, nu = nu+2)
              auxU2    <- MomTrunc::pmvnormt(lower = as.vector(y1),upper = as.vector(LS1), mean = as.vector(muUi),sigma = SigmaUi, nu = nu)
               if(auxU2==0){ auxU2 <- .Machine$double.xmin  }
              MoMT     <- MomTrunc::meanvarTMD(y1,LS1,muUi,SigmaUiA,nu=nu+2, dist="t")
              U0       <- as.numeric(auxU1/auxU2)
              U1       <- auxU1/auxU2*MoMT$mean
              U2       <- auxU1/auxU2*MoMT$EYY
              ty       <- MoMT$mean
              tuy      <- U1
              tuyy     <- U2
              tu       <- U0
              dados1[,,j]<- matrix(ty,p,q)

        } else {
              PsiA <- Vari*nu/(nu+2)
              nu1  <- (nu+length(cc1[cc1==0]))

              muc  <- mu1[cc1==1]+Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
              Sc   <- Vari[cc1==1,cc1==1]-Vari[cc1==1,cc1==0]%*%solve(Vari[cc1==0,cc1==0])%*%Vari[cc1==0,cc1==1]
              Sc   <- (Sc+t(Sc))/2
              ScA  <- nu/(nu+2)*Sc

              Qy1  <- t(y1[cc1==0]-mu1[cc1==0])%*%solve(Vari[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])
              Qy2  <- t(y1[cc1==0]-mu1[cc1==0])%*%solve(PsiA[cc1==0,cc1==0])%*%(y1[cc1==0]-mu1[cc1==0])

              auxcte  <- as.numeric((nu+Qy1)/(nu+length(cc1[cc1==0])))
              auxcte1 <- as.numeric((nu+2+Qy2)/(nu+2+length(cc1[cc1==0])))

              Sc22 <- auxcte*Sc

              muUi    <- muc
              SigmaUi <- Sc22

              SigmaUiA <- auxcte1*ScA


                 auxU1 <- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi),sigma = SigmaUiA, nu = nu1+2)
              auxU2<- MomTrunc::pmvnormt(lower = as.vector(y1[cc1==1]),upper = as.vector(LS1[cc1==1]), mean = as.vector(muUi),sigma = SigmaUi, nu = nu1)
               if(auxU2==0) { auxU2 <- .Machine$double.xmin}
              MoMT <- MomTrunc::meanvarTMD(y1[cc1==1],LS1[cc1==1], muUi, SigmaUiA, nu=nu1+2, dist="t")
              U0 <- as.numeric(auxU1/auxU2)/auxcte
              U1 <- (U0)*( MoMT$mean)
              U2 <- (U0)*(MoMT$EYY)

              Auxtuy <- (matrix(y1,p*q,1))
              ty<- Auxtuy
              ty[cc1==1]<-MoMT$mean
              tuy <- Auxtuy*U0
              tuy[cc1==1]<- U1

              tuyy <- (Auxtuy%*%t(Auxtuy))

              AAx <- tuyy[cc1==0,cc1==0]*U0
              ABx <- Auxtuy[cc1==0]%*%t(U1)
              BAx <- t(ABx)
              BBx <- U2

              tuyy[cc1==0,cc1==0] <- AAx
              tuyy[cc1==0,cc1==1] <- ABx
              tuyy[cc1==1,cc1==0] <- BAx
              tuyy[cc1==1,cc1==1] <- BBx
              tu <- as.numeric(U0)
          dados1[,,j]<- matrix(ty,p,q)
        }

        omega1<- tuyy-(tuy)%*%t(mu1)-(mu1)%*%t(tuy)+(mu1%*%t(mu1))*tu
        omega<- as.matrix(Matrix::nearPD(omega1)$mat)
        omega<- (omega+t(omega))/2
        L1<- t(Matrix::chol(omega))
        aux<- somaL3(L1,Sigma,Psi)
        suma3<-suma3+aux$sSigma

      }

    }
     
    Sigma <- suma3/(q*n)
    Sigma<- (Sigma+t(Sigma))/2
    Vari<- kronecker(Psi,Sigma)
    if(nuI==TRUE) {     
     Enu<-optimize(logLikCenst, c(2.01,150), tol = 0.0000000001, maximum = TRUE, cc, LS, dados, mu, Sigma, Psi)
     nu<-Enu$maximum
     
     loglik[count]<-Enu$objective
     }
    else{
    loglik[count]<- logLikCenst(nu, cc, LS, dados, mu, Sigma, Psi) }

    if (count>2){
      at<- (loglik[count]-loglik[count-1])/(loglik[count-1]-loglik[count-2])
      criterio<-(loglik[count]-loglik[count-1])/(1-at)
    }

    if (count==MaxIter){
      criterio <- 0.000000000000001
    }


  }

  npar <- (p*q)+(p*(p+1)/2)+(q*(q+1)/2)
  BIC   <- -2*loglik[count] + npar*log(n)    # to be minimized

  obj.out <- list(mu=mu,Sigma=Sigma, Psi=Psi, nu=nu, dadosPred=dados1, loglik=loglik[count], BIC = BIC, iter = count)

  class(obj.out) <- "EM.MatrixtC.cens3"
  return(obj.out)

}

################################################################################
#### Maximum Likelihood estimation for complete data: NORMAL CASE ##############
############################### COMPLETE CASE  #################################
################################################################################

EM.MatrixNC2<-function(dados, precision=0.0000001, MaxIter=50){

  p<-dim(dados)[1]
  q<-dim(dados)[2]
  n<-dim(dados)[3]

  loglik<-numeric()
  criterio<-1
  count<-0

  mu<-apply(dados, c(1,2), mean)

  # initial values

  Psi<-diag(q)
  Sigma<-diag(p)

  while(criterio > precision){

    count <- count + 1

    suma2<-matrix(0,q,q)
    suma3<-matrix(0,p,p)

    for (j in 1:n){
      suma2<-suma2+t(dados[,,j]- mu)%*%solve(Sigma)%*%(dados[,,j]- mu)
    }
    Psi <- suma2/det(suma2)^(1/q)

    for (j in 1:n){
      suma3<-suma3+(dados[,,j]- mu)%*%solve(Psi)%*%t(dados[,,j]- mu)
    }
    Sigma <- suma3/(q*n)

    loglik[count]<-loglikel(dados, mu, Sigma, Psi)

    if (count>2){
      at<- (loglik[count]-loglik[count-1])/(loglik[count-1]-loglik[count-2])
      criterio<-(loglik[count]-loglik[count-1])/(1-at)
    }

    if (count==MaxIter){
      criterio <- 0.000000000000001
    }


  }

  npar <- (p*q)+(p*(p+1)/2)+(q*(q+1)/2)-1
  BIC   <- -2*loglik[count] + npar*log(n)    # to be minimized

  obj.out <- list(mu=mu,Sigma=Sigma, Psi=Psi, loglik=loglik[count], BIC=BIC, iter = count)

  class(obj.out) <- "EM.MatrixNC2"
  return(obj.out)

}


################################################################################
### Samples Generator under different settings
################################################################################

GeraMatrix<-function(n=n, cens=cens, Ind=1, M=M, U=U, V=V, nu, dist="Normal"){
################################################################################
## M number of samples
## Cens : percentage of censoring
## Ind: 1 only interval censoring
##      2 Missing values
##      3: Missing and Censoring (50% Missing and 50% Censoring)
## M: Location
## U= Sigma
## V= Psi
## nu= degree of Freedom for Student-t
## dist: "Normal" , "SL" or "t"
################################################################################
p<-ncol(U)
r<-ncol(V)

X.or <- array(data = NA, dim = c(p,r,n))

if(dist=="Normal"){   # simulate matrix normal data
for (i in 1:n) {
  X.or[,,i] <- LaplacesDemon::rmatrixnorm(M=M,U=U,V=V)
  }}

  
  if(dist=="SL"){    # simulatinge matrix Slash
for (i in 1:n) {
     W<-rbeta(1,nu,1)
  X.or[,,i]<-  LaplacesDemon::rmatrixnorm(M=M,U=U/W,V=V)
  
  }}



if(dist=="t"){    # simulating matrix Student-t data
  for (i in 1:n) {
    W<-rgamma(1,nu/2,nu/2)
    X.or[,,i]<-  LaplacesDemon::rmatrixnorm(M=M,U=U/W,V=V)
  }}


if(Ind==1){
#### ALL INTERVAL CENSORED DATA ####
LS <- X.cens <- X.or
cutoff <- quantile(X.cens, prob=cens)
cc<-(X.cens< cutoff)+0

# change the order??
LS[cc==1]<- X.cens[cc==1]+2*sd(X.or[cc==1]) # LS (when censoring) contains the upper limit points (for missing use +Inf)
LS[cc==0]<- 0   ## LS, when observed (cc==0) can be any value for data points
X.cens[cc==1]<- X.cens[cc==1]-2*sd(X.or[cc==1]) # allocate the lower limit in the  data (for missing use -Inf)
}

####  ALL MISSING DATA ####
if(Ind==2){
LS <- X.cens <- X.or
cutoff <- quantile(X.cens, prob=cens)
cc<-(X.cens< cutoff)+0
# change the order??
LS[cc==1]<- +Inf
LS[cc==0]<- 0  ## LS, when observed (cc==0) can be any value for data points
X.cens[cc==1]<- -Inf
         }

if(Ind==3){
####  INTERVAL CENSORED AND MISSING DATA ####
LS <- X.cens <- X.or
cutoff <- quantile(X.cens, prob=cens)
cc<-(X.cens< cutoff)+0

# Consider 50% censored and 50% missing of the total 15%

n1<-round(0.5*sum(cc))
n2<-round(0.5*sum(cc))

# change the order??
LS[cc==1][1:n1]<- +Inf
LS[cc==1][(n1+1):(n1+n2)]<- X.cens[cc==1][(n1+1):(n1+n2)]+2*sd(X.or[cc==1][(n1+1):(n1+n2)])

LS[cc==0]<- 0  ## LS, when observed (cc==0) can be any value for data points
X.cens[cc==1][1:n1]<- -Inf
X.cens[cc==1][(n1+1):(n1+n2)]<- X.cens[cc==1][(n1+1):(n1+n2)]-2*sd(X.or[cc==1][(n1+1):(n1+n2)]) # allocate the lower limit in the  data (for missing use -Inf)

}
return(list(X.cens=X.cens,cc=cc, LS=LS))
}


################################################################################
### AR(1) matrix: for data generation
################################################################################


ar1_cor <- function(n, rho){
  exponent <- abs(outer(1:n,1:n,"-"))
  rho^exponent
}