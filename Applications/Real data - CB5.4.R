library(baytrends)
library(dplyr)

data("dataCensored")
dataOR <- unSurvDF(dataCensored)
`%!in%` = Negate(`%in%`)

stat <- unique(dataOR$station)
vars.cens <- c("chla", "din", "nh4", "no23", "po4", "tdn", "tdp", "tn", "tp")
col.nam <- c("AP","B","BP","S")

stat.sel <- 7 # select the station
sel  <- c(2,3,5) # select the variables in vars.cens (must be sorted)
p <- length(sel) # number of variables 
col.sel <- c(1,2,3,4) # select the columns in col.nam (must be sorted)
r <- length(col.sel) # number of columns

vars.sel <- character(length = length(sel))
for (i in 1:length(sel)) {
  
  vars.sel[i] <- paste(vars.cens[sel[i]],"_lo",sep = "")
  
}

#### Data extraction and arrangement ####

filt1 <- filter(dataOR, station==stat[stat.sel])
times1 <- unique(filt1$date)
elim1 <- numeric(length(times1))

for (i in 1:length(times1)) {
  
  temp1 <- filter(filt1, date==times1[i])
  if(nrow(temp1)!=r) 
    elim1[i] <- 1
  
  temp2 <- select(temp1, vars.sel)
  if(all(is.na(temp2))) 
    elim1[i] <- 1
  
}

times1<- times1[-which(elim1==1)]
n <- length(times1)

X.or <- array(NA, dim = c(p,r,n), dimnames = list(c(vars.cens[sel]),
                                                  c(col.nam[col.sel]),
                                                  c(1:n)))
X.cens <- X.or
cc <- LS <- array(0, dim = c(p,r,n))

for (i in 1:n) {
  
  temp3   <- filter(filt1, date==times1[i])
  temp4   <- select(temp3, vars.sel)
  
  X.cens[,,i] <- X.or[,,i] <- t(temp4)
  
  ## Handling Missing Values (NA) ##
  
  if(any(is.na(X.cens[,,i]))){
    res0 <- which(is.na(X.cens[,,i]))
    cc[,,i][res0] <- 1
    LS[,,i][res0] <- +Inf
    X.cens[,,i][res0] <- -Inf
  }
  
  ## Handling Censoring Values ##
  
  sel2 <- which(colnames(temp3) %in% vars.sel) # positions of the variables in temp3
  
  for (j in 1:p) {
    
    dif <- temp3[,sel2[j]]-temp3[,(sel2[j]+1)]
    if(any(dif!=0, na.rm = T)){
      
      res1 <- which(dif!=0)
      LS[j,res1,i] <- temp3[res1,(sel2[j]+1)]
      cc[j,res1,i] <- 1
      
    }
    
  }
  
}

#### Fitting Part ####

set.seed(1234)

### Fit our model ###

res.cens<- EM.MatrixNC.cens3(X.cens, cc, LS, precision=1e-7, MaxIter=500)
rest.rm<- EM.MatrixtC.cens3(X.cens, cc, LS, nuI=TRUE, precision=1e-7, MaxIter=500) # to estimate nu use nuI=TRUE
  
### Fit Matrix normal after removing na/cens observations ###

ck1 <- numeric(n)

for (i in 1:n) {
  
  if(any(cc[,,i]==1)){
    
    ck1[i] <-1
    
  }
  
}

X.rm <- X.cens[,,which(ck1!=1)]

res.rm<- EM.MatrixNC2(dados = X.rm, precision=1e-7, MaxIter=500)

#### Analysis of the results ####

perc.change <- function(V1,V2){
  ((V2-V1)/abs(V1))*100
}

round(res.cens$mu[c(3,1,2),c(4,1,3,2)],digits = 3)
round(cov2cor(res.cens$Sigma)[c(3,1,2),c(3,1,2)],digits = 3)
round(cov2cor(res.cens$Psi)[c(4,1,3,2),c(4,1,3,2)],digits = 3)

# perc. change 
round(perc.change(res.cens$mu[c(3,1,2),c(4,1,3,2)],res.rm$mu[c(3,1,2),c(4,1,3,2)]),digits = 3)
round(perc.change(cov2cor(res.cens$Sigma)[c(3,1,2),c(3,1,2)],cov2cor(res.rm$Sigma)[c(3,1,2),c(3,1,2)]),digits = 3) 
round(perc.change(cov2cor(res.cens$Psi)[c(4,1,3,2),c(4,1,3,2)],cov2cor(res.rm$Psi)[c(4,1,3,2),c(4,1,3,2)]),digits = 3)

