################################################################################
#### REAL DATA – EE2.1
#### Matrix-variate Student-t (censored), estimating nu
################################################################################

rm(list = ls(all = TRUE))

## ===================== Packages =====================
library(baytrends)
library(dplyr)

library(matrixNormal)
library(MixMatrix)
library(mvtnorm)
library(MomTrunc)
library(mnormt)
library(Matrix)
library(LaplacesDemon)

## ===================== Load EM / ECM functions =====================
codes_dir <- "C:/Users/carlo/Dropbox/Spatial Point Process/MeuProjeto/Artigo1_Hyperbolic/Matrix-t/MVTC-codes"
setwd(codes_dir)
source("Functions File -Student.R")

## ===================== EM controls =====================
precision <- 1e-5
MaxIter   <- 40

## ===================== Load raw data =====================
data("dataCensored")
dataOR <- unSurvDF(dataCensored)

`%!in%` <- Negate(`%in%`)

## ===================== EE2.1 specification =====================
stat      <- unique(dataOR$station)
stat.sel  <- 3                      # EE2.1

vars.cens <- c("chla", "din", "nh4", "no23", "po4", "tdn", "tdp", "tn", "tp")
sel       <- c(5, 6, 7)             # po4, tdn, tdp
p         <- length(sel)

col.nam   <- c("AP","B","BP","S")
col.sel   <- 1:4
r         <- length(col.sel)

vars.sel <- paste0(vars.cens[sel], "_lo")

## ===================== Data extraction =====================
filt1  <- filter(dataOR, station == stat[stat.sel])
times1 <- unique(filt1$date)

## Remove incomplete dates
elim1 <- numeric(length(times1))

for (i in seq_along(times1)) {
  tmp <- filter(filt1, date == times1[i])

  if (nrow(tmp) != r) elim1[i] <- 1
  if (all(is.na(select(tmp, vars.sel)))) elim1[i] <- 1
}

times1 <- times1[elim1 == 0]
n      <- length(times1)

## ===================== Build arrays =====================
X.cens <- array(NA, dim = c(p, r, n),
                dimnames = list(vars.cens[sel], col.nam[col.sel], 1:n))

cc <- array(0, dim = c(p, r, n))
LS <- array(0, dim = c(p, r, n))

for (i in 1:n) {

  tmp   <- filter(filt1, date == times1[i])
  tmp_v <- select(tmp, vars.sel)

  X.cens[,,i] <- t(tmp_v)

  ## ---- Missing values ----
  if (any(is.na(X.cens[,,i]))) {
    idx <- which(is.na(X.cens[,,i]))
    cc[,,i][idx] <- 1
    LS[,,i][idx] <- +Inf
    X.cens[,,i][idx] <- -Inf
  }

  ## ---- Censoring ----
  sel2 <- which(colnames(tmp) %in% vars.sel)

  for (j in 1:p) {
    dif <- tmp[, sel2[j]] - tmp[, sel2[j] + 1]
    if (any(dif != 0, na.rm = TRUE)) {
      idx <- which(dif != 0)
      cc[j, idx, i] <- 1
      LS[j, idx, i] <- tmp[idx, sel2[j] + 1]
    }
  }
}

## ===================== Quick dataset summary (optional) =====================
cat("\n================ EE2.1 DATA SUMMARY ================\n")
cat(sprintf("Array: %d x %d x %d\n", dim(X.cens)[1], dim(X.cens)[2], dim(X.cens)[3]))
cat(sprintf("Censura (cc==1): %d (%.2f%%)\n",
            sum(cc == 1), 100 * sum(cc == 1) / length(cc)))
cat(sprintf("Valores -Inf em X.cens: %d (%.2f%%)\n",
            sum(is.infinite(X.cens) & X.cens < 0, na.rm = TRUE),
            100 * sum(is.infinite(X.cens) & X.cens < 0, na.rm = TRUE) / length(X.cens)))

## ===================== Model fitting =====================
cat("\n>>> Fitting matrix-variate Student-t with censoring (estimating nu)...\n")

ResCt_EE <- EM.MatrixtC.cens3(
  X.cens, cc, LS,
  nuI       = TRUE,
  precision = precision,
  MaxIter   = MaxIter
)

## ===================== Output =====================
cat("\n================ EE2.1 RESULTS ================\n")
cat("\n--- Student-t (censored) ---\n")
print(ResCt_EE$loglik)
print(ResCt_EE$mu)
print(ResCt_EE$Sigma)
print(ResCt_EE$Psi)
print(ResCt_EE$nu)
