##############################################################################################################################
##                                                                                                                          ##
##                                Robust control chart based on mixed effects modeling framework:                           ##
##                                          a case study in NAND flash memory industry                                      ##
##                                                                                                                          ##
##                                                          by DY                                                           ##
##                                                                                                                          ##
## ------------------------------------------------------------------------------------------------------------------------ ##
##                                                                                                                          ##
##  :  R code for simple implementation of the proposed control chart MECCp-EWMA                                            ##
##                                                                                                                          ##
##############################################################################################################################




###############################################################################################################################
# library
###############################################################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
library(mixmeta)
#
#
###############################################################################################################################


###############################################################################################################################
# function for p-value combination
###############################################################################################################################
pvalueCombine <- function(pvals, Gamma, Window){
  
  result <- rep(NA, length(pvals))
  weightvals <- (Gamma)^(Window:1)
  if(length(pvals) >= Window){
    
    for(t in Window:length(pvals)){
      
      tmp.ts <- (t-Window+1):t
      tmp.ws <- weightvals
      tmp.ps <- pvals[tmp.ts]
      tmp.zs <- qnorm(1-tmp.ps)
      tmp.zs[tmp.zs == Inf] <- 12
      tmp.zs[tmp.zs == -Inf] <- -12
      
      result[t] <- sum(tmp.ws*(tmp.zs))/sqrt(sum(tmp.ws^2))
    }
  }
  return(result)
}
#
#
###############################################################################################################################


###############################################################################################################################
# data generation
###############################################################################################################################

### set seed
set.seed(1)
#
##############################################################

### hyperparameters for MECCp-EWMA
lambda <- 0.2     # EWMA
window.c <- 10    # p-value combination
gamma.val <- 0.8  # p-value combination
#
##############################################################

### data generation setting 1
sig <- 1
tau <- 0.3
mu0 <- 0
mu1 <- 1
#
##############################################################

### data generation setting 2
TT0 <- 50
TT1 <- 50
TT2 <- 30
TT3 <- 20
TT <- TT0 + TT1 + TT2 + TT3
nt <- 10000
#
##############################################################

### data generation
mus <- c(rep(mu0,TT0),                         # phase 1 
         rep(mu0,TT1),                         # phase 2 - in control
         mu0 + (mu1 - mu0)/(TT2+1) * (1:TT2),  # phase 2 - out of control (gradual increase)
         rep(mu1,TT3))                         # phase 2 - out of control
mut <- rnorm(TT, mus, tau)
xts <- list()
for(t in 1:TT){ xts[[t]] <- rnorm(nt, mut[t], sig) }
#
#
###############################################################################################################################


###############################################################################################################################
# MECCp-EWMA
###############################################################################################################################

### transformation (sign chart)
gmedian <- xts[1:TT0] %>% unlist %>% median()
yts <- lapply(1:TT, function(t){ as.numeric(xts[[t]]>gmedian) })
ybarts <- lapply(yts, mean) %>% unlist()
#
##############################################################

### (Phase 1) parameter estimation
tmp_mu <- yts[1:TT0] %>% unlist() %>% mean()
hat_sig <- sqrt(tmp_mu * (1-tmp_mu))
mixmeta.fit <- mixmeta(ybarts[1:TT0], rep(hat_sig^2/nt, TT0))
hat_u <- mixmeta.fit$Psi %>% c() %>% sqrt()
hat_mu <- mixmeta.fit$coefficients
#
##############################################################

### EWMA trend
ybarts2 <- ybarts[-(1:TT0)]
wts <- ybarts2[1]; for(t in 2:length(ybarts2)){ wts <- c(wts, lambda * ybarts2[t] + (1-lambda) * wts[t-1]) }
v2ts <- hat_u^2 + hat_sig^2/nt
for(t in 2:length(wts)){
  
  tmp.const <- sum(c(rep(lambda^2,t-1),1) * (1-lambda)^(2*(0:(t-1))))
  tmp.v2 <- tmp.const * (hat_u^2 + hat_sig^2/nt)
  v2ts <- c(v2ts, tmp.v2)
}
#
##############################################################

### p-value combinations
zscores <- (wts - hat_mu)/sqrt(v2ts)
pvalues <- 2*pnorm(-abs(zscores))
zscores.combined <- pvalueCombine(pvalues, gamma.val, window.c)
#
#
###############################################################################################################################


###############################################################################################################################
# Result
###############################################################################################################################

### plot
par(mfrow=c(2,1))

## MECC-EWMA
L <- 3
cl_idx <- 11:(TT1+TT2+TT3)
plot(1:(TT1+TT2+TT3), wts, type='o', 
     ylim=quantile(c(wts,hat_mu-sqrt(v2ts[cl_idx])*L,hat_mu+sqrt(v2ts[cl_idx])*L),c(0,1)),
     main="MECC", xlab="time", ylab="Wt")
lines(cl_idx, hat_mu + sqrt(v2ts[cl_idx])*L, col=2, lwd=2)
lines(cl_idx, hat_mu - sqrt(v2ts[cl_idx])*L, col=2, lwd=2)
abline(v=TT1+1/2, col=4, lwd=2, lty=2)

## MECCp - EWMA
Lq <- 4
plot(1:(TT1+TT2+TT3), zscores.combined, type='o',
     main="MECCp", xlab="time", ylab="Qt")
abline(h=Lq, col=2, lwd=2)
abline(v=TT1+1/2, col=4, lwd=2, lty=2)
#
#
###############################################################################################################################
