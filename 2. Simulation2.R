##############################################################################################################################
##                                                                                                                          ##
##                                Robust control chart based on mixed effects modeling framework:                           ##
##                                          a case study in NAND flash memory industry                                      ##
##                                                                                                                          ##
##                                                          by DY                                                           ##
##                                                                                                                          ##
## ------------------------------------------------------------------------------------------------------------------------ ##
##                                                                                                                          ##
##  :  R code for Simulation study 2 (Motivation 4) in the Supplementary material                                           ##
##                                                                                                                          ##
##     -. Demonstrate the usefulness of the p-value combination method to detect gradual process changes efficiently        ##
##                                                                                                                          ##
##############################################################################################################################




##############################################################################################################################
# library
##############################################################################################################################
library(dplyr)
library(ggplot2)
library(reshape2)
#
#
##############################################################################################################################


##############################################################################################################################
# function for p-value combination
##############################################################################################################################
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
##############################################################################################################################


##############################################################################################################################
# ARL0
##############################################################################################################################

### simulation setting
target.ARL0 <- 100
Llist <- seq(1.5,3,0.01)
Llist2 <- seq(1.8,4,0.01)
sim.num <- 5000
#
##############################################################

### hyperparameters for MECCp-EWMA
window.c <- 10
gamma.val <- 0.8
#
##############################################################

### data generation setting 1
sig <- 1
mu0 <- 0
#
##############################################################

### data generation setting 2
TT0 <- 50
TT1 <- 1000
TT2 <- 50
TT <- TT0 + TT1 + TT2
nt <- 1000
#
##############################################################

### simulation to obtain L when ARL0 = 100
ARL0.SC.list <- ARL0.PC.list <- NULL
for(i in 1:sim.num){
  
  ## set seed
  cat(i, "    ")
  set.seed(i)
  
  ## data generation
  xts <- list()
  for(t in 1:TT){ xts[[t]] <- rnorm(nt, mu0, sig) }
  xbarts <- lapply(xts[-c(1:TT0,((TT0+TT1+1):TT))], mean) %>% unlist()
  
  ## (Phase 1) parameter estimation
  hat_sig <- xts[1:TT0] %>% unlist() %>% sd()
  hat_mu <- xts[1:TT0] %>% unlist() %>% mean()
  
  ## ARL0 for Shewhart xbar chart
  tmp.arl0list <- NULL
  for(tmp.l in Llist){
    
    tmp.lcl <- hat_mu - tmp.l*hat_sig/sqrt(nt)
    tmp.ucl <- hat_mu + tmp.l*hat_sig/sqrt(nt)
    tmp.ooc.idx <- which((xbarts < tmp.lcl) | (xbarts > tmp.ucl))
    if(length(tmp.ooc.idx)>0){
      
      tmp.rl <- tmp.ooc.idx[1]
    }else{ tmp.rl <- TT1+1 }
    tmp.arl0list <- c(tmp.arl0list, tmp.rl)
  }
  ARL0.SC.list <- rbind(ARL0.SC.list, tmp.arl0list)
  
  ## ARL0 for p-value combination
  zscores <- (xbarts - hat_mu)/(hat_sig/sqrt(nt))
  pvalues <- 2*pnorm(-abs(zscores))
  zscores.combined <- pvalueCombine(pvalues, gamma.val, window.c)
  
  tmp.arl0list2 <- NULL
  for(tmp.l in Llist2){
    
    tmp.ucl <- tmp.l
    tmp.ooc.idx <- which((zscores.combined > tmp.ucl))
    if(length(tmp.ooc.idx)>0){
      
      tmp.rl <- tmp.ooc.idx[1]
    }else{ tmp.rl <- TT1+1 }
    tmp.arl0list2 <- c(tmp.arl0list2, tmp.rl)
  }
  ARL0.PC.list <- rbind(ARL0.PC.list, tmp.arl0list2)
}
#
##############################################################


### control limit L
## Shewhart control chart
L.SC.idx <- which(apply(ARL0.SC.list, 2, mean) > target.ARL0)[1]
L.SC <- Llist[L.SC.idx]

## pvalue combination
L.PC.idx <- which(apply(ARL0.PC.list, 2, mean) > target.ARL0)[1]
L.PC <- Llist2[L.PC.idx]
#
#
##############################################################################################################################


##############################################################################################################################
# ARL1
##############################################################################################################################

### data generation setting 1-2
mu1 <- 1
#
##############################################################

### data generation setting 2-2
TT3 <- 50
TT4 <- 50
#
##############################################################

### simulation to obtain the run length distribution when ARL0=100
ARL1.SC.list <- ARL1.PC.list <- NULL
for(i in 1:sim.num){
  
  ## set seed
  cat(i, "    ")
  set.seed(i)
  
  ## data generation
  xts <- list()
  for(t in 1:TT){ xts[[t]] <- rnorm(nt, mu0, sig) }
  for(t in (TT+1):(TT+TT3)){ xts[[t]] <- rnorm(nt, mu0 + (mu1 - mu0)/(TT3+1)*(t-TT), sig) }
  for(t in (TT+TT3+1):(TT+TT3+TT4)){ xts[[t]] <- rnorm(nt, mu1, sig) }
  xbarts1 <- lapply(xts[(TT0+TT1+1):(TT+TT3+TT4)], mean) %>% unlist()
  
  ## (Phase 1) parameter estimation
  hat_sig <- xts[1:TT0] %>% unlist() %>% sd()
  hat_mu <- xts[1:TT0] %>% unlist() %>% mean()
  
  ## zscores of p-value combinations
  zscores <- (xbarts1 - hat_mu)/(hat_sig/sqrt(nt))
  pvalues <- 2*pnorm(-abs(zscores))
  zscores.combined <- pvalueCombine(pvalues, gamma.val, window.c)
  
  ## ARL1 for Shewhart control chart
  tmp.lcl <- hat_mu - L.SC*hat_sig/sqrt(nt)
  tmp.ucl <- hat_mu + L.SC*hat_sig/sqrt(nt)
  tmp.idx <- which((xbarts1[-(1:TT2)] > tmp.ucl) | (xbarts1[-(1:TT2)] < tmp.lcl))[1]
  ARL1.SC.list <- c(ARL1.SC.list, tmp.idx)
  
  ## ARL1 for pvalue combination
  tmp.ucl <- L.PC
  tmp.idx2 <- which(zscores.combined[-(1:TT2)] > tmp.ucl)[1]
  ARL1.PC.list <- c(ARL1.PC.list, tmp.idx2)
}
#
#
##############################################################################################################################
  
##############################################################################################################################
# Result
##############################################################################################################################

### result
## ARL1
ARL1.SC <- mean(ARL1.SC.list)
ARL1.PC <- mean(ARL1.PC.list)
c(ARL1.SC, ARL1.PC)

## SDRL
SDRL.SC <- sd(ARL1.SC.list)
SDRL.PC <- sd(ARL1.PC.list)
c(SDRL.SC, SDRL.PC)
#
#
##############################################################################################################################