##############################################################################################################################
##                                                                                                                          ##
##                                Robust control chart based on mixed effects modeling framework:                           ##
##                                          a case study in NAND flash memory industry                                      ##
##                                                                                                                          ##
##                                                          by DY                                                           ##
##                                                                                                                          ##
## ------------------------------------------------------------------------------------------------------------------------ ##
##                                                                                                                          ##
##  :  R code for Simulation study 1 (Motivation 1) in the Supplementary material                                           ##
##                                                                                                                          ##
##     -. Demonstrate the usefulness of the mixed-effects modeling to include between-week variations                       ##
##                                                                                                                          ##
##############################################################################################################################




##############################################################################################################################
# library
##############################################################################################################################
library(mixmeta)
library(dplyr)
library(ggplot2)
library(reshape2)
library(gridExtra)
#
#
##############################################################################################################################


##############################################################################################################################
# Case 1
##############################################################################################################################

### set seed
set.seed(1000)
#
##############################################################

### data generation setting
sig <- 1
tau <- 1
mu <- 0
TT <- 30
nt <- 10000
mut <- rnorm(TT,mu,tau)
#
##############################################################

### data generation when weekly variation exists
xts <- list()
for(t in 1:TT){ xts[[t]] <- rnorm(nt, mut[t], sig) }
#
##############################################################

### EWMA chart
## phase 1
hat_sig <- xts %>% unlist() %>% sd()
hat_mu <- xts %>% unlist %>% mean()

## plotting statistics & LCL & UCL
xbarts <- lapply(xts, mean) %>% unlist()
ucl <- hat_mu + 3*hat_sig/sqrt(nt)
lcl <- hat_mu - 3*hat_sig/sqrt(nt)
#
##############################################################

### Mixed-effects modeling
## phase 1
xbarts <- lapply(xts, mean) %>% unlist()
xsets <- (lapply(xts, var) %>% unlist())/nt
mixmeta.fit <- mixmeta(xbarts, xsets)

## LCL & UCL
ucl2 <- mixmeta.fit$coefficients + 3*sqrt(xsets+c(mixmeta.fit$Psi))
lcl2 <- mixmeta.fit$coefficients - 3*sqrt(xsets+c(mixmeta.fit$Psi))
#
##############################################################

### plot
plot(1:TT, xbarts, ylim=quantile(c(xbarts,ucl,lcl,ucl2,lcl2),c(0,1)), type='o')
abline(h=ucl, col=2, lwd=2)
abline(h=lcl, col=2, lwd=2)
points(ucl2, col=4, lwd=2, type='l')
points(lcl2, col=4, lwd=2, type='l')
#
##############################################################

### save data
tmp_data1 <- tibble(Case ="Case 1", Week = 1:TT, xbar = xbarts, 
                    Model1_LCL=lcl, Model1_UCL=ucl, Model2_LCL=lcl2, Model2_UCL=ucl2)
#
#
##############################################################################################################################


##############################################################################################################################
# Case 2
##############################################################################################################################

### set seed
set.seed(1000)
#
##############################################################

### data generation setting
sig <- 1
tau <- 1
mu <- 0
TT <- 30
nt <- 100
mut <- rnorm(TT,mu,tau)
#
##############################################################

### data generation when weekly variation exists
xts <- list()
for(t in 1:TT){ xts[[t]] <- rnorm(nt, mut[t], sig) }
#
##############################################################

### EWMA chart
## phase 1
hat_sig <- xts %>% unlist() %>% sd()
hat_mu <- xts %>% unlist %>% mean()

## plotting statistics & LCL & UCL
xbarts <- lapply(xts, mean) %>% unlist()
ucl <- hat_mu + 3*hat_sig/sqrt(nt)
lcl <- hat_mu - 3*hat_sig/sqrt(nt)
#
##############################################################

### Mixed-effects modeling
## phase 1
xbarts <- lapply(xts, mean) %>% unlist()
xsets <- (lapply(xts, var) %>% unlist())/nt
mixmeta.fit <- mixmeta(xbarts, xsets)

## LCL & UCL
ucl2 <- mixmeta.fit$coefficients + 3*sqrt(xsets+c(mixmeta.fit$Psi))
lcl2 <- mixmeta.fit$coefficients - 3*sqrt(xsets+c(mixmeta.fit$Psi))
#
##############################################################

### plot
plot(1:TT, xbarts, ylim=quantile(c(xbarts,ucl,lcl,ucl2,lcl2),c(0,1)), type='o')
abline(h=ucl, col=2, lwd=2)
abline(h=lcl, col=2, lwd=2)
points(ucl2, col=4, lwd=2, type='l')
points(lcl2, col=4, lwd=2, type='l')
#
##############################################################

### save data
tmp_data2 <- tibble(Case ="Case 2", Week = 1:TT, xbar = xbarts, 
                    Model1_LCL=lcl, Model1_UCL=ucl, Model2_LCL=lcl2, Model2_UCL=ucl2)
#
#
##############################################################################################################################


##############################################################################################################################
# Case 3
##############################################################################################################################

### set seed
set.seed(1000)
#
##############################################################

### data generation setting
sig <- 1
tau <- 0.0001
mu <- 0
TT <- 30
nt <- 10000
mut <- rnorm(TT,mu,tau)
#
##############################################################

### data generation when weekly variation exists
xts <- list()
for(t in 1:TT){ xts[[t]] <- rnorm(nt, mut[t], sig) }
#
##############################################################

### EWMA chart
## phase 1
hat_sig <- xts %>% unlist() %>% sd()
hat_mu <- xts %>% unlist %>% mean()

## plotting statistics & LCL & UCL
xbarts <- lapply(xts, mean) %>% unlist()
ucl <- hat_mu + 3*hat_sig/sqrt(nt)
lcl <- hat_mu - 3*hat_sig/sqrt(nt)
#
##############################################################


### Mixed-effects modeling
## phase 1
xbarts <- lapply(xts, mean) %>% unlist()
xsets <- (lapply(xts, var) %>% unlist())/nt
mixmeta.fit <- mixmeta(xbarts, xsets)

## LCL & UCL
ucl2 <- mixmeta.fit$coefficients + 3*sqrt(xsets+c(mixmeta.fit$Psi))
lcl2 <- mixmeta.fit$coefficients - 3*sqrt(xsets+c(mixmeta.fit$Psi))
#
##############################################################

### plot
plot(1:TT, xbarts, ylim=quantile(c(xbarts,ucl,lcl,ucl2,lcl2),c(0,1)), type='o')
abline(h=ucl, col=2, lwd=2)
abline(h=lcl, col=2, lwd=2)
points(ucl2, col=4, lwd=2, type='l')
points(lcl2, col=4, lwd=2, type='l')
#
##############################################################

### save data
tmp_data3 <- tibble(Case ="Case 3", Week = 1:TT, xbar = xbarts, 
                    Model1_LCL=lcl, Model1_UCL=ucl, Model2_LCL=lcl2, Model2_UCL=ucl2)
#
#
##############################################################################################################################


##############################################################################################################################
# Case 4
##############################################################################################################################

### set seed
set.seed(1000)
#
##############################################################

### data generation setting
sig <- 1
tau <- 0.0001
mu <- 0
TT <- 30
nt <- 100
mut <- rnorm(TT,mu,tau)
#
##############################################################

### data generation when weekly variation exists
xts <- list()
for(t in 1:TT){ xts[[t]] <- rnorm(nt, mut[t], sig) }
#
##############################################################

### EWMA chart
## phase 1
hat_sig <- xts %>% unlist() %>% sd()
hat_mu <- xts %>% unlist %>% mean()

## plotting statistics & LCL & UCL
xbarts <- lapply(xts, mean) %>% unlist()
ucl <- hat_mu + 3*hat_sig/sqrt(nt)
lcl <- hat_mu - 3*hat_sig/sqrt(nt)
#
##############################################################

### Mixed-effects modeling
## phase 1
xbarts <- lapply(xts, mean) %>% unlist()
xsets <- (lapply(xts, var) %>% unlist())/nt
mixmeta.fit <- mixmeta(xbarts, xsets)

## LCL & UCL
ucl2 <- mixmeta.fit$coefficients + 3*sqrt(xsets+c(mixmeta.fit$Psi))
lcl2 <- mixmeta.fit$coefficients - 3*sqrt(xsets+c(mixmeta.fit$Psi))
#
##############################################################

# plot
plot(1:TT, xbarts, ylim=quantile(c(xbarts,ucl,lcl,ucl2,lcl2),c(0,1)), type='o')
abline(h=ucl, col=2, lwd=2)
abline(h=lcl, col=2, lwd=2)
points(ucl2, col=4, lwd=2, type='l')
points(lcl2, col=4, lwd=2, type='l')
#
##############################################################

### save data
tmp_data4 <- tibble(Case ="Case 4", Week = 1:TT, xbar = xbarts, 
                    Model1_LCL=lcl, Model1_UCL=ucl, Model2_LCL=lcl2, Model2_UCL=ucl2)
#
#
##############################################################################################################################


##############################################################################################################################
# Result
##############################################################################################################################

### plot
plot_data <- tmp_data1 %>% bind_rows(tmp_data2) %>% bind_rows(tmp_data3) %>% bind_rows(tmp_data4)
plot_data_melted <- melt(plot_data, id.var = c("Case", "Week"))
plot_data_melted %>% head()

ggplot(plot_data_melted, aes(x=Week, y=value)) + 
  geom_line(aes(color=variable)) + geom_point(aes(color=variable)) +
  facet_wrap( ~ Case, ncol=2, scales="free") +
  theme_bw() + ylab("") +
  scale_colour_manual( values=c("black", "#00A9FF", "#00A9FF", "#FF68A1", "#FF68A1"), 
                       name = "Line", 
                       labels = c("Xbar", "Model1 - LCL", "Model1 - UCL", "Model2 - LCL", "Model2 - UCL") )
#
##############################################################

### plot (for report)
## subplot 1
tmp1 <- plot_data_melted %>% filter(Case == "Case 1")
tmp1.p1 <- tmp1 %>% ggplot(aes(x=Week, y=value)) + 
  geom_line(aes(color=variable)) + 
  geom_point(aes(color=variable)) +
  theme_bw() + ylab("") +
  ggtitle("(a)") + theme(plot.title = element_text(hjust = 0.5, size=25)) +
  scale_colour_manual( values=c("black", "#00A9FF", "#00A9FF", "#FF68A1", "#FF68A1"), 
                       name = "Line", 
                       labels = c("Xbar", "Model1 - LCL", "Model1 - UCL", "Model2 - LCL", "Model2 - UCL") )

## subplot 2
tmp2 <- plot_data_melted %>% filter(Case == "Case 2")
tmp2.p2 <- tmp2 %>% ggplot(aes(x=Week, y=value)) + 
  geom_line(aes(color=variable)) + 
  geom_point(aes(color=variable)) +
  theme_bw() + ylab("") +
  ggtitle("(b)") + theme(plot.title = element_text(hjust = 0.5, size=25)) +
  scale_colour_manual( values=c("black", "#00A9FF", "#00A9FF", "#FF68A1", "#FF68A1"), 
                       name = "Line", 
                       labels = c("Xbar", "Model1 - LCL", "Model1 - UCL", "Model2 - LCL", "Model2 - UCL") )

## subplot 3
tmp3 <- plot_data_melted %>% filter(Case == "Case 3")
tmp3.p3 <- tmp3 %>% ggplot(aes(x=Week, y=value)) + 
  geom_line(aes(color=variable)) + 
  geom_point(aes(color=variable)) +
  theme_bw() + ylab("") + 
  ggtitle("(c)") + theme(plot.title = element_text(hjust = 0.5, size=25)) +
  scale_colour_manual( values=c("black", "#00A9FF", "#00A9FF", "#FF68A1", "#FF68A1"), 
                       name = "Line", 
                       labels = c("Xbar", "Model1 - LCL", "Model1 - UCL", "Model2 - LCL", "Model2 - UCL") )

## subplot 4
tmp4 <- plot_data_melted %>% filter(Case == "Case 4")
tmp4.p4 <- tmp4 %>% ggplot(aes(x=Week, y=value)) + 
  geom_line(aes(color=variable)) + 
  geom_point(aes(color=variable)) +
  theme_bw() + ylab("") +
  ggtitle("(d)") + theme(plot.title = element_text(hjust = 0.5, size=25)) +
  scale_colour_manual( values=c("black", "#00A9FF", "#00A9FF", "#FF68A1", "#FF68A1"), 
                       name = "Line", 
                       labels = c("Xbar", "Model1 - LCL", "Model1 - UCL", "Model2 - LCL", "Model2 - UCL") )

## summary
grid.arrange(tmp1.p1, tmp2.p2, tmp3.p3, tmp4.p4, nrow = 2)
#
#
##############################################################################################################################