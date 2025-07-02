#
#  R Code for reproducing results of the article
#
#  Latent class distributional regression
#  for the estimation of smooth reference intervals 
#  from contaminated data sources
#
#  T.Hepp et al.
#
###########################################

### Simulation study

source("lcdr_functions.R")

#-----------------------------------------------------------
require(MASS)
require(gamlss)
require(parallel)
#----------------------------------------------------------

# Data generating process

M <- list(mu1="x+10*sin((x-.5)*sqrt(12)*pi/2)",sigma1="8+5*x",
           mu2="c+(c+1)*x+10*sin((x-.5)*sqrt(12)*pi/2)",sigma2="11+9*x")


#########################
#                       #
#  alpha_h=60%, n=500   #
#                       #
#########################

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c20_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","20",M)),n=500,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c20_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c15_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","15",M)),n=500,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c15_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c10_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","10",M)),n=500,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c10_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c5_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","5",M)),n=500,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c5_n500 complete")


#########################
#                       #
#  alpha_h=70%, n=500   #
#                       #
#########################

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c20_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","20",M)),n=500,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c20_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c15_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","15",M)),n=500,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c15_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c10_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","10",M)),n=500,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c10_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c5_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","5",M)),n=500,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c5_n500 complete")


#########################
#                       #
#  alpha_h=80%, n=500   #
#                       #
#########################

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c20_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","20",M)),n=500,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c20_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c15_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","15",M)),n=500,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c15_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c10_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","10",M)),n=500,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c10_n500 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c5_n500 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","5",M)),n=500,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c5_n500 complete")

#########################
#                       #
#  alpha_h=60%, n=1000  #
#                       #
#########################

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c20_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","20",M)),n=1000,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c20_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c15_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","15",M)),n=1000,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c15_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c10_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","10",M)),n=1000,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c10_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a60_c5_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","5",M)),n=1000,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a60_c5_n1000 complete")


#########################
#                       #
#  alpha_h=70%, n=1000  #
#                       #
#########################

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c20_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","20",M)),n=1000,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c20_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c15_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","15",M)),n=1000,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c15_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c10_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","10",M)),n=1000,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c10_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a70_c5_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","5",M)),n=1000,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a70_c5_n1000 complete")


#########################
#                       #
#  alpha_h=80%, n=1000  #
#                       #
#########################

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c20_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","20",M)),n=1000,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c20_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c15_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","15",M)),n=1000,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c15_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c10_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","10",M)),n=1000,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c10_n1000 complete")

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
M_a80_c5_n1000 <- mclapply(id, function(r){
  DAT <- data_mixture(id=1000+r,param=as.list(gsub("c","5",M)),n=1000,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pnorm(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                          data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("a80_c5_n1000 complete")

#####

