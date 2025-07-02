#
#  Supplementary R Code for the article
#
#  Latent class distributional regression
#  for the estimation of smooth reference intervals 
#  from contaminated data sources
#
#  T.Hepp et al.
#
###########################################

### Gamma Mixture (Simulation study)

source("lcdr_functions.R")

#-----------------------------------------------------------
require(MASS)
require(gamlss)
require(parallel)
#----------------------------------------------------------

# Data generating process

M <- list(mu1="5+x-2.5*(x-1)^2",sigma1="exp(-1+x-x^2-.5*x^3)",
          mu2="9+2*x-2.5*(x-1)^2",sigma2="exp(-1+x-x^2-.5*x^3)")

############
#          #
#  n=500   #
#          #
############

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a60_n500 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=500,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a60_n500 complete")


ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a70_n500 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=500,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a70_n500 complete")


ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a80_n500 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=500,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a80_n500 complete")



#############
#           #
#  n=1000   #
#           #
#############

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a60_n1000 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=1000,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a60_n1000 complete")


ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a70_n1000 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=1000,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a70_n1000 complete")



ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a80_n1000 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=1000,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a80_n1000 complete")


#############
#           #
#  n=2000   #
#           #
#############

id <- 1:1000

ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a60_n2000 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=2000,alpha=c(.6,.4))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a60_n2000 complete")


ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a70_n2000 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=2000,alpha=c(.7,.3))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a70_n2000 complete")



ncr <- detectCores(logical=FALSE)

sink("/dev/null")
GA_a80_n2000 <- mclapply(id, function(r){
  DAT <- data_gamma_mix(id=1000+r,param=as.list(M),n=2000,alpha=c(.8,.2))
  
  ### Fit 
  
  # Naive GAMLSS (no mixture)
  NAIV <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT),trace=F)
  # GAMLSS with full class information (gold standard)
  GOLD <- gamlss(y~bs(x,Boundary.knots=c(0,1)),y~bs(x,Boundary.knots=c(0,1)),
                 family=GA(),data=as.data.frame(DAT)[DAT$z==1,],trace=F)
  
  # Latent Class Distributional Regression
  ### Initial weights for LCDR based on naive fit
  tmp <- pGA(DAT$y,NAIV$mu.fv,NAIV$sigma.fv)
  INIT <- matrix(c(1-tmp,tmp),ncol=2)
  LCDR <- lcdr_gamma(formula=list(mu1=y~bs(x,Boundary.knots=c(0,1)),sigma1=y~bs(x,Boundary.knots=c(0,1)),
                               mu2=y~bs(x,Boundary.knots=c(0,1)),sigma2=y~bs(x,Boundary.knots=c(0,1))),
                  data=DAT,init=INIT,run=r)
  
  return(list(LCDR=rbind(LCDR$mucoefs[[1]],LCDR$sicoefs[[1]]),
              NAIV=rbind(coef(NAIV,parameter="mu"),coef(NAIV,parameter="sigma")),
              GOLD=rbind(coef(GOLD,parameter="mu"),coef(GOLD,parameter="sigma")),
              LCDR.alpha=LCDR$alpha,LCDR.error=LCDR$error))
},mc.cores = ncr)
sink(); print("GA_a80_n2000 complete")

#####