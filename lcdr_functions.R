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

### Custom functions

#----------------------------------------------------------
require(gamlss)
#----------------------------------------------------------

####################################
# Visualization and data generation

data_mixture <- function(id,param,n,alpha=NULL,...){
  
  M <- length(param)/2
  
  if(is.null(alpha)){
    alpha <- rep(1/M,M)
  } else if(sum(alpha)!=1){
    alpha <- alpha/sum(alpha)
  }
  
  set.seed(id)
  X <- Y <- Z <- numeric()
    
  for(m in 1:M){
  x <- sort(runif(floor(n*alpha[m])))
  y <- rnorm(floor(n*alpha[m]),eval(parse(text=param[2*m-1])),eval(parse(text=param[2*m])))
  z <- rep(m,floor(n*alpha[m]))
  
  X <- c(X,x); Y <- c(Y,y); Z <- c(Z,z)
  }
list(x=X,y=Y,z=Z)
}

###################
# Model estimation

lcdr <- function(formula,data=NULL,init=NULL,threshold=0.001,tstart=NULL,
                            run=NULL){ #cluster only
  
  tst <- tstart
  
  data <- data.frame(data)
  M <- length(formula)/2
  
  # initial weights
  if(!is.matrix(init)){
    stop(cat("Please provide weights for initialization in matrix format"))
  }

  # models
  MODS <- lapply(1:M,function(m){
    gamlss(formula = eval(parse(text=formula[2*m-1])),
           sigma.formula = eval(parse(text=formula[2*m])),
           data=data.frame(data),weights=init[,m],trace=F)
  })
  
  # mixt. weights
  A <- apply(init,2,sum)/nrow(data)
  
  # cond. pdf of each comp.
  PMAT <- sapply(1:M,function(m){
    dnorm(data$y,MODS[[m]]$mu.fv,MODS[[m]]$sigma.fv)*A[m]
  })
  
  # New model weights
  WMAT <- PMAT/rowSums(PMAT)
  
  lL <- c(-Inf,sum(log(rowSums(PMAT))))
  lLvec <- lL[2]
  
  error <- tryCatch({
    while((diff(lL))>threshold){
      
      # M-Step
      
      # models
      MODS <- lapply(1:M,function(m){
        gamlss(formula = eval(parse(text=formula[2*m-1])),
               sigma.formula = eval(parse(text=formula[2*m])),
               data=data.frame(data),weights=WMAT[,m],start.from=MODS[[m]],trace=F)
      })
      
      # mixt. weights
      A <- apply(WMAT,2,sum)/nrow(data)
      
      # E-Step
      
      # cond. pdf of each comp.
      PMAT <- sapply(1:M,function(m){
        dnorm(data$y,MODS[[m]]$mu.fv,MODS[[m]]$sigma.fv)*A[m]
      })
      
      # New model weights
      WMAT <- PMAT/rowSums(PMAT)
      
      lL[1] <- lL[2]
      lL[2] <- sum(log(rowSums(PMAT)))
      
      lLvec <- c(lLvec,lL[2])
    }
  },error=function(e){
    cat("ID",run,":",conditionMessage(e), "\n")
    return(conditionMessage(e))
  })

  FITARR <- array(NA,c(nrow(data),2,M))
  FITARR[,1,] <- sapply(1:M,function(m){MODS[[m]]$mu.fv})
  FITARR[,2,] <- sapply(1:M,function(m){MODS[[m]]$sigma.fv})
  
  dimnames(FITARR)[[2]] <- c("mu","sigma")
  dimnames(FITARR)[[3]] <- paste0("Comp",1:M)
  
  list(mucoefs=lapply(1:M,function(m){coef(MODS[[m]],what="mu")}),
       sicoefs=lapply(1:M,function(m){coef(MODS[[m]],what="sigma")}),
       fit=FITARR,lL=lLvec,alpha=A,data=data,error=error)

}


#############################
# Integrated (squared) error

IE.eval <- function(obj, method, p, true, low=0, up=1){
  integrate(function(x){
    qnorm(p,cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][1,],
          exp(cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][2,]))-
      qnorm(p,eval(parse(text=true[1])),eval(parse(text=true[2])))},lower=low,upper=up)$value
}

ISE.eval <- function(obj, method, p, true, low=0, up=1){
  integrate(function(x){
    (qnorm(p,cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][1,],
           exp(cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][2,]))-
       qnorm(p,eval(parse(text=true[1])),eval(parse(text=true[2]))))^2},lower=low,upper=up)$value
}


##########################
#
#  Supplementary Code
#
##########################


### Gamma mixture (custom functions)

#----------------------------------------------------------
require(gamlss)
#----------------------------------------------------------

####################################
# Visualization and data generation

data_gamma_mix <- function(id,param,n,alpha=NULL,...){
  
  M <- length(param)/2
  
  if(is.null(alpha)){
    alpha <- rep(1/M,M)
  } else if(sum(alpha)!=1){
    alpha <- alpha/sum(alpha)
  }
  
  set.seed(id)
  X <- Y <- Z <- numeric()
  
  for(m in 1:M){
    x <- sort(runif(floor(n*alpha[m])))
    y <- rGA(floor(n*alpha[m]),eval(parse(text=param[2*m-1])),eval(parse(text=param[2*m])))
    z <- rep(m,floor(n*alpha[m]))
    
    X <- c(X,x); Y <- c(Y,y); Z <- c(Z,z)
  }
  
  list(x=X,y=Y,z=Z)
}

###################
# Model estimation

lcdr_gamma <- function(formula,data=NULL,init=NULL,threshold=0.001,tstart=NULL,
                       run=NULL){ #cluster only
  
  tst <- tstart
  
  data <- data.frame(data)
  M <- length(formula)/2
  
  # initial weights
  if(!is.matrix(init)){
    stop(cat("Please provide weights for initialization in matrix format"))
  }
  
  # models
  MODS <- lapply(1:M,function(m){
    gamlss(formula = eval(parse(text=formula[2*m-1])),
           sigma.formula = eval(parse(text=formula[2*m])),family=GA(),
           data=data.frame(data),weights=init[,m],trace=F)
  })
  
  # mixt. weights
  A <- apply(init,2,sum)/nrow(data)
  
  # cond. pdf of each comp.
  PMAT <- sapply(1:M,function(m){
    dGA(data$y,MODS[[m]]$mu.fv,MODS[[m]]$sigma.fv)*A[m]
  })
  
  # New model weights
  WMAT <- PMAT/rowSums(PMAT)
  
  lL <- c(-Inf,sum(log(rowSums(PMAT))))
  lLvec <- lL[2]
  
  error <- tryCatch({
    while((diff(lL))>threshold){
      
      # M-Step
      
      # models
      MODS <- lapply(1:M,function(m){
        gamlss(formula = eval(parse(text=formula[2*m-1])),
               sigma.formula = eval(parse(text=formula[2*m])),family=GA(),
               data=data.frame(data),weights=WMAT[,m],start.from=MODS[[m]],trace=F)
      })
      
      # mixt. weights
      A <- apply(WMAT,2,sum)/nrow(data)
      
      # E-Step
      
      # cond. pdf of each comp.
      PMAT <- sapply(1:M,function(m){
        dGA(data$y,MODS[[m]]$mu.fv,MODS[[m]]$sigma.fv)*A[m]
      })
      
      # New model weights
      WMAT <- PMAT/rowSums(PMAT)
      
      lL[1] <- lL[2]
      lL[2] <- sum(log(rowSums(PMAT)))
      
      lLvec <- c(lLvec,lL[2])
    }
  },error=function(e){
    cat("ID",run,":",conditionMessage(e), "\n")
    return(conditionMessage(e))
  })
  
  FITARR <- array(NA,c(nrow(data),2,M))
  FITARR[,1,] <- sapply(1:M,function(m){MODS[[m]]$mu.fv})
  FITARR[,2,] <- sapply(1:M,function(m){MODS[[m]]$sigma.fv})
  
  dimnames(FITARR)[[2]] <- c("mu","sigma")
  dimnames(FITARR)[[3]] <- paste0("Comp",1:M)
  
  list(mucoefs=lapply(1:M,function(m){coef(MODS[[m]],what="mu")}),
       sicoefs=lapply(1:M,function(m){coef(MODS[[m]],what="sigma")}),
       fit=FITARR,lL=lLvec,alpha=A,data=data,error=error)
  
}


#############################
# Integrated (squared) error

IE.eval_gamma <- function(obj, method, p, true, low=0, up=1){
  integrate(function(x){
    qGA(p,exp(cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][1,]),
        exp(cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][2,]))-
      qGA(p,eval(parse(text=true[1])),eval(parse(text=true[2])))},lower=low,upper=up)$value
}

ISE.eval_gamma <- function(obj, method, p, true, low=0, up=1){
  integrate(function(x){
    (qGA(p,exp(cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][1,]),
         exp(cbind(1,bs(x,Boundary.knots=c(0,1)))%*%obj[[method]][2,]))-
       qGA(p,eval(parse(text=true[1])),eval(parse(text=true[2]))))^2},lower=low,upper=up)$value
}
