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

### Results and Figures

source("lcdr_functions.R")

## Figure 1: Example of simulated data

M <- list(mu1="x+10*sin((x-.5)*sqrt(12)*pi/2)",sigma1="8+5*x",
           mu2="c+(c+1)*x+10*sin((x-.5)*sqrt(12)*pi/2)",sigma2="11+9*x")

a <- c(.7,.3)

Expl_dat <- data_mixture(1234,as.list(gsub("c","20",M)),alpha=a,n=500)

pdf("../results/GaussMix_Example.pdf",width = 6.5, height = 5)

layout(matrix(c(1,2,1,3,1,4),2),heights = c(.55,.45))

par(mar=c(4,4,2,4)+.1)
plot(Expl_dat$x,Expl_dat$y,col=c(1,"grey40")[Expl_dat$z],pch=c(20,3)[Expl_dat$z],
     ylim=c(-35,95),xaxs="i",xlim=0:1,las=1,ylab="",xlab="")
lines(seq(0,1,.01),seq(0,1,.01)+10*sin((seq(0,1,.01)-.5)*sqrt(12)*pi/2),lwd=2)
lines(seq(0,1,.01),20+21*seq(0,1,.01)+10*sin((seq(0,1,.01)-.5)*sqrt(12)*pi/2),lwd=2,col="grey40")
mtext("y",2,2.5,font=4,las=1,cex=.8);mtext("x",1,2,font=4,las=1,cex=.8)
abline(v=c(.1,.5,.9),lty=2)

par(mar=c(4,3.5,1,2)+.1)
curve(dnorm(x,.1+10*sin((.1-.5)*sqrt(12)*pi/2),8+5*.1)*a[1]+
        dnorm(x,20+21*.1+10*sin((.1-.5)*sqrt(12)*pi/2),11+9*.1)*a[2],-35,95,
      lwd=2,ylim=c(0,.035),ylab="",xlab="",yaxt="n",frame=T)
curve(dnorm(x,.1+10*sin((.1-.5)*sqrt(12)*pi/2),8+5*.1)*a[1],lty=1,add=T)
curve(dnorm(x,20+21*.1+10*sin((.1-.5)*sqrt(12)*pi/2),11+9*.1)*a[2],lty=1,col="grey40",add=T)
mtext("y",1,2,font=4,las=1,cex=.9)#;legend("topright",legend=expression(pdf(y*"|"*x==.1)),bty="n")
mtext(expression(f(y*"|"*x==.1)),2,.5,font=4,cex=.8)

par(mar=c(4,2.75,1,2.75)+.1)
curve(dnorm(x,.5+10*sin(0),8+5*.5)*a[1]+
        dnorm(x,20+21*.5+10*sin(0),11+9*.5)*a[2],-35,95,
      lwd=2,ylim=c(0,.035),ylab="",xlab="",yaxt="n",frame=T)
curve(dnorm(x,.5+10*sin(0),8+5*.5)*a[1],lty=1,add=T)
curve(dnorm(x,20+21*.5+10*sin(0),11+9*.5)*a[2],lty=1,col="grey40",add=T)
mtext("y",1,2,font=4,las=1,cex=.9)#;legend("topright",legend=expression(pdf(y*"|"*x==.9)),bty="n")
mtext(expression(f(y*"|"*x==.5)),2,.5,font=4,cex=.8)

par(mar=c(4,2,1,3.5)+.1)
curve(dnorm(x,.9+10*sin((.9-.5)*sqrt(12)*pi/2),8+5*.9)*a[1]+
        dnorm(x,20+21*.9+10*sin((.9-.5)*sqrt(12)*pi/2),11+9*.9)*a[2],-35,95,
      lwd=2,ylim=c(0,.035),ylab="",xlab="",yaxt="n",frame=T)
curve(dnorm(x,.9+10*sin((.9-.5)*sqrt(12)*pi/2),8+5*.9)*a[1],lty=1,add=T)
curve(dnorm(x,20+21*.9+10*sin((.9-.5)*sqrt(12)*pi/2),11+9*.9)*a[2],lty=1,col="grey40",add=T)
mtext("y",1,2,font=4,las=1,cex=.9)#;legend("topright",legend=expression(pdf(y*"|"*x==.9)),bty="n")
mtext(expression(f(y*"|"*x==.9)),2,.5,font=4,cex=.8)

dev.off()



### The following section relies on the simulation results
# NOTE: This will take some time to complete
# source("lcdr_simulation.R")

### Use stored data:
load("../data/SimResults.Rdata")


### Aggregate results

sim_names <- grep("M_\\.*",ls(),value=T)
sim_names <- sim_names[c(grep("c5",sim_names),grep("c10",sim_names),
                         grep("c15",sim_names),grep("c20",sim_names))]
sim_names <- sim_names[c(grep("a60",sim_names),grep("a70",sim_names),grep("a80",sim_names))]
sim_names <- sim_names[c(grep("n500",sim_names),grep("n1000",sim_names))]

# save(list=sim_names, file = "../results/SimResults.Rdata")

# Errors/Warnings

errors <- !sapply(sim_names,function(d){
  sapply(get(d), function(x) is.null(x$LCDR.error))
})

mean(errors)
colSums(errors)

# Integrated error:

LCDR.IE <- sapply(sim_names,function(d){
  sapply(get(d), IE.eval,p=.95,"LCDR",M)
})
LCDR.IE[errors] <- NA

GOLD.IE <- sapply(sim_names,function(d){
  sapply(get(d), IE.eval,p=.95,"GOLD",M)
})
GOLD.IE[errors] <- NA

NAIV.IE <- sapply(sim_names,function(d){
  sapply(get(d), IE.eval,p=.95,"NAIV",M)
})
NAIV.IE[errors] <- NA

RES.IE <- cbind(paste0(round(colMeans(LCDR.IE,na.rm=T),3)," (",round(apply(LCDR.IE,2,sd,na.rm=T),3),")"),
      paste0(round(colMeans(GOLD.IE,na.rm=T),3)," (",round(apply(GOLD.IE,2,sd,na.rm=T),3),")"),
      paste0(round(colMeans(NAIV.IE,na.rm=T),3)," (",round(apply(NAIV.IE,2,sd,na.rm=T),3),")"))

row.names(RES.IE) <- substr(sim_names,3,20)
colnames(RES.IE) <- c("LCDR","GOLD","NAIV")

write.csv(RES.IE,file="../results/GaussMix_IntegratedError.csv")


# Integrated squared error:

LCDR.ISE <- sapply(sim_names,function(d){
  sapply(get(d), ISE.eval,p=.95,"LCDR",M)
})
LCDR.ISE[errors] <- NA

GOLD.ISE <- sapply(sim_names,function(d){
  sapply(get(d), ISE.eval,p=.95,"GOLD",M)
})
GOLD.ISE[errors] <- NA

NAIV.ISE <- sapply(sim_names,function(d){
  sapply(get(d), ISE.eval,p=.95,"NAIV",M)
})
NAIV.ISE[errors] <- NA

RES.ISE <- cbind(paste0(round(colMeans(LCDR.ISE,na.rm=T),3)," (",round(apply(LCDR.ISE,2,sd,na.rm=T),3),")"),
                 paste0(round(colMeans(GOLD.ISE,na.rm=T),3)," (",round(apply(GOLD.ISE,2,sd,na.rm=T),3),")"),
                 paste0(round(colMeans(NAIV.ISE,na.rm=T),3)," (",round(apply(NAIV.ISE,2,sd,na.rm=T),3),")"))

row.names(RES.ISE) <- substr(sim_names,3,20)
colnames(RES.ISE) <- c("LCDR","GOLD","NAIV")

write.csv(RES.ISE,file="../results/GaussMix_IntegratedSquaredError.csv")

    

### Figures 2 and 3: Estimated Quantiles (converged only)

################################### 500
    
pdf("../results/GaussMix_n500.pdf",width = 6.5, height = 4, pointsize = 6.5)
    
x <- seq(0,1,.01)

layout(matrix(c(1,2,3,4),nrow=2,byrow=T),heights = c(2.5,2.5))
    
par(mar=c(0.5,3.5,3.5,0.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
axis(3);mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
axis(3);legend("top",legend = "c=5",bty="n",text.font = 2)

matlines(x,qnorm(.95,
                 sapply(M_a60_c5_n500[!errors[,"M_a60_c5_n500"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c5_n500[!errors[,"M_a60_c5_n500"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(0.5,0.5,3.5,3.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(3);axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
axis(3);axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a60_c20_n500[!errors[,"M_a60_c20_n500"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c20_n500[!errors[,"M_a60_c20_n500"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(3.5,3.5,0.5,0.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
legend("top",legend = "c=5",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c5_n500[!errors[,"M_a80_c5_n500"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c5_n500[!errors[,"M_a80_c5_n500"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###


par(mar=c(3.5,0.5,0.5,3.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c20_n500[!errors[,"M_a80_c20_n500"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c20_n500[!errors[,"M_a80_c20_n500"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

dev.off()

################################### 1000

pdf("../results/GaussMix_n1000.pdf",width = 6.5, height = 4, pointsize = 6.5)

x <- seq(0,1,.01)

layout(matrix(c(1,2,3,4),nrow=2,byrow=T),heights = c(2.5,2.5))

par(mar=c(0.5,3.5,3.5,0.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
axis(3);mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
axis(3);legend("top",legend = "c=5",bty="n",text.font = 2)

matlines(x,qnorm(.95,
                 sapply(M_a60_c5_n1000[!errors[,"M_a60_c5_n1000"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c5_n1000[!errors[,"M_a60_c5_n1000"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(0.5,0.5,3.5,3.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(3);axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
axis(3);axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a60_c20_n1000[!errors[,"M_a60_c20_n1000"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c20_n1000[!errors[,"M_a60_c20_n1000"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(3.5,3.5,0.5,0.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
legend("top",legend = "c=5",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c5_n1000[!errors[,"M_a80_c5_n1000"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c5_n1000[!errors[,"M_a80_c5_n1000"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(3.5,0.5,0.5,3.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c20_n1000[!errors[,"M_a80_c20_n1000"]][1:100],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c20_n1000[!errors[,"M_a80_c20_n1000"]][1:100],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

dev.off()

### Supplementary figure for reviewers
    
pdf("../results/GaussMix_not_conv.pdf",width = 6.5, height = 3, pointsize = 7)

x <- seq(0,1,.01)

layout(matrix(c(1,2,5,6,3,4,7,8),nrow=2,byrow=T),heights = c(2.5,2.5))

par(mar=c(0.5,3.5,3.5,0.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
axis(3);mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a60_c5_n500[errors[,"M_a60_c5_n500"]])),
       bty="n",text.font = 2)
axis(3);legend("top",legend = "c=5",bty="n",text.font = 2)

matlines(x,qnorm(.95,
                 sapply(M_a60_c5_n500[errors[,"M_a60_c5_n500"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c5_n500[errors[,"M_a60_c5_n500"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###

par(mar=c(0.5,0.5,3.5,3.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(3);axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a60_c20_n500[errors[,"M_a60_c20_n500"]])),
       bty="n",text.font = 2)
axis(3);axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a60_c20_n500[errors[,"M_a60_c20_n500"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c20_n500[errors[,"M_a60_c20_n500"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###

par(mar=c(3.5,3.5,0.5,0.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
legend("top",legend = "c=5",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a80_c5_n500[errors[,"M_a80_c5_n500"]])),
       bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c5_n500[errors[,"M_a80_c5_n500"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c5_n500[errors[,"M_a80_c5_n500"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###

par(mar=c(3.5,0.5,0.5,3.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a80_c20_n500[errors[,"M_a80_c20_n500"]])),
       bty="n",text.font = 2)
axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c20_n500[errors[,"M_a80_c20_n500"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c20_n500[errors[,"M_a80_c20_n500"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###
    
x <- seq(0,1,.01)

par(mar=c(0.5,3.5,3.5,0.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
axis(3);mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a60_c5_n1000[errors[,"M_a60_c5_n1000"]])),
       bty="n",text.font = 2)
axis(3);legend("top",legend = "c=5",bty="n",text.font = 2)

matlines(x,qnorm(.95,
                 sapply(M_a60_c5_n1000[errors[,"M_a60_c5_n1000"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a60_c5_n1000[errors[,"M_a60_c5_n1000"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###

par(mar=c(0.5,0.5,3.5,3.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(3);axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a60_c20_n1000[errors[,"M_a60_c20_n1000"]])),
       bty="n",text.font = 2)
axis(3);axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
# matlines(x,qnorm(.95,
#                  sapply(M_a60_c20_n1000[errors[,"M_a60_c20_n1000"]],
#                         function(f){
#                           cbind(1,bs(x))%*%f[["LCDR"]][1,]
#                         }),
#                  sapply(M_a60_c20_n1000[errors[,"M_a60_c20_n1000"]],
#                          function(f){
#                            exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
#                         })
#  ),lty=1,col=rgb(.75,0,0,.5)
#  )
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###

par(mar=c(3.5,3.5,0.5,0.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",ylim=c(-40,100),las=1)
mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a80_c5_n1000[errors[,"M_a80_c5_n1000"]])),
       bty="n",text.font = 2)
legend("top",legend = "c=5",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c5_n1000[errors[,"M_a80_c5_n1000"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c5_n1000[errors[,"M_a80_c5_n1000"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

###

par(mar=c(3.5,0.5,0.5,3.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",yaxt="n",ylim=c(-40,100))
axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=1000",bty="n",text.font = 2)
legend("bottom",legend = paste("# non-convergence:",length(M_a80_c20_n1000[errors[,"M_a80_c20_n1000"]])),
       bty="n",text.font = 2)
axis(4,las=1);legend("top",legend = "c=20",bty="n",text.font = 2)
matlines(x,qnorm(.95,
                 sapply(M_a80_c20_n1000[errors[,"M_a80_c20_n1000"]],
                        function(f){
                          cbind(1,bs(x))%*%f[["LCDR"]][1,]
                        }),
                 sapply(M_a80_c20_n1000[errors[,"M_a80_c20_n1000"]],
                        function(f){
                          exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                        })
),lty=1,col=rgb(.75,0,0,.5)
)
lines(x,qnorm(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

dev.off()

    
##########################
#
#  Supplementary R Code
#
##########################

### Gamma mixture

## Figure A1: Example of simulated data

M <- list(mu1="5+x-2.5*(x-1)^2",sigma1="exp(-1+x-x^2-.5*x^3)",
          mu2="9+2*x-2.5*(x-1)^2",sigma2="exp(-1+x-x^2-.5*x^3)")


a <- c(.7,.3)

Expl_dat <- data_gamma_mix(1234,as.list(M),alpha=a,n=500)

pdf("../results/GammaMix_Example.pdf",width = 6.5, height = 5)

layout(matrix(c(1,2,1,3,1,4),2),heights = c(.55,.45))

par(mar=c(4,4,2,4)+.1)
plot(Expl_dat$x,Expl_dat$y,col=c(1,"grey40")[Expl_dat$z],pch=c(20,3)[Expl_dat$z],
     ylim=c(0,27),xaxs="i",xlim=0:1,las=1,ylab="",xlab="")
mtext("y",2,2.5,font=4,las=1,cex=.8);mtext("x",1,2,font=4,las=1,cex=.8)
lines(seq(0,1,.01),5+seq(0,1,.01)-2.5*(seq(0,1,.01)-1)^2,lwd=2)
lines(seq(0,1,.01),9+2*seq(0,1,.01)-2.5*(seq(0,1,.01)-1)^2,lwd=2,col="grey40")
abline(v=c(.1,.5,.9),lty=2)

par(mar=c(4,3.5,1,2)+.1)
curve(dGA(x,5+.1-2.5*(.1-1)^2,exp(-1+.1-.1^2-.5*.1^3))*a[1]+
        dGA(x,9+2*.1-2.5*(.1-1)^2,exp(-1+.1-.1^2-.5*.1^3))*a[2],0,25,
      lwd=2,ylim=c(0,.25),ylab="",xlab="",yaxt="n",frame=T)
curve(dGA(x,5+.1-2.5*(.1-1)^2,exp(-1+.1-.1^2-.5*.1^3))*a[1],0,27,lty=1,add=T)
curve(dGA(x,9+2*.1-2.5*(.1-1)^2,exp(-1+.1-.1^2-.5*.1^3))*a[2],0,27,lty=1,col="grey40",add=T)
mtext("y",1,2,font=4,las=1,cex=.9)#;legend("topright",legend=expression(pdf(y*"|"*x==.1)),bty="n")
mtext(expression(f(y*"|"*x==.1)),2,.5,font=4,cex=.8)

par(mar=c(4,2.75,1,2.75)+.1)
curve(dGA(x,5+.5-2.5*(.5-1)^2,exp(-1+.5-.5^2-.5*.5^3))*a[1]+
        dGA(x,9+2*.5-2.5*(.5-1)^2,exp(-1+.5-.5^2-.5*.5^3))*a[2],0,25,
      lwd=2,ylim=c(0,.25),ylab="",xlab="",yaxt="n",frame=T)
curve(dGA(x,5+.5-2.5*(.5-1)^2,exp(-1+.5-.5^2-.5*.5^3))*a[1],0,27,lty=1,add=T)
curve(dGA(x,9+2*.5-2.5*(.5-1)^2,exp(-1+.5-.5^2-.5*.5^3))*a[2],0,27,lty=1,col="grey40",add=T)
mtext("y",1,2,font=4,las=1,cex=.9)#;legend("topright",legend=expression(pdf(y*"|"*x==.5)),bty="n")
mtext(expression(f(y*"|"*x==.5)),2,.5,font=4,cex=.8)

par(mar=c(4,2,1,3.5)+.1)
curve(dGA(x,5+.9-2.5*(.9-1)^2,exp(-1+.9-.9^2-.5*.9^3))*a[1]+
        dGA(x,9+2*.9-2.5*(.9-1)^2,exp(-1+.9-.9^2-.5*.9^3))*a[2],0,25,
      lwd=2,ylim=c(0,.25),ylab="",xlab="",yaxt="n",frame=T)
curve(dGA(x,5+.9-2.5*(.9-1)^2,exp(-1+.9-.9^2-.5*.9^3))*a[1],0,27,lty=1,add=T)
curve(dGA(x,9+2*.9-2.5*(.9-1)^2,exp(-1+.9-.9^2-.5*.9^3))*a[2],0,27,lty=1,col="grey40",add=T)
mtext("y",1,2,font=4,las=1,cex=.9)#;legend("topright",legend=expression(pdf(y*"|"*x==.9)),bty="n")
mtext(expression(f(y*"|"*x==.9)),2,.5,font=4,cex=.8)

dev.off()


### The following section relies on the simulation results
# NOTE: This will take some time to complete
# source("lcdr_simulation_gamma.R")
    
### Use stored data:    
load("../data/SimResults_gamma.Rdata")
    
### Aggregate results

sim_names <- grep("GA_",ls(),value=T)
sim_names <- sim_names[c(grep("a60",sim_names),grep("a70",sim_names),grep("a80",sim_names))]
sim_names <- sim_names[c(grep("n500",sim_names),grep("n1000",sim_names),grep("n2000",sim_names))]

save(list=sim_names, file = "../results/SimResults_gamma.Rdata")

# Errors/Warnings

errors <- !sapply(sim_names,function(d){
  sapply(get(d), function(x) is.null(x$LCDR.error))
})

mean(errors)
colSums(errors)

# Integrated error:

LCDR.IE <- sapply(sim_names,function(d){
  sapply(get(d), IE.eval_gamma,p=.95,"LCDR",M)
})
LCDR.IE[errors] <- NA

GOLD.IE <- sapply(sim_names,function(d){
  sapply(get(d), IE.eval_gamma,p=.95,"GOLD",M)
})
GOLD.IE[errors] <- NA

NAIV.IE <- sapply(sim_names,function(d){
  sapply(get(d), IE.eval_gamma,p=.95,"NAIV",M)
})
NAIV.IE[errors] <- NA

RES.IE <- cbind(paste0(round(colMeans(LCDR.IE,na.rm=T),3)," (",round(apply(LCDR.IE,2,sd,na.rm=T),3),")"),
                paste0(round(colMeans(GOLD.IE,na.rm=T),3)," (",round(apply(GOLD.IE,2,sd,na.rm=T),3),")"),
                paste0(round(colMeans(NAIV.IE,na.rm=T),3)," (",round(apply(NAIV.IE,2,sd,na.rm=T),3),")"))

row.names(RES.IE) <- substr(sim_names,4,20)
colnames(RES.IE) <- c("LCDR","GOLD","NAIV")

write.csv(RES.IE,file="../results/GammaMix_IntegratedError.csv")


# Integrated squared error:

LCDR.ISE <- sapply(sim_names,function(d){
  sapply(get(d), ISE.eval_gamma,p=.95,"LCDR",M)
})
LCDR.ISE[errors] <- NA

GOLD.ISE <- sapply(sim_names,function(d){
  sapply(get(d), ISE.eval_gamma,p=.95,"GOLD",M)
})
GOLD.ISE[errors] <- NA

NAIV.ISE <- sapply(sim_names,function(d){
  sapply(get(d), ISE.eval_gamma,p=.95,"NAIV",M)
})
NAIV.ISE[errors] <- NA

RES.ISE <- cbind(paste0(round(colMeans(LCDR.ISE,na.rm=T),3)," (",round(apply(LCDR.ISE,2,sd,na.rm=T),3),")"),
                 paste0(round(colMeans(GOLD.ISE,na.rm=T),3)," (",round(apply(GOLD.ISE,2,sd,na.rm=T),3),")"),
                 paste0(round(colMeans(NAIV.ISE,na.rm=T),3)," (",round(apply(NAIV.ISE,2,sd,na.rm=T),3),")"))

row.names(RES.ISE) <- substr(sim_names,4,20)
colnames(RES.ISE) <- c("LCDR","GOLD","NAIV")

write.csv(RES.ISE,file="../results/GammaMix_IntegratedSquaredError.csv")




### Figure A2: Estimated Quantiles (converged only)

################################### 500

pdf("../results/GammaMix_Results.pdf",width = 6.5, height = 4, pointsize = 6.5)

x <- seq(0,1,.01)

layout(matrix(c(1,2,3,4),nrow=2,byrow=T),heights = c(2.5,2.5))

par(mar=c(0.5,3.5,3.5,0.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",ylim=c(0,20),las=1)
axis(3);mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)

matlines(x,qGA(.95,
               sapply(GA_a60_n500[!errors[,"GA_a60_n500"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][1,])
                      }),
               sapply(GA_a60_n500[!errors[,"GA_a60_n500"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                      })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qGA(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(0.5,0.5,3.5,3.5)+.1)
plot(NA,xlab="",xaxt="n",xlim=c(0,1),ylab="",yaxt="n",ylim=c(0,20))
axis(3);axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",3,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.6),bty="n",text.font = 2)
legend("topright",legend = "n=2000",bty="n",text.font = 2)
matlines(x,qGA(.95,
               sapply(GA_a60_n2000[!errors[,"GA_a60_n2000"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][1,])
                      }),
               sapply(GA_a60_n2000[!errors[,"GA_a60_n2000"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                      })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qGA(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(3.5,3.5,0.5,0.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",ylim=c(0,20),las=1)
mtext("y",2,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=500",bty="n",text.font = 2)
matlines(x,qGA(.95,
               sapply(GA_a80_n500[!errors[,"GA_a80_n500"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][1,])
                      }),
               sapply(GA_a80_n500[!errors[,"GA_a80_n500"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                      })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qGA(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)


###

par(mar=c(3.5,0.5,0.5,3.5)+.1)
plot(NA,xlab="",xlim=c(0,1),ylab="",yaxt="n",ylim=c(0,20))
axis(4,las=1);mtext("y",4,2.5,font=4,las=1,cex=1);mtext("x",1,2.5,font=4,las=1,cex=1)
legend("topleft",legend = expression(alpha[1]==0.8),bty="n",text.font = 2)
legend("topright",legend = "n=2000",bty="n",text.font = 2)
matlines(x,qGA(.95,
               sapply(GA_a80_n2000[!errors[,"GA_a80_n2000"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][1,])
                      }),
               sapply(GA_a80_n2000[!errors[,"GA_a80_n2000"]][1:100],
                      function(f){
                        exp(cbind(1,bs(x))%*%f[["LCDR"]][2,])
                      })
),lty=1,col=rgb(.5,.5,.5,.5)
)
lines(x,qGA(.95,eval(parse(text=M$mu1)),eval(parse(text=M$sigma1))),col=1,lty=1,lwd=2)

dev.off()


