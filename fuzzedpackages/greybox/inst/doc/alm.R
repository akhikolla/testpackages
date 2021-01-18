## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  fig.height=4,
  fig.width=6
)

library(greybox)

## ----pdfNormal, echo=FALSE----------------------------------------------------
plot(seq(-5,5,0.1),dnorm(seq(-5,5,0.1)),type="l",
     xlab="y_t",ylab="Density",main="PDF of Normal distribution")
lines(seq(-5,5,0.1),dnorm(seq(-5,5,0.1),0,2), col="blue")
lines(seq(-5,5,0.1),dnorm(seq(-5,5,0.1),1,1), col="red")
legend("topright",legend=c("N(0,1)","N(0,2)","N(1,1)"), lwd=1, col=c("black","blue","red"))

## ----pdfLaplace, echo=FALSE---------------------------------------------------
plot(seq(-5,5,0.1),dlaplace(seq(-5,5,0.1)),type="l",
     xlab="y_t",ylab="Density",main="PDF of Laplace distribution")
lines(seq(-5,5,0.1),dlaplace(seq(-5,5,0.1),0,2), col="blue")
lines(seq(-5,5,0.1),dlaplace(seq(-5,5,0.1),1,1), col="red")
legend("topright",legend=c("Laplace(0,1)","Laplace(0,2)","Laplace(1,1)"), lwd=1, col=c("black","blue","red"))

## ----pdfALaplace, echo=FALSE--------------------------------------------------
plot(seq(-5,5,0.1),dalaplace(seq(-5,5,0.1),0,0.5),type="l",
     xlab="y_t",ylab="Density",main="PDF of Asymmetric Laplace distribution")
lines(seq(-5,5,0.1),dalaplace(seq(-5,5,0.1),0,1,0.25), col="blue")
lines(seq(-5,5,0.1),dalaplace(seq(-5,5,0.1),0,1,0.75), col="red")
legend("topright",legend=c("ALaplace(0,0.5,0.5)","ALaplace(0,1,0.25)","ALaplace(0,1,0.75)"), lwd=1, col=c("black","blue","red"))

## ----pdfS, echo=FALSE---------------------------------------------------------
plot(seq(-5,5,0.01),ds(seq(-5,5,0.01),0,1),type="l",
     xlab="y_t",ylab="Density",main="PDF of S distribution")
lines(seq(-5,5,0.01),ds(seq(-5,5,0.01),0,2), col="blue")
lines(seq(-5,5,0.01),ds(seq(-5,5,0.01),1,1), col="red")
legend("topright",legend=c("S(0,1)","S(0,2)","S(1,1)"), lwd=1, col=c("black","blue","red"))

## ----pdfgnorm, echo=FALSE-----------------------------------------------------
dgnorm <- function(x, mu = 0, alpha = 1, beta = 1,
                   log = FALSE) {
    # A failsafe for NaN / NAs of alpha / beta
    if(any(is.nan(alpha))){
        alpha[is.nan(alpha)] <- 0
    }
    if(any(is.na(alpha))){
        alpha[is.na(alpha)] <- 0
    }
    if(any(alpha<0)){
        alpha[alpha<0] <- 0
    }
    if(any(is.nan(beta))){
        beta[is.nan(beta)] <- 0
    }
    if(any(is.na(beta))){
        beta[is.na(beta)] <- 0
    }
    gnormValues <- (exp(-(abs(x-mu)/ alpha)^beta)* beta/(2*alpha*gamma(1/beta)))
    if(log){
        gnormValues[] <- log(gnormValues)
    }

    return(gnormValues)
}

plot(seq(-5,5,0.1),dgnorm(seq(-5,5,0.1),0,1,2),type="l",
     xlab="y_t",ylab="Density",main="PDF of Generalised Normal distribution")
lines(seq(-5,5,0.1),dgnorm(seq(-5,5,0.1),0,1,1), col="blue")
lines(seq(-5,5,0.01),dgnorm(seq(-5,5,0.01),0,1,0.5), col="red")
lines(seq(-5,5,0.1),dgnorm(seq(-5,5,0.1),0,1,100), col="purple")
legend("topright",legend=c("GN(0,1,2)","GN(0,1,1)","GN(0,1,0.5)","GN(0,1,100)"),
       lwd=1, col=c("black","blue","red","purple"))

## ----pdfLogis, echo=FALSE-----------------------------------------------------
plot(seq(-5,5,0.1),dlogis(seq(-5,5,0.1),0,0.5),type="l",
     xlab="y_t",ylab="Density",main="PDF of Logistic distribution")
lines(seq(-5,5,0.1),dlogis(seq(-5,5,0.1),0,2), col="blue")
lines(seq(-5,5,0.1),dlogis(seq(-5,5,0.1),1,1), col="red")
legend("topright",legend=c("Logis(0,0.5)","Logis(0,1)","Logis(1,0.5)"), lwd=1, col=c("black","blue","red"))

## ----pdfStudent, echo=FALSE---------------------------------------------------
plot(seq(-5,5,0.1),dt(seq(-5,5,0.1),100),type="l",
     xlab="y_t",ylab="Density",main="PDF of Student's t distribution")
lines(seq(-5,5,0.1),dt(seq(-5,5,0.1),10), col="blue")
lines(seq(-5,5,0.1),dt(seq(-5,5,0.1),1), col="red")
legend("topright",legend=c("t(100)","t(10)","t(1)"), lwd=1, col=c("black","blue","red"))

## ----normalDistributionData---------------------------------------------------
xreg <- cbind(rnorm(100,10,3),rnorm(100,50,5))
xreg <- cbind(500+0.5*xreg[,1]-0.75*xreg[,2]+rs(100,0,3),xreg,rnorm(100,300,10))
colnames(xreg) <- c("y","x1","x2","Noise")

inSample <- xreg[1:80,]
outSample <- xreg[-c(1:80),]

## ----normalRegression---------------------------------------------------------
ourModel <- alm(y~x1+x2, data=inSample, distribution="dnorm")
summary(ourModel)
plot(predict(ourModel,outSample,interval="p",level=c(0.9,0.95)))

## ----ALaplaceRegression-------------------------------------------------------
ourModel <- alm(y~x1+x2, data=inSample, distribution="dalaplace",alpha=0.95)
summary(ourModel)
plot(predict(ourModel,outSample))

## ----pdflognorm, echo=FALSE---------------------------------------------------
plot(seq(0,5,0.01),dlnorm(seq(0,5,0.01),0,1),type="l",
     xlab="y_t",ylab="Density",ylim=c(0,1.5),main="PDF of Log Normal distribution")
lines(seq(0,5,0.01),dlnorm(seq(0,5,0.01),1,1), col="blue")
lines(seq(0,5,0.01),dlnorm(seq(0,5,0.01),0,2), col="red")
legend("topright",legend=c("logN(0,1)","logN(1,1)","logN(0,2)"), lwd=1, col=c("black","blue","red"))

## ----pdfBCNorm, echo=FALSE----------------------------------------------------
plot(seq(0,5,0.1),dbcnorm(seq(0,5,0.1),0,1,1),type="l",ylim=c(0,1),
     xlab="y_t",ylab="Density",main="PDF of Box-Cox Normal distribution")
lines(seq(0.01,5,0.01),dbcnorm(seq(0.01,5,0.01),0,1,0.5), col="blue")
lines(seq(0,5,0.1),dbcnorm(seq(0,5,0.1),0,1,2), col="red")
lines(seq(0,5,0.01),dbcnorm(seq(0,5,0.01),0,1,0.01), col="purple")
legend("topright",legend=c("BCN(0,1,1)","BCN(0,1,0.5)","BCN(0,1,2)","BCN(0,1,0.01)"),
       lwd=1, col=c("black","blue","red","purple"))

## ----pdfIG, echo=FALSE--------------------------------------------------------
library(statmod)
plot(seq(0,5,0.01),dinvgauss(seq(0,5,0.01),1,1),type="l",
     xlab="y_t",ylab="Density",main="PDF of Inverse Gaussian distribution")
lines(seq(0.01,5,0.01),dinvgauss(seq(0.01,5,0.01),1,2), col="blue")
lines(seq(0,5,0.01),dinvgauss(seq(0,5,0.01),2,1), col="red")
legend("topright",legend=c("IG(1,1)","IG(1,2)","IG(2,1)"),
       lwd=1, col=c("black","blue","red","purple"))

## ----pdflogLaplace, echo=FALSE------------------------------------------------
plot(seq(0.01,5,0.01),dlaplace(log(seq(0.01,5,0.01)),0,1)/seq(0.01,5,0.01),type="l",ylim=c(0,1.5),
     xlab="y_t",ylab="Density",main="PDF of Log Laplace distribution")
lines(seq(0.01,5,0.01),dlaplace(log(seq(0.01,5,0.01)),0,2)/seq(0.01,5,0.01), col="blue")
lines(seq(0.01,5,0.01),dlaplace(log(seq(0.01,5,0.01)),1,1)/seq(0.01,5,0.01), col="red")
legend("topright",legend=c("logLaplace(0,1)","logLaplace(0,2)","logLaplace(1,1)"),
       lwd=1, col=c("black","blue","red"))

## ----pdflogS, echo=FALSE------------------------------------------------------
plot(seq(0.01,5,0.01),ds(log(seq(0.01,5,0.01)),0,1)/seq(0.01,5,0.01),type="l",ylim=c(0,1.5),
     xlab="y_t",ylab="Density",main="PDF of Log S distribution")
lines(seq(0.01,5,0.01),ds(log(seq(0.01,5,0.01)),0,2)/seq(0.01,5,0.01), col="blue")
lines(seq(0.01,5,0.01),ds(log(seq(0.01,5,0.01)),1,1)/seq(0.01,5,0.01), col="red")
legend("topright",legend=c("logS(0,1)","logS(0,2)","logS(1,1)"),
       lwd=1, col=c("black","blue","red"))

## ----pdflogGN, echo=FALSE-----------------------------------------------------
plot(seq(0.01,5,0.01),dgnorm(log(seq(0.01,5,0.01)),0,1,2)/seq(0.01,5,0.01),type="l",ylim=c(0,1.5),
     xlab="y_t",ylab="Density",main="PDF of Log Generalised Normal distribution")
lines(seq(0.01,5,0.01),dgnorm(log(seq(0.01,5,0.01)),0,1,1)/seq(0.01,5,0.01), col="blue")
lines(seq(0.01,5,0.01),dgnorm(log(seq(0.01,5,0.01)),0,1,0.5)/seq(0.01,5,0.01), col="red")
lines(seq(0.01,5,0.01),dgnorm(log(seq(0.01,5,0.01)),0,1,100)/seq(0.01,5,0.01), col="purple")
legend("topright",legend=c("logGN(0,1,2)","logGN(0,1,1)","logGN(0,1,0.5)","logGN(0,1,100)"),
       lwd=1, col=c("black","blue","red","purple"))

## ----pdfFnorm, echo=FALSE-----------------------------------------------------
plot(seq(0.01,5,0.01),dfnorm(seq(0.01,5,0.01),0,1),type="l",ylim=c(0,1.5),
     xlab="y_t",ylab="Density",main="PDF of Folded Normal distribution")
lines(seq(0.01,5,0.01),dfnorm(seq(0.01,5,0.01),-1,1), col="blue")
lines(seq(0.01,5,0.01),dfnorm(seq(0.01,5,0.01),-2,1), col="red")
legend("topright",legend=c("FN(0,1)","FN(-1,1)","FN(-2,1)"),
       lwd=1, col=c("black","blue","red"))

## ----pdfLogitnorm, echo=FALSE-------------------------------------------------
plot(seq(0.01,0.99,0.01),dlogitnorm(seq(0.01,0.99,0.01),0,1),type="l",ylim=c(0,5),
     xlab="y_t",ylab="Density",main="PDF of Logit-Normal distribution")
lines(seq(0.01,0.99,0.01),dlogitnorm(seq(0.01,0.99,0.01),-1,1), col="blue")
lines(seq(0.01,0.99,0.01),dlogitnorm(seq(0.01,0.99,0.01),1,1), col="purple")
lines(seq(0.01,0.99,0.01),dlogitnorm(seq(0.01,0.99,0.01),0,3), col="red")
legend("topright",legend=c("logitN(0,1)","logitN(-1,1)","logitN(1,1)","logitN(0,3)"),
       lwd=1, col=c("black","blue","purple","red"))

## ----pdfBeta, echo=FALSE------------------------------------------------------
plot(seq(0.01,0.99,0.01),dbeta(seq(0.01,0.99,0.01),1,1),type="l",ylim=c(0,4),
     xlab="y_t",ylab="Density",main="PDF of Beta distribution")
lines(seq(0.01,0.99,0.01),dbeta(seq(0.01,0.99,0.01),0.1,1), col="blue")
lines(seq(0.01,0.99,0.01),dbeta(seq(0.01,0.99,0.01),1,0.1), col="purple")
lines(seq(0.01,0.99,0.01),dbeta(seq(0.01,0.99,0.01),2,2), col="red")
legend("topright",legend=c("Beta(1,1)","Beta(0.1,1)","Beta(1,0.1)","Beta(2,2)"),
       lwd=1, col=c("black","blue","purple","red"))

## ----dataRound----------------------------------------------------------------
xreg[,1] <- round(abs(xreg[,1]))
inSample <- xreg[1:80,]
outSample <- xreg[-c(1:80),]

## ----negBinRegression---------------------------------------------------------
ourModel <- alm(y~x1+x2, data=inSample, distribution="dnbinom")
summary(ourModel)

## ----negBinRegressionWithSize-------------------------------------------------
ourModel <- alm(y~x1+x2, data=inSample, distribution="dnbinom", size=30)
summary(ourModel)

## ----cdfLogis, echo=FALSE-----------------------------------------------------
plot(seq(-5,5,0.01),plogis(seq(-5,5,0.01),0,1),type="l",
     xlab="y_t",ylab="Density",main="CDF of Logistic distribution")
lines(seq(-5,5,0.01),plogis(seq(-5,5,0.01),-1,1), col="blue")
lines(seq(-5,5,0.01),plogis(seq(-5,5,0.01),1,1), col="purple")
lines(seq(-5,5,0.01),plogis(seq(-5,5,0.01),2,2), col="red")
legend("bottomright",legend=c("Logit(0,1)","Logit(-1,1)","Logit(1,1)","Logit(2,2)"),
       lwd=1, col=c("black","blue","purple","red"))

## ----cdfNorm, echo=FALSE------------------------------------------------------
plot(seq(-5,5,0.01),pnorm(seq(-5,5,0.01),0,1),type="l",
     xlab="y_t",ylab="Density",main="CDF of Normal distribution")
lines(seq(-5,5,0.01),pnorm(seq(-5,5,0.01),-1,1), col="blue")
lines(seq(-5,5,0.01),pnorm(seq(-5,5,0.01),1,1), col="purple")
lines(seq(-5,5,0.01),pnorm(seq(-5,5,0.01),2,2), col="red")
legend("bottomright",legend=c("N(0,1)","N(-1,1)","N(1,1)","N(2,2)"),
       lwd=1, col=c("black","blue","purple","red"))

## ----mixtureExampleData-------------------------------------------------------
xreg[,1] <- round(exp(xreg[,1]-400) / (1 + exp(xreg[,1]-400)),0) * xreg[,1]
# Sometimes the generated data contains huge values
xreg[is.nan(xreg[,1]),1] <- 0;
inSample <- xreg[1:80,]
outSample <- xreg[-c(1:80),]

## ----mixtureExampleOccurrence-------------------------------------------------
modelOccurrence <- alm(y~x1+x2+Noise, inSample, distribution="plogis")

## ----mixtureExampleFinal------------------------------------------------------
modelMixture <- alm(y~x1+x2+Noise, inSample, distribution="dlnorm", occurrence=modelOccurrence)

## ----mixtureSummary-----------------------------------------------------------
summary(modelMixture)
summary(modelMixture$occurrence)

## ----mixtureDiagnostics-------------------------------------------------------
par(mfcol=c(3,3))
plot(modelMixture, c(1:9))

## ----mixturePredict-----------------------------------------------------------
predict(modelMixture,outSample,interval="p",level=c(0.8,0.9,0.95))

## ----mixtureExampleFinalAR----------------------------------------------------
modelMixtureAR <- alm(y~x1+x2+Noise, inSample, distribution="dlnorm", occurrence=modelOccurrence, ar=1)
summary(modelMixtureAR)
plot(predict(modelMixtureAR,outSample,interval="p",side="u"))

## ----mixtureExampleFinalARForecast--------------------------------------------
plot(forecast(modelMixtureAR, h=10, interval="p",side="u"))

## -----------------------------------------------------------------------------
lossFunction <- function(actual, fitted, B, xreg){
  return(mean(abs(actual-fitted)^3));
}
modelLossCustom <- alm(y~x1+x2+Noise, inSample, distribution="dnorm", loss=lossFunction)
summary(modelLossCustom)

## -----------------------------------------------------------------------------
summary(modelLossCustom, bootstrap=TRUE, nsim=100)

