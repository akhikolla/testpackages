## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  fig.height=4,
  fig.width=6
)

## ---- echo=FALSE, message=FALSE-----------------------------------------------
library(greybox)

## ----BJxreg1------------------------------------------------------------------
BJxreg <- xregExpander(BJsales.lead,lags=c(-5,-10))

## ----BJxreg2------------------------------------------------------------------
BJxreg <- xregExpander(BJsales.lead,lags=c(7,-5,-10))

## ----BJxreg3------------------------------------------------------------------
BJxreg <- xregExpander(BJsales.lead,lags=c(-10:10))

## ----BJData-------------------------------------------------------------------
BJxreg <- as.data.frame(xregExpander(BJsales.lead,lags=c(-10:10)))
BJxreg <- cbind(as.matrix(BJsales),BJxreg)
colnames(BJxreg)[1] <- "y"
ourModel <- stepwise(BJxreg)

## ----BJStepwise---------------------------------------------------------------
ourModel <- stepwise(BJxreg)

## ----BJStepwiseResult---------------------------------------------------------
ourModel

## ----texregExample, results = 'asis'------------------------------------------
texreg::htmlreg(ourModel)

## ----BJcombine1---------------------------------------------------------------
ourModel <- lmCombine(BJxreg[,-c(3:7,18:22)],bruteforce=TRUE)
summary(ourModel)

## ----BJcombine2---------------------------------------------------------------
ourModel <- lmCombine(BJxreg,bruteforce=FALSE)
summary(ourModel)

## ----BJcombine3---------------------------------------------------------------
BJInsample <- BJxreg[1:130,];
BJHoldout <- BJxreg[-(1:130),];
ourModel <- lmCombine(BJInsample,bruteforce=FALSE)

## ----BJcombinePlot------------------------------------------------------------
summary(ourModel)
plot(ourModel)

## ----BJcombineForecast--------------------------------------------------------
ourForecast <- predict(ourModel,BJHoldout)
plot(ourForecast)

