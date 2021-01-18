## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align="center",
  fig.height=4,
  fig.width=6,
  fig.path='Figs/',
  fig.show='hold',
  warning=FALSE,
  message=FALSE
)

## ----load_libraries, message=FALSE, warning=FALSE-----------------------------
require(greybox)
require(smooth)
require(Mcomp)

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "MMM", lags=c(1,12), distribution="dnorm")
summary(testModel)
plot(forecast(testModel,h=18,interval="parametric"))

## -----------------------------------------------------------------------------
testModel

## -----------------------------------------------------------------------------
plot(forecast(testModel,h=18,interval="simulated"))

## -----------------------------------------------------------------------------
par(mfcol=c(3,4))
plot(testModel,which=c(1:11))
par(mfcol=c(1,1))
plot(testModel,which=12)

## -----------------------------------------------------------------------------
lossFunction <- function(actual, fitted, B){
  return(sum(abs(actual-fitted)^3))
}
testModel <- adam(M3[[1234]], "AAN", silent=FALSE, loss=lossFunction)
testModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[1234]], "MMN", silent=FALSE, distribution="dgnorm", beta=3)

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "ZXZ", lags=c(1,12), silent=FALSE)
testModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "CXC", lags=c(1,12))
testForecast <- forecast(testModel,h=18,interval="semiparametric", level=c(0.9,0.95))
testForecast
plot(testForecast)

## -----------------------------------------------------------------------------
forecast(testModel,h=18,interval="semiparametric", level=c(0.9,0.95,0.99), side="upper")

## -----------------------------------------------------------------------------
testModel <- adam(forecast::taylor, "MMdM", lags=c(1,48,336), silent=FALSE, h=336, holdout=TRUE)
testModel

## -----------------------------------------------------------------------------
testModel <- adam(forecast::taylor, "MMdM", lags=c(1,48,336), silent=FALSE, h=336, holdout=TRUE,
                  maxeval=10000)
testModel

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  testModel$B

## -----------------------------------------------------------------------------
testModel <- adam(forecast::taylor, "MMdM", lags=c(1,48,336), silent=FALSE, h=336, holdout=TRUE,
                  B=testModel$B)
testModel

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  testModel <- adam(forecast::taylor, "MMdM", lags=c(1,48,336), silent=FALSE, h=336, holdout=TRUE,
#                    initial="b")

## -----------------------------------------------------------------------------
testModel <- adam(forecast::taylor, "MMdM", lags=c(1,48,336), silent=TRUE, h=336, holdout=TRUE,
                  initial=list(level=30000, trend=1), persistence=list(beta=0.1))
testModel

## -----------------------------------------------------------------------------
testModel <- adam(rpois(120,0.5), "MNN", silent=FALSE, h=12, holdout=TRUE,
                  occurrence="odds-ratio")
testModel

## -----------------------------------------------------------------------------
adamModel <- adam(M3[[2568]], "CCC")
esModel <- es(M3[[2568]], "CCC")
"adam:"
adamModel
"es():"
esModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[1234]], "NNN", silent=FALSE, orders=c(0,2,2))
testModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "NNN", silent=FALSE, lags=c(1,12),
                  orders=list(ar=c(1,1),i=c(1,1),ma=c(2,2)), distribution="dlnorm")
testModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "ANN", silent=FALSE, lags=c(1,12), persistence=0,
                  orders=list(ar=c(1,1),i=c(1,1),ma=c(2,2)), distribution="dnorm")
testModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "NNN", silent=FALSE, lags=c(1,12),
                  orders=list(ar=c(1,1),i=c(1,1),ma=c(2,2)), distribution="dnorm",
                  arma=list(ar=c(0.1,0.1), ma=c(-0.96, 0.03, -0.12, 0.03)))
testModel

## -----------------------------------------------------------------------------
testModel <- adam(M3[[2568]], "NNN", silent=FALSE, lags=c(1,12),
                  orders=list(ar=c(1,1),i=c(1,1),ma=c(2,0)), distribution="dnorm",
                  initial=list(arima=M3[[2568]]$x[1:24]))
testModel

## -----------------------------------------------------------------------------
BJData <- cbind(BJsales,BJsales.lead)
testModel <- adam(BJData, "AAN", h=18, silent=FALSE)

## -----------------------------------------------------------------------------
BJData <- cbind(as.data.frame(BJsales),as.data.frame(xregExpander(BJsales.lead,c(-7:7))))
colnames(BJData)[1] <- "y"
testModel <- adam(BJData, "ANN", h=18, silent=FALSE, holdout=TRUE, formula=y~xLag1+xLag2+xLag3)
testModel

## -----------------------------------------------------------------------------
testModel <- adam(BJData, "ANN", h=18, silent=FALSE, holdout=TRUE, regressors="select")

## -----------------------------------------------------------------------------
testModel <- adam(BJData, "NNN", h=18, silent=FALSE, holdout=TRUE, regressors="select", orders=c(0,1,1))

## -----------------------------------------------------------------------------
BJData <- BJData[,c("y",names(testModel$initial$xreg))];
testModel <- adam(BJData, "NNN", h=18, silent=TRUE, holdout=TRUE, orders=c(0,1,1),
                  initial=testModel$initial, arma=testModel$arma)
testModel
names(testModel$initial)[1] <- names(testModel$initial)[[1]] <- "level"
testModel2 <- adam(BJData, "ANN", h=18, silent=TRUE, holdout=TRUE,
                   initial=testModel$initial, persistence=testModel$arma$ma+1)
testModel2

## -----------------------------------------------------------------------------
testModel <- adam(BJData, "ANN", h=18, silent=FALSE, holdout=TRUE, regressors="adapt")
testModel$persistence

## -----------------------------------------------------------------------------
testModel <- adam(BJData, "AAN", h=18, silent=FALSE, holdout=TRUE, orders=c(1,0,1))
summary(testModel)

## -----------------------------------------------------------------------------
testModel <- auto.adam(M3[[1234]], "XXX", silent=FALSE,
                       distribution=c("dnorm","dlaplace","ds"))
testModel

## ----eval=FALSE, echo=TRUE----------------------------------------------------
#  testModel <- auto.adam(M3[[1234]], "ZZZ", silent=FALSE, parallel=TRUE)

## -----------------------------------------------------------------------------
testModel <- auto.adam(M3[[1234]], "AAN", orders=list(ar=2,i=2,ma=2), silent=TRUE,
                       distribution=c("dnorm","dlaplace","ds","dgnorm"))
testModel

## -----------------------------------------------------------------------------
testModel <- auto.adam(M3[[1234]], "XXN", orders=list(ar=2,i=2,ma=2,select=TRUE),
                       distribution="default", silent=FALSE)
testModel

## -----------------------------------------------------------------------------
testModel <- auto.adam(Mcomp::M3[[2568]], "PPP", silent=FALSE, outliers="use",
                       distribution="default")
testModel

