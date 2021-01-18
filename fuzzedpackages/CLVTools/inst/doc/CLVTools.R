## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  #fig.path = "figures/WALKTHROUGH-",
  out.width = "100%"
)

## ----install-package-CRAN, eval = FALSE---------------------------------------
#  install.packages("CLVTools")

## ----install-package-GITHUB, eval = FALSE-------------------------------------
#  install.packages("devtools")
#  devtools::install_github("bachmannpatrick/CLVTools", ref = "master")

## ----load-library-------------------------------------------------------------
library("CLVTools")

## ----load-data----------------------------------------------------------------
data("apparelTrans")
apparelTrans

## ----load-CreateObj-----------------------------------------------------------
clv.apparel <- clvdata(apparelTrans,  
                       date.format="ymd", 
                       time.unit = "week",
                       estimation.split = 40,
                       name.id = "Id",
                       name.date = "Date",
                       name.price = "Price")

## ----print-CLVObject----------------------------------------------------------
clv.apparel

## ----summary-CLVObject--------------------------------------------------------
summary(clv.apparel)

## ----estimate-model-----------------------------------------------------------
est.pnbd <- pnbd(clv.data = clv.apparel)
est.pnbd

## ----estimate-model2, eval=FALSE----------------------------------------------
#  est.pnbd <- pnbd(clv.data = clv.apparel,
#                       start.params.model = c(r=1, alpha = 2, s = 1, beta = 2),
#                       optimx.args = list(control=list(trace=5),
#                                         method="Nelder-Mead"
#                                         ))

## ----param-summary------------------------------------------------------------
#Full detailed summary of the parameter estimates
summary(est.pnbd)

#Extract the coefficients only
coef(est.pnbd)
#Alternative: oefficients(est.pnbd.obj)


## ----coef-model---------------------------------------------------------------
#Extract the coefficients only
coef(est.pnbd)
#Alternative: oefficients(est.pnbd.obj)

#Extract the confidence intervals
confint(est.pnbd)


## ----ll-model-----------------------------------------------------------------
# LogLikelihood at maximum
logLik(est.pnbd)

# Variance-Covariance Matrix at maximum
vcov(est.pnbd)


## ----estimate-ggomnbd, eval=FALSE---------------------------------------------
#  est.ggomnbd <- ggomnbd(clv.data = clv.apparel,
#                       start.params.model = c(r=0.7, alpha=5, b=0.005,  s=0.02, beta=0.001),
#                       optimx.args = list(control=list(trace=5),
#                                         method="Nelder-Mead"))

## ----predict-model------------------------------------------------------------
results <- predict(est.pnbd)
print(results)


## ----predict-model2, eval = FALSE---------------------------------------------
#  predict(est.pnbd, prediction.end = 30)

## ----plot-model3, eval = FALSE------------------------------------------------
#  predict(est.pnbd, prediction.end = "2006-05-08")

## ----plot-model, fig.height=4.40, fig.width=9---------------------------------
plot(est.pnbd)


## ----plot-model2, eval = FALSE------------------------------------------------
#  plot(est.pnbd, prediction.end = 30, cumulative = TRUE)

## ----predict-model3, eval = FALSE---------------------------------------------
#  plot(est.pnbd, prediction.end = "2006-05-08", cumulative = TRUE)

## ----Cov-staticData-----------------------------------------------------------
data("apparelStaticCov")
apparelStaticCov

## ----Cov-dynData--------------------------------------------------------------
data("apparelDynCov")
apparelDynCov

## ----Cov-setStatic------------------------------------------------------------
clv.static<- SetStaticCovariates(clv.data = clv.apparel, 
                                      data.cov.life = apparelStaticCov, 
                                      data.cov.trans = apparelStaticCov,
                                      names.cov.life = c("Gender", "Channel"), 
                                      names.cov.trans =c("Gender", "Channel"), 
                                      name.id = "Id")

## ----Cov-setDynamic, eval=FALSE, message=FALSE, warning=TRUE------------------
#  clv.dyn <- SetDynamicCovariates(clv.data = clv.apparel,
#                                       data.cov.life = apparelDynCov,
#                                       data.cov.trans = apparelDynCov,
#                                       names.cov.life = c("Marketing", "Gender", "Channel"),
#                                       names.cov.trans = c("Marketing", "Gender", "Channel"),
#                                       name.id = "Id",
#                                       name.date = "Cov.Date")

## ----Static-cov-estimate, message=TRUE, warning=FALSE-------------------------
est.pnbd.static <- pnbd(clv.static, 
                         start.params.model = c(r=1, alpha = 2, s = 1, beta = 2),
                         start.params.life = c(Gender=0.6, Channel=0.4),
                         start.params.trans = c(Gender=0.6, Channel=0.4))

## ----Dyn-cov-estimate, eval=FALSE---------------------------------------------
#  est.pnbd.dyn <- pnbd(clv.dyn,
#                       start.params.model = c(r=1, alpha = 2, s = 1, beta = 2),
#                       start.params.life = c(Marketing=0.5, Gender=0.6, Channel=0.4),
#                       start.params.trans = c(Marketing=0.5, Gender=0.6, Channel=0.4))

## ----Cov-summary--------------------------------------------------------------
summary(est.pnbd.static)

## ----Cov-cor, eval=FALSE------------------------------------------------------
#  est.pnbd.cor <- pnbd(clv.apparel,
#                       use.cor= TRUE)
#  summary(est.pnbd.cor)

## ----reg-advOptions-----------------------------------------------------------
est.pnbd.reg <- pnbd(clv.static, 
                         start.params.model = c(r=1, alpha = 2, s = 1, beta = 2),
                         reg.lambdas = c(trans=100, life=100))
summary(est.pnbd.reg)

## ----constr-advOptions--------------------------------------------------------
est.pnbd.constr <- pnbd(clv.static, 
                         start.params.model = c(r=1, alpha = 2, s = 1, beta = 2),
                         start.params.constr = c(Gender=0.6),
                         names.cov.constr=c("Gender"))
summary(est.pnbd.constr)

## ----doFuture, eval=FALSE-----------------------------------------------------
#    # disable multithreading for data.table (to avoid nested parallelism)
#    setDTthreads(1)
#  
#    library("doFuture")
#    registerDoFuture()
#    plan("multisession")

## ----spending-load-data and initialize----------------------------------------
data("apparelTrans")
apparelTrans

clv.apparel <- clvdata(apparelTrans,  
                       date.format="ymd", 
                       time.unit = "week",
                       estimation.split = 40,
                       name.id = "Id",
                       name.date = "Date",
                       name.price = "Price")

## ----spending-estimate-model--------------------------------------------------
est.gg<- gg(clv.data = clv.apparel)
est.gg

## ----spendsing-estimate-model2, eval=TRUE-------------------------------------
est.gg<- gg(clv.data = clv.apparel, remove.first.transaction=FALSE)
est.gg

## ----spending-predict-model---------------------------------------------------
results.spending <- predict(est.gg)
print(results.spending)

## ----spending-plot-model4, fig.height=4.40, fig.width=9-----------------------
plot(est.gg)


