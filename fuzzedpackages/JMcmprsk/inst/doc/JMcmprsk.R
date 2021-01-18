## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  require(JMcmprsk)
#  set.seed(123)
#  data(lung)
#  yread <- lung[, c(1,2:11)]
#  cread <- unique(lung[, c(1, 12, 13, 6:10)])

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  jmcfit <- jmc(long_data = yread, surv_data = cread, out = "FVC",
#             FE = c("time", "FVC0", "FIB0", "CYC", "FVC0.CYC",
#                    "FIB0.CYC", "time.CYC"),
#             RE = "linear", ID = "ID",cate = NULL, intcpt = 0,
#              quad.points = 20, quiet = FALSE)
#  jmcfit

## ---- message=FALSE, warning=FALSE, results='asis'----------------------------
require(JMcmprsk)
data(lung)
lungY <- lung[, c(2:11)]
lungC <- unique(lung[, c(1, 12, 13, 6:10)])
lungC <- lungC[, -1]
lungM <- data.frame(table(lung$ID))
lungM <- as.data.frame(lungM[, 2])

## ---- eval=FALSE, message=FALSE, warning=FALSE--------------------------------
#  res1=jmc_0(p=8, lungY, lungC, lungM, point=20, do.trace = FALSE, type_file = FALSE)
#  res1

## ---- echo=TRUE, message=FALSE, warning=FALSE, results='base'-----------------
#because of the long running time, we save the jmo and jmc results within the package
fitfile=system.file("extdata", "runfit.RData", package = "JMcmprsk")
load(fitfile)
jmcfit

## ---- echo=TRUE, message=FALSE, warning=FALSE---------------------------------
coef(jmcfit, coeff = "beta")
coef(jmcfit, coeff = "gamma")

## ---- echo=TRUE, message=FALSE, warning=FALSE---------------------------------
## Linear hypothesis of testing all coefficients of beta's / gamma's equal 0
linearTest(jmcfit,coeff="beta")
linearTest(jmcfit,coeff="gamma")

## ---- echo=TRUE, message=FALSE, warning=FALSE---------------------------------
## Extract the standard errors for the longitudinal portion
summary(jmcfit, coeff = "longitudinal")
## Extract the standard errors for the survival portion
summary(jmcfit, coeff = "survival")

## ---- message=FALSE, warning=FALSE, results='asis'----------------------------
set.seed(123)
require(JMcmprsk)
data(ninds)
yread <- ninds[, c(1, 2:14)]
cread <- ninds[, c(1, 15, 16, 6, 10:14)]
cread <- unique(cread)

## ---- eval=FALSE, message=FALSE, warning=FALSE, results='asis'----------------
#  jmofit <- jmo(yread, cread, out = "Y",
#              FE = c("group", "time3", "time6", "time12", "mrkprior",
#                     "smlves", "lvORcs", "smlves.group", "lvORcs.group"),
#              cate = NULL,RE = "intercept", NP = c("smlves", "lvORcs"),
#              ID = "ID",intcpt = 1, quad.points = 20,
#              max.iter = 1000, quiet = FALSE, do.trace = FALSE)
#  jmofit

## ---- message=FALSE, warning=FALSE, results='asis'----------------------------
require(JMcmprsk)
data(ninds)
yread <- ninds[, c(2:14)]
mread <- as.data.frame(table(ninds$ID))
mread <- as.data.frame(mread[, 2])
cread <- ninds[, c(1, 15, 16, 6, 10:14)]
cread <- unique(cread)
cread <- cread[, -1]

## ---- eval=FALSE, message=FALSE, warning=FALSE, results='asis'----------------
#  jmofit <- jmo_0(p=9,s=2, yread,cread,mread,point=6,do.trace = FALSE, type_file = FALSE)
#  jmofit

## ---- echo=TRUE, message=FALSE, warning=FALSE, results='base'-----------------
#because of the long running time, we save the jmo and jmc results within the package
fitfile=system.file("extdata", "runfit.RData", package = "JMcmprsk")
load(fitfile)
jmofit

## ---- message=FALSE, warning=FALSE, results='base'----------------------------
##extract the parameter estimates of longitudinal proportional odds fixed effects
beta <- coef(jmofit, coeff = "beta")
beta
##extract the parameter estimates of longitudinal non-proportional odds fixed effects
alpha <- coef(jmofit, coeff = "alpha")
alpha
##extract the parameter estimates of longitudinal logit-specific intercept
theta <- coef(jmofit, coeff = "theta")
theta
##extract the parameter estimates of survival fixed effects
gamma <- coef(jmofit, coeff = "gamma")
gamma

## ---- message=FALSE, warning=FALSE, results='base'----------------------------
## Linear hypothesis of testing all coefficients of beta's / gamma's / alpha's equal 0
linearTest(jmofit,coeff="beta")
linearTest(jmofit,coeff="gamma")
linearTest(jmofit,coeff="alpha")

## ---- message=FALSE, warning=FALSE, results='base'----------------------------
## Extract the standard errors of both longitudinal and survival portions
summary(jmofit, coeff = "longitudinal")
summary(jmofit, coeff = "survival")

