## ---- eval = F, echo = T------------------------------------------------------
#  library("devtools")
#  devtools::install_github("rvpanaro/spsurv")

## ---- eval=F, echo = T--------------------------------------------------------
#  spbp.default <-
#    spbp(formula, degree, data,
#              approach = c("mle", "bayes"),
#              model = c("ph", "po", "aft"),
#              priors = list(beta = c("normal(0,4)"),
#                           gamma = "lognormal(0,10)"),
#             scale = TRUE,
#             ...)

## ---- message = F-------------------------------------------------------------
library("KMsurv")
data("larynx")

library(spsurv)
fit <- spsurv::spbp(Surv(time, delta)~ age + factor(stage),
                    approach = "mle",  data = larynx)
summary(fit)                    

## ---- message=F---------------------------------------------------------------
library("KMsurv")
data("larynx")
 
library(spsurv)
fit <- spsurv::spbp(Surv(time, delta)~age + factor(stage),
                     approach = "mle",  data = larynx)


## -----------------------------------------------------------------------------
fit$coefficients
head(model.matrix(fit))
diag(fit$var)

fit$loglik

## -----------------------------------------------------------------------------
fit <- spsurv::spbp(Surv(time, delta)~age + factor(stage),
                     approach = "bayes",  data = larynx,
                     iter = 2000, chains = 1, warmup = 1000, cores = 1)

## -----------------------------------------------------------------------------
summary(fit)

## -----------------------------------------------------------------------------
rstan::traceplot(fit$stanfit, pars = c("beta", "gamma"))
rstan::stan_dens(fit$stanfit, pars = c("beta", "gamma"))

## ---- eval=F------------------------------------------------------------------
#  shinystan::launch_shinystan(fit$stanfit)

## -----------------------------------------------------------------------------
## CoxPH model
fitcoxph <- survival::coxph(Surv(time , delta)~age + factor(stage),
data = larynx)

## Determine the groups of patients
newdata <-  data.frame(age =c(77,77,77,77), stage = c(1,2,3,4))

## survfit Breslow estimator
breslowsurv <- survival::survfit(fitcoxph, newdata = newdata)

## spbp point-wise estimate
spbpsurv <- spsurv::survivor(fit, newdata = newdata)

plot(breslowsurv, bty = "n", lwd = 3, main = "77 years old patient survival per Stage")

points(spbpsurv$time, spbpsurv$survival1, col = 1, pch = 23)
points(spbpsurv$time, spbpsurv$survival2, col = 2, pch = 23)
points(spbpsurv$time, spbpsurv$survival3, col = 3, pch = 23)
points(spbpsurv$time, spbpsurv$survival4, col = 4, pch = 23)

legend("topright", c("Stage I", "Stage II", "Stage III", "Stage IV"), pch = 23, bty = "n", col = 1:4)

