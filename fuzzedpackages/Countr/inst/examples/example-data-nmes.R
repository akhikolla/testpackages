library(Countr)
library(AER)
data("NMES1988")

## this is the model fitted in Cameron (2013) Regression Analysis of Count Data
## Chapter 6.3
## The code used by the authors is available in
## References/Regression Analysis of Count Data (2nd Edition)/racd06p1.R
## The data loaded from pkg AER is slightly different (processed) but the result
## obtained for Poisson /negative binomial is identical
nmes_form <- as.formula(visits ~ health +  chronic  + adl +
                            region + age + afam +  gender  + married  +
                            school + income + employed +  insurance + medicaid)

## Poisson model
poiss <- glm(nmes_form, family = poisson(), data = NMES1988)

## negative binomial from MASS pkg
library(MASS)
nb <- glm.nb(nmes_form, data = NMES1988)

## renewal models
wei <- renewalCount(formula = nmes_form, data = NMES1988, dist = "weibull",
                    computeHessian = FALSE, weiMethod = "conv_dePril",
                    control = renewal.control(trace = 0,
                                              method = c("nlminb", "BFGS")
                                              ),
                    convPars = list(convMethod = "dePril")
                    )

gengam <- renewalCount(formula = nmes_form, data = NMES1988,
                       dist = "gengamma",
                       computeHessian = FALSE,
                       control = renewal.control(trace = 0,
                                                 method = "BFGS"
                                                 ## method = c("nlminb", "spg",
                                                 ##            "Nelder-Mead","BFGS")
                                                 ),
                       convPars = list(convMethod = "dePril")
                       )


## compare models
loglik <- sapply(list(poiss, nb, wei, gengam), logLik)
aic <- sapply(list(poiss, nb, wei, gengam), AIC)
bic <- sapply(list(poiss, nb, wei, gengam), BIC)

## pearson stat
pears <- compareToGLM(poisson_model = poiss,
                      breaks = 0:13,
                      nbinom_model = nb,
                      weibull = wei,
                      gengamma = gengam) 

pearson <- colSums(dplyr::select(pears, contains("_pearson")))
