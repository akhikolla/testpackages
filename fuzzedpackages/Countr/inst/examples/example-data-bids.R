##library(Countr)
library(Ecdat)
library(pscl)
library(MASS)

## load the bids data Bids. This data is described in Cameron(2013) Section 5.2.5
data(Bids)

## Define the model described in Cameron(2013) Section 5.2.5
form <- numbids ~ leglrest + rearest + finrest + whtknght + bidprem + insthold +
    size + I(size^2) + regulatn

## fit the Poisson model
pois <- glm(form, family=poisson(), data = Bids)

## run coef(pois) to see coefficient from Table 5.3 Cameron(2013).
## Note that standard errors are computed differently.

## fit negative binomial model: model does not converge
nb <- glm.nb(form, data = Bids)

## fit the renewal models
wei <- renewalCount(formula = form, data = Bids, dist = "weibull",
                    computeHessian = FALSE, weiMethod = "conv_dePril",
                    control = renewal.control(trace = 0,
                                              method = c("nlminb", 
                                                         "Nelder-Mead","BFGS"))
                    )

gam <- renewalCount(formula = form, data = Bids, dist = "gamma",
                    computeHessian = FALSE,
                    control = renewal.control(trace = 0,
                                              method = "nlminb")
                    )

gengam <- renewalCount(formula = form, data = Bids, dist = "gengamma",
                       computeHessian = FALSE, 
                       control = renewal.control(trace = 0,
                                                 method = "nlminb")
                       )

AIC <- c(poiss = AIC(pois),
         nb = AIC(nb),
         wei = AIC(wei),
         gam = AIC(gam),
         gengam = AIC(gengam)
         )

print(AIC)

## breaks defined at 5 as done in Table 5.6 Cameron (2013)
breaks_ <- 0:5
pears <- compareToGLM(poisson_model = pois,
                      breaks = breaks_,
                      nbinom_model = nb,
                      weibull = wei,
                      gamma = gam,
                      gengamma = gengam)

pearson <- colSums(dplyr::select(pears, contains("_pearson")))

## goofness of fit test
test_wei <- chiSq_gof(wei, breaks_)
test_gam <- chiSq_gof(gam, breaks_)
test_gengam <- chiSq_gof(gengam, breaks_)

## test_poiss$Chisq is exactly the same as the value reported in
## Caemron(2013)[2nd edition] page 195
test_poiss <- chiSq_gof(pois, breaks_)  
test_nb <- chiSq_gof(nb, breaks_)
