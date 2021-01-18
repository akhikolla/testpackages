library(MASS)
data("quine") ## data suggested by referee 2 (over dispersion)

quine_form <- as.formula(Days ~ Eth + Sex + Age + Lrn)

## Poisson and negative binomial
pois <- glm(quine_form, family = poisson(), data = quine)
nb <- glm.nb(quine_form, data = quine)

## various renewal models
wei <- renewalCount(formula = quine_form, data = quine, dist = "weibull",
                          computeHessian = FALSE, weiMethod = "conv_dePril",
                          control = renewal.control(trace = 0,
                              method = c("nlminb", 
                                  "Nelder-Mead","BFGS")),
                          convPars = list(convMethod = "dePril")
                          )

gam <- renewalCount(formula = quine_form, data = quine, dist = "gamma",
                    computeHessian = FALSE, weiMethod = "conv_dePril",
                    control = renewal.control(trace = 0,
                                              method = "nlminb"),
                        convPars = list(convMethod = "dePril")
                    )

gengam <- renewalCount(formula = quine_form, data = quine, dist = "gengamma",
                       computeHessian = FALSE, weiMethod = "conv_dePril",
                       control = renewal.control(trace = 0,
                                                 method = "nlminb"),
                           convPars = list(convMethod = "dePril")
                       )

AIC <- c(poiss = AIC(pois),
         nb = AIC(nb),
         wei = AIC(wei),
         gam = AIC(gam),
         gengam = AIC(gengam)
         )

print(AIC)

breaks_ <- c(0, 1, 3, 5:7, 9, 12, 15, 17, 23, 27, 32)
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
test_poiss <- chiSq_gof(pois, breaks_)
test_nb <- chiSq_gof(nb, breaks_)
