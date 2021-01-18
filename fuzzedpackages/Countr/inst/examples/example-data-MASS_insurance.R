library(Countr)
library(MASS)
data(Insurance)


## Poisson model
pois <- glm(Claims ~ District + Group + Age + log(Holders),
            data = Insurance, family = poisson)

## negative binomial
nb <- glm.nb(Claims ~ District + Group + Age + log(Holders),
             data = Insurance)

mat <- model.matrix(pois)
betas <- coef(pois)
lambda <- exp(mat %*% betas)

## renewal models
wei <- renewalCount(Claims ~ District + Group + Age + log(Holders),
                    data = Insurance, dist = "weibull",
                    weiMethod = "conv_dePril",
                    convPars = list(nsteps = 50, convMethod = "dePril"),
                    control = renewal.control(trace = 0, method = "nlminb"),
                    computeHessian = FALSE)


gam <- renewalCount(Claims ~ District + Group + Age + log(Holders),
                    data = Insurance, dist = "gamma",
                    convPars = list(nsteps = 50, convMethod = "dePril"),
                    control = renewal.control(trace = 0, method = "nlminb"),
                    computeHessian = FALSE)


gengam <- renewalCount(Claims ~ District + Group + Age + log(Holders),
                       data = Insurance, dist = "gengamma",
                       convPars = list(nsteps = 50, convMethod = "dePril"),
                       control = renewal.control(trace = 0, method = "nlminb"),
                       computeHessian = FALSE)


## compare models
loglik <- sapply(list(pois, nb, wei, gam, gengam), logLik)
aic <- sapply(list(pois, nb, wei, gam, gengam), AIC)
bic <- sapply(list(pois, nb, wei, gam, gengam), BIC)

## pearson stat
breaks_ <- c(0, 10, 20, 30, 70)
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
