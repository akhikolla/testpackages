library(Countr)
library("MASS") # for glm.nb()

library(dplyr) 
library(xtable)

data(quine, package = "MASS")

breaks_ <- c(0, 1, 3, 5:7, 9, 12, 15, 17, 23, 27, 32)
freqtable <- 
    count_table(count = quine$Days, breaks = breaks_, formatChar = TRUE)

 print(xtable(freqtable[ , 1:7]), floating = FALSE, only.contents = TRUE)
 cat("\n\\\\[5pt]\n")
 print(xtable(freqtable[ , -(1:7)]), floating = FALSE, only.contents = TRUE)

quine_form <- as.formula(Days ~ Eth + Sex + Age + Lrn)
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

library(lmtest)
pois_nb <- lrtest(pois, nb)
pois_wei <- suppressWarnings(lrtest(pois, wei))
pois_gam <- suppressWarnings(lrtest(pois, gam))
pois_gengam <- suppressWarnings(lrtest(pois, gengam))
pois_res <- data.frame("Alternative model" = c("negative-binomial", "weibull",
                                 "gamma", "generalised-gamma"),
                       Chisq = c(pois_nb$Chisq[2], pois_wei$Chisq[2],
                                 pois_gam$Chisq[2], pois_gengam$Chisq[2]),
                       Df = c(1, 1, 1, 2),
                       Critical_value = c(rep(qchisq(0.95, 1), 3),
                                           qchisq(0.95, 2)),
                       stringsAsFactors = FALSE                    
                       )                  
                         
print(xtable(pois_res, caption = "LR results against Poisson model. Each row compares an alternative model vs the Poisson model. All alternatives are preferable to Poisson.", 
             label = "tab:lr_pois"))

gengam_wei <- lrtest(wei, gengam)
gengam_gam <- lrtest(gam, gengam)

gengam_res <- data.frame(Model = c("weibull", "gamma"),
                         Chisq = c(gengam_wei$Chisq[2], gengam_gam$Chisq[2]),
                         Df = 1,
                         Critical_value = rep(qchisq(0.95, 1), 2),
                         stringsAsFactors = FALSE                    
                         )                  
                         
print(xtable(gengam_res, caption = "LR results against generalised-gamma model", 
             label = "tab:lr_gengam"))

ic <- data.frame(Model = c("weibull", "gamma", "negative-binomial"),
                 AIC = c(AIC(wei), AIC(gam), AIC(nb)),
                 BIC = c(BIC(wei), BIC(gam), BIC(nb)),
                 stringsAsFactors = FALSE                    
                 )                  
                         
print(xtable(ic, caption = "Information criteria results", 
             label = "tab:ic_models"))

  gof <- chiSq_gof(gam, breaks = breaks_)
  print(gof)

save.image()
