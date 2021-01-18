# ----load options---------------------------------------------------------
options(digits = 3, show.signif.stars = FALSE)

# ----load packages--------------------------------------------------------
library("Countr")
library("lmtest")
library("dplyr")
library("xtable")
library("BB")  # needed for optimisation with spg algo only

# ----football data--------------------------------------------------------
data("football", package = "Countr")
table(football$awayTeamGoals)


# ----models for away goals------------------------------------------------
away_poiss <- glm(formula = awayTeamGoals ~ 1,
                  family = poisson,
                  data = football)
away_wei <- renewalCount(formula = awayTeamGoals ~ 1,
                         data = football,
                         dist = "weibull", computeHessian = FALSE,
                         control = renewal.control(trace = 0))

# -------------------------------------------------------------------------
breaks_ <- 0:5
pears <- compareToGLM(poisson_model = away_poiss, breaks = breaks_,
                      weibull = away_wei)


# ----frequency plot-------------------------------------------------------
frequency_plot(pears$Counts, pears$Actual,
               select(pears, contains("_predicted")),
               colours = c("gray", "blue", "green", recursive="black"))


# ----LR test--------------------------------------------------------------
lr <- lrtest(away_poiss, away_wei)
lr


# ----X2 test--------------------------------------------------------------
gof_wei <- chiSq_gof(away_wei, breaks = breaks_)
gof_pois <- chiSq_gof(away_poiss, breaks = breaks_)
rbind(Poisson = gof_pois, "Weibull-count" = gof_wei)


# ----fertility data-------------------------------------------------------
data("fertility", package = "Countr")


# ----head-fertility-------------------------------------------------------
dataCaption <- "First few rows of fertility data."
print(xtable(head(fertility), caption = dataCaption, label = "tbl:data"),
      rotate.colnames = TRUE)


# ----children-table-------------------------------------------------------
freqtable <- count_table(count = fertility$children, breaks = 0:9, formatChar = TRUE)
print(xtable(freqtable, caption = "Fertility data: Frequency distribution of column \\texttt{children}.",
             label = "tbl:freq"))


# ----fertility data processing--------------------------------------------
nam_fac <- sapply(fertility, function(x) !is.numeric(x))
fert_factor <- summary(fertility[ , nam_fac])
fert_num <- t(sapply(fertility[ , !nam_fac], summary)) # summary(fertility[ , !nam_fac])


# ----covariates-table-----------------------------------------------------
print(xtable(fert_factor,
             caption = "Summary of the factor variables",
             label = "tbl:frecfac"))


# ----covariates-table-num-------------------------------------------------
print(xtable(fert_num,
             caption = "Summary of the numeric explanatory variables",
             label = "tbl:frecnum"))


# ----main model-----------------------------------------------------------
regModel <- children ~ german + years_school + voc_train + university +
  religion + rural + year_birth + age_marriage


# ----link functions-------------------------------------------------------
link_weibull <- list(scale = "log", shape = "log")


# ----gamma model--fertility data------------------------------------------
gamModel <- renewalCount(formula = regModel, data = fertility,
                         dist = "gamma",
                         control = renewal.control(trace = 0))


# ----visualize parameters name: weibull-----------------------------------
getParNames("weibull")


# ----visualize parameters name: gamma-------------------------------------
renewalNames(regModel, data = fertility, dist = "gamma")


# ----Poisson initial values-----------------------------------------------
IV <- glm(regModel, family = poisson(), data = fertility)


# ----print initial values-------------------------------------------------
coef(IV)


# ----prepare initial values-----------------------------------------------
startW <- renewalCoef(IV, target = "scale")


# -------------------------------------------------------------------------
startW <- c(startW, "shape_" = log(1))
startW


# ----weibull model with initial values------------------------------------
weiModel <- renewalCount(formula = regModel, data = fertility,
                         dist = "weibull",
                         control = renewal.control(trace = 0, start = startW))


# ----change the optim algo to L-BFGS-B------------------------------------
weiModelA <- renewalCount(formula = regModel, data = fertility,
  dist = "weibull",
  control = renewal.control(trace = 0, method = "L-BFGS-B"))


# ----try few optim algo---------------------------------------------------
weiModel_many <- renewalCount(formula = regModel, data = fertility, dist = "weibull",
		 control = renewal.control(trace = 0, method = c("nlminb", "Nelder-Mead", "BFGS")))

# ----compare their performance--------------------------------------------
t(weiModel_many$optim)


# ----prepare a formula with CountrFormula---------------------------------
CountrFormula(y ~ x1 + x2 + x3, shape = ~x1)


# ----anc formula----------------------------------------------------------
anc <- list(sigma = regModel, Q = regModel)


# ----starting values with anc---------------------------------------------
startA <- renewalCoef(IV, target = "gengamma")
startA[c("Q_", "sigma_")] <- c(1, log(1))
startA


# ----fit the gen gamma model with anc-------------------------------------
gengamModel_ext0 <- renewalCount(formula = regModel, data = fertility,
                                 dist = "gengamma", anc = anc,
                                 control = renewal.control(start = startA, trace = 0),
                                 computeHessian = FALSE)


# ----alternative models for anc parameters--------------------------------
sigmaModel <- ~ german + university + religion + age_marriage
QModel     <- ~ german + religion + age_marriage


# -------------------------------------------------------------------------
anc <- list(sigma = sigmaModel, Q = QModel)


# -------------------------------------------------------------------------
regModelSQ <- Formula::as.Formula(regModel, sigmaModel, QModel)


# -------------------------------------------------------------------------
CountrFormula(regModel, sigma = sigmaModel, Q = QModel)


# ----initial values for new anc models------------------------------------
IV2 <- glm(update(sigmaModel, children ~ .),
           family = poisson(), data = fertility)
IV3 <- glm(update(QModel, children ~ .),
           family = poisson(), data = fertility)


# -------------------------------------------------------------------------
startGG <- c(renewalCoef(IV, target = "mu"),
             renewalCoef(IV2, target = "sigma"),
             renewalCoef(IV3, target = "Q"))
startGG


# ----fit the new models---------------------------------------------------
fm_gengamAnc <- renewalCount(formula = regModel, data = fertility,
                             dist = "gengamma", anc = anc,
                             control = renewal.control(start = startGG, trace = 0),
                             computeHessian = FALSE)


# -------------------------------------------------------------------------
fm_gengam <- renewalCount(formula = regModelSQ, data = fertility,
                          dist = "gengamma",
                          control = renewal.control(start = startGG, trace = 0),
                          computeHessian = FALSE)


# ----refit model with previously obtained IV------------------------------
startBB <- coef(fm_gengam)
fm_gengam_ext <- renewalCount(formula = regModelSQ, data = fertility,
                              dist = "gengamma",
                              control = renewal.control(method = "spg", start = startBB, trace = 0),
                              computeHessian = FALSE)


# ----check convergence status---------------------------------------------
fm_gengam_ext$converged


# ----set link functions---------------------------------------------------
parNames <- c("scale", "shape")
sWei <- function(tt, distP) exp( -distP[["scale"]] * tt ^ distP[["shape"]])
link <- list(scale = "log", shape = "log")


# -------------------------------------------------------------------------
control_custom <- renewal.control(start = startW, trace = 0)


# -------------------------------------------------------------------------
.getExtrapol <- function(distP) {
    c(2, distP[["shape"]])
}


# -------------------------------------------------------------------------
customPars <- list(parNames = parNames, survivalFct = sWei,
                   extrapolFct = .getExtrapol)


# ----user passed inter-arrival time model---------------------------------
weiModelCust <- renewalCount(formula = regModel, data = fertility,
                             dist = "custom", link = link,
                             control = control_custom,
                             customPars = customPars,
                             computeHessian = FALSE)


# -------------------------------------------------------------------------
summary(gamModel)

# -------------------------------------------------------------------------
summary(weiModel)


# ----bootstrap------------------------------------------------------------
se_boot <- se.coef(object = weiModel, type =  "boot", R = 5)
confint_boot <- confint(object = weiModel, type = "boot", R = 5)


# ----prepare data for prediction------------------------------------------
newData <- head(fertility)


# ----perform predictions--------------------------------------------------
predNew.response <- predict(weiModel, newdata = newData,
                            type = "response", se.fit = TRUE)
predNew.prob <- predict(weiModel, newdata = newData, type = "prob", se.fit = TRUE)

# -------------------------------------------------------------------------
options(digits = 5)


# -------------------------------------------------------------------------
predtable <- data.frame(newData$children, predNew.prob$values, predNew.response$values)
names(predtable) <- c("Y", "P(Y=y|x)", "E(Y|x)")
predtable


# -------------------------------------------------------------------------
options(digits = 3)


# -------------------------------------------------------------------------
cbind(builtIn = coef(weiModel), user = coef(weiModelCust))


# ----load-data------------------------------------------------------------
data("quine", package = "MASS")


# ----quine-table----------------------------------------------------------
breaks_ <- c(0, 1, 3, 5:7, 9, 12, 15, 17, 23, 27, 32)
freqtable <- count_table(count = quine$Days, breaks = breaks_, formatChar = TRUE)


# -------------------------------------------------------------------------
 print(xtable(freqtable[ , 1:7]), floating = FALSE, only.contents = TRUE)
 cat("\n\\\\[5pt]\n")
 print(xtable(freqtable[ , -(1:7)]), floating = FALSE, only.contents = TRUE)



# ----quine data models----------------------------------------------------
quine_form <- as.formula(Days ~ Eth + Sex + Age + Lrn)
pois <- glm(quine_form, family = poisson(), data = quine)
nb <- MASS::glm.nb(quine_form, data = quine)

## various renewal models
wei <- renewalCount(formula = quine_form, data = quine, dist = "weibull",
                    computeHessian = FALSE, weiMethod = "conv_dePril",
                    control = renewal.control(trace = 0))

gam <- renewalCount(formula = quine_form, data = quine, dist = "gamma",
                    computeHessian = FALSE,
                    control = renewal.control(trace = 0))

gengam <- renewalCount(formula = quine_form, data = quine, dist = "gengamma",
                       computeHessian = FALSE,
                       control = renewal.control(trace = 0))


# ----lr-test-poisson------------------------------------------------------
library("lmtest")
pois_nb <- lrtest(pois, nb)
pois_wei <- suppressWarnings(lrtest(pois, wei))
pois_gam <- suppressWarnings(lrtest(pois, gam))
pois_gengam <- suppressWarnings(lrtest(pois, gengam))
pois_res <- data.frame("Alternative model" =
                           c("Negative-binomial", "Weibull", "Gamma",
                             "Generalized-gamma"),
                       Chisq = c(pois_nb$Chisq[2], pois_wei$Chisq[2],
                                 pois_gam$Chisq[2], pois_gengam$Chisq[2]),
                       Df = c(1, 1, 1, 2),
                       Critical_value = c(rep(qchisq(0.95, 1), 3), qchisq(0.95, 2)),
                       stringsAsFactors = FALSE)


# -------------------------------------------------------------------------
print(xtable(pois_res, caption = "LR results against Poisson model. Each row compares an alternative model vs the Poisson model. All alternatives are preferable to Poisson.
The critical value corresponds to a significance level of 5\\%",
             label = "tab:lr_pois"))


# ----lr-test-renewal------------------------------------------------------
gengam_wei <- lrtest(wei, gengam)
gengam_gam <- lrtest(gam, gengam)
gengam_res <- data.frame(Model = c("Weibull", "Gamma"),
  Chisq = c(gengam_wei$Chisq[2], gengam_gam$Chisq[2]), Df = 1,
  Critical_value = rep(qchisq(0.95, 1), 2), stringsAsFactors = FALSE)


# -------------------------------------------------------------------------
print(xtable(gengam_res, caption = "LR results against generalized-gamma model",
             label = "tab:lr_gengam"))


# ----AIC-models-----------------------------------------------------------
ic <- data.frame(Model = c("Gamma", "Weibull", "Negative-binomial",
                           "Generalized-gamma", "Poisson"),
                 AIC = c(AIC(gam), AIC(wei), AIC(nb), AIC(gengam), AIC(pois)),
                 BIC = c(BIC(gam), BIC(wei), BIC(nb), BIC(gengam), BIC(pois)),
		 stringsAsFactors = FALSE)


# ----table of AIC results-------------------------------------------------
print(xtable(ic, caption = "Information criteria results", label = "tab:ic_models"),
      hline.after = c(0, 3), type = "latex")


# ----goodness of fit------------------------------------------------------
gof <- chiSq_gof(gam, breaks = breaks_)
gof
