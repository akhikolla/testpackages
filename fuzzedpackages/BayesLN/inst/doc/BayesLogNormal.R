## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

library(BayesLN)

## -----------------------------------------------------------------------------
# Load dataset
data("EPA09")

# Bayes estimator under relative quadratic loss and optimal prior setting
LN_Mean(x = EPA09, x_transf = FALSE, method = "optimal", CI = FALSE)


## -----------------------------------------------------------------------------
LN_Mean(x = EPA09, x_transf = FALSE, method = "weak_inf", alpha_CI = 0.05, type_CI = "UCL")

## -----------------------------------------------------------------------------
LN_Quant(x = EPA09, quant = 0.95, method = "optimal", CI = FALSE)

## -----------------------------------------------------------------------------
# Load data
data("fatigue")

# Design matrices
Xtot <- cbind(1, log(fatigue$stress), log(fatigue$stress)^2)
X <- Xtot[-c(1,13,22),]
y <- fatigue$cycle[-c(1,13,22)]
Xtilde <- Xtot[c(1,13,22),] # units to predict

#Estimation

LN_MeanReg(y = y,
           X = X, Xtilde = Xtilde,
           method = "weak_inf", y_transf = FALSE)

## -----------------------------------------------------------------------------
# Load the dataset included in the package 
data("ReadingTime")

# Define data.frame containing the covariate patterns to investigate
data_pred_new <- expand.grid(so=c(-1,1), subj=factor(12), item=factor(8))

# Model estimation 
Mod_est_RT <- LN_hierarchical(formula_lme = log_rt ~ so +(1|subj)+(1|item), 
                              data_lme = ReadingTime, data_pred = data_pred_new, 
                              functional = c("Marginal", "Subject"), 
                              nsamp = 25000, burnin = 5000, n_thin = 5)

## -----------------------------------------------------------------------------
# Prior parameters
Mod_est_RT$par_prior

## ---- fig.width = 6.5---------------------------------------------------------
# coda package
library(coda)
# Traceplots model parameters
oldpar <- par(mfrow=c(2,3))
traceplot(Mod_est_RT$samples$par[, 1:6])
par(oldpar)

## -----------------------------------------------------------------------------
# Posterior summaries
Mod_est_RT$summaries$marg
Mod_est_RT$summaries$subj

