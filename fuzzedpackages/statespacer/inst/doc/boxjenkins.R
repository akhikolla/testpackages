## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- setup-------------------------------------------------------------------
# Load statespacer
library(statespacer)

# Load the dataset
library(datasets)
Data <- matrix(window(sunspot.year, start = 1770, end = 1869))

# Estimate the ARIMA model
fit <- statespacer(y = Data,
                   H_format = matrix(0),
                   local_level_ind = TRUE,
                   arima_list = list(c(3, 0, 0)),
                   format_level = matrix(0),
                   initial = c(0.5*log(var(Data)), 0, 0, 0),
                   verbose = TRUE)

## -----------------------------------------------------------------------------
# Coefficients of the ARMA component
arma_coeff <- rbind(
   fit$system_matrices$AR$ARIMA1,
   fit$standard_errors$AR$ARIMA1
)
arma_coeff <- cbind(
   arma_coeff,
   c(fit$smoothed$level[1],
     sqrt(fit$system_matrices$Z_padded$level %*%
          fit$smoothed$V[,,1] %*%
          t(fit$system_matrices$Z_padded$level))
   )
)
rownames(arma_coeff) <- c("coefficient", "std_error")
colnames(arma_coeff) <- c("ar1", "ar2", "ar3", "intercept")
arma_coeff

goodness_fit <- rbind(
   fit$system_matrices$Q$ARIMA1,
   fit$diagnostics$loglik,
   fit$diagnostics$AIC
)
rownames(goodness_fit) <- c("Variance", "Loglikelihood", "AIC")
goodness_fit

## ---- warning = FALSE---------------------------------------------------------
# Load the dataset
Data <- matrix(log(AirPassengers))

# The SARIMA specification, must be a list containing lists!
sarima_list <- list(list(s = c(12, 1), ar = c(0, 0), i = c(1, 1), ma = c(1, 1)))

# Fit the SARIMA model
fit <- statespacer(y = Data,
                   H_format = matrix(0),
                   sarima_list = sarima_list,
                   initial = c(0.5*log(var(diff(Data))), 0, 0),
                   verbose = TRUE)

## -----------------------------------------------------------------------------
# Coefficients of the ARMA component
arma_coeff <- rbind(
   c(fit$system_matrices$SMA$SARIMA1$S1, fit$system_matrices$SMA$SARIMA1$S12),
   c(fit$standard_errors$SMA$SARIMA1$S1, fit$standard_errors$SMA$SARIMA1$S12)
)

rownames(arma_coeff) <- c("coefficient", "std_error")
colnames(arma_coeff) <- c("ma1 s = 1", "ma1 s = 12")
arma_coeff

goodness_fit <- rbind(
   fit$system_matrices$Q$SARIMA1,
   fit$diagnostics$loglik,
   fit$diagnostics$AIC
)
rownames(goodness_fit) <- c("Variance", "Loglikelihood", "AIC")
goodness_fit

