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
Data <- Seatbelts

# Data preparation
# The log of "drivers", "front", and "rear" will be used as dependent variables
# The log of "PetrolPrice" and "kms" will be used as explanatory variables
Data[, c("drivers", "front", "rear", "PetrolPrice", "kms")] <- log(Data[, c("drivers", "front", "rear", "PetrolPrice", "kms")])

## -----------------------------------------------------------------------------
# Dependent variable
y <- as.matrix(Data[, "drivers"])

# Period of the seasonal component
BSM_vec <- 12

# Explanatory variables
# Note: Must be a list of matrices! 
#       Each dependent gets its own matrix of explanatory variables.
addvar_list <- list(as.matrix(Data[, c("PetrolPrice", "law")]))

# Format of the variance - covariance matrix of the level component
# By setting the elements of this matrix to 0, 
# the component becomes deterministic.
format_level <- matrix(0)

# Format of the variance - covariance matrix of the seasonal component
# Note: This format must be a list of matrices, because multiple 
#       seasonalities can be specified!
format_BSM_list <- list(matrix(0))

# Fitting the model
fit <- statespacer(y = y,
                   local_level_ind = TRUE,
                   BSM_vec = BSM_vec,
                   addvar_list = addvar_list,
                   format_level = format_level, 
                   format_BSM_list = format_BSM_list,
                   method = "BFGS",
                   initial = 0.5 * log(var(y)),
                   verbose = TRUE)

## ---- fig.height = 4.5, fig.width = 7-----------------------------------------
# The estimated variance of the observation disturbance
fit$system_matrices$H$H

# Smoothed estimate of the level
fit$smoothed$level[1,]

# Smoothed estimate of the coefficient of log "PetrolPrice"
fit$smoothed$addvar_coeff[1, 1]

# Smoothed estimate of the coefficient of the "law" dummy
fit$smoothed$addvar_coeff[1, 2]

# Plot the data next to the smoothed level + explanatory variables components
plot(Data[, c("drivers")], type = "l", ylim = c(6.95, 8.1),
     xlab = "year", ylab = "logarithm of drivers")
lines(seq(tsp(Data)[1], tsp(Data)[2], 1/tsp(Data)[3]), 
      fit$smoothed$level + fit$smoothed$addvar, type = 'l', col = "red")
legend(1978, 8.09, c("log(drivers)", "level + regression effects"), 
       lty = c(1,1), lwd=c(2.5, 2.5), col = c("black", "red"))

## ---- fig.height = 4.5, fig.width = 7, warning = FALSE------------------------
# By setting the entries in the format to 1, the component becomes stochastic
format_level <- matrix(1)
format_BSM_list <- list(matrix(1))
fit <- statespacer(y = y,
                   local_level_ind = TRUE,
                   BSM_vec = BSM_vec,
                   addvar_list = addvar_list,
                   format_level = format_level,
                   format_BSM_list = format_BSM_list,
                   method = "BFGS",
                   initial = log(var(y)),
                   verbose = TRUE)

# The estimated variance of the observation disturbance
fit$system_matrices$H$H

# The estimated variance of the level disturbance
fit$system_matrices$Q$level

# The estimated variance of the seasonal disturbance
fit$system_matrices$Q$BSM12

# Smoothed estimate of the coefficient of log "PetrolPrice"
fit$smoothed$addvar_coeff[1, 1]

# Smoothed estimate of the coefficient of the "law" dummy
fit$smoothed$addvar_coeff[1, 2]

# Plot the data next to the smoothed level + explanatory variables components
plot(Data[, c("drivers")], type = "l", ylim = c(6.95, 8.1),
     xlab = "year", ylab = "logarithm of drivers")
lines(seq(tsp(Data)[1], tsp(Data)[2], 1/tsp(Data)[3]), 
      fit$smoothed$level + fit$smoothed$addvar, type = 'l', col = "red")
legend(1978, 8.09, c("log(drivers)", "level + regression effects"), 
       lty = c(1,1), lwd=c(2.5, 2.5), col = c("black", "red"))

## ---- fig.height = 4.5, fig.width = 7-----------------------------------------
# Plot the stochastic seasonal
plot(seq(tsp(Data)[1], tsp(Data)[2], 1/tsp(Data)[3]), 
     fit$smoothed$BSM12,
     type = "l", ylim = c(-0.2, 0.3),
     xlab = "year", ylab = "stochastic seasonal")
abline(h = 0)

## ---- fig.height = 4.5, fig.width = 7-----------------------------------------
plot(seq(tsp(Data)[1], tsp(Data)[2], 1/tsp(Data)[3]), 
     fit$smoothed$epsilon,
     type = "l", ylim = c(-0.15, 0.15),
     xlab = "year", ylab = "irregular component")
abline(h = 0)

## ---- warning = FALSE---------------------------------------------------------
# Dependent variable
y <- as.matrix(Data[, c("front", "rear")])

# Explanatory variables
# Note: Must be a list of matrices! 
#       Each dependent gets its own matrix of explanatory variables.
X <- as.matrix(Data[, c("PetrolPrice", "kms", "law")])
addvar_list <- list(X, X)

# Format of the variance - covariance matrix of the level component
# Note: Only the lower triangular part of the format is used.
#       The format specifies which elements in the matrices L and D should be 
#       non-zero, where L and D are the matrices of the Cholesky LDL decomposition.
#       The diagonal is used to specify which elements of the Diagonal matrix
#       should be non-zero. The lower triangular part excluding the diagonal 
#       specifies which elements in the Loading matrix should be non-zero.
format_level <- matrix(1, 2, 2)

# Format of the variance - covariance matrix of the seasonal component
# Note: This format must be a list of matrices, because multiple seasonalities
#       can be specified!
format_BSM_list <- list(matrix(0, 2, 2))

# Format of the variance - covariance matrix of the observation disturbances
H_format <- matrix(1, 2, 2)

# Fitting the model
fit <- statespacer(y = y,
                   H_format = H_format,
                   local_level_ind = TRUE,
                   BSM_vec = BSM_vec,
                   addvar_list = addvar_list,
                   format_level = format_level, 
                   format_BSM_list = format_BSM_list,
                   method = "BFGS",
                   initial = 0.5 * log(diag(var(y))),
                   verbose = TRUE)

## ---- fig.height = 4.5, fig.width = 7-----------------------------------------
# The estimated variance - covariance matrix of the observation disturbance
fit$system_matrices$H$H

# The estimated variance - covariance matrix of the level disturbance
fit$system_matrices$Q$level

# Coefficients + T-stats
coeff <- cbind(
  c("front PetrolPrice", "front kms", "front law", 
    "rear PetrolPrice", "rear kms", "rear law"),
  fit$smoothed$addvar_coeff[1,],
  fit$smoothed$addvar_coeff[1,] / fit$smoothed$addvar_coeff_se[1,]
)
colnames(coeff) <- c("Variable", "Coefficient", "T-Stat")
coeff

# plot of level for "front"
plot(seq(tsp(Data)[1], tsp(Data)[2], 1/tsp(Data)[3]), 
     fit$smoothed$level[, 1], type = "l",
     xlab = "year", ylab = "level of front passengers")

# plot of level for "rear"
plot(seq(tsp(Data)[1], tsp(Data)[2], 1/tsp(Data)[3]), 
     fit$smoothed$level[, 2], type = "l",
     xlab = "year", ylab = "level of rear passengers")

## -----------------------------------------------------------------------------
fit$system_matrices$Q_correlation_matrix$level

## -----------------------------------------------------------------------------
format_level <- matrix(0, 2, 2)
format_level[, 1] <- 1

