## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- setup-------------------------------------------------------------------
# Load statespacer
library(statespacer)

# Load the dataset
library(YieldCurve)
data("FedYieldCurve")
y <- window(FedYieldCurve, start = "1984-12-31", end = "2000-12-01")
years <- index(y) # Used for plots later on
y <- as.matrix(y)

## -----------------------------------------------------------------------------
# Specifying the list
self_spec_list <- list()

# We want to specify the H matrix ourselves
self_spec_list$H_spec <- TRUE

# We have got 6 state parameters: 3 factors and 3 fixed means
self_spec_list$state_num <- 6

# In total we need 20 parameters:
#   1 for lambda
#   1 for sigma2 (H)
#   6 for the variance - covariance matrix Sigma_eta (Q)
#   9 for the vector autoregressive coefficient matrix Phi
#   3 for the means mu
self_spec_list$param_num <- 20

# R is a fixed diagonal matrix
self_spec_list$R <- diag(1, 6, 6)

# P_inf is a matrix of zeroes, as all state parameters are stationary
self_spec_list$P_inf <- matrix(0, 6, 6)

# Needed because we want to use collapse = TRUE
# The fixed means only appear in the state equations, 
# not in the observation equations. So the 4th, 5th, and 6th state parameters
# are state_only.
self_spec_list$state_only <- 4:6

## -----------------------------------------------------------------------------
self_spec_list$sys_mat_fun <- function(param) {
  
  # Maturities of the interest rates
  maturity <- c(3, 6, 12, 24, 36, 60, 84, 120)
  
  # The constant lambda
  lambda <- exp(2 * param[1])
  
  # The variance of the observation errors
  sigma2 <- exp(2 * param[2])
  H <- sigma2 * diag(1, 8, 8)

  # Z matrix corresponding to the factors
  lambda_maturity <- lambda * maturity
  z <- exp(-lambda_maturity)
  Z <- matrix(1, 8, 3)
  Z[, 2] <- (1 - z) / lambda_maturity
  Z[, 3] <- Z[, 2] - z

  # Variance of the state disturbances
  Q <- Cholesky(param = param[3:8], decompositions = FALSE, format = matrix(1, 3, 3))
  
  # Vector autoregressive coefficient matrix, enforcing stationarity
  Tmat <- CoeffARMA(A = array(param[9:17], dim = c(3, 3, 1)),
                 variance = Q,
                 ar = 1, ma = 0)$ar[,,1]

  # Initial uncertainty of the factors
  T_kronecker <- kronecker(Tmat, Tmat)
  Tinv <- solve(diag(1, dim(T_kronecker)[1], dim(T_kronecker)[2]) - T_kronecker)
  vecQ <- matrix(Q)
  vecPstar <- Tinv %*% vecQ
  P_star <- matrix(vecPstar, dim(Tmat)[1], dim(Tmat)[2])

  # Adding parts corresponding to the fixed means to the system matrices
  Z <- cbind(Z, matrix(0, 8, 3)) # Not used in the observation equation
  Q <- BlockMatrix(Q, matrix(0, 3, 3)) # Fixed, so no variance in its errors
  a1 <- matrix(param[18:20], 6, 1) # Fixed means go into the initial guess
  Tmat <- cbind(Tmat, diag(1, 3, 3) - Tmat)
  Tmat <- rbind(Tmat, cbind(matrix(0, 3, 3), diag(1, 3, 3)))
  P_star <- BlockMatrix(P_star, matrix(0, 3, 3))

  # Return the system matrices
  return(list(H = H, Z = Z, Tmat = Tmat, Q = Q, a1 = a1, P_star = P_star))
}

## -----------------------------------------------------------------------------
self_spec_list$transform_fun <- function(param) {
  lambda <- exp(2 * param[1])
  sigma2 <- exp(2 * param[2])
  means <- param[18:20]
  return(c(lambda, sigma2, means))
}

## -----------------------------------------------------------------------------
initial <- c(-1, -2, -1, -1, 1, 0, 0, 0, 4, 0, 0, 0, 3, 0, 0, 0, 2, 0, 0, 0)

## -----------------------------------------------------------------------------
# Optimal parameters
initial <- c(-1.27000018, -2.82637721, -1.35280609, -1.74569078,
             -0.89761151, -0.77767132, 1.01894143, 0.42965982,
             13.69072496, 3.46050373, -10.36767668, -0.07334641,
             6.68658053, 0.76975206, 0.02852844, 0.50448668,
             2.99984132, 8.16851107, -2.28360681, -0.45333494)

## -----------------------------------------------------------------------------
fit <- statespacer(y = y,
                   self_spec_list = self_spec_list,
                   collapse = TRUE,
                   initial = initial,
                   method = "BFGS",
                   verbose = TRUE)

## -----------------------------------------------------------------------------
# The level beta_1
plot(years, fit$smoothed$a[, 1], type = 'l', 
     xlab = "year", ylab = "level of yield curve")

# The slope beta_2
plot(years, fit$smoothed$a[, 2], type = 'l',
     xlab = "year", ylab = "slope of yield curve")

# The shape beta_3
plot(years, fit$smoothed$a[, 3], type = 'l',
     xlab = "year", ylab = "shape of yield curve")

## -----------------------------------------------------------------------------
parameters <- cbind(
  c("lambda", "sigma2", "mu1", "mu2", "mu3"), 
  fit$system_matrices$self_spec,
  fit$standard_errors$self_spec
)
colnames(parameters) <- c("Parameter", "Value", "Standard Error")
parameters

# Vector autoregressive coefficient matrix
fit$system_matrices$T$self_spec[1:3, 1:3]

# Variance of the state disturbances
fit$system_matrices$Q$self_spec[1:3, 1:3]

