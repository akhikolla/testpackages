## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=6
)

## -----------------------------------------------------------------------------
# Load the hdme package
library(hdme)

## -----------------------------------------------------------------------------
create_example_data <- function(n, p, s = 5, sdX = 1, sdU = 0.5, 
                                sdEpsilon = 0.1, family = "gaussian") {
  # Independent true covariates with mean zero and standard deviation sdX
  X <- matrix(rnorm(n * p, sd = sdX), nrow = n, ncol = p)
  # if Gaussian, response has standard deviation sdEpsilon, and zero intercept
  # if binomial, response is binomial with mean (1 + exp(-X %*% beta))^(-1)
  beta <- c(-2, -1, 0.5, 1, 2, rep(0, p - s))
  
  if(family == "gaussian") {
    # True coefficient vector has s non-zero elements and p-s zero elements
    y <- X %*% beta + rnorm(n, sd = sdEpsilon)  
  } else if (family == "binomial") {
    # Need an amplification in the binomial case
    beta <- beta * 3
    y <- rbinom(n, size = 1, prob = (1 + exp(-X %*% beta))**(-1))
  }
  
  # The measurements W have mean X and standard deviation sdU. 
  # We assume uncorrelated measurement errors
  W <- X + matrix(rnorm(n * p, sd = sdU), nrow = n, ncol = p)
  
  return(list(X = X, W = W, y = y, beta = beta, sigmaUU = diag(p) * sdU))  
}

## ---- message=FALSE-----------------------------------------------------------
n <- 100
p <- 500
set.seed(1000)
ll <- create_example_data(n, p)

## ---- message=FALSE-----------------------------------------------------------
library(glmnet)
library(dplyr)
# Lasso with cross-validation on data without measurement error
fit1 <- cv.glmnet(ll$X, ll$y)
# Lasso with cross-validation on data with measurement error
fit2 <- cv.glmnet(ll$W, ll$y)
# Create a data frame with results ([-1] because we drop the intercept)
lassoEstimates <- tibble(
  index = rep(1:p, times = 3),
  beta = c(ll$beta, as.numeric(coef(fit1)[-1]), coef(fit2)[-1]),
  label = c(rep("True values", p), rep("No measurement error", p), rep("Measurement error", p))
  )


## -----------------------------------------------------------------------------
library(ggplot2)
ggplot(lassoEstimates, aes(x = index, y = beta, color = label)) +
  geom_point() +
  xlab("p") +
  theme(legend.title=element_blank()) + 
  ggtitle("Measurement error leading to false positives")

## ---- message=FALSE, warning=FALSE--------------------------------------------
library(tidyr) 
estimatesOfNonzero <- lassoEstimates %>% 
  spread(key = label, value = beta) %>% 
  filter(`True values` != 0) %>% 
  gather(key = label, value = beta, -index)

ggplot(estimatesOfNonzero, aes(x = index, y = beta, color = label)) +
  geom_point() +
  xlab("p") +
  theme(legend.title=element_blank()) + 
  ggtitle("Measurement error leading to attenuation")

## -----------------------------------------------------------------------------
# Number of samples
n <- 1000
# Number of covariates
p <- 50
# Create example data
ll <- create_example_data(n, p, family = "binomial")

## -----------------------------------------------------------------------------
args(gds)

## -----------------------------------------------------------------------------
# Fit the Generalized Dantzig Selector
gds_estimate <- gds(ll$X, ll$y, family = "binomial")

## -----------------------------------------------------------------------------
class(gds_estimate)

## -----------------------------------------------------------------------------
str(gds_estimate)

## -----------------------------------------------------------------------------
set.seed(1000)
# Generate example data
ll <- create_example_data(n, p)
# Fit the corrected lasso
corrected_fit <- corrected_lasso(W = ll$W, y = ll$y, sigmaUU = ll$sigmaUU)

## -----------------------------------------------------------------------------
# Class of the object
class(corrected_fit)
# The coef() method prints the number of nonzero estimates as a function of the radius
coef(corrected_fit)

## -----------------------------------------------------------------------------
args(corrected_lasso)

## -----------------------------------------------------------------------------
plot(corrected_fit)

## -----------------------------------------------------------------------------
plot(corrected_fit, type = "path")

## -----------------------------------------------------------------------------
set.seed(323)
n <- 100
p <- 50
ll <- create_example_data(n, p, sdU = 0.2, family = "binomial")

## -----------------------------------------------------------------------------
corrected_fit <- corrected_lasso(ll$W, ll$y, ll$sigmaUU, family = "binomial")

## -----------------------------------------------------------------------------
plot(corrected_fit)

## -----------------------------------------------------------------------------
plot(corrected_fit, type = "path")

## -----------------------------------------------------------------------------
set.seed(1000)
# Generate example data
ll <- create_example_data(n, p)
# Run lasso with cross-validation
cv_corrected_fit <- cv_corrected_lasso(W = ll$W, y = ll$y, sigmaUU = ll$sigmaUU)

## -----------------------------------------------------------------------------
class(cv_corrected_fit)

## -----------------------------------------------------------------------------
str(cv_corrected_fit)

## -----------------------------------------------------------------------------
plot(cv_corrected_fit)

## -----------------------------------------------------------------------------
corrected_fit <- corrected_lasso(ll$W, ll$y, ll$sigmaUU, radii = cv_corrected_fit$radius_1se)

## -----------------------------------------------------------------------------
str(corrected_fit)

## -----------------------------------------------------------------------------
set.seed(1)
# Number of samples
n <- 1000
# Number of covariates
p <- 50
# Generate data
ll <- create_example_data(n, p, sdU = 0.2)

## -----------------------------------------------------------------------------
mus_fit <- mus(ll$W, ll$y)
class(mus_fit)

## -----------------------------------------------------------------------------
coef(mus_fit)

## -----------------------------------------------------------------------------
plot(mus_fit)

## -----------------------------------------------------------------------------
mus_fit <- mus(ll$W, ll$y, delta = 0.1)

## -----------------------------------------------------------------------------
plot(mus_fit)

## -----------------------------------------------------------------------------
set.seed(323)
n <- 100
p <- 50
ll <- create_example_data(n, p, sdU = 0.2, family = "binomial")
gmus_fit <- gmus(ll$W, ll$y, family = "binomial")

## -----------------------------------------------------------------------------
class(gmus_fit)
str(gmus_fit)

## -----------------------------------------------------------------------------
plot(gmus_fit)

## -----------------------------------------------------------------------------
gmus_fit <- gmus(ll$W, ll$y, delta = 0.1, family = "binomial")

## -----------------------------------------------------------------------------
plot(gmus_fit)

## -----------------------------------------------------------------------------
set.seed(323)
n <- 100
p <- 50
ll <- create_example_data(n, p, sdU = 0.2, family = "binomial")
gmu_lasso_fit <- gmu_lasso(ll$W, ll$y, family = "binomial")

## -----------------------------------------------------------------------------
class(gmu_lasso_fit)
str(gmu_lasso_fit)

## -----------------------------------------------------------------------------
plot(gmu_lasso_fit)

