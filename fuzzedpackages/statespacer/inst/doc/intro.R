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
y <- matrix(Nile)

## -----------------------------------------------------------------------------
fit <- statespacer(y = y,
                   local_level_ind = TRUE,
                   initial = 0.5*log(var(y)),
                   verbose = TRUE)

## -----------------------------------------------------------------------------
c(fit$system_matrices$H$H, fit$system_matrices$Q$level)

## ---- fig.height = 4.5, fig.width = 7-----------------------------------------
plot(1871:1970, fit$function_call$y, type = 'p', ylim = c(500, 1400), 
     xlab = NA, ylab = NA, 
     sub = "The filtered level with 90% confidence intervals, 
            and the observed data points"
)
lines(1871:1970, fit$filtered$level, type = 'l')
lines(1871:1970, fit$filtered$level + qnorm(0.95) * sqrt(fit$filtered$P[1,1,]), 
      type = 'l', col = 'gray'
)
lines(1871:1970, fit$filtered$level - qnorm(0.95) * sqrt(fit$filtered$P[1,1,]), 
      type = 'l', col = 'gray'
)

## ---- fig.height = 4.5, fig.width = 7-----------------------------------------
plot(1871:1970, fit$function_call$y, type = 'p', ylim = c(500, 1400), 
     xlab = NA, ylab = NA, 
     sub = "The smoothed level with 90% confidence intervals, 
            and the observed data points")
lines(1871:1970, fit$smoothed$level, type = 'l')
lines(1871:1970, fit$smoothed$level + qnorm(0.95) * sqrt(fit$smoothed$V[1,1,]),
      type = 'l', col = 'gray'
)
lines(1871:1970, fit$smoothed$level - qnorm(0.95) * sqrt(fit$smoothed$V[1,1,]),
      type = 'l', col = 'gray'
)

