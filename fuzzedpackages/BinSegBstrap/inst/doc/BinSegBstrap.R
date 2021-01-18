## ----setup, echo = FALSE, message = FALSE--------------------------------
library(knitr)
library(BinSegBstrap)

## ----estimateSingleCpFixedBandwidth, fig.cap = 'Observations (grey points), underlying signal (black line) and estimated signal (red line).'----
set.seed(1)
n <- 100
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5

y <- rnorm(n) + signal

# call of estimateSingleCp with fixed bandwidth 0.1
est <- estimateSingleCp(y = y, bandwidth = 0.1)

# estimated location
est$cp

# estimated jump size
est$size

# plot of observations, true and estimated signal
plot(y, pch = 16, col = "grey30")
lines(signal)
lines(est$est, col = "red")

## ----estimateSingleCp, fig.cap = 'Observations (grey points), underlying signal (black line) and estimated signal (red line).'----
set.seed(1)
n <- 100
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5

y <- rnorm(n) + signal

# call of estimateSingleCp with crossvalidated bandwidth
est <- estimateSingleCp(y = y)

# crossvalidated bandwidth
est$bandwidth

# estimated location
est$cp

# estimated jump size
est$size

# plot of observations, true and estimated signal
plot(y, pch = 16, col = "grey30")
lines(signal)
lines(est$est, col = "red")

## ----BstrapTest----------------------------------------------------------
set.seed(1)
n <- 100
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5

y <- rnorm(n) + signal

test <- BstrapTest(y = y)

# whether the test rejected
test$outcome

# p-Value
test$pValue

## ----BinSegBstrap, fig.cap = 'Observations (grey points), underlying signal (black line) and estimated signal (red line).'----
set.seed(1)
n <- 200
signal <- sin(2 * pi * 1:n / n)
signal[51:100] <- signal[51:100] + 5
signal[151:200] <- signal[151:200] + 5

y <- rnorm(n) + signal

est <- BinSegBstrap(y = y)

# estimated change-points
est$cps

# plot of observations, true and estimated signal
plot(y, pch = 16, col = "grey30")
lines(signal)
lines(est$est, col = "red")

