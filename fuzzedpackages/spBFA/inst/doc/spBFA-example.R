## ---- echo = FALSE-------------------------------------------------------
# ###Start with a clean space
# rm(list = ls())
# 
# ###Take care of some stuff that I don't want the user to see...
# path.package <- "/Users/sam/Documents/Postdoc/Software/spBFA/"
# suppressMessages(devtools::load_all(path.package)) #loads scripts
# suppressMessages(devtools::document(path.package)) #creates documentation
###Make sure to remove devtools from Suggests line in DESCRIPTION before submission

## ------------------------------------------------------------------------
library(womblR)
library(spBFA)

## ------------------------------------------------------------------------
head(VFSeries)

## ---- fig.align="center", fig.width = 5.5, fig.height = 5.5--------------
PlotVfTimeSeries(Y = VFSeries$DLS,
                 Location = VFSeries$Location,
                 Time = VFSeries$Time,
                 main = "Visual field sensitivity time series \n at each location",
                 xlab = "Days from baseline visit",
                 ylab = "Differential light sensitivity (dB)",
                 line.col = 1, line.type = 1, line.reg = FALSE)

## ------------------------------------------------------------------------
blind_spot <- c(26, 35) # define blind spot
VFSeries <- VFSeries[order(VFSeries$Location), ] # sort by location
VFSeries <- VFSeries[order(VFSeries$Visit), ] # sort by visit
VFSeries <- VFSeries[!VFSeries$Location %in% blind_spot, ] # remove blind spot locations
dat <- data.frame(Y = VFSeries$DLS / 10) # create data frame with scaled data

## ------------------------------------------------------------------------
Time <- unique(VFSeries$Time) / 365 # years since baseline visit
print(Time)

## ------------------------------------------------------------------------
W <- HFAII_Queen[-blind_spot, -blind_spot] # visual field adjacency matrix
M <- dim(W)[1] # number of locations

## ------------------------------------------------------------------------
TimeDist <- as.matrix(dist(Time))
BPsi <- log(0.025) / -min(TimeDist[TimeDist > 0])
APsi <- log(0.975) / -max(TimeDist)

## ------------------------------------------------------------------------
K <- 10
O <- 1
Hypers <- list(Sigma2 = list(A = 0.001, B = 0.001),
               Kappa = list(SmallUpsilon = O + 1, BigTheta = diag(O)),
               Delta = list(A1 = 1, A2 = 20),
               Psi = list(APsi = APsi, BPsi = BPsi),
               Upsilon = list(Zeta = K + 1, Omega = diag(K)))

## ------------------------------------------------------------------------
Starting <- list(Sigma2 = 1,
                 Kappa = diag(O),
                 Delta = 2 * (1:K),
                 Psi = (APsi + BPsi) / 2,
                 Upsilon = diag(K))

## ------------------------------------------------------------------------
Tuning <- list(Psi = 1)

## ------------------------------------------------------------------------
MCMC <- list(NBurn = 1000, NSims = 1000, NThin = 2, NPilot = 5)

## ---- include = FALSE----------------------------------------------------
data(reg.bfa_sp)

## ---- eval = FALSE-------------------------------------------------------
#  reg.bfa_sp <- bfa_sp(Y ~ 0, data = dat, dist = W, time = Time,  K = 10,
#                       starting = Starting, hypers = Hypers, tuning = Tuning, mcmc = MCMC,
#                       L = Inf,
#                       family = "tobit",
#                       trials = NULL,
#                       temporal.structure = "exponential",
#                       spatial.structure = "discrete",
#                       seed = 54,
#                       gamma.shrinkage = TRUE,
#                       include.space = TRUE,
#                       clustering = TRUE)
#  
#  ## Burn-in progress:  |*************************************************|
#  ## Sampler progress:  0%..  10%..  20%..  30%..  40%..  50%..  60%..  70%..  80%..  90%..  100%..

## ------------------------------------------------------------------------
names(reg.bfa_sp)

## ------------------------------------------------------------------------
library(coda)

## ------------------------------------------------------------------------
Sigma2_1 <- as.mcmc(reg.bfa_sp$sigma2[, 1])

## ---- fig.width = 5.2, fig.height = 5.2, echo = FALSE--------------------
par(mfrow = c(1, 1))
traceplot(Sigma2_1, ylab = expression(paste(sigma^2 ~ "(" ~ s[1]~ ")")), main = expression(paste("Posterior" ~ sigma^2 ~ "(" ~ s[1]~ ")")))

## ---- echo = FALSE-------------------------------------------------------
geweke.diag(Sigma2_1)$z

## ------------------------------------------------------------------------
Diags <- spBFA::diagnostics(reg.bfa_sp, diags = c("dic", "dinf", "waic"), keepDeviance = TRUE)

## ---- fig.align = 'center', fig.width = 4, fig.height = 3.3--------------
Deviance <- as.mcmc(Diags$deviance)
traceplot(Deviance, ylab = "Deviance", main = "Posterior Deviance")

## ---- eval = FALSE-------------------------------------------------------
#  print(Diags)

## ---- echo = FALSE-------------------------------------------------------
unlist(Diags$dic)
unlist(Diags$dinf)
unlist(Diags$waic)

## ------------------------------------------------------------------------
NewTimes <- 3

## ------------------------------------------------------------------------
Predictions <- predict(reg.bfa_sp, NewTimes)

## ------------------------------------------------------------------------
names(Predictions)

## ---- fig.align = 'center', fig.width = 4.5, fig.height = 4.5------------
PlotSensitivity(Y = apply(Predictions$Y$Y10, 2, mean) * 10,
                main = "Posterior mean prediction\n at 3 years",
                legend.lab = "Posterior Mean", legend.round = 2,
                bins = 250, border = FALSE, zlim = c(0, 40))

## ---- fig.align = 'center', fig.width = 4.5, fig.height = 4.5------------
PlotSensitivity(Y = apply(Predictions$Y$Y10 * 10, 2, sd),
                main = "Posterior standard deviation\n (SD) at 3 years",
                legend.lab = "Posterior SD", legend.round = 2,
                bins = 250, border = FALSE, zlim = c(0, 40))

