## ----setup, echo=FALSE, results='hide', message=TRUE---------------------
rm(list = ls())

do.calc <- FALSE # if TRUE, all calculations are conducted. if FALSE, calculations are omitted and the vignette output is built from pre-caculated data.
do.save <- FALSE # if TRUE, calculation parts will be saved to disk. This works only correctly, if building from source. Therefore, if you are not building the vignette from source: do.save <- FALSE

if (do.calc == FALSE) message("For performance reasons, this vignette builds from pre-calculated data.\n\n To run all calculations, set 'do.calc <- TRUE' in the vignette's first code chunk. \n Building the vignette will take a while.")

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center", 
  dev = "png",
  dpi=150, fig.height=7, fig.width=7,
  dev.args = list(),
  out.width = "90%"
)

op <- par(no.readonly = TRUE)

## ----load-guts-----------------------------------------------------------
library(GUTS)
packageVersion("GUTS")

## ----load diazinon data--------------------------------------------------
data(diazinon)
str(diazinon)

## ----setup-GUTS-OBJ-A_SD-------------------------------------------------
guts_object <- list( 
  C1 = guts_setup(
    C = diazinon[["C1"]], Ct = diazinon[["Ct1"]],
    y = diazinon[["y1"]], yt = diazinon[["yt1"]],
    model = "Proper", dist = "loglogistic"
    ),
  C2 = guts_setup(
    C = diazinon[["C2"]], Ct = diazinon[["Ct2"]],
    y = diazinon[["y2"]], yt = diazinon[["yt2"]],
    model = "Proper", dist = "loglogistic"
    ),
  C3 = guts_setup(
    C = diazinon[["C3"]], Ct = diazinon[["Ct3"]],
    y = diazinon[["y3"]], yt = diazinon[["yt3"]],
    model = "Proper", dist = "loglogistic"
    )
)

## ----load optimization routines, echo = TRUE, results = 'hide'-----------
library('adaptMCMC') # Function `MCMC()`, Monte Carlo Markov Chain.

## ----define-log-posterior------------------------------------------------
logposterior <- function( pars, guts_objects, isOutOfBoundsFun) {
	if ( isOutOfBoundsFun(pars) ) return(-Inf)
  return(
	  sum(sapply( guts_objects, function(obj) guts_calc_loglikelihood(obj, pars) )))
}

## ----define-out-of-bounds-fun-Proper-------------------------------------
is_out_of_bounds_fun_Proper <- function(p) any( 
  is.na(p), 
  is.infinite(p), 
  p < 0, 
  p[3] > 30 , 
  p[5] <= 1, 
  exp(8/p[5]) * p[4] > 1e200
)

## ----guess-initial-values-Proper, echo = TRUE, eval = do.calc------------
#  
#  pars_start_Proper <- c(0.05, 0.5, 1, 10, 5)
#  names(pars_start_Proper) <- c("hb", "ke", "kk", "mn", "beta")
#  
#  # to avoid conflicts with boundaries of L-BFGS-B, the minimum logposteriror value is limited to -1e16
#  optim_fun <- function(pars, guts_objects, isOutOfBoundsFun) {
#    return(
#      max(-1e16,logposterior(pars, guts_objects, isOutOfBoundsFun) )
#        )
#  }
#  
#  optim_result_Proper <- optim(pars_start_Proper, optim_fun, lower = c(1e-6, 1e-6, 1e-6, 1e-6, 1e-6), upper = c(1, 1, 30, 40, 20), method = "L-BFGS-B", control = list(trace = 0, fnscale = -1), guts_objects = guts_object, isOutOfBoundsFun = is_out_of_bounds_fun_Proper)

## ----save-initial-values-Proper, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(optim_result_Proper, file = file.path("..", "inst", "extdata", "vignetteGUTS-Proper-initialValues.Rdata"))

## ----load-initial-values-Proper, echo = FALSE, results = 'hide', eval = !do.calc----
load(system.file("extdata", "vignetteGUTS-Proper-initialValues.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----show-initial-values-Proper------------------------------------------
if (optim_result_Proper$convergence != 0) {
  warning("Optimizing initial values has not converged. Using non-optimized initial values.")
  optim_result_Proper$par <- pars_start_Proper
}

print(optim_result_Proper)

## ----run-MCMC-Proper, echo = TRUE, results = 'hide', eval = do.calc------
#  mcmc_pars_Proper <- optim_result_Proper$par
#  mcmc_sigma_Proper <- diag( (mcmc_pars_Proper/10)^2 + .Machine$double.eps )
#  mcmc_result_Proper <- MCMC(p = logposterior,
#    init = mcmc_pars_Proper, scale = mcmc_sigma_Proper, adapt = 20000, acc.rate = 0.4, n = 150000,
#    guts_objects = guts_object, isOutOfBoundsFun = is_out_of_bounds_fun_Proper
#  )
#  
#  #exclude burnin and thin by 20
#  mcmc_result_Proper$samples <- mcmc_result_Proper$samples[seq(50001, 150000, by = 20),]
#  mcmc_result_Proper$log.p <- mcmc_result_Proper$log.p[seq(50001, 150000, by = 20)]

## ----save-MCMC-results-Proper, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(mcmc_result_Proper, file = file.path("..", "inst", "extdata", "vignetteGUTS-Proper-MCMCresults.Rdata"))

## ----load-MCMC-results-Proper, echo = FALSE, results = 'hide', eval = !do.calc----
load(system.file("extdata", "vignetteGUTS-Proper-MCMCresults.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----display-MCMC-Proper-------------------------------------------------

if (all(is.finite(mcmc_result_Proper$log.p))) {
  par( mfrow = c(dim(mcmc_result_Proper$samples)[2] + 1, 2) , mar = c(5,4,1,0.5))
  plot(as.mcmc(cbind(mcmc_result_Proper$samples, LL = mcmc_result_Proper$log.p)), auto.layout = FALSE)
  par(op)
} else {
  par( mfrow = c(dim(mcmc_result_Proper$samples)[2], 2) , mar = c(5,4,1,0.5))
  plot(as.mcmc(mcmc_result_Proper$samples), auto.layout = FALSE)
  par(op)
}

panel.cor <- function(x, y, ...)
{
par(usr = c(0, 1, 0, 1))
txt <- as.character(format(cor(x, y), digits=2))
text(0.5, 0.5, txt, cex = max(0.1, 4 * abs(cor(x, y)) + 1))
}
pairs(mcmc_result_Proper$samples, upper.panel=panel.cor)

## ----evalMCMC-fun, echo = TRUE, results = 'hide'-------------------------
eval_MCMC <- function(sampMCMC, plot = TRUE) {
  bestFit <- sampMCMC$samples[which.max(sampMCMC$log.p),]
  qu <- apply(sampMCMC$samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
  if (plot) {
    plot(seq(dim(sampMCMC$samples)[2]), bestFit, pch = 20, cex = 2, ylim = range(qu), 
      xaxt = "n", xlab = "Model parameter", ylab = "Parameter value")
    arrows(x0 = seq(dim(sampMCMC$samples)[2]), y0 = qu[1,], y1 = qu[3,], angle = 90, length = 0.1, code = 3)
    axis(side = 1, at = seq(dim(sampMCMC$samples)[2]), dimnames(sampMCMC$samples)[[2]])
  }
  res <- rbind(bestFit, qu)
  rownames(res)[1] <- "best"
  return(res)
}

## ----evaluate-MCMC-SD-dat-Proper-----------------------------------------
eval_MCMC(mcmc_result_Proper)

## ----forecast-setup-GUTS-object------------------------------------------
guts_obj_forecast <-
  guts_setup(
    C = c(60, 40, 6, 0, 0, 60, 40, 6, 0, 0, 60, 40, 6, 0),
    Ct = c(0, 2.2, 4, 6, 9.9, 10, 12.2, 14, 16, 19.9, 20, 22.2, 24, 26),
    y = c(100, rep(0,26)),
    yt = seq(0,26),
    model = "Proper", dist = "loglogistic", N = 1000, M = 10000
  )

## ----forecast, echo = TRUE, results = 'hide', eval = do.calc-------------
#  mcmc_forecasts_paras <- mcmc_result_Proper$samples
#  
#  forec <- apply(mcmc_forecasts_paras, 1,
#          function(par) list(
#            guts_calc_survivalprobs(gobj = guts_obj_forecast, par = par),
#            guts_report_damage(gobj = guts_obj_forecast)$damage
#          )
#        )

## ----extract-forecast-results-Proper, echo = TRUE, results = 'hide', eval = do.calc----
#  # extract the damage matrix
#  damage <- sapply(forec, function(x) x[[2]])
#  survProb <- sapply(forec, function(x) x[[1]])
#  rm(forec)
#  
#  # analyse damage
#  damage.qu <- apply(damage, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
#  
#  
#  survProb.qu <- apply(survProb, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
#  
#  hazard <- (- apply(survProb, 2, diff, 1) ) / survProb[seq(2,dim(survProb)[1]),]
#  
#  hazard.qu <- apply(hazard, 1, function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
#  

## ----save-forecast-Proper, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(survProb.qu, damage.qu, hazard.qu, file = file.path("..", "inst", "extdata", "vignetteGUTS-Proper-forecast.Rdata"))

## ----load-forecast-Proper, echo = FALSE, results = 'hide', eval = !do.calc----
load(system.file("extdata", "vignetteGUTS-Proper-forecast.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----plot-forecast-------------------------------------------------------
par(mfrow = c(3,1), mar = c(0.5,4,0.5,0.5), oma = c(2.5,0,0,0), las = 1, cex = 1)
plot(guts_obj_forecast$Ct, guts_obj_forecast$C, 
  xlim = range(guts_obj_forecast$yt), 
  ylim = range(c(guts_obj_forecast$C, damage.qu)), 
  xlab = "", ylab = "exposure dose", 
  pch = 20, cex = 2, lwd = 2, col = "red",
  type = "b", xaxt = "n"
)
damage.ti <- seq(0, max(guts_obj_forecast$yt), 
  length.out = dim(damage.qu)[2]) 
lines(damage.ti, damage.qu[2,], type = "l", ylim = range(damage.qu), xlab = "time", ylab = "damage D", lwd = 2)
lines(damage.ti, damage.qu[1,], lty = 2, lwd = 2)
lines(damage.ti, damage.qu[3,], lty = 2, lwd = 2)
legend("topright", legend = c("external", "proj. internal", "  (damage)"), fill = c("red", "black", NA), border = c("red", "black", NA), cex = 0.8, bty = "n")

prob.ti <- seq(0, max(guts_obj_forecast$yt))
plot(prob.ti, survProb.qu[2,],
  ylim = range(survProb.qu),
  xlab = "", ylab = "proj. survival", xaxt = "n",
  pch = 20, cex = 2)
arrows(prob.ti[-1], survProb.qu[1,-1], prob.ti[-1], survProb.qu[3,-1], angle = 90, code = 3, length = 0.1, lwd = 2)

plot(prob.ti[-1] - 0.5, hazard.qu[2,],
  xlim = range(prob.ti), ylim = range(hazard.qu),
  xlab = "", ylab = "proj. daily hazard", pch = 20, cex = 2, 
  col = "blue")
arrows(prob.ti[-1] - 0.5, hazard.qu[1,], prob.ti[-1] - 0.5, hazard.qu[3,], angle = 90, code = 3, length = 0.1, col = "blue", lwd = 2)
mtext(text = "time", side = 1, line = 1.4, outer = TRUE)

