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

## ----prepare-xlsx-table-local, include = FALSE, eval = TRUE--------------
JA_data_file_name <- system.file("extdata", "Data_for_GUTS_software_ring_test_A_v05.xlsx", package = "GUTS", mustWork = TRUE)

## ----setup-A-par---------------------------------------------------------
par_A <- data.frame(symbols = c("hb", "ke", "kk", "mn", "beta"), JAsymbols = c("hb", "kd", "kk", "mw", "beta"), SD = c(0.01, 0.8, 0.6, 3, NA), IT = c(0.02, 0.8, NA, 5, 5.3))
par_A

## ----A-SD-read-all-comments, message = FALSE-----------------------------
library(xlsx)
read.xlsx(
  file = paste0(JA_data_file_name),
  sheetName = "Data A",
  rowIndex = c(1:11),
  colIndex = seq(which(LETTERS == "A"), which(LETTERS == "G")),
  header = TRUE
  )

## ----A-SD-read-----------------------------------------------------------
data_A_SD <- read.xlsx(
  file = paste0(JA_data_file_name),
  sheetName = "Data A",
  rowIndex = c(4:11),
  colIndex = seq(which(LETTERS == "B"), which(LETTERS == "G")),
  header = FALSE
  )
con_A_SD <- as.numeric(data_A_SD[1,])
data_A_SD <- data_A_SD[-1,] 
#the first row contains the concentrations
# all subsequent rows: number of alive organisms at specific day.

day_A_SD <-
  as.numeric(
    t(
      read.xlsx(
        file = paste0(JA_data_file_name),
        sheetName = "Data A",
        rowIndex = c(5:11),
        colIndex = seq(which(LETTERS == "A")),
        header = FALSE
      )
    )
  )

# name the data.frame
names(data_A_SD) <- paste0("c", con_A_SD)
rownames(data_A_SD) <- paste0("d", day_A_SD)

## ----setup-GUTS-OBJ-A_SD-------------------------------------------------
GUTS_A_SD <- list( 
  C0 = guts_setup(
    C = rep_len(con_A_SD[1], length(day_A_SD)), Ct = day_A_SD,
    y = data_A_SD$c0, yt = day_A_SD,
    model = "SD"
    ),
  C2 = guts_setup(
    C = rep_len(con_A_SD[2], length(day_A_SD)), Ct = day_A_SD,
    y = data_A_SD$c2, yt = day_A_SD,
    model = "SD"
    ),
  C4 = guts_setup(
    C = rep_len(con_A_SD[3], length(day_A_SD)), Ct = day_A_SD,
    y = data_A_SD$c4, yt = day_A_SD,
    model = "SD"
    ),
  C6 = guts_setup(
    C = rep_len(con_A_SD[4], length(day_A_SD)), Ct = day_A_SD,
    y = data_A_SD$c6, yt = day_A_SD,
    model = "SD"
    ),
  C8 = guts_setup(
    C = rep_len(con_A_SD[5], length(day_A_SD)), Ct = day_A_SD,
    y = data_A_SD$c8, yt = day_A_SD,
    model = "SD"
    ),
  C16 = guts_setup(
    C = rep_len(con_A_SD[6], length(day_A_SD)), Ct = day_A_SD,
    y = data_A_SD$c16, yt = day_A_SD,
    model = "SD"
    )
)

## ----load optimization routines------------------------------------------
library('adaptMCMC') # Function `MCMC()`, Monte Carlo Markov Chain.

## ----define-log-posterior------------------------------------------------
logposterior <- function( pars, guts_objects, 
  isOutOfBoundsFun = function(p) any( is.na(p), is.infinite(p) )  ) {
	if ( isOutOfBoundsFun(pars) ) return(-Inf)
  return(
	  sum(sapply( guts_objects, function(obj) guts_calc_loglikelihood(obj, pars) ))
  )
}

## ----define out of bounds-fun-SD-----------------------------------------
is_out_of_bounds_fun_SD <- function(p) any( is.na(p), is.infinite(p), p < 0, p["kk"] > 30 )

## ----run-MCMC-SD, echo = TRUE, results = 'hide', eval = do.calc----------
#  pars_start_SD <- rep_len (0.5, 4)
#  names(pars_start_SD) <- par_A$JAsymbols[-which(is.na(par_A$SD))]
#  
#  mcmc_result_SD <- MCMC(p = logposterior,
#    init = pars_start_SD, adapt = 5000, acc.rate = 0.4, n = 150000,
#    guts_objects = GUTS_A_SD,
#    isOutOfBoundsFun = is_out_of_bounds_fun_SD
#  )
#  
#  #exclude burnin and thin by 3
#  mcmc_result_SD$samples <- mcmc_result_SD$samples[seq(50001, 150000, by = 20),]
#  mcmc_result_SD$log.p <- mcmc_result_SD$log.p[seq(50001, 150000, by = 20)]

## ----save-MCMC-results-SD, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(mcmc_result_SD, file = file.path("..", "inst", "extdata", "vignetteGUTS-ringTest-SD-MCMCresults.Rdata"))

## ----load-MCMC-results-SD, echo = FALSE, results = 'hide', eval = !do.calc----
load(system.file("extdata", "vignetteGUTS-ringTest-SD-MCMCresults.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----display-MCMC-SD-----------------------------------------------------

if (all(is.finite(mcmc_result_SD$log.p))) {
  par( mfrow = c(dim(mcmc_result_SD$samples)[2] + 1, 2) , mar = c(5,4,1,0.5))
  plot(as.mcmc(cbind(mcmc_result_SD$samples, LL = mcmc_result_SD$log.p)), auto.layout = FALSE)
  par(op)
} else {
  par( mfrow = c(dim(mcmc_result_SD$samples)[2], 2) , mar = c(5,4,1,0.5))
  plot(as.mcmc(mcmc_result_SD$samples), auto.layout = FALSE)
  par(op)
}

## ----evalMCMC-fun, echo = TRUE, results = 'hide'-------------------------
eval_MCMC <- function(sampMCMC, expectedVal = NULL, plot = TRUE) {
  bestFit <- sampMCMC$samples[which.max(sampMCMC$log.p),]
  qu <- apply(sampMCMC$samples, 2, quantile, probs = c(0.025, 0.5, 0.975))
  if (plot) {
    if(is.null(expectedVal)) expectedVal <- rep(NA, dim(sampMCMC$samples)[2])
    plot(seq(dim(sampMCMC$samples)[2]), expectedVal, pch = 20, col = "darkgrey", cex = 2, ylim = range(qu), 
      xaxt = "n", xlab = "Model parameter", ylab = "Parameter value")
    arrows(x0 = seq(dim(sampMCMC$samples)[2]), y0 = qu[1,], y1 = qu[3,], angle = 90, length = 0.1, code = 3)
    points(x = seq(dim(sampMCMC$samples)[2]), y = bestFit, pch = "-", cex = 4)
    axis(side = 1, at = seq(dim(sampMCMC$samples)[2]), dimnames(sampMCMC$samples)[[2]])
  }
  res <- rbind(bestFit, qu)
  rownames(res)[1] <- "best"
  if (!all(is.na(expectedVal))) {
    res <- rbind(res, expectedVal)
    rownames(res)[dim(res)[1]] <- "expect"
  }
  return(res)
}

## ----evaluate-MCMC-SD----------------------------------------------------
eval_MCMC(mcmc_result_SD, expectedVal = par_A$SD[-which(is.na(par_A$SD))])

## ----A-IT-read-----------------------------------------------------------
data_A_IT <- read.xlsx(
  file = paste0(JA_data_file_name),
  sheetName = "Data A",
  rowIndex = c(17:24),
  colIndex = seq(which(LETTERS == "B"), which(LETTERS == "G")),
  header = FALSE
  )
con_A_IT <- as.numeric(data_A_IT[1,])
data_A_IT <- data_A_IT[-1,] 
#the first row contains the concentrations
# all subsequent rows: number of alive organisms at specific day.

day_A_IT <-
  as.numeric(
    t(
      read.xlsx(
        file = paste0(JA_data_file_name),
        sheetName = "Data A",
        rowIndex = c(18:24),
        colIndex = seq(which(LETTERS == "A")),
        header = FALSE
      )
    )
  )

# name the data.frame
names(data_A_IT) <- paste0("c", con_A_IT)
rownames(data_A_IT) <- paste0("d", day_A_IT)

## ----setup-GUTS-OBJ-A-IT-------------------------------------------------
GUTS_A_IT <- lapply(seq(length(con_A_IT)), 
  function(i, dat, days, con) guts_setup(
    C = rep_len(con[i], length(days)), Ct = days,
    y = dat[,i], yt = days,
    model = "IT", dist = "loglogistic"
  ), dat = data_A_IT, days = day_A_IT, con = con_A_IT
) 
names(GUTS_A_IT) <- paste0("c", con_A_IT)

## ----define out of bounds-fun-IT-----------------------------------------
is_out_of_bounds_fun_IT <- function(p) any( is.na(p), is.infinite(p), p < 0, p[4] <= 1, exp(8/p[4]) * p[3] > 1e200)

## ----run-MCMC-IT, echo = TRUE, results = 'hide', eval = do.calc----------
#  pars_start_IT <- rep_len(0.5, 4)
#  names(pars_start_IT) <- par_A$JAsymbols[-which(is.na(par_A$IT))]
#  mcmc_result_IT <- MCMC(p = logposterior,
#    init = pars_start_IT, adapt = 5000, acc.rate = 0.4, n = 150000,
#    guts_objects = GUTS_A_IT, isOutOfBoundsFun = is_out_of_bounds_fun_IT
#  )
#  
#  #exclude burnin and thin by 3
#  mcmc_result_IT$samples <- mcmc_result_IT$samples[seq(50001, 150000, by = 20), ]
#  mcmc_result_IT$log.p <- mcmc_result_IT$log.p[seq(50001, 150000, by = 20)]

## ----save-MCMC-results-IT, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(mcmc_result_IT, file = file.path("..", "inst", "extdata", "vignetteGUTS-ringTest-IT-MCMCresults.Rdata"))

## ----load-MCMC-results-IT, echo = FALSE, results = 'hide', eval = !do.calc----
load(system.file("extdata", "vignetteGUTS-ringTest-IT-MCMCresults.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----display-MCMC-IT-----------------------------------------------------
if (all(is.finite(mcmc_result_IT$log.p))) {
  par( mfrow = c(dim(mcmc_result_IT$samples)[2] + 1, 2) , mar = c(5,4,1,0.5))
  plot(as.mcmc(cbind(mcmc_result_IT$samples, LL = mcmc_result_IT$log.p)), auto.layout = FALSE)
  par(op)
} else {
  par( mfrow = c(dim(mcmc_result_IT$samples)[2], 2) , mar = c(5,4,1,0.5))
  plot(as.mcmc(mcmc_result_IT$samples), auto.layout = FALSE)
  par(op)
}

## ----evaluate-MCMC-IT----------------------------------------------------
eval_MCMC(mcmc_result_IT, expectedVal = par_A$IT[-which(is.na(par_A$IT))])

## ----forecast-setup-GUTS-object------------------------------------------
conc <- seq(0, 16)
guts_obj_forecast <- lapply(conc,
  function(concentration) guts_setup(
    C = rep(concentration, 10),
    Ct = seq(0,9),
    y = c(100, rep(0,9)),
    yt = seq(0,9),
    model = "IT", dist = "loglogistic", N = 1000
  )
)

## ----forecast-paras------------------------------------------------------
mcmc_forecasts_paras <- mcmc_result_IT$samples
mcmc_forecasts_paras[,1] <- 0

## ----forecast, echo = TRUE, results = 'hide', eval = do.calc-------------
#  forec <- lapply(guts_obj_forecast,
#    function(gobj, mcmc_res)
#      rbind(
#        rep(gobj$C[1], dim(mcmc_forecasts_paras)[1]),
#        apply(mcmc_res, 1,
#          function(par) guts_calc_survivalprobs(gobj = gobj, par, external_dist = NULL)
#        )
#      ),
#    mcmc_res = mcmc_forecasts_paras
#  )
#  
#  forec <- do.call("cbind", forec)
#  forec <- as.data.frame(t(forec))
#  names(forec) <- c("conc", paste0("day", guts_obj_forecast[[1]]$yt))

## ----save-forecast, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(forec, file = file.path("..", "inst", "extdata", "vignetteGUTS-ringTest-forecast.Rdata"))

## ----load-forecast, echo = FALSE, results = 'hide', eval = !do.calc------
load(system.file("extdata", "vignetteGUTS-ringTest-forecast.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----plot-forecast-------------------------------------------------------
par(mfrow = c(3,3), mar = c(5,4,3, 0.5))
invisible(
  sapply(seq(9),
    function(day)
      plot(as.factor(forec$conc), forec[,paste0("day",day)], 
        ylim = c(0,1),
        xlab = "concentration (micromol/l)", ylab = "probability of survival", 
        main = paste("day",day))
  )
)
par(op)

## ----estimate-4d-LC50, echo = TRUE, results = 'hide', eval = do.calc-----
#  library("drc")
#  logLC50s <- sapply(seq_len(dim(mcmc_forecasts_paras)[1]),
#    function(indParaset, forDat, concentrations) {
#      dat <- forDat[dim(mcmc_forecasts_paras)[1] * (seq_along(concentrations) - 1) + indParaset,"day4"]
#      return(
#        coefficients(
#          drm(data.frame(dat, concentrations), fct = LL2.3(names = c("Slope", "upper", "logLC50")))
#        )[3]
#      )
#    },
#    forDat <- forec,
#    concentrations = conc
#  )

## ----save-LC50, echo = FALSE, results = 'hide', eval = do.calc & do.save----
#  save(logLC50s, file = file.path("..", "inst", "extdata", "vignetteGUTS-ringTest-logLC50.Rdata"))

## ----load-LC50, echo = FALSE, results = 'hide', eval = !do.calc----------
load(system.file("extdata", "vignetteGUTS-ringTest-logLC50.Rdata", 
  package = "GUTS", mustWork = TRUE)
)

## ----calc-LC50-stats-----------------------------------------------------
LC50 <- quantile(exp(logLC50s), c(0.025, 0.5, 0.975))

