## ----setup, include=FALSE-------------------------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  error = TRUE,
  comment = "#>"
)
op <- options(width = 100, digits = 4)

## ----load-pkg-data, eval=TRUE---------------------------------------------------------------------
library("fddm")
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0), ]

## ----log-likelihood, eval=TRUE--------------------------------------------------------------------
ll_fun <- function(pars, rt, resp, truth, err_tol) {
  v <- numeric(length(rt))

  # the truth is "upper" so use vu
  v[truth == "upper"] <- pars[[1]]
  # the truth is "lower" so use vl
  v[truth == "lower"] <- pars[[2]]

  dens <- dfddm(rt = rt, response = resp, a = pars[[3]], v = v, t0 = pars[[4]],
                w = pars[[5]], sv = pars[[6]], log = TRUE, err_tol = 1e-6)

  return( ifelse(any(!is.finite(dens)), 1e6, -sum(dens)) )
}

## ----prep-simp-data, eval=TRUE--------------------------------------------------------------------
onep <- med_dec[ med_dec[["id"]] == "2" & med_dec[["group"]] == "experienced", ]
onep[["resp"]] <- ifelse(onep[["response"]] == "blast", "upper", "lower")
onep[["truth"]] <- ifelse(onep[["classification"]] == "blast", "upper", "lower")
str(onep)

## ----show-simp-fit, eval=TRUE, warning=FALSE------------------------------------------------------
fit <- nlminb(c(0, 0, 1, 0, 0.5, 0), objective = ll_fun,
              rt = onep[["rt"]], resp = onep[["resp"]], truth = onep[["truth"]],
              # limits:   vu,   vl,   a,  t0, w,  sv
              lower = c(-Inf, -Inf, .01,   0, 0,   0),
              upper = c( Inf,  Inf, Inf, Inf, 1, Inf))
fit

## ----fitting-fun, eval=TRUE-----------------------------------------------------------------------
rt_fit <- function(data, id_idx = NULL, rt_idx = NULL, response_idx = NULL,
                   truth_idx = NULL, response_upper = NULL) {

  # Format data for fitting
  if (all(is.null(id_idx), is.null(rt_idx), is.null(response_idx),
      is.null(truth_idx), is.null(response_upper))) {
    df <- data # assume input data is already formatted
  } else {
    if(any(data[,rt_idx] < 0)) {
      stop("Input data contains negative response times; fit will not be run.")
    }
    if(any(is.na(data[,response_idx]))) {
      stop("Input data contains invalid responses (NA); fit will not be run.")
    }

    nr <- nrow(data)
    df <- data.frame(id = character(nr),
                     rt = double(nr),
                     response = character(nr),
                     truth = character(nr),
                     stringsAsFactors = FALSE)

    if (!is.null(id_idx)) { # relabel identification tags
      for (i in 1:length(id_idx)) {
        idi <- unique(data[,id_idx[i]])
        for (j in 1:length(idi)) {
          df[["id"]][data[,id_idx[i]] == idi[j]] <- paste(
            df[["id"]][data[,id_idx[i]] == idi[j]], idi[j], sep = " ")
        }
      }
      df[["id"]] <- trimws(df[["id"]], which = "left")
    }

    df[["rt"]] <- as.double(data[,rt_idx])

    df[["response"]] <- "lower"
    df[["response"]][data[,response_idx] == response_upper] <- "upper"

    df[["truth"]] <- "lower"
    df[["truth"]][data[,truth_idx] == response_upper] <- "upper"
  }

  # Preliminaries
  ids <- unique(df[["id"]])
  nids <- max(length(ids), 1) # if inds is null, there is only one individual
  ninit_vals <- 5

  # Initilize the output dataframe
  cnames <- c("ID", "Convergence", "Objective",
              "vu_fit", "vl_fit", "a_fit", "t0_fit", "w_fit", "sv_fit")
  out <- data.frame(matrix(ncol = length(cnames), nrow = nids))
  colnames(out) <- cnames
  temp <- data.frame(matrix(ncol = length(cnames)-1, nrow = ninit_vals))
  colnames(temp) <- cnames[-1]

  # Loop through each individual and starting values
  for (i in 1:nids) {
    out[["ID"]][i] <- ids[i]

    # extract data for id i
    dfi <- df[df[["id"]] == ids[i],]
    rti <- dfi[["rt"]]
    respi <- dfi[["response"]]
    truthi <- dfi[["truh"]]

    # starting value for t0 must be smaller than the smallest rt
    min_rti <- min(rti)

    # create initial values for this individual
    init_vals <- data.frame(vu = rnorm(n = ninit_vals, mean = 4, sd = 2),
                            vl = rnorm(n = ninit_vals, mean = -4, sd = 2),
                            a  = runif(n = ninit_vals, min = 0.5, max = 5),
                            t0 = runif(n = ninit_vals, min = 0, max = min_rti),
                            w  = runif(n = ninit_vals, min = 0, max = 1),
                            sv = runif(n = ninit_vals, min = 0, max = 5))

    # loop through all of the starting values
    for (j in 1:ninit_vals) {
      mres <- nlminb(init_vals[j,], ll_fun,
                     rt = rti, resp = respi, truth = truthi,
                     # limits:   vu,   vl,   a,  t0, w,  sv
                     lower = c(-Inf, -Inf, .01,   0, 0,   0),
                     upper = c( Inf,  Inf, Inf, Inf, 1, Inf))
      temp[["Convergence"]][j] <- mres[["convergence"]]
      temp[["Objective"]][j] <- mres[["objective"]]
      temp[j, -c(1, 2)] <- mres[["par"]]
    }

    # determine best fit for the individual
    min_idx <- which.min(temp[["Objective"]])
    out[i, -1] <- temp[min_idx,]
  }
  return(out)
}

## ----fitting-run, eval=TRUE, warning=FALSE--------------------------------------------------------
data(med_dec, package = "fddm")
med_dec <- med_dec[which(med_dec[["rt"]] >= 0),]
fit <- rt_fit(med_dec, id_idx = c(2,1), rt_idx = 8, response_idx = 7,
              truth_idx = 5, response_upper = "blast")
fit

## ----plot, eval=TRUE------------------------------------------------------------------------------
library("reshape2")
library("ggplot2")

fitp <- data.frame(fit[, c(1, 4, 5)]) # make a copy to manipulate for plotting
colnames(fitp)[-1] <- c("vu", "vl")

for (i in 1:length(unique(fitp[["ID"]]))) {
  first <- substr(fitp[["ID"]][i], 1, 1)
  if (first == "n") {
    fitp[["ID"]][i] <- "novice"
  } else if (first == "i") {
    fitp[["ID"]][i] <- "inexperienced"
  } else {
    fitp[["ID"]][i] <- "experienced"
  }
}

fitp <- melt(fitp, id.vars = "ID", measure.vars = c("vu", "vl"),
             variable.name = "vuvl", value.name = "estimate")

ggplot(fitp, aes(x = factor(ID, levels = c("novice", "inexperienced", "experienced")),
                 y = estimate,
                 color = factor(vuvl, levels = c("vu", "vl")))) +
  geom_point(alpha = 0.4, size = 4) +
  labs(title = "Parameter Estimates for vu and vl",
       x = "Experience Level", y = "Parameter Estimate",
       color = "Drift Rate") +
  theme_bw() +
  theme(panel.border = element_blank(),
        plot.title = element_text(size = 23),
        plot.subtitle = element_text(size = 16),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        axis.title.x = element_text(size = 20,
                                    margin = margin(10, 5, 5, 5, "pt")),
        axis.title.y = element_text(size = 20),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 16))

## ----session-info, collapse=TRUE------------------------------------------------------------------
sessionInfo()

## ----reset-options, include=FALSE---------------------------------------------
options(op)  # reset options

