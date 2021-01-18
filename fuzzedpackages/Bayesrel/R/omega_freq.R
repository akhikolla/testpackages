# gives freq omega, and loadings and errors
#

omegaFreqData <- function(data, interval, omega.int.analytic, pairwise, n.boot = 1e3, callback = function(){},
                          parametric = F
                          ){
  p <- ncol(data)
  file <- lavOneFile(data)
  colnames(data) <- file$names

  lam_names <- paste("l", 1:p, sep = "")
  err_names <- paste("e", 1:p, sep = "")
  model <- paste0("f1 =~ ")
  loadings <- paste(paste(lam_names, "*", file$names, sep = ""),
                       collapse = " + ")
  factors <- "f1 ~~ 1*f1\n"
  errors <- paste(paste(file$names, " ~~ ", err_names, "*",
                           file$names, sep = ""), collapse = "\n")
  sum_loads <- paste("loading :=", paste(lam_names, collapse = " + "),
                      "\n")
  sum_errs <- paste("error :=", paste(err_names, collapse = " + "),
                    "\n")
  omega <- "omega := (loading^2) / ((loading^2) + error) \n"
  mod <- paste(model, loadings, "\n", factors, errors,
                 "\n", sum_loads, sum_errs, omega)

  if (pairwise) {
    fit <- fitmodel_mis(mod, data)
  } else {
    fit <- fitmodel(mod, data)
  }
  if (is.null(fit)) {
    return(list(omega = NA, fit.object = NULL))
  } else {
    params <- lavaan::parameterestimates(fit, level = interval)
    omega <- params$est[params$lhs=="omega"]
    if (omega.int.analytic) {
      om_low <- params$ci.lower[params$lhs=="omega"]
      om_up <- params$ci.upper[params$lhs=="omega"]
      om_obj <- NA
    } else {
      if (parametric) {
        bb <- lavaan::bootstrapLavaan(fit, type = "parametric", R = n.boot)
        if (dim(bb)[1] < 2) {
          om_low <- NA
          om_up <- NA
          om_obj <- NA
        } else {
          llow <- apply(bb[, 1:p], 2, quantile, prob = (1-interval)/2)
          elow <- apply(bb[, (p+1):(p*2)], 2, quantile, prob = (1-interval)/2)
          lup <- apply(bb[, 1:p], 2, quantile, prob = interval+(1-interval)/2)
          eup <- apply(bb[, (p+1):(p*2)], 2, quantile, prob = interval+(1-interval)/2)
          om_low <- omegaBasic(llow, elow)
          om_up <- omegaBasic(lup, eup)
          suml <- apply(bb[, 1:p], 1, sum)
          sume <- apply(bb[, (p+1):(p*2)], 1, sum)
          om_obj <- suml^2 / (suml^2 + sume)
        }

      } else {
        bb <- lavaan::bootstrapLavaan(fit, type = "nonparametric", R = n.boot)
        if (dim(bb)[1] < 2) {
          om_low <- NA
          om_up <- NA
          om_obj <- NA
        } else {
          llow <- apply(bb[, 1:p], 2, quantile, prob = (1-interval)/2)
          elow <- apply(bb[, (p+1):(p*2)], 2, quantile, prob = (1-interval)/2)
          lup <- apply(bb[, 1:p], 2, quantile, prob = interval+(1-interval)/2)
          eup <- apply(bb[, (p+1):(p*2)], 2, quantile, prob = interval+(1-interval)/2)
          om_low <- omegaBasic(llow, elow)
          om_up <- omegaBasic(lup, eup)
          suml <- apply(bb[, 1:p], 1, sum)
          sume <- apply(bb[, (p+1):(p*2)], 1, sum)
          om_obj <- suml^2 / (suml^2 + sume)
        }
      }
    }

    fit_tmp <- lavaan::fitMeasures(fit)
    indic <- c(fit_tmp["chisq"], fit_tmp["df"], fit_tmp["pvalue"],
               fit_tmp["rmsea"], fit_tmp["rmsea.ci.lower"], fit_tmp["rmsea.ci.upper"],
               fit_tmp["srmr"])
  }
  callback()
  return(list(omega = omega, omega_lower = om_low, omega_upper = om_up, indices = indic, fit.object = fit,
              omega_boot = om_obj))
}


fitmodel <- function(mod, data) {
  out <- tryCatch(
    {
      lavaan::cfa(mod, data, std.lv = T)
    },
    error = function(cond) {
      return(NULL)
    },
    warning = function(cond) {
      return(NULL)
    },
    finally = {}
  )
  return(out)
}

fitmodel_mis <- function(mod, data) {
  out <- tryCatch(
    {
      lavaan::cfa(mod, data, std.lv = T, missing = "ML")
    },
    error = function(cond) {
      return(NULL)
    },
    warning = function(cond) {
      return(NULL)
    },
    finally = {}
  )
  return(out)
}
