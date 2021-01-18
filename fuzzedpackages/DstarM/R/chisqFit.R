#' Calculate model fit
#'
#' @param resObserved either output from \code{\link{estObserved}} or a matrix containing custom densities to calculate the fitness for.
#' @param data A dataframe containing data.
#' @param DstarM Logical. Should the DstarM fit measure be calculated or the traditional fit measure?
#' @param tt time grid custom densities where calculated on. Should only be supplied if \code{resOberved} is a matrix containing custom densities
#' @param formula Optional formula argument, for when columns names in the data are different from those used to obtain the results.
#'
#' @details This function allows a user to manually calculate a chi-square goodness of fit measure for model densities.
#' This is useful for comparing a traditional analysis and a D*M analysis. For completion, this function can also calculate a
#' D*M fit measure. We do not recommend usage of the D*M measure. While the chi-square fit measure is
#' identical to the value of the optimizer when fitting, the DstarM fit measure is not equal to that of a DstarM analysis.
#' This is because this function calculates the DstarM fit measure on the complete distribution, not on the
#' model distributions, as is done during the optimization.
#'
#' @examples
#' tt = seq(0, 5, .1)
#' pars = c(.8, 2, .5, .5, .5, # condition 1
#'          .8, 3, .5, .5, .5,  # condition 2
#'          .8, 4, .5, .5, .5)  # condition 3
#' pdfND = dbeta(tt, 10, 30)
#'
#' # simulate data
#' allDat = simData(n = 3e3, pars = pars, tt = tt, pdfND = pdfND, return.pdf = TRUE)
#' truePdf = allDat$pdfUnnormalized
#' dat = allDat$dat
#' chisqFit(resObserved = truePdf, data = dat, tt = tt)
#' \dontrun{
#'# estimate it
#' define restriction matrix
#' restr = matrix(1:5, 5, 3)
#' restr[2, 2:3] = 6:7 # allow drift rates to differ
#' # fix parameters for speed up
#' fixed = matrix(c('z1', 'a1 / 2', 'sz1', .5, 'sv1', .5), 2, 3)
#' resD = estDstarM(data = dat, tt = tt, restr = restr, fixed = fixed,
#'                  Optim = list(parallelType = 1))
#' resN = estND(resD, Optim = list(parallelType = 1))
#'
#' resO = estObserved(resD, resN, data = dat)
#' resO$fit # proper fit
#'}

#' @export
chisqFit <- function(resObserved, data, DstarM = FALSE, tt = NULL, formula = NULL) {
  
  if (is.DstarM.fitObs(resObserved)) {
    
    if (is.null(tt)) {
      
      tt <- resObserved$tt  # time grid
      m <- resObserved$obs  # model implied densities
      
    } else {
      
      m <- getPdfs(resObserved$resDecision, tt)
      
    }
    
    formula <- resObserved[["resDecision"]][["formula"]]
    
  } else {
    # assuming custom densities
    
    m <- resObserved
    
  }
  
  by <- unique(zapsmall(diff(tt)))  # stepsize of time grid
  
  if (length(by) != 1) {
    stop("Time grid tt must be equally spaced and length(unique(zapsmall(diff(tt)))) == 1 must be TRUE.", 
      call. = FALSE)
  }
  
  # calculate data based densities
  if (!is.data.frame(data)) {
    stop(sprintf("Argument data should be dataframe. The supplied object has mode %s", 
      mode(data)))
  } else {
    
    if (!is.matrix(resObserved)) 
      formula <- resObserved$resDecision$formula
    
    data <- getData(formula, data)
    rtime <- data[["rtime"]]
    response <- data[["response"]]
    condition <- data[["condition"]]
    hasConditions <- data[["hasConditions"]]
    data <- data[["data"]]
    ncondition <- length(unique(data[[condition]]))  # get number of conditions
    
    # sanity check: do number of conditions in objects data and resObserved
    # match?
    if (is.DstarM.fitObs(resObserved) && !is.null(resObserved$ncondition) && 
      ncondition != resObserved$ncondition) {
      
      stop(sprintf("Number of conditions in resObserved (%d) does not match number of conditions in the data (%d)", 
        as.integer(resObserved$ncondition), as.integer(ncondition)))
      
      
    } else if (!is.DstarM.fitObs(resObserved) && 2 * ncondition != NCOL(resObserved)) {
      
      stop(sprintf("Number of conditions in resObserved (%d) does not match number of conditions in the data (%d)", 
        as.integer(NCOL(resObserved)), as.integer(2 * ncondition)))
      
    }
    
    # helper matrix
    mm <- matrix(0, ncondition * 2, ncondition)
    mm[1:dim(mm)[1L] + dim(mm)[1L] * rep(1:dim(mm)[2L] - 1, each = 2)] <- 1
    
    rt <- split(data[[rtime]], list(data[[response]], data[[condition]]))
    ncr <- lengths(rt)
    g <- getGhat(rt = rt, tt = tt, ncondition = ncondition, mm = mm, 
      by = by)
    
  }
  
  # calculate the chi square goodness of fit
  if (DstarM) {
    tmp <- 1:(NCOL(m) - 1)
    ii <- rep(tmp, times = rev(tmp))
    jj <- unlist(lapply(tmp, function(x, m) (x + 1):m, m = ncol(m)))
    out <- numeric(length(ii))
    for (l in 1:length(ii)) {
      a <- customConvolveO(g[, ii[l]], by * rev(m[, jj[l]]))[seq_along(tt)]
      b <- customConvolveO(g[, jj[l]], by * rev(m[, ii[l]]))[seq_along(tt)]
      out[i] <- chisq(tt = tt, a = a, b = b) * 100 * (ncr[ii[l]] + 
        ncr[jj[l]])/sum(ncr)
    }
    
  } else {
    out <- numeric(NCOL(m))
    for (i in 1:NCOL(m)) {
      out[i] <- chisq(tt = tt, a = m[, i], b = g[, i]) * 100 * ncr[i]/sum(ncr)
    }
  }
  
  return(list(sum = sum(out), chisq = out))
  
}


