#' Calculate Mean of the nondecision distribution.
#'
#' @param res An object of class D*M.
#' @param data The data object used to create \code{res}.
#' @param formula Optional formula argument, for when columns names in the data are different from those used to obtain the results.
#'
#' @return A vector containing estimates for the mean of the nondecision densities.

#' @details The object \code{res} can either be output from \code{estDstarM} or output from \code{estND}.
#' If the former is supplied it is also necessary to supply the data used for the estimation.
#' The mean will then be estimated by subtracting the mean of the model densities from the mean of the data density.
#' If the latter is supplied than this is not required; the mean will be calculated by
#' integrating the nondecision distribution.

# calculate Ter according to splits used in analyses; returns scalar |
# groups
#' @export
getTer <- function(res, data, formula = NULL) {

  if (!(is.DstarM.fitD(res) || is.DstarM.fitND(res)))
    stop("res should be output from either estDstarM or estND.")

  # nondecision output is easy
  if (is.DstarM.fitND(res)) {
    return(apply(res$r.hat, 2, nth.momentS, x = res$tt))
  }
  if (dim(data)[1L] != res$n) {
    warning(sprintf("Number of observations used in analysis (%g) does not match number of observations in data provided (%g).",
      res$n, dim(data)[1L]), call. = FALSE, immediate. = TRUE)
  }

  if (is.null(res$splits) & !is.null(res$split))
    res$splits <- res$split  # backward compatibility

  data <- getData(res[["formula"]], data)
  rtime <- data[["rtime"]]
  response <- data[["response"]]
  condition <- data[["condition"]]
  hasConditions <- data[["hasConditions"]]
  data <- data[["data"]]

  ncondition <- res$ncondition
  splits <- res$splits
  m <- res$modelDist

  group <- groups(ncondition, splits)
  mm2 <- matrix(0, 2 * ncondition, dim(group)[2L])
  for (i in 1:dim(group)[2L]) {
    mm2[group[, i], i] <- 1
  }

  m <- m %*% mm2
  m <- m %*% (diag(dim(m)[2L])/(colSums(mm2)/2))
  uniq <- unique(data[[condition]])
  group <- groups(ncondition, splits, TRUE)
  for (i in 1:length(group)) {
    group[i] <- uniq[i]
  }
  muDat <- rep.int(0, dim(group)[2L])
  for (i in dim(group)[2L]) {
    muDat[i] <- mean(data[[rtime]][data[[condition]] %in% group[, i]])
  }

  muMod <- apply(m, 2, nth.momentS, x = res$tt)
  return(muDat - muMod)
}

