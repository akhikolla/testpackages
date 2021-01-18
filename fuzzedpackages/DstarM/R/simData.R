#' Simulate data from a given density function via multinomial sampling
#'
#' @param n Number of observations to be sampled
#' @param pars Parameter values for the density function to be evaluated with. \code{length(pars)} must be a multiple of npars.
#' @param tt time grid on which the density function will be evaluated. Responses not in this time grid cannot appear.
#' @param fun.density Density function to use.
#' @param args.density Additional arguments to be passed to \code{fun.density}, aside from \code{tt}, \code{pars}, and a boundary argument ('upper' or 'lower')
#' @param npars Number of parameters \code{fun.density} must be evaluated with. If \code{length(pars) > npars} each \code{npars} values in \code{pars} will be seen as the parameter values of a condition.
#' @param return.pdf Logical, if TRUE \code{genData} returns a list containing the probability density function used and the data, if FALSE \code{genData} returns a dataframe with simulated data.
#' @param pdfND either a vector of length tt specifying the nondecision density for all condition-response pairs,
#' or a matrix where columns corresponds to the nondecision densities of condition-response pairs. Supplying \code{NULL} implies no nondecision distribution.
#' @param normalizePdfs Logical, should the pdf of the nondecision distribution be normalized?
#'
#' @details Simulate data via multinomial sampling.
#' The response options to sample from should be provided in \code{tt}.
#' The number of conditions is defined as \code{length(pars) / npars}.
#'
#' @return A sorted dataframe where rows represent trials. It contains: a column named rt
#' containing reaction times in seconds, a column named response containing either
#' response option lower or upper, and a column named condition indicating which
#' condition a trials belongs to. If \code{return.pdf} is TRUE it returns a list where the
#' first element is the sorted dataframe, the second through the fifth elements are lists
#' that contain densities used for simulating data.
#'
#'
#' @examples
#' tt = seq(0, 5, .01)
#' pdfND = dbeta(tt, 10, 30)
#' n = 100
#' pars = c(1, 2, .5, .5, .5)
#' dat = simData(n, pars, tt, pdfND)
#' head(dat)

#' @export
# generate data via base rmultinom and drd; returns dataframe
simData <- function(n, pars, tt, pdfND, fun.density = Voss.density, args.density = list(prec = 3),
  npars = 5, return.pdf = FALSE, normalizePdfs = TRUE) {
  by <- unique(zapsmall(diff(tt)))
  if (length(by) != 1) {
    stop("Time grid tt must be equally spaced and length(unique(zapsmall(diff(tt)))) == 1 must be TRUE.",
      call. = FALSE)
  }
  if (length(pars)/npars != ceiling(length(pars)/npars))
    stop("Wrong length pars")
  ncondition <- length(pars)/npars
  if (is.null(dim(pdfND)) & length(pdfND) == length(tt)) {
    pdfND <- matrix(pdfND, length(tt), 2 * ncondition)
  } else if (any(dim(pdfND) != c(length(tt), 2 * ncondition))) {
    stop("pdfND must either be a vector of length tt, or a matrix where every columns represents a nondecision distribution with dimensions length(tt) x (ncondition * 2). Every column represents the nondecision distribution of a condition-response pair, hence the *2.",
      call. = FALSE)
  }
  hasND <- !is.null(pdfND)
  # normalize all ND
  if (hasND && normalizePdfs) {
    pdfND <- pdfND %*% (diag(dim(pdfND)[2L])/apply(pdfND, 2, simpson, x = tt))
  } #else if (!hasND && npars == 5L) {
  #  warning("If fun.density")
  #}

  if (n/ncondition != ceiling(n/ncondition))
    stop("n is not a multiple of the number of conditions")
  pars.list <- lapply(1:ncondition, function(x, pars, npars) pars[(1 +
    npars * (x - 1)):(npars * x)], pars = pars, npars = npars)
  # mm is a helper matrix
  mm <- matrix(0, ncondition * 2, ncondition)
  mm[1:dim(mm)[1L] + dim(mm)[1L] * rep(1:dim(mm)[2L] - 1, each = 2)] <- 1
  # get all pdfs
  if (!is.null(args.density$DstarM)) {
    DstarM <- args.density$DstarM
    args.density <- args.density[names(args.density) != "DstarM"]
  } else if (is.null(args.density$DstarM) && npars == 7L && identical(fun.density, Voss.density)) {
    message("NOTE: args.density$DstarM has been set to FALSE. Set args.density$DstarM explicitly to avoid this behavior.")
    DstarM <- FALSE
  } else {
    DstarM <- TRUE
  }
  pdfD <- getPdf(pars.list = pars.list, tt = tt, DstarM = DstarM,
                 mm = mm, oscPdf = FALSE, fun.density = fun.density, args.density = args.density)
  pdf <- pdfD
  # convolve all pdfs with the non-decision distribution
  if (hasND) {
    for (i in 1:dim(pdfD)[2L]) {
      pdf[, i] <- customConvolveO(pdf[, i], by * rev(pdfND[, i]))[seq_along(tt)]
    }
  }
  # omit negative numerical artefacts
  pdf[pdf < 0] <- 0

  pdftot <- pdf

  dim(pdf) <- c(length(tt) * 2, ncondition)
  pdf <- pdf %*% (diag(ncondition)/colSums(pdf))  # normalize so the colSums are 1
  freq <- 0 * pdf
  rts <- rep.int(0, n)
  condition <- rts
  response <- rts
  stopifnot(n/ncondition == round(n/ncondition))
  for (i in 1:dim(pdf)[2L]) {
    freq[, i] <- rowSums(stats::rmultinom(1, n/ncondition, pdf[, i]))
    for (j in 1:2) {
      vals <- rep(tt, freq[(1 + (j - 1) * length(tt)):(j * length(tt)),
        i])
      ind <- which(rts == 0)[1]
      if (length(vals)) {
        # if there are no observations nothing needs to be filled in
        rts[ind:(ind + length(vals) - 1)] <- vals
        response[ind:(ind + length(vals) - 1)] <- ifelse(j == 1,
          "lower", "upper")
      }
    }
    condition[which(rts != 0 & condition == 0)] <- i
  }
  response <- factor(response, levels = c("lower", "upper"))
  dat <- data.frame(rt = rts, response = response, condition = condition)
  if (return.pdf) {
    colnames(pdftot) <- colnames(pdfD) <- colnames(pdfD) <- paste(rep(1:ncondition,
      each = 2), c("lower", "upper"))
    cor <- apply(pdftot, 2, simpson, x = tt)
    pdfNormalized <- pdftot %*% (diag(dim(pdftot)[2L])/cor)
    return(list(dat = dat, pdfNormalized = pdfNormalized, pdfUnnormalized = pdftot,
      pdfDecision = pdfD, pdfSimulate = pdf))
  } else {
    return(dat)
  }
}

