#' Class SummaryHMB
#'
#' Class \code{SummaryHMB} defines summary information for HMB object.
#'
#' @name SummaryHMB-class
#' @rdname SummaryHMB-class
#' @exportClass SummaryHMB
setClass(
  'SummaryHMB',
  slots = c(
    method = 'character',
    estimation = 'matrix',
    samples = 'matrix',
    betacoef = 'matrix',
    alphacoef = 'matrix',
    gammacoef = 'list'
  )
)

# #' Method show
# #'
# #' Display model outputs
# #'
# #' @name show
# #' @rdname show-methods
# #' @exportMethod show
# #' @param obj Object of class SummaryHMB
# #' @examples
# #' pop_U  = sample(nrow(HMB_data), 20000)
# #' pop_Sa = sample(pop_U, 5000)
# #' pop_S  = sample(pop_U, 300)
# #'
# #' y_S    = HMB_data[pop_S, "GSV"]
# #' X_S    = HMB_data[pop_S, c("hMAX", "h80", "CRR", "pVeg")]
# #' X_Sa   = HMB_data[pop_Sa, c("hMAX", "h80", "CRR", "pVeg")]
# #' Z_Sa   = HMB_data[pop_Sa, c("B20", "B30", "B50")]
# #' Z_U    = HMB_data[pop_U, c("B20", "B30", "B50")]
# #'
# #' hmb_model = hmb(y_S, X_S, X_Sa, Z_Sa, Z_U)
# #' show(hmb_model)
# setGeneric(
#   name = "show",
#   def = function(obj) standardGeneric("show")
# )

#' Method show
#'
#' Display model summary properties
#'
#' @name show
#' @rdname show-methods
#' @exportMethod show
#' @aliases show,SummaryHMB-method
#' @examples
#' pop_U  = sample(nrow(HMB_data), 20000)
#' pop_Sa = sample(pop_U, 5000)
#' pop_S  = sample(pop_U, 300)
#'
#' y_S    = HMB_data[pop_S, "GSV"]
#' X_S    = HMB_data[pop_S, c("hMAX", "h80", "CRR", "pVeg")]
#' X_Sa   = HMB_data[pop_Sa, c("hMAX", "h80", "CRR", "pVeg")]
#' Z_Sa   = HMB_data[pop_Sa, c("B20", "B30", "B50")]
#' Z_U    = HMB_data[pop_U, c("B20", "B30", "B50")]
#'
#' hmb_model = hmb(y_S, X_S, X_Sa, Z_Sa, Z_U)
#' show(summary(hmb_model))
setMethod(
  'show',
  'SummaryHMB',
  definition = function(object) {
    digits = max(3L, getOption('digits') - 3L)
    cat('Summary for ', object@method, '\n', sep = '')

    cat('\nEstimate:\n', sep = '')
    print(object@estimation, digits = digits)

    cat('\nSample sizes:\n', sep = '')
    print(object@samples, digits = digits)

    cat('\nBeta coefficients:\n', sep = '')
    printCoefmat(object@betacoef, digits = digits, signif.stars = TRUE)

    cat('\nAlpha coefficients:\n', sep = '')
    printCoefmat(object@alphacoef, digits = digits, signif.stars = TRUE)

    if (object@method %in% c('TSMB', 'GTSMB')) {
      cat('\nGamma coefficients:\n', sep='')
      print(object@gammacoef$gamma, digits = digits)
    }
  }
)

