#' Class HMB
#'
#' Class \code{HMB} is the base class for the HMB-package
#'
#' @name HMB-class
#' @rdname HMB-class
#' @exportClass HMB
#' @seealso \code{\link{hmb}}, \code{\link{ghmb}}, \code{\link{tsmb}}, \code{\link{gtsmb}}
setClass(
  'HMB',
  slots = c(
    method = 'character',
    n = 'list',
    data = 'list',
    modelArgs = 'list',

    Alpha = 'matrix',
    Beta = 'matrix',
    Gamma = 'matrix',
    AlphaCov = 'matrix',
    BetaCov = 'matrix',
    mu = 'numeric',
    muVar = 'numeric',
    resids = 'list'
  )
)

#### Set validitiy ####
setValidity(
  'HMB',
  function(object) {
    errors = character()
    method = toupper(object@method)

    # Check if method is correct
    if (!(method %in% c(
      'HMB',
      'TSMB',
      'GHMB',
      'GTSMB'
    ))) {
      msg = "Method is not of correct type, either SMB or TSMB."
      errors = c(errors, msg)
    }

    if (length(errors) == 0) {
      return(TRUE)
    } else {
      return(errors)
    }
  }
)




#' Method getSpec
#'
#' Get model specifications of HMB-class object
#'
#' @name getSpec
#' @rdname getSpec-methods
#' @exportMethod getSpec
#' @param obj Object of class HMB
#' @return A list containing the estimated parameters, together with model arguments
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
#' getSpec(hmb_model)
setGeneric(
  name = "getSpec",
  def = function(obj) standardGeneric("getSpec")
)

#' @rdname getSpec-methods
setMethod(
  "getSpec",
  "HMB",
  definition = function(obj) {
    validObject(obj)

    retlist = list(
      Alpha = obj@Alpha,
      AlphaCov = obj@AlphaCov,
      Beta = obj@Beta,
      BetaCov = obj@BetaCov
    )

    if (length(obj@modelArgs) > 0) {
      retlist$modelArgs = obj@modelArgs
    }

    if (length(obj@resids) > 0) {
      retlist$resids = obj@resids
    }

    if (obj@method %in% c('TSMB', 'GTSMB')) {
      retlist$Gamma = obj@Gamma
    }

    return(retlist)
  }
)


## ' Method show
## '
## ' Display model outputs
## '
## ' @name show
## ' @rdname show-methods
## ' @exportMethod show
## ' @param object Object of class HMB
## ' @examples
## ' pop_U  = sample(nrow(HMB_data), 20000)
## ' pop_Sa = sample(pop_U, 5000)
## ' pop_S  = sample(pop_U, 300)
## '
## ' y_S    = HMB_data[pop_S, "GSV"]
## ' X_S    = HMB_data[pop_S, c("hMAX", "h80", "CRR", "pVeg")]
## ' X_Sa   = HMB_data[pop_Sa, c("hMAX", "h80", "CRR", "pVeg")]
## ' Z_Sa   = HMB_data[pop_Sa, c("B20", "B30", "B50")]
## ' Z_U    = HMB_data[pop_U, c("B20", "B30", "B50")]
## '
## ' hmb_model = hmb(y_S, X_S, X_Sa, Z_Sa, Z_U)
## ' show(hmb_model)
##
## setGeneric(
##   name = "show",
##   def = function(object) standardGeneric("show")
## )

#' Method show
#'
#' Display model outputs
#'
#' @name show
#' @rdname show-methods
#' @exportMethod show
#' @param object Object of class HMB
#' @aliases show,HMB-method
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
#' show(hmb_model)
setMethod(
  "show",
  "HMB",
  definition = function(object) {
    cat('Estimated population mean:', object@mu, '\n')
    cat('Estimated variance:', object@muVar, '\n')
  }
)


#' Method summary
#' 
#'Summary of HMB model
#'
#' @name summary
#' @rdname summary-methods
#' @exportMethod summary
#' @param obj Object of class HMB
#' @return Summary of HMB model.
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
#' S_Sa_map = matrix(pop_S, nrow = nrow(X_S), ncol = nrow(X_Sa))
#' S_Sa_map = t(apply(S_Sa_map, 1, function(x) {
#'   return(x == pop_Sa)
#' })) * 1
#'
#' hmb_model = hmb(y_S, X_S, X_Sa, Z_Sa, Z_U)
#' summary(hmb_model)
setGeneric(
  name = "summary",
  def = function(obj) standardGeneric("summary")
)

#' @rdname summary-methods
setMethod(
  "summary",
  "HMB",
  definition = function(obj) {
    validObject(obj)

    res = new('SummaryHMB')
    res@method = obj@method

    res@samples = matrix(
      c(obj@n$S, obj@n$Sa, obj@n$U),
      1L, 3L, dimnames = list('')
    )
    colnames(res@samples) = c('S', 'Sa', 'U')

    res@estimation = matrix(
      c(obj@mu, obj@muVar, obj@mu + qnorm(c(.025, .975), 0, sqrt(obj@muVar))),
      1L, 4L, dimnames = list('')
    )
    colnames(res@estimation) = c('Mean', 'Variance', 'Lower 95 % conf', 'Upper 95 % conf')

    betastd = sqrt(diag(obj@BetaCov))
    res@betacoef = cbind(
      obj@Beta, betastd, obj@Beta / betastd,
      2 * pnorm(abs(obj@Beta), 0, betastd, lower.tail = FALSE)
    )

    dimnames(res@betacoef) = list(
      c('(Intercept)', colnames(obj@data$X_S)[-1]),
      c('Estimate', 'Std. error', 'Z value', 'Pr(>|Z|)')
    )

    alphastd = sqrt(diag(obj@AlphaCov))
    res@alphacoef = cbind(
      obj@Alpha, alphastd, obj@Alpha / alphastd,
      2 * pnorm(abs(obj@Alpha), 0, alphastd, lower.tail = FALSE)
    )
    dimnames(res@alphacoef) = list(
      c('(Intercept)', colnames(obj@data$Z_Sa)[-1]),
      c('Estimate', 'Std. error', 'Z value', 'Pr(>|Z|)')
    )

    if (obj@method %in% c('TSMB', 'GTSMB')) {
      res@gammacoef$gamma = obj@Gamma
      rownames(res@gammacoef$gamma) = colnames(obj@data$Z_Sa)
      colnames(res@gammacoef$gamma) = colnames(obj@data$X_Sa)
      rownames(res@gammacoef$gamma)[1] = '(Intercept)'
      colnames(res@gammacoef$gamma)[1] = '(Intercept)'
    }

    return(res)
  }
)
