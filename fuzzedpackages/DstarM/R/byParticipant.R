#' #' Helper function to do DstarM analyses by participant.
#' #'
#' #' @description A helper function to do run DstarM functions per participant.
#' #' Supported functions are: \code{\link{estDstarM}}, \code{\link{estND}}, and
#' #' \code{\link{estObserved}}.
#' #'
#' #' @param data data.frame with a $id column. For other requirements, see
#' #' \code{\link{estDstarM}}.
#' #' @param resD list of decision models (all of class DstarM).
#' #' @param resND list of nondecision models (all of class DstarM).
#' #' @param argsEstDstarM Additional arguments for \code{\link{estDstarM}}.
#' #' @param argsEstND Additional arguments for \code{\link{estND}}.
#' #' @param argsEstObserved Additional arguments for \code{\link{estObserved}}.
#' #' @param uniqueArgs Logical, are the additional arguments specified using one list
#' #' for all participants or does the list of additional arguments contain lists with
#' #' participant specific arguments?
#' #'
#' #' @examples
#' #' \dontrun{
#' #' library(DstarM)
#' #' tt = seq(0, 5, .1)
#' #' pars = c(.8, 2, .5, .5, .5, # condition 1
#' #'          .8, 3, .5, .5, .5, # condition 2
#' #'          .8, 4, .5, .5, .5) # condition 3
#' #' pdfND = dbeta(tt, 10, 30)
#' #' # simulate data for 2 participants
#' #' n = 1.2e3
#' #' data = rbind(simData(n = n, pars = pars, tt = tt, pdfND = pdfND),
#' #'              simData(n = n, pars = pars, tt = tt, pdfND = pdfND))
#' #' # add participant column
#' #' data$pp = rep(1:2, each = n)
#' #' # define restriction matrix
#' #' restr = matrix(1:5, 5, 3)
#' #' restr[2, -1] = 6:7 # allow drift rates to differ
#' #' # fix variance parameters
#' #' fixed = matrix(c('sz1', .5, 'sv1', .5), 2, 2)
#' #' ## Run D*M analysis by pp
#' #' resDpp = byParticipant(data = data,
#' #'                        argsEstDstarM = list(tt = tt,
#' #'                                             restr = restr,
#' #'                                             fixed = fixed,
#' #'                                             Optim = list(parallelType = 1)))
#' #' coef(resDpp)
#' #' plot(resDpp)
#' #' resNDpp = byParticipant(resD = resDpp,
#' #'                         argsEstND = list(Optim = list(parallelType = 1)))
#' #' plot(resNDpp, xlim = 0:1)
#' #' lines(tt, pdfND, col = 3, type = 'b')
#' #' resObsPp = byParticipant(resD = resDpp,
#' #'                          resND = resNDpp)
#' #' plot(resObsPp)
#' #' }
#' #' @export
#' byParticipant = function(data, resD, resND,
#'                          argsEstDstarM = list(),
#'                          argsEstND = list(),
#'                          argsEstObserved = list(),
#'                          uniqueArgs = FALSE) {
#'   if (!missing(data) & missing(resD) & missing(resND)) {
#'     # estDstarM
#'     if (!('pp' %in% colnames(data))) {
#'       stop('Argument data must have a column indicating participant id.')
#'     }
#'     x = lapply(split(data, data$pp), list)
#'     what = 'estDstarM'
#'     args = argsEstDstarM
#'   } else if (missing(data) & !missing(resD) & missing(resND)) {
#'     # estND
#'     idx1 = which(names(resD) == 'byPp')
#'     if (!all(sapply(resD[-idx1], is.DstarM))) {
#'       stop('resD must be a list where every output is of class D*M.')
#'     }
#'     x = lapply(resD[-idx1], list)
#'     what = 'estND'
#'     args = argsEstND
#'
#'   } else if (missing(data) & !missing(resD) & !missing(resND)) {
#'     # estObserved
#'     idx1 = which(names(resD) == 'byPp')
#'     idx2 = which(names(resND) == 'byPp')
#'     if (!all(sapply(resD[-idx1], is.DstarM),
#'              sapply(resND[-idx2], is.DstarM))) {
#'       stop('resD must be a list where every output is of class D*M.')
#'     }
#'     x = Map(function(x, y) list(x, y),
#'             x = resD[-idx1],
#'             y = resND[-idx2])
#'     what = 'estObserved'
#'     args = argsEstObserved
#'
#'   } else {
#'     stop(sprintf('Unsupported input.'),
#'          call. = FALSE)
#'   }
#'   # function to loop over -- to be declared elsewhere
#'   fun = function(x, what, args) {
#'     args[names(formals(what))[seq_along(x)]] = x
#'     do.call(what, args)
#'   }
#'   if (uniqueArgs) {
#'     out = Map(fun, x = x, what = what, args = args)
#'   } else {
#'     out = lapply(x, fun, what = what, args = args)
#'   }
#'   out$byPp = what
#'   class(out) = 'DstarM'
#'   return(out)
#' }
#'
#'
