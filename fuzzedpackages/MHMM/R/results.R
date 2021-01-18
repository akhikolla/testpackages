########################################################################################################################
## Classe S4 mhmmresults
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{mhmmresults}} class
##'
##'  
##' \describe{
##'   \item{param}{\linkS4class{mhmmparam}. MLE of the model parameters}
##'   \item{data}{\linkS4class{mhmmdata}. Data}
##'   \item{partitions}{list. Elements of the latent variables (partition among subjects and latent states)}
##'   \item{probabilities}{list. Posterior probabilities of the latent variables.}
##'   \item{meantimesperstates}{matrix. Mean time spent by each subject (row) into each activity level (column).}
##'   \item{meanvalueperstates}{numeric. Summary statistics of the states.}
##'   \item{loglike}{numeric. Loglikelihood.}
##'   \item{bic}{numeric. BIC.}
##' }
##'
#' @examples
#'   getSlots("mhmmresults")
#'
#' @name mhmmresults-class
#' @rdname mhmmresults-class
#' @exportClass mhmmresults
setClass(
  Class = "mhmmresults",
  representation = representation(
    param="mhmmparam",
    data="mhmmdata",
    partitions="list",
    probabilities="list",
    meantimesperstates = "matrix",
    meanvalueperstates = "numeric",
    loglike="numeric",
    bic="numeric"
  ),
  prototype = prototype(
    param=new("mhmmparam"),
    data=new("mhmmdata"),
    partitions=list(),
    probabilities=list(),
    meantimesperstates = matrix(),
    meanvalueperstates = numeric(),
    loglike=numeric(),
    bic=numeric()
  )
)
