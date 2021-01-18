########################################################################################################################
## Classe S4 mhmmdata
########################################################################################################################
###################################################################################
##' Constructor of \code{\linkS4class{mhmmdata}} class
##'
##'  
##' \describe{
##'   \item{nobs}{numeric. number of subjects}
##'   \item{yi}{list. each element corresponds to the sequences of a single subject.}
##'   \item{nbseq}{nbseq. number of sequences for each subject.}
##'   \item{nbtimeobs}{list. length of each sequence.}
##'   \item{tstart}{list. starting time of each sequence.}
##' }
##'
#' @examples
#'   getSlots("mhmmdata")
#'
#' @name mhmmdata-class
#' @rdname mhmmdata-class
#' @exportClass mhmmdata
setClass(
  Class = "mhmmdata",
  representation = representation(
    nobs="numeric",
    yi="list",
    nbseq="numeric",
    nbtimeobs="list",
    tstart = "list"
  ),
  prototype = prototype(
    nobs=numeric(),
    yi=list(),
    nbseq=numeric(),
    nbtimeobs=list(),
    tstart = list()
  )
)

mhmmData <- function(y){
  nobs <- length(y)
  yi <- list()
  tstart <- list()
  for (i in 1:nobs){
    beforenotna <- TRUE
    tstart[[i]] <- 1
    yi[[i]] <- list()
    block <- 1
    loc <- 1
    yi[[i]][[block]] <- numeric(0)
    for (t in 1:length(y[[i]])){
      if (!is.na(y[[i]][t])){
        yi[[i]][[block]][loc] <- y[[i]][t]
        loc <- loc + 1
        if (!beforenotna) tstart[[i]] <- c(tstart[[i]], t)
        beforenotna <- TRUE
      }else{
        if (beforenotna){
          block <- block + 1
          loc <- 1
          yi[[i]][[block]] <- numeric(0)
        }
        beforenotna <- FALSE
      }
    }
  }
  nbseq <- sapply(yi, function(u) length(u))
  nbtimeobs <- list()
  for (i in 1:nobs){
    nbtimeobs[[i]] <- sapply(yi[[i]], function(u) length(u))
  }
  new("mhmmdata", nobs=nobs, yi=yi, nbseq=nbseq, nbtimeobs=nbtimeobs, tstart=tstart)
}

