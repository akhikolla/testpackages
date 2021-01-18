###Function to get posteior adjacency values for each location
#'
#' PosteriorAdj
#'
#' Calculates the posterior mean and standard deviation for the neighborhood adjacencies
#'  from the \code{\link{STBDwDM}} model.
#'
#' @param object a \code{\link{STBDwDM}} model object for which predictions
#'  are desired from.
#'
#' @details The function \code{PosteriorAdj} calculates the posterior mean and standard
#'  deviation of the neighborhood adjacencies for each pairwise location. The neighborhood
#'  structure used to do this comes from Berchuck et al. 2017.
#'
#' @return \code{PosteriorAdj} returns a matrix containing the following columns.
#'
#'   \describe{
#'
#'   \item{\code{i}}{Location \code{i} (i.e. which row/column on the adjacency matrix W).}
#'
#'   \item{\code{j}}{Location \code{j} (i.e. which row/column on the adjacency matrix W).}
#'
#'   \item{\code{DM}}{The dissimilarity metric between locations \code{i} and \code{j}.}
#'
#'   \item{\code{meant}}{The posterior mean of the neighborhood adjacency between location
#'    \code{i} and \code{j} at time \code{t, t = 1, ... , Nu}.}
#'
#'   \item{\code{sdt}}{The posterior mean of the neighborhood adjacency between location
#'    \code{i} and \code{j} at time \code{t, t = 1, ... , Nu}.}
#'
#'   }
#'
#' @author Samuel I. Berchuck
#' @export
PosteriorAdj <- function(object) {

  ###Check Inputs
  if (missing(object)) stop('"object" is missing')
  if (!is.STBDwDM(object)) stop('"object" must be of class STBDwDM')

  ###Set data objects
  DatObj <- object$datobj
  Z <- DatObj$Z
  # AdjacentEdgesBoolean <- DatObj$AdjacentEdgesBoolean
  Nu <- DatObj$Nu
  ScaleDM <- DatObj$ScaleDM
  W <- DatObj$W
  AdjacentEdgesBoolean <- (W == 1) & (!lower.tri(W))

  ###Set parameter objects
  NKeep <- dim(object$mu)[1]
  Alpha <- object$alpha

  ###Get posterior Wij
  ZMat <- as.matrix(Z, ncol = 1)
  Wij <- cbind(which(AdjacentEdgesBoolean, arr.ind = TRUE), ZMat * ScaleDM)
  WijMat <- matrix(0, ncol = length(Z), nrow = NKeep)
  for (t in 1 : Nu) {
    for (s in 1 : NKeep) {
      WijMat[s, ] <- exp( -Alpha[s, t] * Z)
    }
    Wij <- cbind(Wij, apply(WijMat, 2, mean), apply(WijMat, 2, sd))
  }
  colnames(Wij) <- c("i", "j", "DM", paste(c("mean", "sd"), rep(1 : Nu, each = 2), sep=""))
  Wij <- structure(Wij, class = "PosteriorAdj")
  return(Wij)
}



###Verify the class of our regression object------------------------------------------------------------------------
#' is.PosteriorAdj
#'
#' \code{is.PosteriorAdj} is a general test of an object being interpretable as a
#' \code{\link{PosteriorAdj}} object.
#'
#' @param x object to be tested.
#'
#' @details The \code{\link{PosteriorAdj}} class is defined as the posterior adjacency
#'  object that results from the \code{\link{PosteriorAdj}} function.
#'
#' @export
is.PosteriorAdj <- function(x) {
  identical(attributes(x)$class, "PosteriorAdj")
}

