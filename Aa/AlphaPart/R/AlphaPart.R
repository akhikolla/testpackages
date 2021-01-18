#' AlphaPart.R
#'
#' A function to partition breeding values by a path variable. The partition method is
#' described in García-Cortés et al., 2008: Partition of the genetic trend to validate multiple selection decisions.
#' Animal : an international journal of animal bioscience. DOI: 10.1017/S175173110800205X
#'
#' @details
#' Pedigree in \code{x} must be valid in a sense that there are:\itemize{
#' \item{no directed loops (the simplest example is that the individual identification is equal to the identification of a father or mother)}
#' \item{no bisexuality, e.g., fathers most not appear as mothers}
#' \item{father and/or mother can be unknown (missing) - defined with any "code" that is dif ferent from existing identif ications}
#' }
#'
#' Unknown (missing) values for breeding values are propagated down the pedigree
#' to provide all available values from genetic evaluation. Another option is
#' to cut pedigree links - set parents to unknown and remove them from pedigree
#' prior to using this function - see \code{\link[AlphaPart]{pedSetBase}} function.
#' Warning is issued in the case of unknown (missing) values.
#'
#' In animal breeding/genetics literature the model with the underlying pedigree type
#' \code{"IPP"} is often called animal model, while the model for pedigree type \code{"IPG"}
#' is often called sire - maternal grandsire model. With a combination of \code{colFid}
#' and \code{colMid} mother - paternal grandsire model can be accomodated as well.
#'
#' Argument \code{colBy} can be used to directly perform a summary analysis by
#' group, i.e., \code{summary(AlphaPart(...), by="group")}. See \code{\link[AlphaPart]{summary.AlphaPart}}
#' for more. This can save some CPU time by skipping intermediate steps. However,
#' only means can be obtained, while \code{summary} method gives more flexibility.
#'
#' @seealso
#' \code{\link[AlphaPart]{summary.AlphaPart}} for summary method that works on output of \code{AlphaPart},
#' \code{\link[AlphaPart]{pedSetBase}} for setting base population,
#' \code{\link[AlphaPart]{pedFixBirthYear}} for imputing unknown (missing) birth years,
#' \code{\link[pedigree]{orderPed}} in \pkg{pedigree} package for sorting pedigree
#'
#' @references Garcia-Cortes, L. A. et al. (2008) Partition of the genetic trend to validate multiple selection
#' decisions. Animal, 2(6):821-824. \url{http://dx.doi.org/10.1017/S175173110800205X}
#'
#' @param x data.frame , with (at least) the following columns: individual, father, and mother identif ication,
#' and year of birth; see arguments \code{colId},
#' \code{colFid}, \code{colMid}, \code{colPath}, and \code{colBV}; see also details about the validity of pedigree.
#' @param pathNA Logical, set dummy path (to "XXX") where path information is unknown (missing).
#' @param recode Logical, internally recode individual, father and, mother identification to
#' \code{1:n} codes, while missing parents are defined with \code{0}; this option
#' must be used if  identif ications in \code{x} are not already given as \code{1:n}
#' codes, see also argument \code{sort}.
#' @param unknown Value(s) used for representing unknown (missing) parent in \code{x}; this options
#' has an effect only when \code{recode=FALSE} as it is only needed in that situation.
#' @param sort Logical, initially sort \code{x} using \code{orderPed()} so that children follow
#' parents in order to make imputation as optimal as possible (imputation is performed
#' within a loop from the first to the last unknown birth year); at the end original
#' order is restored.
#' @param verbose Numeric, print additional information: \code{0} - print nothing, \code{1} - print
#' some summaries about the data.
#' @param profile Logical, collect timings and size of objects.
#' @param printProfile Character, print profile info on the fly (\code{"fly"}) or at the end (\code{"end"}).
#' @param pedType Character, pedigree type: the most common form is \code{"IPP"} for Individual, Parent
#' 1 (say father), and Parent 2 (say mother) data; the second form is \code{"IPG"} for
#' Individual, Parent 1 (say father), and one of Grandparents of Parent 2 (say maternal
#' grandfather).
#' @param colId Numeric or character, position or name of a column holding individual identif ication.
#' @param colFid Numeric or character, position or name of a column holding father identif ication.
#' @param colMid Numeric or character, position or name of a column holding mother identif ication or
#' maternal grandparent identif ication if  \code{pedType="IPG"} .
#' @param colPath Numeric or character, position or name of a column holding path information.
#' @param colBV Numeric or character, position(s) or name(s) of column(s) holding breeding Values.
#' @param colBy Numeric or character, position or name of a column holding group information (see details).
#'
#' @example inst/examples/examples_AlphaPart.R
#' @return An object of class \code{AlphaPart}, which can be used in further analyses - there is a handy summary
#' method (\code{\link[AlphaPart]{summary.AlphaPart}} works on objects of \code{AlphaPart} class) and a plot method
#' for its output (\code{\link[AlphaPart]{plot.summaryAlphaPart}} works on objects of \code{summaryAlphaPart} class).
#' Class \code{AlphaPart} is a list. The first \code{length(colBV)} components (one for each trait and named with
#' trait label, say trt) are data frames. Each data.frame contains:
#'   \item{\code{x}}{columns from initial data \code{x}}
#'   \item{trt_pa}{parent average}
#'   \item{trt_w}{Mendelian sampling term}
#'   \item{trt_path1, trt_path2, ...}{breeding value partitions}
#'
#' The last component of returned object is also a list named \code{info} with the following components holding
#' meta information about the analysis:
#'   \item{path}{column name holding path information}
#'   \item{nP}{number of paths}
#'   \item{lP}{path labels}
#'   \item{nT}{number of traits}
#'   \item{lT}{trait labels}
#'   \item{warn}{potential warning messages associated with this object}
#'
#' If  \code{colBy!=NULL} the resulting object is of a class \code{summaryAlphaPart},
#' see \code{\link[AlphaPart]{summary.AlphaPart}} for details.
#'
#' If  \code{profile=TRUE}, profiling info is printed on screen to spot any computational bottlenecks.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export
#' @importFrom utils str
#' @importFrom pedigree orderPed
#' @importFrom gdata NAToUnknown
#' @importFrom gdata NAToUnknown
#' @importFrom gdata unknownToNA
#' @importFrom gdata object.size
#' @importFrom stats aggregate




AlphaPart <- function (x, pathNA=FALSE, recode=TRUE, unknown=NA, sort=TRUE, verbose=1, profile=FALSE,
  printProfile="end", pedType="IPP", colId=1, colFid=2, colMid=3, colPath=4, colBV=5:ncol(x),
  colBy=NULL) {

  # TODO: move BV to another object (to simplif y work with McMC or some other
  # TODO: sortPedigree: A rabimo tole nujno za to funkcijo ali samo za summarizing? Hmm, za sortiranje, kajne? Vidis, to nisem lepo sprogramiral – ena funkcija naj bi pocela samo eno stvar na enkrat – to poenostavi kodo. Pusti za sedaj. Future work. Lahko das v TODO file v paketu;)

  ## --- Setup ---

  test <- (length(colId) > 1 | length(colFid) > 1 | length(colMid) > 1 | length(colPath) > 1 | length(colBy) > 1)
  if (test) {
    stop("arguments 'colId', 'colFid', 'colMid', 'colPath', and 'colBy' must be of length 1")
  }

  if (is.null(colBy)) {
    groupSummary <- FALSE
  } else {
    groupSummary <- TRUE
  }

  test <- pedType %in% c("IPP", "IPG")
  if (any(!test)) {
    stop("'pedType' must be either 'IPP' or 'IPG'")
  }

  if (profile) {
    time0 <- Sys.time()
    cat("\nStart:", format(time0), "\n")
    timeRet <- data.frame(task="Start", timeP=time0, time=0, timeCum=0, memory=0, memoryCum=0,
                          stringsAsFactors=FALSE)
    .profilePrint <- function(x, task, printProfile, time, mem, update=FALSE)
    {
      i <- nrow(x)
      x[i + 1, "task"]        <- task
      x[i + 1, "timeP"]       <- time
      x[i + 1, "time"]        <- timeTMP1 <- round(time - x[i, "timeP"], digits=1L)
      x[i + 1, "timeCum"]     <- timeTMP2 <- round(time - x[1, "timeP"], digits=1L)
      if (!update) {
        x[i + 1, "memory"]    <- memTMP1  <- round(mem/1024^2, digits=1L)
        x[i + 1, "memoryCum"] <- memTMP2  <- round(mem/1024^2, digits=1L) + x[i, "memoryCum"]
      } else {
        x[i + 1, "memory"]    <- memTMP1  <- round(mem/1024^2, digits=1L)
        x[i + 1, "memory"]    <- abs(x[i + 1, "memory"] - x[i, "memory"])
        x[i + 1, "memoryCum"] <- memTMP2  <- x[i + 1, "memory"] + x[i, "memoryCum"]
      }
      if (printProfile == "fly") {
        cat("\n", task, ":\n", sep="")
        cat(" - time (this task):",     format(timeTMP1),     "\n")
        cat(" - time (all tasks):",     format(timeTMP2),     "\n")
        cat(" - memory (this object):", paste(memTMP1, "Mb"), "\n")
        cat(" - memory (all objects):", paste(memTMP2, "Mb"), "\n")
      }
      x
    }
  }

  ## --- Sort and recode pedigree ---

  ## Make sure that identifications are numeric if  recode=FALSE
  test <- !sapply(x[, c(colId, colFid, colMid)], is.numeric) & !recode
  if (any(test)) {
    stop("argument 'recode' must be 'TRUE' when identif ications in 'x' are not numeric")
  }

  ## Make sure that colBV columns are numeric
  test <- !sapply(x[, c(colBV)], is.numeric)
  if (any(test)) {
    stop("colBV columns must be numeric!")
    str(x)
  }

  ## Sort so that parents preceede children
  if (sort) {
    recode <- TRUE
    x <- x[order(orderPed(ped=x[, c(colId, colFid, colMid)])), ]
  }

  ## Recode all ids to 1:n
  if (recode) {
    y <- cbind( id=1:nrow(x),
               fid=match(x[, colFid], x[, colId], nomatch=0),
               mid=match(x[, colMid], x[, colId], nomatch=0))
  } else {
    y <- as.matrix(x[, c(colId, colFid, colMid)])
    ## Make sure we have 0 when recoded data is provided
    if (is.na(unknown)) {
      y[, c(colFid, colMid)] <- NAToUnknown(x=y[, c(colFid, colMid)], unknown=0)
    } else {
      if (unknown != 0)  {
        y[, c(colFid, colMid)] <- NAToUnknown(x=unknownToNA(x=y[, c(colFid, colMid)], unknown=unknown), unknown=0)
      }
    }
  }
  y <- cbind(y, as.matrix(x[, colBV]))

  ## Test if  father and mother codes preceede children code - computational engine needs this
  test <- y[, 2] >= y[, 1]
  if (any(test)) {
    print(x[test, ])
    print(sum(test))
    stop("sorting/recoding problem: parent (father in this case) code must preceede children code - use arguments 'sort' and/or 'recode'")
  }

  test <- y[, 3] >= y[, 1]
  if (any(test)) {
    print(x[test, ])
    print(sum(test))
    stop("sorting/recoding problem: parent (mother in this case) code must preceede children code - use arguments 'sort' and/or 'recode'")
  }

  if (profile) {
    timeRet <- .profilePrint(x=timeRet, task="Sort and/or recode pedigree", printProfile=printProfile,
                             time=Sys.time(), mem=(object.size(x) + object.size(y)))
  }

  ## --- Dimensions and Paths ---

  ## Pedigree size
  nI <- nrow(x)

  ## Traits
  lT <- colnames(x[, colBV, drop=FALSE])
  nT <- length(lT) # number of traits
  colnames(y)[4:ncol(y)] <- lT

  ## Missing values
  nNA <- sapply(x[, colBV, drop=FALSE], function(z) sum(is.na(z)))
  names(nNA) <- lT

  ## Paths - P matrix
  test <- is.na(x[, colPath])
  if (any(test)) {
    if (pathNA) {
      x[, colPath] <- as.character(x[, colPath])
      x[test, colPath] <- "XXX"
    } else {
      stop("unknown (missing) value for path not allowed; use 'pathNA=TRUE'")
    }
  }
  if (!is.factor(x[, colPath])) x[, colPath] <- factor(x[, colPath])
  lP <- levels(x[, colPath])
  nP <- length(lP) # number of paths
  P <- as.integer(x[, colPath]) - 1

  ## Groups
  if (groupSummary) {
    test <- is.na(x[, colBy])
    if (any(test)) {
      if (pathNA) {
        x[, colBy] <- as.character(x[, colBy])
        x[test, colBy] <- "XXX"
      } else {
        stop("unknown (missing) value for group not allowed; use 'pathNA=TRUE'")
      }
    }
    if (!is.factor(x[, colBy])) x[, colBy] <- factor(x[, colBy])
    lG <- levels(x[, colBy])
    nG <- length(lG)
    g <- as.integer(x[, colBy])
  }

  if (verbose > 0) {
    cat("\nSize:\n")
    cat(" - individuals:", nI, "\n")
    cat(" - traits: ", nT, " (", paste(lT, collapse=", "), ")", "\n", sep="")
    cat(" - paths: ",  nP, " (", paste(lP, collapse=", "), ")", "\n", sep="")
    if (groupSummary) {
      cat(" - groups: ", nG, " (", paste(lG, collapse=", "), ")", "\n", sep="")
    }
    cat(" - unknown (missing) values:\n")
    print(nNA)
  }

  if (any(nNA > 0)) stop("unknown (missing) values are propagated through the pedigree and therefore not allowed")

  if (profile) {
    timeRet <- .profilePrint(x=timeRet, task="Dimensions and Matrices P", printProfile=printProfile,
                             time=Sys.time(), mem=object.size(P))
  }

  ## --- Compute ---

  ## Prepare stuff for C++
  c1 <- c2 <- 0.5
  if (pedType == "IPG") c2 <- 0.25

  ## Add "zero" row (simplif ies computations with missing parents!)
  y <- rbind(y[1, ], y)
  y[1, ] <- 0
  P <- c(0, P)
  if (groupSummary) g <- c(0, g)

  ## Compute
  if (!groupSummary) {
    tmp <- .Call("AlphaPartDrop",
                 c1_=c1, c2_=c2,
                 nI_=nI, nP_=nP, nT_=nT,
                 y_=y, P_=P, Px_=cumsum(c(0, rep(nP, nT-1))),
                 PACKAGE="AlphaPart")
  } else {
    N <- aggregate(x=y[-1, -c(1:3)], by=list(by=x[, colBy]), FUN=length)
    tmp <- vector(mode="list", length=3)
    names(tmp) <- c("pa", "w", "xa")
    tmp$pa <- tmp$w <- matrix(data=0, nrow=nG+1, ncol=nT)
    tmp$xa <- .Call("AlphaPartDropGroup",
                 c1_=c1, c2_=c2,
                 nI_=nI, nP_=nP, nT_=nT, nG_=nG,
                 y_=y, P_=P, Px_=cumsum(c(0, rep(nP, nT-1))), g_=g,
                 PACKAGE="AlphaPart")
  }

  ## Assign nice column names
  colnames(tmp$pa) <- paste(lT, "_pa", sep="")
  colnames(tmp$w)  <- paste(lT, "_w", sep="")
  colnames(tmp$xa) <- c(t(outer(lT, lP, paste, sep="_")))

  if (profile) {
    timeRet <- .profilePrint(x=timeRet, task="Computing", printProfile=printProfile,
                             time=Sys.time(), mem=object.size(tmp))
  }

  ## --- Massage results ---

  ## Put partitions for one trait in one object (-1 is for removal of the "zero" row)
  ret <- vector(mode="list", length=nT+1)
  t <- 0
  colP <- colnames(tmp$pa)
  colW <- colnames(tmp$w)
  colX <- colnames(tmp$xa)
  for (j in 1:nT) { ## j <- 1
    Py <- seq(t+1, t+nP)
    ret[[j]] <- cbind(tmp$pa[-1, j], tmp$w[-1, j], tmp$xa[-1, Py])
    colnames(ret[[j]]) <- c(colP[j], colW[j], colX[Py])
    t <- max(Py)
  }

  if (profile) {
    timeRet <- .profilePrint(x=timeRet, task="Massage results", printProfile=printProfile,
                             time=Sys.time(), mem=object.size(ret))
  }

  ## Add initial data
  if (!groupSummary) {
    for (i in 1:nT) {
      ## Hassle in order to get all columns and to be able to work with numeric
      ##   or character column "names"
      colX <- colX2 <- colnames(x)
      names(colX) <- colX; names(colX2) <- colX2
      ## ... put current agv in the last column in original data
      colX <- c(colX[!(colX %in% colX[colBV[i]])], colX[colBV[i]])
      ## ... remove other traits
      colX <- colX[!(colX %in% colX2[(colX2 %in% colX2[colBV]) & !(colX2 %in% colX2[colBV[i]])])]
      ret[[i]] <- cbind(x[, colX], as.data.frame(ret[[i]]))
      rownames(ret[[i]]) <- NULL
    }
  }

  ## Additional (meta) info. on number of traits and paths for other methods
  tmp <- colnames(x); names(tmp) <- tmp
  ret[[nT+1]] <- list(path=tmp[colPath], nP=nP, lP=lP, nT=nT, lT=lT, warn=c())
  ## names(ret)[nT+1] <- "info"
  names(ret) <- c(lT, "info")

  if (profile) {
    timeRet <- .profilePrint(x=timeRet, task="Finalizing returned object + adding initial data", printProfile=printProfile,
                             time=Sys.time(), mem=object.size(ret), update=TRUE)
  }

  if  (profile & printProfile == "end") {
    print(timeRet)
  }

  ## --- Return ---

  class(ret) <- c("AlphaPart", class(ret))
  if (groupSummary) {
    ret$by <- colBy
    ret$N <- N
    summary(object=ret, sums=TRUE)
  } else {
    ret
  }


}



