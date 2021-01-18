#' pedFixBirthYear.R
#'
#' A function to fix (impute) missing birth years in pedigree.
#'
#' @details
#' Warnings are issued when there is no information to use to impute birth years or missing
#' values (\code{NA}) are propagated.
#'
#' Arguments \code{down} and \code{na.rm} allow for repeated use of this function, i.e., with
#' \code{down=FALSE} and with \code{down=TRUE} (both in combination with \code{na.rm=TRUE}) in order to
#' propagate information over the pedigree until "convergence".
#'
#' This function can be very slow on large pedigrees with extensive missingness of birth years.
#'
#' @seealso
#' \code{\link[pedigree]{orderPed}} in \pkg{pedigree} package
#'
#' @param x data.frame , with (at least) the following columns: individual, father, and mother identification,
#' and year of birth; see arguments \code{colId},
#' \code{colFid}, \code{colMid}, and \code{colBY}
#' @param interval Numeric, a value for generation interval in years.
#' @param down Logical, the default is to impute birth years based on the birth year of children
#' starting from the youngest to the oldest individuals, while with \code{down=TRUE}
#' birth year is imputed based on the birth year of parents in the opposite order.
#' @param na.rm Logical, remove \code{NA} values when searching for the minimal (maximal) year of birth
#' in children (parents); setting this to \code{FALSE} can lead to decreased success of
#' imputation
#' @param sort Logical, initially sort \code{x} using \code{orderPed()} so that children follow
#' parents in order to make imputation as optimal as possible (imputation is performed
#' within a loop from the first to the last unknown birth year); at the end original
#' order is restored.
#' @param direct Logical, insert inferred birth years immediately so they can be used for successive
#' individuals within the loop.
#' @param report Logical, report success.
#' @param colId Numeric or character, position or name of a column holding individual identification.
#' @param colFid Numeric or character, position or name of a column holding father identification.
#' @param colMid Numeric or character, position or name of a column holding mother identification.
#' @param colBY Numeric or character, position or name of a column holding birth year.
#'
#' @example inst/examples/examples_pedFixBirthYear.R
#'
#' @return Object \code{x} with imputed birth years based on the birth year of children or parents.
#' If \code{report=TRUE} success is printed on the screen as the number of initially, fixed,
#' and left unknown birth years is printed.
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#'
#' @export
#' @importFrom pedigree orderPed

pedFixBirthYear <- function (x, interval, down=FALSE, na.rm=TRUE, sort=TRUE, direct=TRUE,
  report=TRUE, colId=1, colFid=2, colMid=3, colBY=4) {



  ## TODO: Profile function with Rprof() - done one round of optimization - can it be done better? 
  ##
  ## TODO: Can we add another loop (say via argument automatic=TRUE) to repeatedly call this function?
  ##       The idea is to start say with down=FALSE then in next run with down=TRUE and so on. 
  ##       Would have to check if there is no progress between iterations to avoid infinite loop.
  ## 
  ## TODO: subset data from the start to avoid excessive search and update the table at the end?
  ## 
  ## TODO: data.table for "binary (indexed) search"?
  ##        - http://datatable.r-forge.r-project.org/
  ##        - https://docs.google.com/viewer?url=http%3A%2F%2Fdatatable.r-forge.r-project.org%2Fdatatable-faq.pdf
  ##
  ## TODO: sqldf (index?)
  
  ## --- Setup ---

  test <- (length(colId) > 1 | length(colFid) > 1 | length(colMid) > 1 | length(colBY) > 1 | length(interval) > 1)
  if (test) stop("arguments 'interval', 'colId', 'colFid', 'colMid', and 'colBY'  must be of length 1")

  ## Sort so that children follow parents
  if (sort) {
    ordOrig <- rownames(x) 
    x <- x[order(orderPed(ped=x[, c(colId, colFid, colMid)])), ]
  }

  ## --- Core ---

  ## Individuals with missing birth year
  miss <- x[is.na(x[, colBY]), colId]
  ## miss <- miss[1:100] ## handy for Rprof()
  nM <- length(miss)

  ## Setup list and functions regarding the source of information
  if (!down) { ## from children - start from the youngest individual
    miss <- rev(miss)
    funF <- min
    interval <- -interval
  } else {    ## from parents  - start from the oldest individual
    funF <- max
  }

  ## Rprof()

  ## Loop
  if (direct) {

    k <- 0
    for (i in miss) {
      if (!down) { ## birth years of children
        tmp <- x[x[, colMid] %in% i | x[, colFid] %in% i,                                colBY]
      } else {    ## collect birth years of parents
        tmp <- x[x[, colId] %in% Reduce(f=c, x=x[x[, colId] %in% i, c(colFid, colMid)]), colBY]
      }
      if (length(tmp) > 0) {
        tmp2 <- suppressWarnings(funF(tmp, na.rm=na.rm) + interval)
        ## suppressWarnings to avoid warnings about Inf when tmp <- c(NA) is passed to funF with na.rm=TRUE
        if (is.na(tmp2) | is.infinite(tmp2)) {
          warning(paste("missing value for id=", i, sep=""))
        } else {
          x[x[, colId] %in% i, colBY] <- tmp2
          k <- k + 1
        }
      } else {
        warning(paste("no information found for id=", i, sep=""))
      }
    }

  } else {

    ret <- vector(length=nM)
    for (i in 1:nM) {
      if (!down) { ## birth years of children
        tmp <- x[x[, colMid] %in% miss[i] | x[, colFid] %in% miss[i],                          colBY]
      } else {    ## collect birth years of parents
        tmp <- x[x[, colId] %in% Reduce(f=c, x=x[x[, colId] %in% miss[i], c(colFid, colMid)]), colBY]
      }
      if (length(tmp) > 0) {
        ret[i] <- suppressWarnings(funF(tmp, na.rm=na.rm) + interval)
        ## suppressWarnings to avoid warnings about Inf when tmp <- c(NA) is passed to funF with na.rm=TRUE  
      } else {
        warning(paste("no information found for id=", miss[i], sep=""))
      }
    }

    ret[is.infinite(ret)] <- NA
    tmp <- is.na(ret)
    if (any(tmp)) {
      warning(paste("missing value for id:", paste(miss[tmp], collate=", "), sep=""))    
    }

    x <- merge(x=x, y=data.frame(.id=miss, .by=ret), by.x=colnames(x)[colId], by.y=".id", all.x=TRUE)
    tmp <- is.na(x[, colBY]) & !is.na(x[, ".by"]) & x[, ".by"] != 0
    k <- sum(tmp)
    x[tmp, colBY] <- x[tmp, ".by"] 
    x[, ".by"] <- NULL

  }

  ## Rprof(NULL)
  ## return(summaryRprof())  
  
  ## --- Return ---

  if (report) {
    cat("Summary:\n")
    cat(" - initially:",     nM, "\n")
    cat(" - fixed:",          k, "\n")
    cat(" - left:",      nM - k, "\n")
  }
  if (sort) x <- x[ordOrig, ]
  x



}

