#' write.csv.R
#'
#' Save summaries of partitioned breeding values to CSV files on disk for further
#' analyses of processing with other software or just for saving (backing up)
#' results.
#'
#' Function \code{\link[utils]{write.csv}} from the \pkg{utils} package works
#' when exported object is a \code{\link[base]{data.frame}} or a
#' \code{\link[base]{matrix}}. This is an attempt to make this function generic
#' so that one can define \code{write.csv} methods for other objects.
#'
#' @seealso \code{\link[utils]{write.csv}} help page on the default \code{write.csv} and \code{write.csv2}
#' methods in the \pkg{utils} package;
#' \code{\link[AlphaPart]{summary.AlphaPart}} and \code{\link[AlphaPart]{AlphaPart}}
#' help pages on the objects of \code{summaryAlphaPart} and \code{AlphaPart} classes.
#'
#' @param x AlphaPart, object returned from \code{\link[AlphaPart]{AlphaPart}} function or
#' summaryAlphaPart, object returned from \code{\link[AlphaPart]{summary.AlphaPart}} function.
#' @param file Character, file name with or without .csv extension, e.g., both "file" and "file.csv" are valid.
#' @param traitsAsDir Logical, should results be saved within trait folders;
#' the construction is \code{file.path(dirname(file), trait, basename(file))};
#' folders are created if they do not exist.
#' @param csv2 Logical, export using \code{\link[utils]{write.csv2}} or \code{\link[utils]{write.csv}}.
#' @param row.names Logical, export row names as well?
#' @param ... Other options passed to \code{\link[utils]{write.csv2}} or \code{\link[utils]{write.csv}}.
#'
#' @example inst/examples/examples_write.csv.R
#'
#' @return  \item{write.csv}{See \code{\link[utils]{write.csv}} for details.}
#'          \item{write.csv.AlphaPart}{For each trait (list component in \code{x}) a file is saved on disk with name
#' "AlphaPart_trait.csv", where the file will hold original data and breeding value partitions.
#' With \code{traitsAsDir=TRUE} files are saved as "trait/file_trait.csv".
#' File names are printed on screen during the process of export and at the end invisibly returned.}
#'          |item{write.csv.summaryAlphaPart}{For each trait (list component in \code{x}) a file partitions named
# "file_trait.csv" is saved on disk.
#' With \code{traitsAsDir=TRUE} files are saved as "trait/file_trait_*.csv". File names
#' are printed on screen during the process of export and at the end invisibly returned.}
#'
#' @useDynLib AlphaPart, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom utils write.csv2


#' @export
write.csv  <- function (...) {
  UseMethod("write.csv")
}

#' @describeIn write.csv Default \code{write.csv} method.
write.csv.default <- function (...) {
  utils::write.csv(...)
}

#' @describeIn write.csv  Save partitioned breeding values to CSV files on disk on disk for further
#' analyses or processing with other software or just for saving (backing up)
#' results.
#'
#' @export
write.csv.AlphaPart <- function (x, file, traitsAsDir=FALSE, csv2=TRUE, row.names=FALSE, ...) {
  ## --- Setup ---

  if(length(file) > 1) stop("'file' argument must be of length one")
  if(!("AlphaPart" %in% class(x))) stop("'x' must be of a 'AlphaPart' class")
  fileOrig <- sub(pattern=".csv$", replacement="", x=file)
  ret <- c()

  ## --- Code ---

  for(i in 1:(length(x)-1)) { ## loop over traits
    if(traitsAsDir) {
      dir.create(path=file.path(dirname(fileOrig), x$info$lT[i]), recursive=TRUE, showWarnings=FALSE)
      file <- file.path(dirname(fileOrig), x$info$lT[i], basename(fileOrig))
    }
    fileA <- paste(file, "_", x$info$lT[i], ".csv", sep="")
    ret   <- c(ret, fileA)
    cat(fileA, "\n")
    
    if(csv2) {
      write.csv2(x=x[[i]], file=fileA, row.names=row.names, ...)
    } else {
      write.csv(x=x[[i]], file=fileA, row.names=row.names, ...)
    }
  }
  
  ## --- Return ---
  
  invisible(ret)
}


#' @describeIn write.csv Save summaries of partitioned breeding values to CSV files on disk for further
#' analyses of processing with other software or just for saving (backing up)
#' results.
#'
#' @export
write.csv.summaryAlphaPart <- function (x, file, traitsAsDir=FALSE, csv2=TRUE, row.names=FALSE, ...) {
  ## --- Setup ---

  if(length(file) > 1) stop("'file' argument must be of length one")
  if(!("summaryAlphaPart" %in% class(x))) stop("'x' must be of a 'summaryAlphaPart' class")
  fileOrig <- sub(pattern=".csv$", replacement="", x=file)
  ret <- c()

  ## --- Code ---

  for(i in 1:(length(x)-1)) { ## loop over traits
    if(traitsAsDir) {
      dir.create(path=file.path(dirname(fileOrig), x$info$lT[i]), recursive=TRUE, showWarnings=FALSE)
      file <- file.path(dirname(fileOrig), x$info$lT[i], basename(fileOrig))
    }
    fileA <- paste(file, x$info$lT[i], ".csv", sep="_")
    ret  <- c(ret, fileA)
    cat(fileA, "\n")

    
    if(csv2) {
      write.csv2(x=x[[i]], file=fileA, row.names=row.names, ...)
    } else {
      write.csv(x=x[[i]], file=fileA, row.names=row.names, ...)
    }
  }
  
  ## --- Return ---
  invisible(ret)
}

