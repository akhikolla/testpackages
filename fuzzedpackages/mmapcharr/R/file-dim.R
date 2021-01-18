################################################################################

#' Number of lines
#' 
#' Number of lines of a file.
#'
#' @param file Path to file.
#'
#' @return The number of lines of the file.
#' @export
#'
#' @examples
#' tmpfile <- tempfile()
#' write(1:5, tmpfile, ncolumns = 1)
#' nline(tmpfile)
nline <- function(file) {
  nline_cpp(charSepXPtr(file, n = file.size(file), m = 1, r = 0))
}

################################################################################

#' Size of line
#' 
#' Number of elements of each line of a file.
#'
#' @param file Path to file.
#'
#' @return The number of elements of each line of a file.
#' @export
#'
#' @examples
#' tmpfile <- tempfile()
#' write(1:10, tmpfile, ncolumns = 2)
#' nline(tmpfile)
nelem <- function(file) {
  firstline <- strsplit(readLines(file, n = 1), split = "", fixed = TRUE)[[1]]
  length(firstline)
}

################################################################################

#' File dimensions
#' 
#' Number of lines and columns of file (and extra 'return' characters).
#'
#' @param file Path to file.
#'
#' @return The number of lines and columns of file 
#'   (and extra 'return' characters).
#' @export
#'
#' @examples
#' tmpfile <- tempfile()
#' write(0:9, tmpfile, ncolumns = 2)
#' dim_file(tmpfile)
dim_file <- function(file) {
  
  s <- file.size(file)
  n <- nline(file)
  cat(paste(n, "lines detected.\n"))
  p <- nelem(file)
  if (!(is_int(s / n) && is_int(m <- (p + 1) / 2))) stop2(ERROR_FILE)
  cat(paste(m, "columns detected.\n"))
  
  # s = n * (p + nchar_newline)
  r <- s / n - p  ## should be 2 on Windows and 1 on Unix-like systems 
  is_win <- (.Platform$OS.type == "windows")
  if (r == 1) {
    if (is_win) warning2("Only one 'return' characters detected, yet Windows")
  } else if (r == 2) {
    if (!is_win) warning2("Two 'return' characters detected, yet not Windows")
  } else {
    stop2(ERROR_DIM)
  }
  
  res <- c(nrow = n, ncol = m, nextra = r)
  storage.mode(res) <- "integer"
  
  res
}

################################################################################