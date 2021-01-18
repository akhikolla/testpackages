#' @useDynLib hipread, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
NULL

#' Get path to hipread's example datasets
#'
#' Get access to example extracts.
#'
#' @param path Name of file. If `NULL`, the example files will be listed.
#' @return The filepath to an example file, or if path is empty, a vector of all
#'   available files.
#' @export
#' @examples
#' hipread_example() # Lists all available examples
#' hipread_example("test-basic.dat") # Gives filepath for a basic example
hipread_example <- function(path = NULL) {
  if (is.null(path)) {
    file <- dir(system.file("extdata", package = "hipread"))
  } else {
    file <- system.file("extdata", path, package = "hipread")
    if (!file.exists(file)) {
      all_files <- paste(dir(system.file("extdata", package = "hipread")), collapse = ", ")
      stop(paste0(
        "Could not find file '", path, "' in hipread examples. Available files are:\n",
        all_files
      ))
    }
  }
  file
}


show_progress <- function() {
  # If readr show_progress option isn't set, this means readr hasn't been loaded
  # I think this means we can assume that the user has not set the progress indicator
  # so it will be teh default of TRUE
  readr_progress <- getOption("readr.show_progress")
  if (is.null(readr_progress)) readr_progress <- TRUE

  isTRUE(readr_progress) && # user has not disabled progress bar
    interactive() && # an interactive session
    is.null(getOption("knitr.in.progress")) # Not actively knitting a document
}
