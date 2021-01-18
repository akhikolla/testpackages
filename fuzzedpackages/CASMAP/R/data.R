#' Get the path to the example data file for regionGWAS mode
#'
#' Path to \code{CASMAP_example_data_1.txt} in \code{inst/extdata}.
#' A dataset containing binary samples for the regionGWAS method.
#' There are accompanying labels and covariates dataset.
#'
#' @format A matrix of \code{0}s and \code{1}s, with 1000 rows (features) 
#'         and 100 columns
#'         (samples). In other words, each column is a sample, and each sample
#'         has 1000 binary features.
#' 
#' @details Path to the file containing the data, for reading in to
#'          CASMAP object using the \code{readFiles} function.
#'          Note that the significant region is \code{[99, 102]}.
#'
#'
#' @seealso \code{getExampleLabelsFilename}, 
#'          \code{getExampleCovariatesFilename}
#'
#' @export
#' @examples
#' datafile <- getExampleDataFilename()
getExampleDataFilename <- function(){
    filename <-  system.file("extdata", "CASMAP_example_data_1.txt", 
        package = "CASMAP", mustWork = TRUE)
    return(filename)
}


#' Get the path to the example labels file for regionGWAS mode
#'
#' Path to \code{CASMAP_example_labels_1.txt} in \code{inst/extdata}.
#' A dataset containing the binary labels for the data in the file 
#' \code{CASMAP_example_data_1.txt}, the path to which is given by 
#' \code{getExampleDataFilename}.
#'
#' @format A single column of 100 labels, each of which is either \code{0}
#'         or \code{1}.
#'
#' @details Path to the file containing the labels, for reading in to
#'          CASMAP object using the \code{readFiles} function.
#'
#' 
#' @seealso \code{getExampleDataFilename}, 
#'          \code{getExampleCovariatesFilename}
#'
#' @export
#' @examples
#' labelsfile <- getExampleLabelsFilename()
getExampleLabelsFilename <- function(){
    filename <-  system.file("extdata", "CASMAP_example_labels_1.txt", 
                             package = "CASMAP", mustWork = TRUE)
    return(filename)
}



#' Get the path to the example covariates file for regionGWAS mode
#'
#' Path to \code{CASMAP_example_covariates_1.txt} in \code{inst/extdata}.
#' The covariates categories for the data set  
#' \code{CASMAP_example_data_1.txt}, the path to which is given by 
#' \code{getExampleDataFilename}.
#'
#' @format A single column vector of 100 labels, each of which
#'         is \code{0} or \code{1} (same format as labels file).
#' 
#' @details Path to the file containing the labels, for reading in to
#'          CASMAP object using the \code{readFiles} function.
#'
#' @seealso \code{getExampleDataFilename}, 
#'          \code{getExampleLabelsFilename}
#'
#' @export
#' @examples
#' covfile <- getExampleCovariatesFilename()
getExampleCovariatesFilename <- function(){
    filename <-  system.file("extdata", "CASMAP_example_covariates_1.txt", 
                             package = "CASMAP", mustWork = TRUE)
    return(filename)
}


#' Get the path to the example significant intervals file
#'
#' Path to \code{CASMAP_example_covariates_1.txt} in 
#' \code{inst/extdata}.
#'
#' @keywords internal
#' @examples
#' sigregfile <- getExampleSignificantRegionsFilename()
getExampleSignificantRegionsFilename <- function(){
    filename <-  system.file("extdata", "CASMAP_example_sig_regions_1.txt", 
                             package = "CASMAP", mustWork = TRUE)
    return(filename)
}



