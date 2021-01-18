#' @useDynLib fastcmh
#' @importFrom Rcpp sourceCpp
#' @import bindata
#' @importFrom stats rbinom
#' @importFrom utils str write.table



fastcmhInfo <- function(){
    cat("For fastcmh documentation, type \"?fastcmh\". Note that using the compiler flag -O3 in ~/.R/Makevars could improve the speed of the program.\n")
}
#the above function exists to avoid roxygen2 messing up NAMESPACE


#' Run the fastcmh algorithm 
#'
#' This function runs the FastCMH algorithm on a particular data set. 
#'
#' @param folder  The folder in which the data is saved. If the any of 
#' \code{data}, \code{label} and \code{pvalue} arguments are not specified, 
#' then filenames must have following a naming convention inside the folder: 
#' the data file is \code{"data.txt"} 
#' (i.e. the full path is \code{"folder/data.txt"}), the phenotype label file 
#' is \code{label.txt}, and covariate label file is 
#' \code{cov.txt}. More details on the structure of these files is given 
#' below, or the user can use the \code{\link{makefastcmhdata}} function to 
#' see an example of the correct data formats. If \code{folder="/data/"}, the 
#' data in \code{fastcmh/inst/extdata} is used.
#'
#' @param data The filename for the data file. Default is \code{NULL}. The 
#' data file must be an \code{L x n} txt file containing only \code{0}s and 
#' \code{1}s, which are space-separated in each row, while each row is on a 
#' separate newline.
#'
#' @param label The filename for the phenotype label file. Default is 
#' \code{NULL}. The label file should consist of a single column (i.e. each 
#' row is on a separate line) of \code{0}s and \code{1}s.
#'
#' @param cov The filename for the covariate label file. Default is 
#' \code{NULL}. The \code{cov} file contains a single column of positive 
#' integers. The first row, containing value \eqn{n_1}, specifies that the 
#' first \eqn{n_1} columns have covariate value \code{1}; the second row, 
#' containing \eqn{n_2}, specifies that the next \eqn{n_2} rows have 
#' covariate value \code{2}, etc.
#'
#' @param alpha The value of the FWER; must be a number between 0 and 1. 
#' Default is \code{alpha=0.05}.
#'
#' @param Lmax The maximum length of significant intervals which is 
#' considered. Must be a non-negative integer. For example, \code{Lmax=10} 
#' searches for significant intervals up to length 10. Setting \code{Lmax=0}
#' will search for significant intervals up to any length (with algorithm 
#' pruning appropriately). Default is \code{Lmax=0}.
#'
#' @param showProcessing A flag which will turn printing to screen on/off. 
#' Default is \code{FALSE} (which is \dQuote{off}).
#'
#' @param saveAllPvals A flag which controls whether or not all the intervals
#' (less than minimum attainable pvalue) will be returned. Default is 
#' \code{FALSE} (which is \dQuote{no, do notreturn all intervals}).
#'
#' @param doFDR A flag which controls whether or not Gilbert's Tarone FDR 
#' procedure (while accounting for positive regression dependence) is 
#' performed. Default is \code{FALSE} (which is \dQuote{no, do not do FDR}).
#'
#' @param useDependenceFDR A flag which controls whether or not Gilbert's 
#' Tarone FDR procedure uses the dependent formulation by Benjamini and 
#' Yekutieli (2001), which further adjusts alpha by dividing by the harmonic 
#' mean. This flag is only used if \code{doFDR==TRUE}. Default is \code{FALSE}.
#'
#' @param saveToFile A flag which controls whether or not the results are 
#' saved to file. By default, \code{saveToFile=FALSE}, and the data frame is 
#' returned in R. See the examples below.
#'
#' @param saveFilename A string which gives the filename to which the output 
#' is saved (needs to have \code{saveToFile=TRUE}) as an RData file. Default 
#' is \code{"fastcmhresults.RData"}.
#'
#' @param saveFolder A string which gives the path to which the output will 
#' be saved (needs to have \code{saveToFile=TRUE}). Default is \code{"./"}.
#'
#' 
#' @section Details:
#' This function runs the FastCMH algorithm on a particular data set in 
#' order to discover intervals that are statistically significantly 
#' associated with a particular label, while accounting for categorical 
#' covariates. 
#'
#' The user must either supply the folder, which contains files named 
#' \code{"data.txt"}, \code{"label.txt"} and \code{"cov.txt"}, or the 
#' non-default filenames must be specified individually. See the descriptions of arguments \code{data}, \code{label} and \code{cov} to see the format of 
#' the input files, or make a small sample data file using the 
#' \code{\link{makefastcmhdata}} function.
    
#' By default, filtered results are provided. The user also has the option 
#' of using an FDR procedure rather than the standard FWER-preserving 
#' procedure. 
#'
#' 
#' @section Value:
#' \code{runfastcmh} will return a list if \code{saveToFile=FALSE} (default 
#' setting), otherwise it will save the list in an .RData file. The fields 
#' of the list are:
#' 
#' \describe{
#' \item{\code{sig}}{a dataframe listing the significant intervals, after 
#' filterting. Columns \code{start}, \code{end} and \code{pvalue} indicate 
#' the start and end points of the interval (inclusive), and the 
#' \emph{p}-value for that interval.}
#' 
#' \item{\code{unfiltered}}{a dataframe listing all the significant intervals
#' before filtering. The filtering compares the overlapping intervals and 
#' returns the interval with the smallest p-value in each cluster of 
#' overlapping intervals. Dataframe has has structure as \code{sig}.}
#'
#' \item{\code{fdr}}{(if doFDR==TRUE) significant intervals using Gilbert's 
#' FDR-Tarone procedure, after filtering. Dataframe has same structure as 
#' \code{sig}.}
#'
#' \item{\code{unfilteredFdr}}{(if doFDR==TRUE) a dataframe listing all the significant intervals before filtering. See description of \code{unfiltered}.}
#' 
#' \item{\code{allTestablle}}{(if saveAllPvals==TRUE) a dataframe listing all
#' the testable intervals, many of which will not be significant. Dataframe 
#' has same structure as \code{sig}.}
#' 
#' \item{\code{histObs}}{Together with histFreq gives a histogram of maximum 
#' attainable CMH statistics.}
#' 
#' \item{\code{histFreq}}{Histogram of maximum attainable CMH statistics (only
#' reliable in the testable range).}
#' 
#' \item{\code{summary}}{a character string summarising the results. Use 
#' \code{cat(...$summary)} to print the results with the correct 
#' indentation/new lines.}
#' 
#' \item{\code{timing}}{a list containing (i) \code{details}, a character 
#' string summarising the runtime values for the experiment - use 
#' \code{cat(...$timing$details)} for correct indentation, etc. 
#' (ii) \code{exec}, the total execution time. (iii) \code{init}, the time 
#' to initialise the objects. (iv) \code{fileIO}, the time to read the input 
#' files. (v) \code{compSigThresh}, the time to compute the significance 
#' threshold. (vi) \code{compSigInt}, the time to compute the significant 
#' intervals.} 
#'} 
#' 
#' 
#' 
#' @section Author(s):
#' Felipe Llinares Lopez, Dean Bodenham
#' 
#'
#' @section See Also:
#' \code{\link{makefastcmhdata}}
#'
#'
#' @section References:
#'
#' Gilbert, P. B. (2005) \emph{A modified false discovery rate 
#' multipl-comparisons procedure for discrete data, applied to human 
#' immunodeficiency virus genetics}. Journal of the Royal Statistical 
#' Society: Series C (Applied Statistics), 54(1), 143-158.
#'
#' Benjamini, Y., Yekutieli, D. (2001). \emph{The control of the false 
#' discovery rate in multiple testing under dependency}. 
#' Annals of Statistics, 29(4), 1165-1188.
#'
#'
#' @examples
#' #Example with default naming convention used for data, label and cov files
#' # Note: using "/data/" as the argument for folder 
#' #       accesses the data/ directory in the fastcmh package folder
#' mylist <- runfastcmh("/data/") 
#'
#' #Example where the progress will be shown
#' mylist <- runfastcmh(folder="/data/", showProcessing=TRUE) 
#'
#' #Example where many parameters are specified
#' mylist <- runfastcmh(folder="/data/", data="data2.txt", alpha=0.01, Lmax=7)
#'
#' #Example where Gilbert's Tarone-FDR procedure is used
#' mylist <- runfastcmh("/data/", doFDR=TRUE) 
#'
#' #Example where FDR procedure takes some dependence structures into account
#' mylist <- runfastcmh("/data/", doFDR=TRUE, useDependenceFDR=TRUE) 
#'
#'
#' @export
runfastcmh <- function(folder=NULL, data=NULL, label=NULL, cov=NULL, alpha=0.05, Lmax=0, showProcessing=FALSE, saveAllPvals=FALSE, doFDR=FALSE, useDependenceFDR=FALSE, saveToFile=FALSE, saveFilename="fastcmhresults.RData", saveFolder=NULL){

    #check if any of the file names are NULL; not specified
    if ( (is.null(data)) || (is.null(label)) || (is.null(cov)) )
        notAllFileNamesProvided <- TRUE
    else
        notAllFileNamesProvided <- FALSE

    #if folder not specified AND not all filenames are provided, stop and throw error.
    if ( (is.null(folder)) && (notAllFileNamesProvided)){
    filenameErrorString <- paste0("Folder is not provided and not all files specified. runfastcmh aborted. \nPlease either specify the folder, and then use 'data.txt', 'label.txt', and 'cov.txt' as the filenames for the data, phenotype labels and covariate labels, respectively, or specify all filenames, including paths.\n")

      stop(filenameErrorString, call. = FALSE)
    }


    #get default filenames, in case some are missing
    if (!is.null(folder)){

        #special case:
        #now make it access data folder in R package
        if (folder=="/data/"){
            #for some reason, automatically goes to fastcmh/inst/
            folder <- file.path(system.file(package="fastcmh"), "extdata")
            folder <- paste0(folder, .Platform$file.sep)
        }

        #normalise path for Windows
#        folder <- normalizePath(folder)
        folder <- fixSeparator(folder)
        #add / to end of folder, if not present
        folder <- checkdir(folder)

        #get defaults
        df <- loadDefaultFileNames(folder) 

        #if data filename is null, use default - later will check existence    
        if ( is.null(data) ) {
            data <- df$xfilename
        } else {
            data <- paste0(folder, data)
        }

        #if label filename is null, use default - later will check existence    
        if ( is.null(label) ) {
            label <- df$yfilename
        } else {
            label <- paste0(folder, label)
        }

        #if cov filename is null, use default - later will check existence    
        if ( is.null(cov) ) {
            cov <- df$cfilename
        } else {
            cov <- paste0(folder, cov)
        }
    }

    #check that data file exists
    if (!file.exists(data)){
        dataErrorString <- paste0("data file ", data, " does not exist. runfastcmh aborted\n")
        stop(dataErrorString, call.=FALSE)
    } else {
        #expand path in case of tilde ~
        data <- path.expand(data)
    }

    #check that label file exists
    if (!file.exists(label)){
        labelErrorString <- paste0("label file ", label, " does not exist. runfastcmh aborted\n")
        stop(labelErrorString, call. = FALSE)
    } else {
        #expand path in case of tilde ~
        label <- path.expand(label)
    }

    #check that cov file exists
    if (!file.exists(cov)){
        covErrorString <- paste0("cov file ", cov, " does not exist. runfastcmh aborted\n")
        stop(covErrorString, call. = FALSE)
    } else {
        #expand path in case of tilde ~
        cov <- path.expand(cov)
    }


    #check alpha
    if (!is.numeric(alpha)){
        stop("alpha is not numeric. runfastcmh aborted. Please ensure alpha is numeric and in interval(0, 1).\n")
    }
    if (  (is.numeric(alpha)) && ( (alpha <= 0) || (alpha >= 1) )  ){
        stop("alpha is not in the interval (0, 1). runfastcmh aborted. Please ensure alpha is numeric and in interval(0, 1).\n")
    }

    #check Lmax is an integer
    tol <- 1e-5
    #this code checks float part - no easy function in R to check an integer 
    LMaxRemainder <- min(abs(c(Lmax%%1, Lmax%%1-1)))
    if (is.numeric(Lmax)){
        if ((LMaxRemainder > tol) || (LMaxRemainder < 0) ){
            stop("Lmax is not a positive integer. runfastcmh aborted.\n")
        }
    } else { 
        stop("Lmax is not a positive integer. runfastcmh aborted.\n")
    }



    #check showProcessing:
    #TODO
    #should there be error messages if the booleans are messed up?
    #seems that there is no base is.boolean function; easy way to force non-true
    #to be FALSE
    if (is.null(showProcessing)){
        showProcessing <- FALSE
    } 
    if(isTRUE(showProcessing)){
        showProcessing <- TRUE
    } else {
        showProcessing <- FALSE
    }

    #check saveAllPvals:
    if (is.null(saveAllPvals)){
        saveAllPvals <- FALSE
    } 
    if(isTRUE(saveAllPvals)){
        saveAllPvals <- TRUE
    } else {
        saveAllPvals <- FALSE
    }

    #check doFDR:
    if (is.null(doFDR)){
        doFDR <- FALSE
    } 
    if(isTRUE(doFDR)){
        doFDR <- TRUE
        #if doFDR is true, need saveAllPvals to be true...just a quick fix for now
        saveAllPvals <- TRUE
    } else {
        doFDR <- FALSE
    }

    #check useDependenceFDR:
    if (is.null(useDependenceFDR)){
        useDependenceFDR <- FALSE
    } 
    if(isTRUE(useDependenceFDR)){
        useDependenceFDR <- TRUE
    } else {
        useDependenceFDR <- FALSE
    }

    #filtered intervals filename

    #check saveToFile:
    if (is.null(saveToFile)){
        saveToFile <- FALSE
    } 
    if(isTRUE(saveToFile)){
        saveToFile <- TRUE
    } else {
        saveToFile <- FALSE
    }

    #check saveFilename:
    if (is.null(saveFilename)){
        saveFilename <- "fastcmhresults.RData" 
    } 

    #check saveFolder:
    if (is.null(saveFolder)){
        saveFolder <- "./" 
#        saveFolder <- normalizePath(saveFolder)
        saveFolder <- fixSeparator(saveFolder)
        saveFolder <- checkIsFolder(saveFolder)
    } 


    #now return df
    fastcmhresults <- main_fastcmh_cpp(data, label, cov, alpha, Lmax, showProcessing, saveAllPvals, doFDR, useDependenceFDR)

    #if save to file:
    if (saveToFile){
        saveFolder <- checkIsFolder(saveFolder)
        saveFolder <- checkFolderExists(saveFolder)
        if (saveFolder != "./"){
            saveFilename <- paste0(saveFolder, saveFilename)
        }
        save(fastcmhresults, file=saveFilename)
        cat("saving fastcmh results to: ", saveFilename, "\n", sep="")
    } else {
        #return results:
        return(fastcmhresults)
    }

}







#make sure it is a directory (last character is "/")
checkIsFolder <- function(str){
    returnString <- str
    mysep <- .Platform$file.sep
    if( substr(str, nchar(str), nchar(str)) == mysep ){
        #do nothing
    } else{
        returnString <- paste0(returnString, mysep)
    }
    return(returnString)
}


#make sure it is a directory (last character is "/")
#also, if folder does not exist, create it.
checkFolderExists <- function(str){
    returnString <- str
    mysep <- .Platform$file.sep
    if( substr(str, nchar(str), nchar(str)) == mysep ){
        #do nothing
    } else{
        returnString <- paste0(returnString, mysep)
    }

    #now check if folder exists; if not, create it
    if (file.exists(returnString)==FALSE){
        cat("folder ", returnString, " does not exist. Creating folder now...\n", sep="")
        dir.create(returnString, showWarnings=FALSE, recursive=TRUE)
    }
    return(returnString)
}




#get a data.frame of default names, with folder prefix
loadDefaultFileNames <- function(folder){
    #make sure last character is "/", otherwise add it
    folder <- checkIsFolder(folder)

    #expand path in case of tilde ~
    folder <- path.expand(folder)

    xfilename <- paste0(folder, "data.txt")
    yfilename <- paste0(folder, "label.txt")
    cfilename <- paste0(folder, "cov.txt")
    Lmax <- 10
    alpha <- 0.05
    basefilename <- paste0(folder, "result")
    pvalfilename <- paste0(folder, "pval.txt")
    filteredfilename <- paste0(folder, "filtered.txt")
    dlist <-list(xfilename=xfilename, yfilename=yfilename, cfilename=cfilename, Lmax=Lmax, alpha=alpha, basefilename=basefilename, pvalfilename=pvalfilename, filteredfilename=filteredfilename)
    return(dlist)
}



#the main function
main_fastcmh_cpp <- 
function(xfilename, yfilename, cfilename, alpha, Lmax, showProcessing, saveAllPvals, doFDR, useDependenceFDR){
    returnList <- tryCatch(main_fastcmh2(xfilename, yfilename, cfilename, alpha, Lmax, showProcessing, saveAllPvals, doFDR, useDependenceFDR),    error = print)
    return(returnList)
}


