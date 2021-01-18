
#just wait for a key before continuing
readkey <- function()
{
    cat ("Press [enter] to continue")
    line <- readline()
}

#' Demo of fastcmh
#'
#' This function runs a demo for fastcmh, by first creating a sample data set
#' and then running fastcmh on this data set.
#'
#' @param saveToFolder A flag indicating whether or not the data files created
#' for the demo should be saved to file. The default is \code{FALSE}, i.e. no 
#' files are saved to the folder. The only reason to save demo data to a 
#' folder is for the user to be able to have a look at the files after the 
#' demo. 
#' 
#' @param folder The folder in which the data for the demo will be saved. 
#' Default is the current directory, \code{"./"}. The demo data will created 
#' in \code{folder/data} and the results will be saved in 
#' \code{folder/results} as an RData file.
#'
#' @section Details:
#' This function will first create a sample data set in \code{folder/data}, 
#' and will then run \code{runfastcmh} on this data set, before saving the 
# results in \code{folder/results}. The method runs in several steps, with 
#' each step showing the R code that can be used to do the step, then running 
#' that R code, and then waiting for the user to press enter before moving 
#' onto the next step. If \code{saveToFolder=FALSE}, (default) then no files 
#' are saved and all the results are kept in memory.
#'
#' @section See Also:
#' \code{\link{runfastcmh}}
#'
#' @examples
#' demofastcmh() 
#' @export
demofastcmh <- function(saveToFolder=FALSE, folder=NULL){

    #the amount of delay, when a sleep delay was used instead of readkey()
    sleepDelay <- 2


    #check saveAllPvals:
    if (is.null(saveToFolder)){
        saveToFolder <- FALSE
    } 
    if(isTRUE(saveToFolder)){
        saveToFolder <- TRUE
    } else {
        saveToFolder <- FALSE
    }

    #check saveFolder:
    if (is.null(folder)){
        folder <- "./fastcmhdemo/" 
    } 

    #normalise the path
#    folder <- normalizePath(folder)
    folder <- fixSeparator(folder)

    #check it is a folder
    folder <- checkIsFolder(folder)

    cat("\n")
    cat("====================================================\n")
    cat("fastcmh Demo:\n")
    cat("====================================================\n\n")

    #creating folder
    if (saveToFolder==TRUE){

        cat("Creating folders for demo: ", folder, "\n", sep="")
        folder <- checkFolderExists(folder)
        #create data folder and results folder
        datafolder <- paste0(folder, "data")
        datafolder <- checkFolderExists(datafolder)
        resultsfolder <- paste0(folder, "results")
        resultsfolder <- checkFolderExists(resultsfolder)
        #add a delay
        cat("\n")
#    Sys.sleep(sleepDelay)
        readkey();


        cat("----------------------------------------------------\n\n")
        cat("Step 1: Creating data files for demo \n", sep="")
        cat("**(Run command)**\n")
        cat("> makefastcmhdata(\"", datafolder, "\", showOutput=TRUE)\n\n", sep="")
        readkey();

        makefastcmhdata(datafolder, showOutput=TRUE)
        cat("\n\n")
        cat("Sample data created in ", datafolder, "\n\n", sep="")
        readkey();

        cat("\n\n----------------------------------------------------\n")
        cat("\n")
        cat("Step 2: Now running fastcmh, and saving results to file\n")
        cat("**(Run command)**\n")
        cat("> runfastcmh(folder=\"", datafolder, ", saveToFile=TRUE, saveFolder=\"", resultsfolder, "\")\n\n", sep="")
        readkey();


        #Do fastcmh
        runfastcmh(folder=datafolder, saveToFile=TRUE, saveFolder=resultsfolder)

        resultsfile <- paste0(resultsfolder, "fastcmhresults.RData")
        cat("\n\n----------------------------------------------------\n")
        cat("\n")
        cat("Step 3: Inspect the results\n")
        cat("**(Run commands)**\n")
        cat("> load(\"", resultsfile, "\")\n", sep="")
        cat("> str(fastcmhresults)\n\n", sep="")
        fastcmhresults <- runfastcmh(folder=datafolder, saveToFile=FALSE)
        print(str(fastcmhresults))
        readkey();


    } else {

        cat("----------------------------------------------------\n\n")
        cat("Step 1: Creating data files for demo \n", sep="")
        cat("**(Run command)**\n")
        cat("> makefastcmhdata(folder=/data/, showOutput=TRUE)\n\n", sep="")

        cat("( files can be found in fastcmh/inst/extdata/ )\n\n")
        readkey();

        mylist <- makefastcmhdata(showOutput=TRUE, saveToList=TRUE)
        readkey();


        cat("\n\n----------------------------------------------------\n")
        cat("\n")
        cat("Step 2: Now running fastcmh, and saving results to list\n")
        cat("**(Run command)**\n")
        cat("> fastcmhresults <- runfastcmh(folder=/data/, saveToFile=FALSE)\n\n", sep="")
        fastcmhresults <- runfastcmh(folder="/data/", saveToFile=FALSE)
        readkey();

        cat("\n\n----------------------------------------------------\n")
        cat("\n")
        cat("Step 3: Inspect the results\n")
        cat("**(Run commands)**\n")
        cat("> str(fastcmhresults)\n\n", sep="")
        print(str(fastcmhresults))
        readkey();
    
    }

    cat("\n\n----------------------------------------------------\n")
    cat("Step 4: Or, access the results directly\n")
    cat("\n\n")
    cat("The data frame of significant intervals (filtered) can be found using:\n")
    cat("**(Run command)**\n")
    cat("> fastcmhresults$sig\n", sep="")
    cat("\n")
    print(fastcmhresults$sig)

    cat("\n")
    cat("In other words, ", sep="") 
    if (nrow(fastcmhresults$sig)==1){
        cat("the only ", sep="")
    } else { 
        cat("the first", sep="")
    }
    cat("significant interval is [", as.numeric(fastcmhresults$sig[1]), ", ", as.numeric(fastcmhresults$sig[2]), "],\n", sep="")
    cat("and the confounded interval [200, 203] is not found/is ignored.\n")


    cat("\n\n")
    cat("====================================================\n")
    cat("End of Demo\n")
    cat("====================================================\n\n")

}

