#' Gets the regionGWAS string
#'
#' A getter for the global \code{regionGWAS} value, a string
#' for the mode parameter.
#'
#' @keywords internal 
getRegionGWASString  <- function(){
    return( get("regionGWASString", envir=CASMAPenv) )
}



#' Gets the higherOrderEpistasis string
#'
#' A getter for the global \code{higherOrderEpistasis} value, a string
#' for the mode parameter.
#'
#' @keywords internal 
getHigherOrderEpistasisString  <- function(){
    return( get("higherOrderEpistasisString", envir=CASMAPenv) )
}


#' Gets the minModeLength
#'
#' A getter for the global \code{minModeLength} value, a string
#' for the mode parameter.
#'
#' @keywords internal 
getMinModeLength <- function(){
    return( get("minModeLength", envir=CASMAPenv) )
}





#' Checks if substring is part of regionGWAS
#'
#' Using\code{grepl} to compare strings, ignoring case.
#'
#' @param x The string which will be compared to 'regionGWAS'
#'
#' @details
#' Uses \code{grepl} to search for exact match. Case will be ignored.
#'
#' @return \code{TRUE} if the string is a substring of 'regionGWAS',
#'         otherwise returns \code{FALSE}.
#'
#' @keywords internal 
isRegionGWASString <- function(x){
    isMatch <- grepl(x, getRegionGWASString(), ignore.case=T)
    return(isMatch)
}


#' Checks if substring is part of higherOrderEpistasis
#'
#' Using grep to search through vector of strings
#'
#' @param x The string which will be compared to 'higherOrderEpistasis'
#'
#' @details
#' Uses \code{grep} to search for exact match.
#'
#'
#' @return \code{TRUE} if the string is a substring of 'higherOrderEpistasis',
#'         otherwise returns \code{FALSE}.
#'
#' @keywords internal 
isHigherOrderEpistasisString <- function(x){
    isMatch <- grepl(x, getHigherOrderEpistasisString(), ignore.case=T)
    return(isMatch)
}


#' Error message for mode
#'
#' Return the appropriate error message for incorrect mode input
#'
#' @keywords internal 
modeErrorMessage <-function(){
    message <- paste0("'mode' needs to be specified as a character string, ")
    message <- paste0(message, "either '", getRegionGWASString(), "' or ")
    message <- paste0(message, "'", getHigherOrderEpistasisString(), "'.")
    return(message)
}


#' Minimum length of modeb
#'
#' Gets the minimum mode character length (should be 3)
#'
#' @keywords internal 
getMinModeLength <-function(){
    return( get("minModeLength", envir=CASMAPenv) )
}


#' Checks mode string is long enough
#'
#' Checks mode string is at least minimum length
#'
#' @keywords internal 
modeNeedsMoreChars <-function(mode){
    return( nchar(mode) < getMinModeLength() )
}


#' Error message for mode, if too short
#'
#' Return the appropriate error message for incorrect mode input
#'
#' @keywords internal 
modeLengthErrorMessage <-function(){
    message <- paste0("'mode' needs be specified as a character string, ")
    message <- paste0(message, "either '", getRegionGWASString(), "' or ")
    message <- paste0(message, "'", getHigherOrderEpistasisString(), "',")
    mc <- getMinModeLength()
    message <- paste0(message, "and should be at least ") 
    message <- paste0(message, getMinModeLength(), "long.")
    return(message)
}


#' Get the function name
#'
#' Uses \code{match.call} and \code{as.character}.
#' @keywords internal
getParentFunctionName <- function(){
    name <- sys.call(-1)
    return(as.character(name))
}
