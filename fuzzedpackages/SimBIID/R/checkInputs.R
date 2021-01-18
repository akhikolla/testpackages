## Checks input arguments for consistency
checkInput <- function(input, type = NULL, length = NA, 
        nrow = NA, ncol = NA, int = FALSE, naAllow = FALSE, 
        gt = NA, gte = NA, lt = NA, lte = NA, inSet = NULL, 
        uni = FALSE, ident = NULL
) {
    ## input is input object
    ## type is set of characters denoting is.*() type functions to test input
    ## length is numeric testing length of input (similarly for nrow and ncol)
    ## int is logical testing for integer
    ## naAllow determines whether NA values are allowed or not
    ## gt, gte, lt, lte define greater than, greater than or equal to
    ##          less than and less than or equal to respectively
    ## inSet specifies a set of values that all values in input must
    ##          be part of
    ## uni defines whether elements of input need to be unique
     
    ## extract name of original object
    inputName <- as.character(sys.call())[2]
    
    ## check type of input
    if(!is.null(type)) {
        type <- paste0("is.", type)
        for(i in 1:length(type)){
            output <- do.call(type[i], list(input))
            if(!output){
                stop(paste0(inputName, " object is not: ", gsub("is.", "", type), "\n"))
            }
        }
    }
    if(!is.na(length)) {
        output <- length(input) == length
        if(!output){
            stop(paste0(inputName, " object is not of length ", length))
        }
    }
    if(!is.na(nrow)) {
        output <- nrow(input) == nrow
        if(!output){
            stop(paste0(inputName, " object does not have ", nrow, " rows"))
        }
    }
    if(!is.na(ncol)) {
        output <- ncol(input) == ncol
        if(!output){
            stop(paste0(inputName, " object does not have ", ncol, " cols"))
        }
    }
    if(int) {
        output <- all(input %% 1 == 0)
        if(!output){
            stop(paste0(inputName, " object is not a vector of integers"))
        }
    }
    if(!naAllow){
        if(!is.function(input)) {
            output <- all(!is.na(input))
            if(!output){
                stop(paste0(inputName, " object contains NAs"))
            }
        }
    }
    if(!is.na(gt)) {
        output <- min(input[!is.na(input)]) > gt
        if(!output){
            stop(paste0(inputName, " elements must be greater than ", gt))
        }
    }
    if(!is.na(gte)) {
        output <- min(input[!is.na(input)]) >= gte
        if(!output){
            stop(paste0(inputName, " elements must be greater than or equal to ", gte))
        }
    }
    if(!is.na(lt)) {
        output <- max(input[!is.na(input)]) < lt
        if(!output){
            stop(paste0(inputName, " elements must be less than ", lt))
        }
    }
    if(!is.na(lte)) {
        output <- max(input[!is.na(input)]) <= lte
        if(!output){
            stop(paste0(inputName, " elements must be less than or equal to ", lte))
        }
    }
    if(!is.null(inSet)) {
        output <- all(input[!is.na(input)] %in% inSet)
        if(!output){
            stop(paste0(inputName, " elements must be in (", paste(inSet, collapse = ", "), ")"))
        }
    }
    if(!is.null(ident)) {
        output <- identical(input, ident)
        if(!output){
            stop(paste0(inputName, " elements must be identical to (", paste(ident, collapse = ", "), ")"))
        }
    }
    if(uni) {
        output <- !any(duplicated(input[!is.na(input)]))
        if(!output){
            stop(paste0(inputName, " elements must be unique"))
        }
    }
}
