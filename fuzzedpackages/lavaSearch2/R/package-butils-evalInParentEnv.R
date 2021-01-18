### evalInParentEnv.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt  5 2017 (10:48) 
## Version: 
## last-updated: sep 28 2018 (12:05) 
##           By: Brice Ozenne
##     Update #: 12
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

#' @title Find Object in the Parent Environments
#' 
#' @description Search an object in the parent environments. For internal use.
#' 
#' @param name [character] the name of the object to get.
#'
#' @concept extractor
#' @keywords internal
evalInParentEnv <- function(name){
  
  ### ** find parent environments
  frames <- sys.status()
  all.frames <- sapply(1:length(frames$sys.frames), function(x){identical(parent.frame(x),globalenv())})
  index.parents <- which(all.frames==FALSE)
  n.parents <- length(index.parents)
  
  ### ** look in parent environments
  iParent <- 1
  res <- NULL
  cv <- FALSE
  while((iParent <= n.parents) && (cv == FALSE)){ # iParent <- 1
    
    res <- try(eval(name, envir = parent.frame(iParent)), silent = TRUE)
      
      if("try-error" %in% class(res)){
        iParent <- iParent + 1     
      }else{
        cv <- TRUE
      }
    
  }
  
  return(res)
}


#----------------------------------------------------------------------
### evalInParentEnv.R ends here
