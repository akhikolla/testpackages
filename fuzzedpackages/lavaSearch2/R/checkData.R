### checkData.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: okt 26 2017 (14:25) 
## Version: 
## last-updated: aug  7 2018 (10:54) 
##           By: Brice Ozenne
##     Update #: 38
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * documentation - checkData
#' @title Check that Validity of the Dataset
#' @description Check whether the dataset can be used to fit the \code{lvm} object.
#' 
#' @name checkData
#' 
#' @param object a \code{lvm} object.
#' @param data [data.frame] the dataset used to obtain the object.
#' @param trace [logical] when \code{TRUE}, the outcome of the check will be displayed. 
#'
#' @return Invisible  \code{TRUE} or \code{FALSE}.
#' 
#' @examples 
#' m <- lvm()
#' regression(m) <- c(y1,y2,y3)~u
#' regression(m) <- u~x
#' latent(m) <- ~u
#'
#' d <- lava::sim(m,1e2)
#'
#' try(checkData(m, data = d)) # return an error
#' 
#' checkData(m, data = d[,-4])
#' 
#' try(checkData(m, data = d[,-(3:4)])) # return an error
#'
#' @concept diagnostic
#' @export
`checkData` <-
  function(object, data, trace) UseMethod("checkData")

## * checkData.lvm
#' @rdname checkData
#' @export
checkData.lvm <- function(object, data, trace = TRUE){

    ## ** normalize arguments
    data <- as.data.frame(data)
        
    ## ** check missing names
    vars <- lava::vars(object)
    latent <- lava::latent(object)    
    missingVars <- vars[vars %in% names(data) == FALSE]

    if(length(latent) == 0){
        
        if(length(missingVars)>0){
            if(trace){
                cat("Missing variable in data: ",paste(missingVars, collapse = " "),"\n", sep ="")
            }
            return(invisible(FALSE))
        }else{
            if(trace){
                cat("No issue detected \n")
            }
            return(invisible(TRUE))
        }
        
    }else{
        
        if(!identical(sort(latent),sort(missingVars))){
            if(trace){
                cat("Wrong specification of the latent variables \n",
                    "latent variables according to the LVM: ",paste(latent, collapse = " "),"\n",
                    "missing variables in data: ",paste(missingVars, collapse = " "),"\n", sep = "")
            }
            return(invisible(FALSE))
        }else{
            if(trace){
                cat("No issue detected \n")
            }
            return(invisible(TRUE))
        }
        
    }
}


#----------------------------------------------------------------------
### checkData.R ends here
