### summary.glht2.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: maj  2 2018 (09:20) 
## Version: 
## Last-Updated: maj  2 2018 (16:58) 
##           By: Brice Ozenne
##     Update #: 34
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

summary.glht2 <- function(object, ...){

    class(object) <- setdiff(class(object), "glht2")
    output <- summary(object, ...)
    class(output) <- append("summary.glht2",class(output))
    
    return(output)
}

print.summary.glht2 <- function(object, ...){

    class(object) <- setdiff(class(object), "summary.glht2")
    output <- utils::capture.output(print(object))
    
    txt.robust <- switch(as.character(object$robust),
                         "TRUE" = "Robust standard errors",
                         "FALSE" = "Model-based standard errors"
                         )
    txt.correction <- switch(as.character(object$bias.correct),
                             "TRUE" = " corrected for small sample bias",
                             "FALSE" = ""
                             )
    output[length(output)] <- paste0("(",txt.robust,txt.correction,")\n")

    cat(paste0(output,collapse = "\n"))
    
    return(invisible(object))
}


######################################################################
### summary.glht2.R ends here
