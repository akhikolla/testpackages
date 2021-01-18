### summary.modelsearch2.R --- 
#----------------------------------------------------------------------
## author: Brice Ozenne
## created: aug 30 2017 (10:46) 
## Version: 
## last-updated: jun 27 2019 (14:21) 
##           By: Brice Ozenne
##     Update #: 115
#----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
#----------------------------------------------------------------------
## 
### Code:

## * method summary.modelsearch2
#' @title summary Method for modelsearch2 Objects
#' @description summary method for modelsearch2 objects.
#'
#' @param object output of the \code{modelsearch2} function.
#' @param print should the summary be printed in the terminal.
#' @param ... [internal] only used by the generic method.
#' 
#' @method summary modelsearch2
#'
#' @details
#' The column \code{dp.Info} contains the percentage of extended models (i.e. model with one additional link)
#' for which the information matrix evaluated at the value of the parameters of the initial model is non positive definie.
#' 
#' @export
summary.modelsearch2 <- function(object, print = TRUE, ...){

    p.value <- NULL # [:for CRAN check] subset
    
    ## ** extract data from object
    n.step <- nStep(object)
    tableTest <- do.call(rbind,lapply(object$sequenceTest, function(iTest){
        out <- iTest[which.max(iTest$statistic), c("link","statistic","p.value","adjusted.p.value","dp.Info","selected","nTests")]
        out[,"dp.Info"] <- mean(iTest[,"dp.Info"], na.rm = TRUE)
        return(out)
    }))
    n.selected <- sum(tableTest$selected)

    keep.cols <- c("link","nTests","statistic","adjusted.p.value")
    ## ** output
    out <- list()
    out$message.pre <- paste0("Sequential search for local dependence using the score statistic \n")
    if(n.selected==0){
        out$message.pre <- c(out$message.pre,
                             "The variable selection procedure did not retain any variable \n")
    }else{
        out$message.pre <- c(out$message.pre,
                             paste0("The variable selection procedure retained ",n.selected," variable",
                                    if(n.selected>1){"s"},":\n")
                             )     
    }

    out$table <- tableTest
    rownames(out$table) <- NULL
    out$message.post <- paste0("Confidence level: ",1-object$alpha," (two sided, adjustement: ",object$method.p.adjust,")\n")  

    ## ** display
    if(print){
        cat(out$message.pre,sep="")
        print(out$table)
        cat(out$message.post,sep="")
        if(any(na.omit(out$table[,"dp.Info"])<1)){
            cat("WARNING: some of the score tests could not be correctly computed, probably because extended models are not all identifiable\n",
                "        consider using the argument \'link\' to specify only identifiable models \n")
        }
    }
    
    ## ** export
    return(invisible(out))
}



#----------------------------------------------------------------------
### summary.modelsearch2.R ends here
