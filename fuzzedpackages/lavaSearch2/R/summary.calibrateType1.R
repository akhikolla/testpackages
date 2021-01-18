### summary.calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr 23 2018 (12:58) 
## Version: 
## Last-Updated: jul 16 2018 (16:37) 
##           By: Brice Ozenne
##     Update #: 71
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * method summary.modelsearch2
#' @title Display the Type 1 Error Rate 
#' @description Display the type 1 error rate from the simulation results.
#'
#' @param object output of the \code{calibrateType1} function.
#' @param alpha [numeric, 0-1] the confidence levels.
#' @param robust [character] should the results be displayed for both model-based and robust standard errors (\code{TRUE}),
#' only  model-based standard error (\code{FALSE}), or only robust standard error (\code{"only"})?
#' @param type [character] should the type 1 error rate be diplayed (\code{"type1error"}) or the bias (\code{"bias")}.
#' @param digits [integer >0] the number of decimal places to use when displaying the summary.
#' @param log.transform [logical] should the confidence intervals be computed on the logit scale.
#' @param print should the summary be printed in the terminal.
#' @param ... [internal] only used by the generic method.
#' 
#' @method summary calibrateType1
#' @export
summary.calibrateType1 <- function(object, robust = FALSE, type = "type1error",
                                   alpha = 0.05, log.transform = TRUE,
                                   digits = 5, print = TRUE, ...){


    type <- match.arg(type, c("type1error","bias"))
    if(is.logical(robust)){robust <- as.character(robust)}
    robust <- match.arg(robust, c("TRUE","FALSE","only"))
    
    ## ** compute bias
    if(type == "bias"){
        stop("not implemented yet!")
    }
    
    ## ** compute type 1 error
    if(type == "type1error"){
        dfLong <- melt(object$p.value,
                       measure.vars = grep("^p.",names(object$p.value),value = TRUE),
                       value.name = "p.value",
                       variable.name = "method")
        dfS <- stats::aggregate(dfLong$p.value,
                                by = list(n = dfLong$n, method = dfLong$method, link = dfLong$link),
                                FUN = function(x){c(n.rep = length(x), type1error = mean(x<=alpha, na.rm = TRUE))},
                                simplify = FALSE)
        dfS <- cbind(dfS[,c("n","method","link")],
                     do.call(rbind,dfS[,"x"]))

        ## binom::binom.confint(x = round(dfS$type1error*dfS$n.rep), dfS$n.rep, method = "logit")
        
        logit.p <- log(dfS$type1error/(1-dfS$type1error))

        if(log.transform){
            ## delta method: ln(x/(1-x))' = ln(x)' - ln(1-x)' = 1/x + 1/(1-x) = 1/(x(1-x))
            var.logit.p <- 1/(dfS$n.rep*dfS$type1error*(1-dfS$type1error))

            ## on logit scale
            dfS$ci.inf <- logit.p + stats::qnorm(0.025) * sqrt(var.logit.p)
            dfS$ci.sup <- logit.p + stats::qnorm(0.975) * sqrt(var.logit.p)

            ## on original scale
            dfS$ci.inf <- 1/(1+exp(-dfS$ci.inf))
            dfS$ci.sup <- 1/(1+exp(-dfS$ci.sup))
        }else{
            dfS$ci.inf <- dfS$type1error + stats::qnorm(0.025) * sqrt(dfS$type1error*(1-dfS$type1error)/dfS$n.rep)
            dfS$ci.sup <- dfS$type1error + stats::qnorm(0.975) * sqrt(dfS$type1error*(1-dfS$type1error)/dfS$n.rep)
        }        

        ## names
        rownames(dfS) <- NULL
        dfS$correction <- sapply(as.character(dfS$method),switch,
                                 p.Ztest = "Gaussian approx.",
                                 p.Satt = "Satterthwaite approx.",
                                 p.SSC = "small sample correction",
                                 p.KR = "Satterthwaite approx. with small sample correction",
                                 p.robustZtest = "Gaussian approx.",
                                 p.robustSatt = "Satterthwaite approx.",
                                 p.robustSSC = "small sample correction",
                                 p.robustKR = "Satterthwaite approx. with small sample correction")
        
        dfS$statistic <- sapply(as.character(dfS$method),switch,
                                p.Ztest = "Wald",
                                p.Satt = "Wald",
                                p.SSC = "Wald",
                                p.KR = "Wald",
                                p.robustZtest = "robust Wald",
                                p.robustSatt = "robust Wald",
                                p.robustSSC = "robust Wald",
                                p.robustKR = "robust Wald")

        ## display
        if(print){
            seqN <- unique(dfS$n)
            seqRep <- setNames(dfS$n.rep[duplicated(dfS$n) == FALSE],seqN)

            vecTrans <- c("Gaus" = "Gaussian approx.",
                          "Satt" = "Satterthwaite approx.",
                          "SSC" = "small sample correction",
                          "SSC + Satt" = "Satterthwaite approx. with small sample correction")

            ls.print <- lapply(seqN, function(iN){ # iN <- 73
                df.tempo <- dfS[dfS$n==iN, c("link","statistic","correction","type1error","ci.inf","ci.sup")]
                if(!is.null(digits)){
                    df.tempo$type1error <- round(df.tempo$type1error, digits = digits)
                    df.tempo$ci.inf <- round(df.tempo$ci.inf, digits = digits)
                    df.tempo$ci.sup <- round(df.tempo$ci.sup, digits = digits)
                }
                df.tempo$CI <- paste0("[",df.tempo$ci.inf," ; ",df.tempo$ci.sup,"]")
                df.tempo$ci.inf <- NULL
                df.tempo$ci.sup <- NULL

                if(robust == "FALSE"){
                    df.tempo <- df.tempo[df.tempo$statistic=="Wald",,drop=FALSE]
                }else if(robust == "only"){
                    df.tempo <- df.tempo[df.tempo$statistic=="robust Wald",,drop=FALSE]
                }
                
                df.tempo <- df.tempo[order(df.tempo$link,df.tempo$statistic),,drop=FALSE]
                df.tempo$statistic[duplicated(cbind(df.tempo$link,df.tempo$statistic))] <- ""
                df.tempo$link[duplicated(df.tempo$link)] <- ""

                rownames(df.tempo) <- NULL
                df.tempo$correction <- factor(df.tempo$correction,
                                              levels = as.character(vecTrans),
                                              labels = names(vecTrans))
                ## df.tempo$correction <- NULL
                
                return(df.tempo)
            })
            names(ls.print) <- as.character(seqN)
       
            cat("Estimated type 1 error rate [95% confidence interval] \n")
            lapply(seqN, function(iN){
                cat("  > sample size: ",iN," | number of simulations: ",seqRep[as.character(iN)],"\n",sep="")
                print(ls.print[[as.character(iN)]],row.names = FALSE)
            })
            cat("\n")
            cat("Corrections: Gaus = Gaussian approximation \n",
                "             SSC  = small sample correction \n",
                "             Satt = Satterthwaite approximation \n",sep="")
            
        }        
    }

    return(invisible(dfS))
}
######################################################################
### summary.calibrateType1.R ends here
