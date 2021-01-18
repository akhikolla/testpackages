### autoplot.calibrateType1.R --- 
##----------------------------------------------------------------------
## Author: Brice Ozenne
## Created: apr  5 2018 (13:20) 
## Version: 
## Last-Updated: maj 23 2018 (09:53) 
##           By: Brice Ozenne
##     Update #: 29
##----------------------------------------------------------------------
## 
### Commentary: 
## 
### Change Log:
##----------------------------------------------------------------------
## 
### Code:

## * Documentation
##' @title Graphical Display of the Bias or Type 1 Error
##' @description Graphical display of the bias or type 1 error
##' for the output of \code{\link{calibrateType1}}.
##' @name autoplot_calibrateType1
##' 
##' @param object output of \code{\link{calibrateType1}}.
##' @param type [character] if type equals \code{"bias"} the bias will be displayed.
##' Otherwise if it equals \code{"type1error"} the type 1 error will be displayed.
##' @param plot [logical] should the plot be displayed?
##' @param color.threshold [character] the color for the line representing the expected value(s).
##' @param type.bias [character] if type.bias equals \code{"absolute"} the absolute bias will be used.
##' Otherwise if it equals \code{"relative"} the relative bias will be used.
##' Only relevant when type equals \code{"bias"}.
##' @param alpha [numeric, 0-1] the significance threshold to consider.
##' Only relevant when type equals \code{"type1error"}.
##' @param nrow.legend [integer, >0] the number of rows for the legend.
##' Only relevant when type equals \code{"type1error"}.
##' @param name2label [named character vector] the label for the legend.
##' The vector should contain the method names (see details).
##' Only relevant when type equals \code{"type1error"}.
##' @param color [character vector] a vector of colours to be used to color the lines.
##' Only relevant when type equals \code{"type1error"}.
##' @param keep.method [character vector] the methods names for which the type 1 error should be displayed.
##' Only relevant when type equals \code{"type1error"}.
##' @param ... [internal] Only used by the generic method.
##' 
##' @details Method names:
##' \itemize{
##' \item \code{p.Ztest}
##' \item \code{p.Satt}
##' \item \code{p.KR}
##' \item \code{p.robustZtest}
##' \item \code{p.robustSatt}
##' \item \code{p.robustKR}
##' }
##' @return An list containing:
##' \itemize{
##' \item plot: a ggplot object.
##' \item data: the dataset used to generate the ggplot object.
##' }
##' 


## * plotBias.calibrateType
##' @rdname autoplot_calibrateType1
##' @method autoplot calibrateType1
##' @export
autoplot.calibrateType1 <- function(object, type = "bias", plot = TRUE, color.threshold = "red",
                                    type.bias = "absolute",
                                    alpha = 0.05, nrow.legend = NULL, name2label = NULL, color = NULL, keep.method = NULL,
                                    ...){

    type <- match.arg(type, choices = c("bias","type1error"))

    ## ** display bias
    if(type == "bias"){
        ## *** data
        type.bias <- match.arg(type.bias, choices = c("absolute","relative"))

        if(type.bias == "absolute"){
            object$bias$bias.ML <- object$bias$estimate.truth-object$bias$estimate.ML
            object$bias$bias.MLcorrect <- object$bias$estimate.truth-object$bias$estimate.MLcorrect
        }else if(type.bias == "relative"){
            object$bias$bias.ML <- object$bias$estiamte.truth-object$bias$estiamte.ML
            object$bias$bias.MLcorrect <- object$bias$estimate.truth-object$bias$estimate.MLcorrect
        }

        df1 <- cbind(object$bias[,c("n","rep","name","type","bias.ML")], corrected = FALSE)
        names(df1)[5] <- "bias"
        df2 <- cbind(object$bias[,c("n","rep","name","type","bias.MLcorrect")], corrected = TRUE)
        names(df2)[5] <- "bias"

        df.gg <- rbind(df1,df2)

        ## *** display
        gg <- ggplot2::ggplot() + ggplot2::geom_boxplot(data = df.gg, aes_string(x = "type", y = "bias")) 
        gg <- gg + ggplot2::facet_wrap(corrected ~ n, labeller = ggplot2::label_both)
        df.line <- data.frame(x = c(-Inf,Inf), y = 0)
        gg <- gg + ggplot2::geom_line(data = df.line, aes_string(x = "x", y = "y"), color = color.threshold)
        
    }

    ## ** display type 1 error
    if(type == "type1error"){
        df.gg <- summary(object, alpha = alpha, type = type, display = FALSE)
        
        ## *** display
        if(is.null(keep.method)){
            keep.method <- as.character(unique(df.gg$method))
        }
        if(is.null(name2label)){
            name2label <- setNames(unique(paste0(df.gg$statistic,", ",df.gg$correction)),unique(df.gg$method))
        }
        if(is.null(color)){
            ## from ggthemes::colorblind_pal()(8)
            color <- c("#000000","#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7")
        }

        label.method <- name2label[keep.method]
        n.method <- length(keep.method)

       y.range <- range(c(alpha,df.gg$type1error))
        gg <- ggplot2::ggplot(df.gg, aes_string(x = "n", y = "type1error", group = "method", color = "method", shape = "method"))
        gg <- gg + ggplot2::geom_point(size = 3) + ggplot2::geom_line(size = 2)
        gg <- gg + ggplot2::facet_grid(~link, labeller = ggplot2::label_parsed)
        gg <- gg + ggplot2::geom_abline(intercept = alpha, slope = 0, color = color.threshold)
        gg <- gg + ggplot2::xlab("sample size")
        gg <- gg + ggplot2::ylab("type 1 error rate")
        gg <- gg + ggplot2::theme(legend.position = "bottom")
        gg <- gg + ggplot2::coord_cartesian(ylim = y.range)

        if(!is.null(nrow.legend)){
            gg <- gg + ggplot2::guides(color = ggplot2::guide_legend(nrow=nrow.legend,byrow=TRUE))
        }
        
        gg <- gg + ggplot2::scale_color_manual("",
                                               breaks = keep.method,
                                               label = label.method,
                                               values = color[1:n.method])
        gg <- gg + ggplot2::scale_shape_manual("",
                                               breaks = keep.method,
                                               label = label.method,
                                               values = seq(15,by=1,length=n.method))
        
         
    }
    
    ## ** export
    if(plot){
        print(gg)
    }
    rownames(df.gg) <- NULL
    return(invisible(list(plot = gg,
                          data = df.gg)))
}


######################################################################
### autoplot.calibrateType1.R ends here
