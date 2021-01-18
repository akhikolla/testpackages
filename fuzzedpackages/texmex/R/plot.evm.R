#' Plots for evmOpt objects
#'
#' Various plots for \code{evmOpt} objects. These differ depending on
#' whether or not there are covariates in the model.  If there are no
#' covariates then the diagnostic plots are PP- and QQ-plots, a return
#' level plot (produced by \code{plotrl.evmSim}) and a histogram of
#' the data with superimposed density estimate.  These are all
#' calculated using the data on the original scale. If there are
#' covariates in the model then the diagnostics consist of PP- and QQ-
#' plots calculated by using the model residuals (which will be
#' standard exponential devaiates under the GPD model and standard
#' Gumbel deviates under the GEV model), and plots of residuals versus
#' fitted model parameters.
#' 
#' The PP- and QQ-plots show simulated pointwise tolerance intervals.
#' The region is a \eqn{100(1 - \alpha)\%}{100(1-alpha)\%} region based
#' on \code{nsim} simulated samples.
#'
#' @param x an object of class \code{evmOpt}
#' @param main titles for diagnostic plots. Should be a vector of
#'     length 4, with values corresponding to the character strings to
#'     appear on the titles of the pp, qq, return level, and density
#'     estimate plots respectively.
#' @param xlab As for \code{main} but labels for x-axes rather than
#'     titles.
#' @param nsim The number of replicates to be simulated to produce the
#'     simulated tolerance intervals.
#' @param alpha A \eqn{100(1 - \alpha)\%}{100(1 - alpha)\%} simulation
#'     envelope is produced.
#' @param ... FIXME
#' @seealso \code{\link{evm}}
#' @export
plot.evmOpt <-
function(x, main=rep(NULL,4), xlab=rep(NULL,4), nsim=1000, alpha=.05, ...){
    if (!missing(main)){
        if (length(main) != 1 & length(main) != 4){
            stop("main should have length 1 or 4")
        }
        else if (length(main) == 1){ main <- rep(main, 4) }
    }

    plot(ppevm(x, nsim=nsim, alpha=alpha), main=main[1], xlab=xlab[1])
    plot(qqevm(x, nsim=nsim, alpha=alpha), main=main[2], xlab=xlab[2])
        
    if (all(sapply(x$data$D,ncol) == 1)){
        plotrl.evmOpt(x, main=main[3], xlab=xlab[3], smooth=FALSE, ...)
        plot(hist.evmOpt(x), main=main[4], xlab=xlab[4])
    }
    else { # Covariates in the model
        np <- length(x$data$D)
        lp <- predict(x,type="lp", unique.=FALSE)$obj$link
        Which <- as.logical(apply(lp[,1:np],2,var)) # identifies which cols have covariates

        for(i in (1:length(x$data$D))[Which]){
          ParName <- names(x$data$D[i])
          plot(lp[,i],resid(x),main=paste("Residuals vs fitted",ParName),xlab=paste("Fitted",ParName),ylab="Residuals")
          panel.smooth(lp[,i], resid(x), col.smooth=2)
        }
    }
    
    invisible()
}

