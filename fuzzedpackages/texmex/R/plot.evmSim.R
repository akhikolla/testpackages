#' Plots for evmSim objects
#'
#' This function produces diagnostic plots for the Markov chains used
#' to simulate from the posterior distributions for the model
#' parameters. If the chains have converged on the posterior
#' distributions, the trace plots should look like
#' "fat hairy caterpillars" and their cumulative means should converge
#' rapidly. Moreover, the autocorrelation functions should converge
#' quickly to zero.
#'
#' @param x an object of class \code{evmSim}
#' @param which.plots Which plots to produce. Option 1 gives kernel
#'     density estimates, 2 gives traces of the Markov chains with
#'     superimposed cumulative means, 3 gives autocorrelation
#'     functions.
#' @param chain Which chain to use in the trace and ACF plots. Only used if
#'   more than one chain was run. Defaults to \code{chain = 1}. If you ran
#'   multiple chains, you'll want to look at them all.
#' @param density.adjust In \code{plot} method for class
#'     \code{evmSim}.  Passed into \code{density}. Controls the amount
#'     of smoothing of the kernel density estimate.
#' @param print.seed Whether or not to print the seed used in the
#'     simulations, or to annotate the plots with it.
#' @param ... ignored
#' @seealso \code{\link{evm}}
#' @seealso \code{\link[stats]{density}}
#' @export
plot.evmSim <-
  function(x, which.plots=1:3, chain = 1, density.adjust=2, print.seed = FALSE , ...){

    if (length(x$chains) > 1){
      msg <- paste0("Trace and ACF plots for chain ", chain, " only.")
      if (chain == 1){
        msg <- paste(msg, "Use the 'chain' argument to specify another chain.")
      }
      message(msg)
    }

    varnames <- names(coef(x))

    if (1 %in% which.plots){
      kdes <- list()

      for (i in 1:length(varnames)){
        kdes[[i]] <- density(x$param[, i], n=200, adjust=density.adjust)
      }

      for(i in 1:length(kdes)){
        plot(kdes[[i]], type = "l" , xlab = varnames[i] , ylab="Density" ,
             main = paste("Posterior for", varnames[i]))
        polygon(c(rev(kdes[[i]]$x), kdes[[i]]$x),
                c(rep(0, length(kdes[[i]]$y)) , kdes[[i]]$y),
                col=6)
        m <- mean(x$param[, i])
        ym <- kdes[[i]]$y[abs(kdes[[i]]$x - m) == min(abs(kdes[[i]]$x - m))]
        segments(x0=m, y0=0, y1=ym , x1=m)
        if (print.seed)
          title(sub = paste(c("Seed: ", x$seed) , collapse = " "),
                adj=0)
      }
    } # Close if (1 %in% which plots

    if(2 %in% which.plots){
      # Trace plots and moving averages
      for (i in 1:length(varnames)){
        plot(x$chains[[chain]] [, i], type = "l" , xlab = "Step number", ylab = paste(varnames[i], "& cumulative mean"))
        lines(cumsum(x$chains[[chain]][ , i]) / (1:nrow(x$chains[[chain]])), lwd=3, col=8)
        axis(4)
        abline(v = x$burn, lty=4, lwd=3, col=2)
        if (print.seed)
          title(sub = paste(c("Seed: ", x$seed) , collapse = " "), adj=0)
      }
    }

    if(3 %in% which.plots){
      # Plot ACFs

      x$chains <- x$chains[chain]
      x <- thinAndBurn(x)


      for(i in 1:length(varnames)){
        acf(x$param[ , i] , main = paste("ACF for", varnames[i], "\n(thin =",x$thin,")"))
        if (print.seed)
          title(sub = paste(c("Seed: ", x$seed) , collapse = " "), adj=0)
      }
    }
    invisible()
  }

