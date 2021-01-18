#' @title Calculate multiple-test cut-off
#' @description Calculate the log-likelihood, AIC, or AICc cut-off necessary for type-one error to reach acceptable levels
#' @param phy An object of class \code{phylo} (see \pkg{ape}).
#' @param n Number of simulations
#' @param mc.cores Number of cores for parallel processing for linux-type systems (not applicable to Windows)
#' @param model Evolutionary model, typically "tm1", "tm2", or "timeSlice", which is used to test the empirical data
#' @param measure Measure used to summarise the model. One of "lnL" (log-likelihood), "AIC", or "AICc"
#' @param alpha.error Target for the desired type-one error rate for the model (default 0.05)
#' @param ... Arguments to be passed to \code{\link{transformPhylo.ML}} which should be identical to the model applied to empirical data
#' @return The cut-off requred to produce an type-one error rate equal to quantile.cut.off (default = 0.05) when data are simulated under Brownian motion, and these data are analysed under the appropriate model.
#' @author Mark Puttick
#' @seealso \code{\link{transformPhylo.ML}}, \code{\link{transformPhylo.ll}}, \code{\link{transformPhylo}}, \code{\link{transformPhylo.MCMC}}
#' @examples
#' data(anolis.tree)
#' set.seed(393)
#' # calculated necessary AICc cut-off to reduce type-one error to 5% 
#' # for a timeSlice model with a split at 30Ma (only 5 simulations used,
#' # it's recommend to use 1000 for analyses)
#' calcCutOff(anolis.tree, n=5, model="timeSlice", splitTime=30)
#' @export

calcCutOff <- function(phy, n=1000, mc.cores=1, model, measure="AICc", alpha.error=0.05, ...) {

	quiet <- function(x) { 
	  sink(tempfile()) 
	  on.exit(sink()) 
	  invisible(force(x)) 
	} 
	
	sim.data <- transformPhylo.sim(phy, n=n, model="bm")
	y.out <- quiet(
				parallel::mclapply(1:n, mc.cores=mc.cores, function(xx) {
					model.fit <- transformPhylo.ML(as.matrix(sim.data[, xx]), phy, model = model, ...)
					if(model == "timeSlice") {
						model.fit <- model.fit[[1]]
						model.fit[, measure]
					} else {
						model.fit <- c(transformPhylo.ML(as.matrix(sim.data[, xx]), phy, model = "BM")[[measure]], model.fit[[measure]])
						model.fit
            		}
            	}	
			)
		)
	
	y.out <- matrix(unlist(y.out), nrow=2)
	aicc.add <- quantile(-apply(y.out, 2, diff), (1-alpha.error))
	return(aicc.add)
	}