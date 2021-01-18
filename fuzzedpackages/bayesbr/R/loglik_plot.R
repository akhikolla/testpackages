#' @title Plot Chain of Loglik Using GGplot
#' @aliases loglikPlot
#' @name loglikPlot
#' @description The function receives a vector containing the model's loglik chain and displays it through a graph using the GGplot package, through this graph it is possible to see if the model has converged.
#' @usage loglikPlot(loglik)
#' @param loglik A vector with the estimated model loglik chain
#' @seealso \code{\link{logLik.bayesbr}},\code{\link{envelope}}
#' @examples
#' data("CarTask", package = "bayesbr")
#'
#' bbr = bayesbr(probability~task + NFCCscale,data = CarTask,
#'              iter = 100, mean_betas = c(1, 0.5,1.2))
#'loglik = bbr$loglik
#'
#'loglikPlot(loglik)
#'@export
loglikPlot = function(loglik){
  n = length(loglik)
  ggplot() + geom_line(aes(y = loglik,x=1:n)) + ylab("Log-Likelihood") +
    xlab("Iterations without warmups") + ggtitle("Loglik's chains of the estimated model")
}
