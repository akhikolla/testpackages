#'@name envelope
#'@aliases envelope
#'@title Envelope Graph: Residuals vs Half-Normal Values
#'@description A graph showing the absolute values of the residuals ordered against the quantiles of simulations of the half-normal distribution.
#'@usage
#'envelope(x, sim = 1000, conf = 0.95, resid.type = c("",
#' "quantile", "sweighted","pearson","ordinary"))
#'@param x an object of the class \emph{bayesbr}, containing the list returned from the \code{\link{bayesbr}} function.
#'@param sim a positive integer containing the number of simulations of the half-normal distribution.
#'@param conf a probability containing the confidence level for the quantiles made under the half-normal samples.
#'@param resid.type the residual type that will be used in the graph
#'@details Atkinson (1985) proposed to use quantiles from a simulated population of the halfnormal distribution, this is used because (blablabla read the book, right). From the distribution of the absolute values of the residual in the graph, it is possible to measure the quality of the model estimation.
#'@return  A graph showing the absolute values of the residuals ordered against the quantiles of simulations of the half-normal distribution.
#'@references
#' Atkinson, A. C. (1985). Plots, transformations, and regression: An introduction to graphical methods of diagnostic regression analysis. \emph{Oxford: Clarendon Press}.
#'@seealso \code{\link{residuals.bayesbr}}, \code{\link{loglikPlot}}, \code{\link{bayesbr}}
#'@examples
#'data("CarTask", package = "bayesbr")
#'
#'bbr = bayesbr(probability~task + NFCCscale, iter = 100,
#'             data=CarTask, mean_betas = c(1, 0.5,1.2),variance_betas=10)
#'
#'envelope(bbr,sim = 100, conf=0.9, resid.type="quantile")
#'\donttest{
#'envelope(bbr,sim = 1000, conf=0.99, resid.type="ordinary")}
#' @export
envelope = function(x,sim = 1000,conf = 0.95,resid.type  = c("","quantile","sweighted", "pearson","ordinary")){
  resid.type = match.arg(resid.type)
  if(resid.type==""){
    resid.type = x$residuals.type
  }
  res = residuals.bayesbr(x,resid.type)
  if(!is.numeric(sim)){
    warning("The parameter sim must be a number, sim will be set equal to its default value 1000",call.=TRUE)
    sim = 1000
  }
  else{
    if(sim<10){
      warning("The minimum value that the parameter sim can take is 10, sim will be set equal to 10",call.=TRUE)
      sim = 10
    }
    if(sim%%1>0){
      warning("Parameter sim must be an integer, sim will be defined ignoring decimal places",call.=TRUE)
      sim = as.integer(sim)
    }
  }
  if(!is.numeric(conf)){
    warning("The conf parameter must be a probability, conf will be set equal to its default value 0.95",call.=TRUE)
    conf = 0.95
  }
  else{
    if(conf>1 | conf<0){
      warning("The parameter conf must be a probability, that is, be contained between 0 and 1, conf will be set equal to its default value 0.95",call.=TRUE)
      conf = 0.95
    }
    if(conf==0){
      warning("It is not possible to perform a confidence interval with parameter conf equal to 0, conf will be set equal to its default 0.95",call.=TRUE)
      conf = 0.95
    }
    if(conf==1){
      warning("It is not possible to perform a confidence interval with parameter conf equal to 1, conf will be set equal to its default 0.95",call.=TRUE)
      conf = 0.95
    }
  }
  n = length(res)
  samples = matrix(nrow = sim, ncol=n)
  q1 = c()
  q2 = c()
  q3 = c()
  confstar = (1 - conf)/2
  for(i in 1:sim){
    sample =  rhalfnorm(n)
    sample = sort(sample)
    samples[i,] = sample
  }
  for(i in 1:n){
    sample = samples[,i]
    q1 = c(q1,as.numeric(quantile(sample,0+confstar)))
    q2 = c(q2,as.numeric(quantile(sample,0.500)))
    q3 = c(q3,as.numeric(quantile(sample,1-confstar)))
  }
  res = sort(abs(res))
  index = 1:n
  p = (index - 0.5)/n
  q = qhalfnorm(p)
  ggplot() + geom_point(mapping = aes(y = res,x = q)) +
    geom_line(aes(y=q1,x=q)) + geom_line(aes(y=q2,x=q),linetype = "dashed")+
    geom_line(aes(y=q3,x=q)) +
    scale_y_continuous(breaks = seq(0,6,0.5)) +
    ylab("Residuals abs. ord.") +
    xlab("Half-normal quantiles") +
    ggtitle("Half-normal residual plots with simulated envelope",subtitle = paste0("With ",conf*100,"% confindence"))
}
