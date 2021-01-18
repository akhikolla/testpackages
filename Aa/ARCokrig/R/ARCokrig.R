#' @title Fit the AR-Cokriging model and make predictions
#' @description This is a simple and high-level funciton to fit autoregressive
#' cokriging models to multifidelity computer model outputs. 
#'  
#' @param formula a list of \eqn{s} elements, each of which contains the formula to specify fixed basis functions or regressors.
#' @param output a list of \eqn{s} elements, each of which contains a matrix of computer model outputs.
#' @param input a list of \eqn{s} elements, each of which contains a matrix of inputs.
#' 
#' @param cov.model a string indicating the type of covariance
#' function in AR-cokriging models. Current covariance functions include
#' \describe{
#' \item{exp}{product form of exponential covariance functions.} 
#' \item{matern_3_2}{product form of Matern covariance functions with 
#' smoothness parameter 3/2.}
#' \item{matern_5_2}{product form of Matern covariance functions with
#' smoothness parameter 5/2.}
#' \item{Gaussian}{product form of Gaussian covariance functions.}
#' \item{powexp}{product form of power-exponential covariance functions with roughness parameter fixed at 1.9.}
#' }
#' 
#' @param nugget.est a logical value indicating whether nugget parameter is included or not. Default value is \code{FALSE}.
#' @param input.new a matrix including new inputs for making prediction
#' @param prior a list of arguments to setup the prior distributions
#' \describe{
#'  \item{name}{the name of the prior. Current implementation includes 
#'  \code{JR}, \code{Reference}, \code{Jeffreys}, \code{Ind_Jeffreys}}
#'  \item{hyperparam}{hyperparameters in the priors. 
#'  For jointly robust (JR) prior, three parameters are included: 
#' \eqn{a} refers to the polynomial penalty to avoid singular correlation 
#'   matrix with a default value 0.2; \eqn{b} refers to the exponenetial penalty to avoid 
#'   diagonal correlation matrix with a default value 1; nugget.UB is the upper
#' bound of the nugget variance with default value 1, which indicates that the
#' nugget variance has support (0, 1).}
#'
#'}
#'
#' @param opt a list of arguments to setup the \code{\link{optim}} routine.
#' @param NestDesign a logical value indicating whether the 
#' experimental design is hierarchically nested within each level
#' of the code.
#' 
#' @param tuning a list of arguments to control the MCEM algorithm for non-nested
#' design. It includes the arguments 
#' \describe{
#'    \item{maxit}{the maximum number of MCEM iterations.}
#'    \item{tol}{a tolerance to stop the MCEM algorithm. If the parameter 
#'    difference between any two consecutive MCEM algorithm is less than 
#'    this tolerance, the MCEM algorithm is stopped.}
#'    \item{n.sample}{the number of Monte Carlo samples in the 
#'    MCEM algorithm.}
#'    \item{verbose}{a logical value to show the MCEM iterations if it is true.}
#'}
#'
#'
#' @param info a list that contains 
#' \describe{
#'  \item{iter}{number of iterations used in the MCEM algorithm}
#'  \item{eps}{parameter difference after the MCEM algorithm stops}
#'}
#' @return The main call inside \code{\link{ARCokrig}} consists of 
#' \code{\link{cokm}}, \code{\link{cokm.fit}}, and \code{\link{cokm.predict}}.
#' Thus, the function returns the \code{\link{cokm}} object and predictions 
#' over new inputs.
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' @export
#' 
#' 
#' @seealso \code{\link{cokm}}, \code{\link{cokm.param}}, \code{\link{cokm.fit}}, \code{\link{cokm.predict}}
#' 
#'
#' @references {
#' \itemize{
#' \item{Ma, P. (2019). ``Objective Bayesian Analysis of a Cokriging Model for Hierarchical Multifidelity Codes." arXiv:1910.10225. \url{https://arxiv.org/abs/1910.10225}. }
#' \item{Ma, P., Karagiannis, G., Konomi, B., Asher, T., Toro, G., and Cox, A. (2019) ``Multifidelity Computer Model Emulation with High-Dimensional Output: An Application to Storm Surge."
#'      arXiv:1909.01836. \url{https://arxiv.org/abs/1909.01836}.}
#' }
#' }
#'
#' @examples 
#' 
#' ##############################################################
#' ##############################################################
#' ############ Example
#' Funcc = function(x){
#'   return(0.5*(6*x-2)^2*sin(12*x-4)+10*(x-0.5)-5)
#' }
#' 
#' Funcf = function(x){
#'   z1 = Funcc(x)
#'   z2 = 2*z1-20*x+20 + sin(10*cos(5*x))
#'   return(z2)
#' }
#' 
#' #####################################################################
#' ###### Nested design 
#' #####################################################################
#' Dc <- seq(-1,1,0.1)
#' indDf <- c(1, 3, 6, 8, 10, 13, 17, 21)
#' zc <- Funcc(Dc)
#' Df <- Dc[indDf]
#' zf <- Funcf(Df)
#' 
#' input.new = as.matrix(seq(-1,1,length.out=200))
#' 
#' 
#' ## fit and predict with the AR-Cokriging model
#' 
#' out = ARCokrig(formula=list(~1,~1+x1), output=list(c(zc), c(zf)),
#'               input=list(as.matrix(Dc), as.matrix(Df)),
#'               cov.model="matern_5_2",
#'               input.new=input.new)
#' 
#' ## plot results
#' \donttest{
#' library(ggplot2)
#' cokrig = out$cokrig
#' df.l1 = data.frame(x=c(Dc), y=c(zc))
#' df.l2 = data.frame(x=c(Df), y=c(zf))
#' CI.lower = cokrig$lower95[[2]]
#' CI.upper = cokrig$upper95[[2]]
#' df.CI = data.frame(x=c(input.new),lower=CI.lower, upper=CI.upper)
#' df.pred = data.frame(x=c(input.new), y=cokrig$mu[[2]])
#'
#' g = ggplot(data.frame(x=c(-1,1)), aes(x)) + 
#'  stat_function(fun=Funcc, geom="line", aes(colour="level 1"), n=500) +
#'  stat_function(fun=Funcf, geom="line", aes(colour="level 2"), n=500) 
#'
#' g = g + geom_point(data=df.l1, mapping=aes(x=x, y=y), shape=16, size=2, color="black") + 
#'  geom_point(data=df.l2, mapping=aes(x=x, y=y), shape=17, size=2, color="black")
#'
#' g = g + geom_line(data=df.pred, aes(x=x, y=y, colour="cokriging"), inherit.aes=FALSE) +
#'  geom_ribbon(data=df.CI, mapping=aes(x=x,ymin=lower, ymax=upper), fill="gray40",
#'              alpha=0.3, inherit.aes=FALSE)
#' g = g + scale_colour_manual(name=NULL, values=c("red","blue", "turquoise3"), 
#'                            breaks=c("cokriging","level 1", "level 2"))  
#'
#' g = g + ggtitle("A Two-Level Example") + 
#'  theme(plot.title=element_text(size=14),
#'        axis.title.x=element_text(size=14),
#'        axis.text.x=element_text(size=14),
#'        axis.title.y=element_text(size=14),
#'        axis.text.y=element_text(size=14),
#'        legend.text = element_text(size=12),
#'        legend.direction = "horizontal",
#'        legend.position = c(0.6, 0.1)) + xlab("") + ylab("")
#' print(g)
#' 
#' }
#' 
#' 
#' 



ARCokrig <- function(formula=list(~1,~1), output, input, cov.model="matern_5_2", nugget.est=FALSE,
	input.new, prior=list(), opt=list(), NestDesign=TRUE, tuning=list(), info=list()){
  
  
  
  ## check the arguments 
  .check.arg.ARCokrig(formula=formula, output=output, input=input,
                      input.new=input.new, prior=prior, opt=opt, 
                      NestDesign=NestDesign, tuning=tuning)



  #cat("\n Constructing the AR-Cokriging model object.\n\n")
  obj = cokm(formula=formula, output=output, input=input, cov.model=cov.model,
    nugget.est=nugget.est, prior=prior, opt=opt, NestDesign=NestDesign, 
    tuning=tuning, info=info)

  ## fit the AR-Cokriging model
  #cat("\n Fit the AR-Cokriging model.\n\n")
  obj = cokm.fit(obj)
  
  ## predict with the AR-Cokriging model
	#cat("\n Predict with the AR-Cokriging model\n\n")
	pred = cokm.predict(obj=obj, input.new=input.new)

	return(list(obj=obj, cokrig=pred))

}


#####################################################################
#####################################################################
.check.arg.ARCokrig <- function(formula, output, input, input.new, 
  prior, opt, NestDesign, tuning){
  
  if(!is(formula, "list")){
    stop("\n\n formula should be a list contaning the regressors at each code level.\n\n")
  }
  if(!is(output, "list")){
    stop("\n\noutput should be a list of responses. Each element in a list should 
         contain output from a code level. The first level should contain
         output from the code with the lowest fidelity.\n\n")
  }
  
  s = length(output)
  
  
  if(!is(input, "list")){
    stop("\n\ninput should be a list of inputs in computer models.\n\n")
  }

  for(t in 1:s){
    if(!is(input[[t]], "matrix")){
      message("\n\n coerce input to a matrix format.\n\n")
      input[[t]] = as.matrix(input[[t]])
    }
  }
  
  
  # if(!is(param, "list")){
  #   stop("\n\nparam should be a list with each element containing initial values for  
  #        correlation parameters and nugget variance parameter (if needed).\n\n")
  # }

  if(!is(input.new, "matrix")){
    stop("\n\ninput.new should be a matrix.\n\n")
  }

  if(!is(prior, "list")){
    stop("\n\nhyperparam should be a list containing parameters 
         to setup the prior distributions.\n\n")
  }

  if(!is(opt, "list")){
    stop("\n\nopt should be a list with each element containing optimization arguments 
         at each code level.\n\n")
  }  

  if(!is(NestDesign, "logical")){
    stop("NestDesign should be a logical value indicating whether the design is 
         hierarchically nested.")
  }

  if(!is(tuning, "list")){
    stop("\n\n tuning should be a list containing tuning parameters to setup the 
          MCEM algorithm in non-nested design.\n")
  }
  
  
}

