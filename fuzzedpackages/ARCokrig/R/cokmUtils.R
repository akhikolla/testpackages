#' @title Construct the cokm object
#' @description This function constructs the cokm object in
#'  autogressive cokriging models
#'  
#' @param formula a list of \eqn{s} elements, each of which contains the formula to specify fixed basis functions or regressors.
#' @param output a list of \eqn{s} elements, each of which contains a matrix of computer model outputs.
#' @param input a list of \eqn{s} elements, each of which contains a matrix of inputs.
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
#' @param nugget.est a logical value indicating whether the nugget is included or not. Default value is \code{FALSE}.
#' @param prior a list of arguments to setup the prior distributions with the reference prior as default.
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
#' @author Pulong Ma <mpulong@gmail.com>
#'  
#' @seealso \code{\link{ARCokrig}}, \code{\link{cokm.fit}}, \code{\link{cokm.predict}}
#' @export


cokm <- function(formula=list(~1,~1), output, input, cov.model="matern_5_2", nugget.est=FALSE,
	prior=list(), opt=list(), NestDesign=TRUE, tuning=list(), info=list()){

## check the arguments 
  .check.arg.cokm(formula=formula, output=output, input=input,
                  prior=prior, opt=opt, NestDesign=NestDesign, tuning=tuning,
                  info=info)

  S = length(output)  # number of code levels

  Dim = dim(input[[1]])[2]
  # if(length(param[[1]])==Dim){
  #   is.nugget=FALSE
  # }else{
  #   is.nugget=TRUE
  # }
  is.nugget = nugget.est
  param = list()
  for(i in 1:S){
    param.max = apply(input[[i]], 2, max)
    param.min = apply(input[[i]], 2, min)
    param[[i]] = (param.max-param.min)/2
  }
  
  phi = do.call(cbind, param)
  
  if(length(opt)==0){
    opt$maxit = 1000
    if(dim(phi)[1]==1){
      opt$method = "Brent"
      opt$lower = -10
      opt$upper = 20
    }else{
      opt$method = "Nelder-Mead"
      opt$lower = -Inf
      opt$upper = Inf
    }
  }else{
    if(dim(phi)[1]==1){
      if(!exists("method", where=opt)){
        opt$method = "Brent"
      }
      if(!exists("lower", where=opt)){
        opt$lower = -10
      }
      if(!exists("upper", where=opt)){
        opt$upper = 20
      }
    }else{
      opt$method = "Nelder-Mead"
      opt$lower = -Inf
      opt$upper = Inf
    }
  }


  if(NestDesign){
    
    #cat("\n Constructing cokm object for nested design.\n\n")
    tuning = list()

  }else{
    #cat("\n Constructing cokm object for non-nested design.\n\n")

    if(!exists("maxit", where=tuning)){
      tuning$maxit = 30
    }

    if(!exists("tol", where=tuning)){
      tuning$tol = 1e-3
    }

    if(!exists("n.sample", where=tuning)){
      tuning$n.sample = 30
    }
    if(!exists("verbose", where=tuning)){
      tuning$verbose = TRUE
    }
  }
  
  if(!exists("name", where=prior)){
    prior$name = "Reference"
  }
  
  if(!exists("hyperparam", where=prior)){
    prior$hyperparam = list()
    for(i in 1:S){
      prior$hyperparam[[i]] = list(a=0.2, b=1, nugget.UB=1)
    }
  }else{
    if(length(prior$hyperparam)==1){
      for(i in 2:S){
        hyperparam[[i]] = hyperparam[[1]]
      }
      prior$hyperparam = hyperparam
    }
  }



## construct the cokm object
 new("cokm",
   formula = formula,
   output = output,
   input = input,
   param = param,
   cov.model = cov.model,
   nugget.est = is.nugget,
   prior = prior,
   opt = opt,
   NestDesign = NestDesign,
   tuning = tuning,
   info = info
  )
  
}


#####################################################################
#####################################################################





#####################################################################
#####################################################################
.check.arg.cokm <- function(formula, output, input, prior, opt, 
  NestDesign, tuning, info){
  
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

  if(!is(prior, "list")){
    stop("\n\nprior should be a list containing arguments to
         setup the prior distributions.\n\n")
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
  



setMethod("summary", signature(object="cokm"),
          function(object, ...){
            message("cokm object\n")
            message("\n")
            message(paste0("Code levels:", length(object@data)))
            message("\n")
            message(paste0("Is nugget included:", object@nugget.est))
            message("\n")
            message(paste0("Nested Design:", object@NestDesign))
            message("\n\n")
          })

