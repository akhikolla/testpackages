#' @title Prediction at new inputs in autoregressive cokriging models for multivarite output
#' @description This function makes prediction in
#'  the parallel partial cokriging model. If a nested design is used, the predictive mean and predictive variance are
#'  computed exactly; otherwise, Monte Carlo simulation from the predictive distribution is used to approximate
#' the predictive mean and predictive variance. 
#' @param obj a \code{\link{mvcokm}} object construted via the function \code{\link{mvcokm}} in 
#' this package
#' @param input.new a matrix including new inputs for making prediction
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @seealso \code{\link{mvcokm}}, \code{\link{mvcokm.fit}}, \code{\link{mvcokm.condsim}}, \code{\link{ARCokrig}}
#' 

mvcokm.predict <- function(obj, input.new){

formula = obj@formula
output = obj@output 
input = obj@input 
param = obj@param
cov.model = obj@cov.model 
NestDesign = obj@NestDesign

phi = do.call(cbind, param)

if(!is.matrix(input.new)){
  stop("input.new should be a matrix.")
}

if(NestDesign){
  pred.list = predict.ND(formula=formula, output=output, input=input, 
              input.new=input.new, phi=phi, cov.model=cov.model)
}else{
  n.sample = obj@tuning$n.sample
  pred.list = predict.NN(formula=formula, output=output, input=input,
              input.new=input.new, phi=phi, cov.model=cov.model,
              nsample=n.sample)
}



return(pred.list)



}
