#' @title Conditional simulation at new inputs in autoregressive cokriging models for multivarite output
#' @description This function makes prediction based on conditional simulation in
#'  autogressive cokriging models for multivariate output
#'  
#' @param obj a \code{\link{mvcokm}} object construted via the function \code{\link{mvcokm}} in 
#' this package
#' @param input.new a matrix including new inputs for making prediction
#' @param nsample a numerical value indicating the number of samples
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @seealso \code{\link{mvcokm}}, \code{\link{mvcokm.fit}}, \code{\link{mvcokm.predict}}, \code{\link{ARCokrig}}
#' 

mvcokm.condsim <- function(obj, input.new, nsample=30){

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
  pred.list = condsim.ND(formula=formula, output=output, input=input, 
              input.new=input.new, phi=phi, cov.model=cov.model,
              nsample=nsample)
}else{
  pred.list = condsim.NN(formula=formula, output=output, input=input,
              input.new=input.new, phi=phi, cov.model=cov.model,
              nsample=nsample)
}



return(pred.list)



}
