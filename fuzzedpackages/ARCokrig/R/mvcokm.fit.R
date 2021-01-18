#' @title fit the autoregressive cokriging model for multivariate output
#' @description This function estimates parameters in
#'  the parallel partial cokriging model
#'
#' @param obj a \code{\link{mvcokm}} object construted via the function \code{\link{mvcokm}} in 
#' this package
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' @seealso \code{\link{mvcokm}}, \code{\link{mvcokm.predict}}, \code{\link{mvcokm.condsim}}, \code{\link{ARCokrig}}


mvcokm.fit = function(obj){


formula = obj@formula
output = obj@output 
input = obj@input 
param = obj@param
cov.model = obj@cov.model 
nugget.est = obj@nugget.est
prior = obj@prior
NestDesign = obj@NestDesign
opt = obj@opt 
tuning = obj@tuning

hyperparam = prior$hyperparam


phi = do.call(cbind, param)


Dim = dim(input[[1]])[2]
p.x = Dim
# if(dim(phi)[1]==Dim){
# 	is.nugget=FALSE
# }else{
# 	is.nugget=TRUE
# }
is.nugget = nugget.est



###################################################################
#### begin parameter esimation
###################################################################

if(NestDesign){

fit = fit.ND(formula=formula, output=output, input=input, phi=phi,
		cov.model=cov.model, prior=prior, opt=opt)

phi.new = fit$par 
phi.new = split(phi.new, col(phi.new, as.factor = TRUE))


}else{

fit = fit.NN(formula=formula, output=output, input=input, phi=phi,
		cov.model=cov.model, prior=prior, opt=opt, MCEM=tuning)

phi.new = fit$par 
phi.new = split(phi.new, col(phi.new, as.factor = TRUE))

obj@info = list(iter=fit$iter, eps=fit$eps)

}

# colnames(phi.new) = paste0("Level", seq(1:S), "")

obj@param = phi.new 


return(obj)




}