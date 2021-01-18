#' @title Get model parameters in the autoregressive cokriging model
#' @description This function compute estimates for regression and variance
#'  parameters given the correlation parameters are known. It is used to show
#' all model parameters in one place.
#' @param obj a \code{\link{cokm}} object construted via the function \code{\link{cokm}} in 
#' this package
#' 
#' @return a list of model parameters including regression coefficients \eqn{\beta}, 
#' scale discrepancy \eqn{\gamma}, variance parameters 
#' \eqn{\sigma^2}, and correlation parameters \eqn{\phi} in covariance functions.
#' If nugget parameters are included in the model, then nugget parameters are shown in \eqn{\phi}.
#'
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @seealso \code{\link{cokm}}, \code{\link{cokm.fit}}, \code{\link{cokm.condsim}}, \code{\link{ARCokrig}}


cokm.param <- function(obj){



formula = obj@formula
output = obj@output 
input = obj@input 
param = obj@param
cov.model = obj@cov.model 

phi = do.call(cbind, param)

Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
	is.nugget=FALSE
}else{
	is.nugget=TRUE
}

S = length(output)
y = output

###################################################################
#### create covariates
###################################################################
H = list()
for(t in 1:S){
	colnames(input[[t]]) = paste0("x", 1:p.x)
	df = data.frame(input[[t]])
	H[[t]] = model.matrix(formula[[t]], df)

}


###################################################################
#### compute intermediate quantities
###################################################################

distlist = list()
for(t in 1:S){
	distlist[[t]] = compute_distance(input[[t]], input[[t]])
}


###################################################################
#### get b and sigma
###################################################################
X = list()
b = list()
sigma = rep(NA, S)

for(t in 1:S){

	R = buildcov(phi[ ,t], distlist[[t]], cov.model, is.nugget)
	n = dim(R)[1]
	U = chol(R)
	RInv = chol2inv(R)
	if(t==1){
		X[[t]] = H[[t]]
	}else{
		IB = match.input(input[[t]], input[[t-1]])$IB
		y_t1 = y[[t-1]][IB]
		X[[t]] = cbind(H[[t]], y_t1)
	}
	q = dim(X[[t]])[2]
	RInvX = RInv%*%X[[t]]
	XRXInv = solve(t(X[[t]])%*%RInvX)

	b[[t]] = c(XRXInv%*%(t(X[[t]])%*%(RInv%*%y[[t]])))

	Q = RInv - RInvX%*%XRXInv%*%t(RInvX)
	sigma[t] = t(y[[t]])%*%Q%*%y[[t]] / (n-q+2)

}

names(b) = paste0("Level", seq(1:S), "")
names(sigma) = paste0("Level", seq(1:S), "")

out = list(corr=param, coef=b, var=sigma)
return(out)




}
