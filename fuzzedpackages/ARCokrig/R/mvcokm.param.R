#' @title Get model parameters in autoregressive cokriging models for multivarite output
#' @description This function computes estimates for regression and variance parameters
#'  given the correlation parameters are known. It is used to show all model
#' parameters in one place.
#'  
#' @param obj a \code{\link{mvcokm}} object construted via the function \code{\link{mvcokm}} in 
#' this package
#' @return a list of model parameters including regression coefficients \eqn{\beta}, 
#' scale discrepancy \eqn{\gamma}, variance parameters 
#' \eqn{\sigma^2}, and correlation parameters \eqn{\phi} in covariance functions.
#' If nugget parameters are included in the model, then nugget parameters are shown in \eqn{\phi}.
#'
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @seealso \code{\link{mvcokm}}, \code{\link{mvcokm.fit}}, \code{\link{mvcokm.predict}}, \code{\link{ARCokrig}}
#' 

mvcokm.param <- function(obj){

formula = obj@formula
output = obj@output 
input = obj@input 
param = obj@param
cov.model = obj@cov.model 
NestDesign = obj@NestDesign

phi = do.call(cbind, param)

Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
	is.nugget=FALSE
}else{
	is.nugget=TRUE
}


###################################################################
#### augment input
###################################################################
S = length(output)   # number of code
out = augment.input(input)
input.union = out$union
input.miss = out$miss
input.list = list(input=input, input.miss=input.miss)

Cl = list()
for(t in 1:S){
	input.max = apply(input.list$input[[t]], 2, max)
	input.min = apply(input.list$input[[t]], 2, min)
	Cl[[t]] =  abs(input.max-input.min)
}

y = output
###################################################################
#### create covariates
###################################################################
H = list()
Hm = list()
for(t in 1:S){
	colnames(input[[t]]) = paste0("x", 1:p.x)
	df = data.frame(input[[t]])
	H[[t]] = model.matrix(formula[[t]], df)

	if(t<S){
		colnames(input.miss[[t]]) = paste0("x", 1:p.x)
		df = data.frame(input.miss[[t]])
		Hm[[t]] = model.matrix(formula[[t]], df)
	}

}


###################################################################
#### compute intermediate quantities
###################################################################
dist.o = list()
dist.m = list()
dist.mo = list()
distlist = list()

for(t in 1:S){
	dist.o[[t]] = compute_distance(input[[t]], input[[t]])

	if(t<S){
		dist.m[[t]] = compute_distance(input.miss[[t]], input.miss[[t]])
		dist.mo[[t]] = compute_distance(input.miss[[t]], input[[t]])
	}

	distlist[[t]] = compute_distance(input.union[[t]], input.union[[t]])
}


n.aug = rep(NA, S)
q = rep(NA, S)
for(t in 1:S){
	n.aug[t] = dim(y[[t]])[1]
	q[t] = dim(H[[t]])[2]
}



if(NestDesign){
	stop("Not implemented yet.")
}else{
	###################################################################
	#### estimate sigma and beta with marginal posteriors 
	###################################################################
	sigma2.hat = list()
	beta.hat = list()
	ym.hat = list()
	### estimate ym based on [ym|y,phi] for t=1

	t=1

	R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
	U = chol(R)
	RInv = chol2inv(U)
	HRHInv = solve(t(H[[t]])%*%RInv%*%H[[t]])

	Rm = buildcov(phi[ ,t], dist.m[[t]], covmodel=cov.model, nugget=is.nugget)
	Rmo = buildcov(phi[ ,t], dist.mo[[t]], covmodel=cov.model, nugget=FALSE)

	# compute conditional mean
	RmoRInv = Rmo%*%RInv
	KW = (Hm[[t]]-RmoRInv%*%H[[t]]) %*% HRHInv %*% t(H[[t]]) %*% RInv + RmoRInv
	ym.hat[[t]] = KW %*% y[[t]]  # n.m-by-N matrix

	### estimate sigma2 based on [sigma2|ym,y,phi] for t=1
	Q = RInv - RInv%*%H[[t]]%*%HRHInv%*%t(H[[t]])%*%RInv
	sigma2.hat[[t]] = compute_Svec(output=y[[t]], Q=Q) / (n.aug[t]-q[t]+2)

	### estimate beta based on [beta|ym,y,phi] for t=1
	beta.hat[[t]] = HRHInv%*%t(H[[t]])%*%RInv%*%y[[t]]

	### estimate sigma and beta for t>1

	for(t in 2:S){
		R = buildcov(phi[ ,t], dist.o[[t]], covmodel=cov.model, nugget=is.nugget)
		U = chol(R)
		RInv = chol2inv(U)

		y_t1 = create.w(t=t, input=input, input.miss=input.miss[[t-1]], 
			y=y[[t-1]], ym=ym.hat[[t-1]])

		out = compute_param(y_t=y[[t]], Ht=H[[t]], y_t1=y_t1, RInv)
		sigma2.hat[[t]] = out$sigma2
		beta.hat[[t]] = out$beta 
	}


	names(beta.hat) = paste0("Level", seq(1:S), "")
	names(sigma2.hat) = paste0("Level", seq(1:S), "")

}



return(list(corr=param, coeff=beta.hat, var=sigma2.hat))



}
