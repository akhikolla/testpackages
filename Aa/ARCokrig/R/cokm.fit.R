#' @title fit the autoregressive cokriging model
#' @description This function estimates parameters in
#'  autogressive cokriging models
#'
#' @param obj a \code{\link{cokm}} object construted via the function \code{\link{cokm}} in 
#' this package
#' 
#' @author Pulong Ma <mpulong@gmail.com>
#' 
#' @export
#' 
#' @seealso \code{\link{cokm}}, \code{\link{cokm.param}}, \code{\link{cokm.predict}}, \code{\link{ARCokrig}}
#' 
#' @examples 
#' 
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
#' ## create the cokm object
#' prior = list(name="JR")
#' obj = cokm(formula=list(~1,~1+x1), output=list(c(zc), c(zf)),
#'               input=list(as.matrix(Dc), as.matrix(Df)),
#'               prior=prior, cov.model="matern_5_2")
#'
#' ## update model parameters in the cokm object
#' 
#' obj = cokm.fit(obj)
#' 
#' 



cokm.fit <- function(obj){

formula = obj@formula
output = obj@output 
input = obj@input 
param = obj@param
cov.model = obj@cov.model 
prior = obj@prior
hyperparam = prior$hyperparam
opt = obj@opt 
NestDesign = obj@NestDesign


phi = do.call(cbind, param)

Dim = dim(input[[1]])[2]
p.x = Dim
if(dim(phi)[1]==Dim){
	is.nugget=FALSE
}else{
	is.nugget=TRUE
}


S = length(output)

###################################################################
#### create covariates
###################################################################
H = list()
for(t in 1:S){
	colnames(input[[t]]) = paste0("x", 1:p.x)
	df = data.frame(input[[t]])
	H[[t]] = model.matrix(formula[[t]], df)

}



Cl = list()
for(t in 1:S){
  input.max = apply(input[[t]], 2, max)
  input.min = apply(input[[t]], 2, min)
  Cl[[t]] =  abs(input.max-input.min)
}

for(t in 1:S){
  hyperparam[[t]]$Cl = Cl[[t]]
}






###################################################################
#### begin optimization algorithm for Nested Design
###################################################################
if(NestDesign){

	###################################################################
	#### compute intermediate quantities
	###################################################################

	distlist = list()
	for(t in 1:S){
		distlist[[t]] = compute_distance(input[[t]], input[[t]])
	}

	phi.new = phi 

	for(t in 1:S){

		if(is.nugget){
			nu = log(phi[p.x+1, t]) - log(hyperparam[[t]]$nugget.UB-phi[p.x+1, t])   # logit of nugget
			init.val = c(-log(phi[1:p.x, t]), nu)	
		}else{
			init.val = -log(phi[ ,t])
		}

		fit = try(optim(init.val, margin.posterior, input=input, output=output, level=t,
				        H=H, dist=distlist[[t]], cov.model=cov.model,
				        prior=prior$name, hyperparam=hyperparam[[t]], 
				        control=list(fnscale=-1, maxit=opt$maxit),
				        method=opt$method, lower=opt$lower, upper=opt$upper),
				silent=T)

		if(inherits(fit, "try-error")){
			phi.new[ ,t] = phi[ ,t]
			print(paste0("optimization error, skip t=", as.character(t)))
			print(fit)
		}else{
			if(is.nugget){
				nugget = hyperparam[[t]]$nugget.UB*exp(fit$par[p.x+1]) / (1+exp(fit$par[p.x+1]))
				phi.new[ ,t] = c(exp(-fit$par[1:p.x]), nugget)
			}else{
				phi.new[ ,t] = exp(-fit$par)
			}	
		}

	}

}else{

	for(t in 1:S){
	  if(!is.matrix(output[[t]])){
	    output[[t]] = as.matrix(output[[t]])
	  }
	}
  

	###################################################################
	#### augment input
	###################################################################
	out = augment.input(input)
	input.union = out$union
	input.miss = out$miss
	input.list = list(input=input, input.miss=input.miss)



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

	tuning = obj@tuning
	n.sample = tuning$n.sample
	maxit = tuning$maxit
	tol = tuning$tol
	verbose = tuning$verbose 


###################################################################
#### begin MCEM algorithm for Non-Nested Design
###################################################################
	phi.new = phi
	conv = FALSE
	iter = 1
	while(!conv){

		###############################################################
		#### generate M Monte Carlo samples for missing data
		###############################################################
		# y.m = list()
		#system.time(
		# for(k in 1:n.sample){
		# 	y.m[[k]] = sample.ym(y=output,input=input,param=phi,Ho=H,Hm=Hm,dist.o=dist.o,
		# 		      dist.m=dist.m,dist.mo=dist.mo,cov.model=cov.model)
		# }
		y.m = sample.ym(y=output,input=input,param=phi,Ho=H,Hm=Hm,dist.o=dist.o,
				      dist.m=dist.m,dist.mo=dist.mo,cov.model=cov.model, 
				      nsample=n.sample)
		#)
		###############################################################
		#### compute and maximize the Q function at each fidelity
		###############################################################
		
		for(t in 1:S){

			if(is.nugget){
				nu = log(phi[p.x+1, t]) - log(hyperparam[[t]]$nugget.UB-phi[p.x+1, t])   # logit of nugget
				init.val = c(-log(phi[1:p.x, t]), nu)	
			}else{
				init.val = -log(phi[ ,t])
			}

			fit = try(optim(init.val, compute.g.univ, input.list=input.list, level=t, y=output, H=H, ym=y.m, Hm=Hm,
					        dist=distlist[[t]], hyper=hyperparam[[t]], cov.model=cov.model,
				            control=list(fnscale=-1, maxit=opt$maxit),
				            method=opt$method, lower=opt$lower, upper=opt$upper),
					silent=T)

			# fit = try(optimr(init.val, compute.Q.default, input=input, level=t, y=output, H=H, y.m=y.m, Hm=Hm,
			# 	        distlist=distlist, cov.model=cov.model,
			#             control=list(fnscale=-1, maxit=opt$maxit),
			#             method=opt$method, lower=opt$lower, upper=opt$upper),
			# 	silent=T)

			if(inherits(fit, "try-error")){
				phi.new[ ,t] = phi[ ,t]
				print(paste0("optimization error, skip t=", as.character(t)))
				print(fit)
			}else{
				if(is.nugget){
					nugget = hyperparam[[t]]$nugget.UB*exp(fit$par[p.x+1]) / (1+exp(fit$par[p.x+1]))
					phi.new[ ,t] = c(exp(-fit$par[1:p.x]), nugget)
				}else{
					phi.new[ ,t] = exp(-fit$par)
				}	
			}

		}

		###############################################################
		#### check convergence
	    if(inherits(fit, "try-error")){
	    	diff = tol + 1
	    }else{
	    	diff = mean((phi.new - phi)^2)
	    }

		if(verbose){
			print(paste0("iter=", as.character(iter)))
		}
		

		iter = iter + 1

		if(iter>maxit){
			conv = TRUE
		}else{
			if(diff<tol){
				conv = TRUE
			}
		}

		phi = phi.new

	}

obj@info=list(iter=iter, eps=diff)

}

colnames(phi.new) = paste0("Level", seq(1:S), "")

phi.new = split(phi.new, col(phi.new, as.factor = TRUE))


obj@param = phi.new 

return(obj)

}
