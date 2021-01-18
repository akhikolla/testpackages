#' Matrix-normal maximum likelihood estimator
#' @param X r x p x N array of training data
#' @param class vector of class labels for N training data observations
#' @param max.iter maximum number of iterations for flip-flop algorithm
#' @param tol covergence tolerance for flip flop algorithm
#' @param quiet should print obj fnc values at each iteration?
#' @return a list of class "MN" containing
#' @return Mean the r x p x C estimate of class means
#' @return U the estimated U
#' @return V the estimated V
#' @return pi.list marginal class probabilites based on training data
#' @export
MN_MLE = function(X, class, max.iter = 1000, tol=1e-6, quiet=TRUE){

	r = dim(X)[1]
	p = dim(X)[2]
	N = dim(X)[3]
	C = length(unique(class))
	K = C*(C-1)/2
	nc = count(class)[,2]


	mean.array = array(0, dim=c(r, p, C))
	
	for(jj in 1:C){
		mean.array[,,jj] = apply(X[,,class==jj], c(1,2), mean)
	}


	## mean center observations to avoid redundant calculations
	X.c = CENTER_X(X, class, mean.array)

	orig.objfnc = EVAL_OBJ_FUNC(X, class, mean.array, diag(1,r), diag(1,p), matrix(0, nrow=K*r, ncol=p), 0,0, S.u= NULL, S.v = NULL)$out

	UV.iter = 1
	UVconverged = FALSE

	# initialize V
	U = diag(1, r)
	
	old.objfnc = orig.objfnc

	while(!UVconverged){

		Vinv.S = V_SAMPLE(X.c, U)
		tempV = solve(Vinv.S)

		Uinv.S = U_SAMPLE(X.c, tempV)
		tempU = solve(Uinv.S)


		new.objfnc = EVAL_OBJ_FUNC(X, class, mean.array, tempU, tempV, matrix(0, nrow=K*r, ncol=p), 0,0, S.u= NULL, S.v = NULL)$out

		r1 = (old.objfnc - new.objfnc)/orig.objfnc

		if(r1 < tol){
			UVconverged=TRUE
		} else {
			if(UV.iter >= max.iter){
				UVconverged=TRUE
			}
		}

		UV.iter = UV.iter + 1
		V = tempV
		U = tempU
		old.objfnc = new.objfnc
		}

		out = list("U" = U, "V" = V, "Mean" = mean.array, "pi.list" = nc/N)
		class(out) = "MN"

	return(out)
}

