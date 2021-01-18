


#' MATLDA: function fit the penalized matrix normal
#' model for a single pair of tuning parameters (lam1,lam2)
#' @param X an r x p x N array of training data
#' @param          class N-vector of trainig class labels {1,...,C}
#' @param	   lambda1 nonneg tp for mean penalty
#' @param	   lambda2 nonneg tp for Kronecker penalty
#' @param	   quiet should print iterations and obj fnc vals?
#' @param	   Xval r x p x Ntest array of test data
#' @param		   classval Ntest-vector of test class data labels
#' @param		   k.iter max iterations for full algorithm
#' @param	   cov.tol tolerance for glasso algorithms
#' @param	   m.tol tolerance for AMA algorithms
#' @param		   full.tol tolerance for full algorithm
#' @return 	Val.mat matrix of misclassification rates on test set for all pairs lam1, lam2
#' @return	G.mat matrix of number of pairwise mean differences that are zero for all pairs lam1, lam2
#' @return	U.mat number of zeroes in U for all lam1, lam2
#' @return	V.mat number of zeroes in V for all lam1, lam2
#' @export

MatLDA_Grid = function(X, class, lambda1, lambda2, quiet=TRUE,
	Xval=NULL, classval= NULL,
	k.iter = 100, cov.tol=1e-5, m.tol=1e-5, full.tol=1e-6){

# ----- Preliminaries; order n1 > .. > nC
	nc.orig = count(class)[,2]
	N = sum(nc.orig)
	class.orig = class
	class.val.orig =  classval

	out = rank(-nc.orig, ties.method = "first")
	class = rep(0, N)
	for (gg in 1:N) {
		class[gg] = out[class.orig[gg]]
	}

	nc = count(class)[,2]
	C = length(nc)
	K = (C*(C-1))/2
	N = sum(nc)
	p = dim(X)[2]
	r = dim(X)[1]
	pihat = nc/N

	if (r < p) {
		warning("Row dimension is greater than column dimension; consider transposing data.")
	}

	classval = rep(0, dim(Xval)[3])
	for (gg in 1:dim(Xval)[3]) {
		classval[gg] = out[class.val.orig[gg]]
	}


# ---- initialize U and V at diagonal from MLE
	out = MN_MLE(X, class)
	V.unscaled = diag(diag(out$V))
	U.unscaled =  diag(diag(out$U))

# ---- normalize so ||U||_1 = r
	out1 = sum(abs(U.unscaled))
	UInit = (r/out1)*U.unscaled
	VInit = (out1/r)*V.unscaled

	Vinv = solve(VInit)
	Uinv = solve(UInit)

# --- compute rho
	rho.max = (4*min(diag(UInit))*min(diag(VInit))*min(nc))/(C*N)
	rho = rho.max/10

# ---- make sample mean array
	mean.array = array(0, dim=c(r,p,C))
	for (k in 1:C) {
		mean.array[,,k] = apply(X[,,class==k], c(1,2), mean)
	}

# ----- compute weight matrix
	weightmat = matrix(0, nrow=K*r, ncol=p)
	for (c in 1:(C-1)) {
		for(cc in (c+1):C){
		k = (C - 1)*(c-1) + (cc-1) - (c*(c-1)/2)
		weightmat[((k-1)*r + 1):(k*r), ] = 1/abs(mean.array[,,c] - mean.array[,,cc])
		}
	}

# ----- set upper bound weight matrix
	weightmat[which(weightmat>1e6)] = 1e6

# ---- allocate memory
	G.mat = matrix(0, nrow=length(lambda1),  ncol=length(lambda2))
	U.mat = matrix(0, nrow=length(lambda1),  ncol=length(lambda2))
	V.mat = matrix(0, nrow=length(lambda1),  ncol=length(lambda2))

	if(!is.null(Xval)){
		Val = matrix(0, nrow=length(lambda1), ncol=length(lambda2))
	}

	outer.converge=FALSE
	meanstacked = STACKED(mean.array)
	M.update = mean.array

# -------------------------------------
# Iterating 
# -------------------------------------
	for(j in 1:length(lambda1)){

	# --- update means
		D = matrix(0, nrow=K*r, ncol=p)
		out = M_Update_Rcpp(Mean=meanstacked,
			 D=D, Uinv=Uinv, Vinv=Vinv, nc = nc, rho=rho, weightmat = weightmat,
		  lambda=lambda1[j], C=C, r=r, mtol = m.tol)

	# --- store D for warm starts
		D = out$D

	# ---- threshold mean update to impose exact equality
		M = THRESH_MEAN(out$M, out$G, C)
		M.update = array(0, dim=c(r, p, C))
		for(c in 1:C){
			M.update[,,c] = M[((c-1)*r +1):(c*r),]
		}

		# ---- center and compute path of V for all lambda2
		Xc = CENTER_X(X, class, M.update)
       	S.v = V_SAMPLE(X.c=Xc, U = UInit)
  		V.temp = glassopath(S.v, rholist=lambda2, trace=0, thr=cov.tol, penalize.diagonal=TRUE)$wi

  		for(i in 1:length(lambda2)){
  			init.objfncval = EVAL_OBJ_FUNC(X, class, mean.array, UInit, VInit,  weightmat, lambda1[j], lambda2[i])$out
			if(!quiet){
				cat("obj.func after V update=", EVAL_OBJ_FUNC(X, class, M.update, UInit, V.temp[,,i], weightmat, lambda1[j], lambda2[i], S.v = S.v)$out, "\n")
			}

			lambda.u = (lambda2[i]*sum(abs(V.temp[,,i]))/p)
			S.u = U_SAMPLE(X.c = Xc, V = V.temp[,,i])
			U.temp = glasso(S.u, rho = lambda.u, thr=cov.tol, trace=0, penalize.diagonal=TRUE)$wi

			if(!quiet){
			cat("obj.func after U update=", EVAL_OBJ_FUNC(X, class, M.update, U.temp, V.temp[,,i], weightmat, lambda1[j], lambda2[i], S.u = S.u)$out, "\n")
			}

			out.inner = GENERAL_ITER(X, class, lambda1 = lambda1[j], lambda2 = lambda2[i],
			quiet=quiet, UInit = UInit, VInit = VInit, U=U.temp, V=V.temp[,,i],
					 meanstacked = meanstacked, mean.array= mean.array, M.update = M.update, D=D, weightmat=weightmat,
					  cov.tol=cov.tol, m.tol=m.tol, full.tol = full.tol, init.objfncval=init.objfncval,k.iter=k.iter)

  			# --- if Valset is non-NULL, compute misclassification rate
  			if(!is.null(Xval)){
				Val[j,i] = MN_CLASSIFICATION(X=Xval, class=classval, out.inner$Mean, out.inner$U, out.inner$V, pi.list=pihat, C=C)$misclass
			} else {
				Val = NULL
			}

			G.mat[j,i] = 0
			for (jj in 1:(C-1)) {
				for (kk in (jj+1):C) {
					G.mat[j,i] = G.mat[j,i] + sum(abs(out.inner$Mean[,,jj] - out.inner$Mean[,,kk])==0)
				}
			}

			U.mat[j,i] = sum(out.inner$U==0)
			V.mat[j,i] = sum(out.inner$V==0)

		}
	}

	return(list("Val.mat" = Val, "G.mat" = G.mat, "U.mat" = U.mat, "V.mat" = V.mat))

}
