

GENERAL_ITER = function(X, class, lambda1, lambda2,  quiet, UInit, VInit, U, V, 
	meanstacked, mean.array, M.update, D, weightmat, cov.tol, m.tol, full.tol, init.objfncval, k.iter){

# ----- preliminaries
	nc = count(class)[,2]
	C = length(nc)
	K = (C*(C-1))/2
	N = sum(nc)
	p = dim(X)[2]
	r = dim(X)[1]
	
# ------ normalize covariance matrices
	out1 = sum(abs(U))
	U.inner = (r/out1)*U
	V.inner = (out1/r)*V
	Uinv.inner = solve(U.inner)
	Vinv.inner = solve(V.inner)

# ----- set initial step size for AMA
	out1 = min(eigen(U.inner)$val)
	out2 = min(eigen(V.inner)$val)
	rho = ((4/C)*(min(nc)/N)*out1*out2)/10

# ------ initialize
	k.outer = 1
	outer.converge = FALSE

# ------ compute obj func values
	old.objfncval = EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.inner,  weightmat, lambda1, lambda2)$out
	
# ----- print first objective function value 
	if(!quiet){
		cat("k.outer=", k.outer,"objfncval = ", old.objfncval,"\n")
	}
	old.objfncval = init.objfncval

# -------------------------------------
# iterating 
# -------------------------------------

	while(!outer.converge){
		if(lambda1!=0){
		# ---- Mean update using Rcpp/AMA
			out = M_Update_Rcpp(Mean=meanstacked, 
				 D=D, Uinv=Uinv.inner, Vinv=Vinv.inner, nc = nc, rho = rho, weightmat = weightmat, 
				 lambda=lambda1, C=C, r=r, mtol=m.tol)
			D = out$D
			
		# ---- threshold means based on theta value 
			M = THRESH_MEAN(out$M, out$G, C)
			M.update = array(0, dim=c(r, p, C))
			for(c in 1:C){
				M.update[,,c] = M[((c-1)*r +1):(c*r),]
			}

		# --- Print progress 
			if(!quiet){
	 			cat("obj.func after M update=", temp <- EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.inner, weightmat, lambda1, lambda2)$out, "\n")
	 		}
	 		
	 	} else {
	 		M.update = mean.array
	 	}

 	# ----- center X and update V and U
		Xc = CENTER_X(X, class, M.update)
   		S.v = V_SAMPLE(Xc, U=U.inner)
		V.t = glasso(S.v, rho=lambda2, trace=0, thr=cov.tol, penalize.diagonal=TRUE)$wi
		if(!quiet){
			cat("obj.func after V update=", EVAL_OBJ_FUNC(X, class, M.update, U.inner, V.t, weightmat, lambda1, lambda2, S.v = S.v)$out, "\n")
		}

		S.u = U_SAMPLE(Xc, V=V.t)
		lambda.u = (lambda2/p)*sum(abs(V.t))
		U.t = glasso(S.u, rho=lambda.u, trace=0, thr=cov.tol, penalize.diagonal=TRUE)$wi
		new.objfncval = EVAL_OBJ_FUNC(X, class, M.update, U.t, V.t, weightmat, lambda1, lambda2, S.u = S.u)$out
		if(!quiet){
			cat("obj.func after U update=", new.objfncval, "\n")
		} 

	# --- update counter
		k.outer = k.outer + 1

	# ------ check obj.func value
		resid = (old.objfncval - new.objfncval)/abs(init.objfncval)

	# ----- print if not quiet
		if(!quiet){
			cat("k.outer=", k.outer, "resid=", resid,"\n")
		}

		if(resid < full.tol){
		 		outer.converge=TRUE
		}

		if (k.outer > k.iter) {
		 	outer.converge=TRUE
		} 

	# ----- replace old obj.func value 
		old.objfncval = new.objfncval

	# ------ normalize U and V
		out1 = sum(abs(U.t))
	  	V.inner = (V.t*out1)/r
	  	U.inner = (U.t*r)/out1
	  	Uinv.inner = solve(U.inner)
		Vinv.inner = solve(V.inner)

	# -- if not converged, update step size parameter
		if(!outer.converge){
			out1 = sqrt(sum(Uinv.inner^2)) # upper bound largest eigenvalue
			out2 = sqrt(sum(Vinv.inner^2))
			rho = (4/C)*(min(nc)/N)*(out1^(-1)*out2^(-1))/10
		}

	}

	return(list("Mean" = M.update, "U" = U.inner, "V" = V.inner))

}

