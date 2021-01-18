

#########################################################
#### Center inputs given class labels and means #########
#########################################################
### Input: X; an r x p x N ARRAY_CLASS_MEAN           ###
##  	   class; a vector of class labels {1,..., C} ###
## 		   mean; an r x p x C array of class means    ###
#########################################################
## Output: X.c; X centered by class mean              ###
#########################################################

CENTER_X = function(X, class, mean){
	
	p = dim(X)[2]
	r = dim(X)[1]
	N = dim(X)[3]

	X.c = array(0, dim=c(r, p, N))
	
	for (jjj in 1:N) {
		dd = class[jjj]
		X.c[,,jjj] = X[,,jjj] - mean[,,dd]
	}
	
	return(X.c)
}



######################################################
#### Compute U given Vinv and centered Y        ######
######################################################
### Input: Y.c; see above
###		   Vinv; an estimate of  
######################################################
## Output: Y.c; Y centered by class mean 
######################################################

U_SAMPLE = function(X.c, V){

	p = dim(X.c)[2]
	r = dim(X.c)[1]
	N = dim(X.c)[3]

	S = matrix(0, nrow=r, ncol=r)

	for (kkk in 1:N) {
		S = S + tcrossprod(X.c[,,kkk], tcrossprod(X.c[,,kkk], V))
	}

	U = S/(N*p)
	return(U)

}


## given centered data, compute V sample covariance
V_SAMPLE = function(X.c, U){

	p = dim(X.c)[2]
	r = dim(X.c)[1]
	N = dim(X.c)[3]

	S = matrix(0, nrow=p, ncol=p)
	
	for (k in 1:N) {
		S = S + tcrossprod(crossprod(X.c[,,k], U), t(X.c[,,k]))
	}

	V = S/(N*r)

	return(V)
}


STACKED = function(X){ 
	t(apply(apply(X, 2, rbind), 1, cbind)) 
}

ARRAY_CLASS_MEAN = function(X, class){

	r = dim(X)[1]
	p = dim(X)[2]
	N = dim(X)[3]
	C = length(unique(class))

	mean.array = array(0, dim=c(r,p,C))
	for (l in 1:C) {
		mean.array[,,l] = apply(X[,,class==l], c(1,2), mean)
	}

	return(mean.array)
}



##################################################################
###### Algorithm used at the convergence of the AMA  #############
###### to enforce zeros in mean differences based on #############
#####  constraint variables ######################################
##################################################################
### inputs: hatM, hatG -- final iterates from M udpate     #######
###		  : C: number of classes   						   #######
##################################################################
##################################################################
### output: Cr x p matrix of thresholded means             #######
###		  					  							   #######
##################################################################
THRESH_MEAN = function(M, G, C){
	
	K = C*(C-1)/2
	r = dim(M)[1]/C
	p = dim(M)[2]
	
	outM = M
	outG = matrix(0, nrow=r, ncol=p)

	for (q in 1:(C-1)) {
		for (d in (q+1):C) {
			kk = (C-1)*(q-1) + d - 1 - q*(q-1)/2
			tempG = G[((kk-1)*r +1):(kk*r), ]
			for (j in 1:r) {
				for (k in 1:p) {
					if (tempG[j,k] == 0) {
						outM[(d-1)*r + j, k] <- outM[(q-1)*r +j, k]
					}
				}
			}
		}
	}

	return(outM)
}



################################################
## Evaluate objective function #################
################################################

EVAL_OBJ_FUNC = function(X, class, M.update, U.hat, V.hat, weightmat, lambda1, lambda2, 
	S.u = NULL, S.v = NULL){

	#  ------ preliminaries
	r = dim(X)[1]
	p = dim(X)[2]
	N = dim(X)[3]
	C = length(unique(class))
	nc = count(class)[,2]

	# ----- compute determinants
	detU = - p*determinant(U.hat, logarithm=TRUE)$modulus[1] 
	detV = - r*determinant(V.hat, logarithm=TRUE)$modulus[1] 

	# -----  compute penalties
	mus = 0 
	for (q in 1:(C-1)) {
		for (d in (q+1):C) {
			kk = (C-1)*(q-1) + d - 1 - q*(q-1)/2
			mus = mus + sum(weightmat[((kk-1)*r + 1):(kk*r), ]*abs(M.update[,,q] - M.update[,,d]))
		}
	}

	means = lambda1*mus
	inv.covs = lambda2*(sum(abs(U.hat))*sum(abs(V.hat)))

	# -----  evaluate negative log likelihood
	if (is.null(S.u) && is.null(S.v)) {
		L = rep(0, N)
		for (jjj in 1:N) {
			new = X[,,jjj] - M.update[,,class[jjj]]; 
			L[jjj] = sum(diag(tcrossprod(crossprod(U.hat, new), tcrossprod(new, V.hat))))
		}
		L.total = sum(L)
	} else {
		if (!is.null(S.u)) {
			L.total = sum(diag(U.hat%*%S.u*(p*N)))
		} else {
			L.total = sum(diag(S.v%*%V.hat*(r*N)))
		}
	}

	out.val = sum(c(N^(-1)*L.total, detU, detV, means, inv.covs))

	return(list("out" = out.val, "MeanPen" = means, "CovPen" = inv.covs))
}








