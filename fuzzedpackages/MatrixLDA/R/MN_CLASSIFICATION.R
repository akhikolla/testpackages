# ----------------------------------------------------------
# Function to classify new observations via MN linear discrimiant 
# ----------------------------------------------------------
MN_CLASSIFICATION = function(X, class=NULL, M, U, V, pi.list, C=NULL){

	N = dim(X)[3]
	if (is.null(C)) {
		C = dim(M)[3]
	}

	class.mat = matrix(0, nrow=N, ncol=C)
	pred.class = rep(0, N)

	for (j in 1:N) {
		for (k in 1:C) {
			new = X[,,j] - M[,,k]
			class.mat[j,k] = .5*sum(rowSums(crossprod(U, new)*t(tcrossprod(V, new))))  - log(pi.list[k])
		}
		pred.class[j] = which.min(class.mat[j,])
	}

	if(!is.null(class)){
		misclass = sum(pred.class != class)/N
	} else {
			misclass = NA
	}

	return(list("pred.class" = pred.class, "misclass" = misclass, "prob.mat" = class.mat))

}


#' Predict_MN
#' @param object An object of type "MN"; output MatLDA or MN_MLE
#' @param newdata New data to be classified; array of r x p x Ntest
#' @param newclass Class labels for new data if available
#' @return pred.class a vector of length N_test of predicted class memberships
#' @return misclass if newclass is specified, the misclassification rate
#' @return prob.mat a matrix of N_test x C with value of discriminant function evaluted at each test data point
#' @export
#'
PredictMN = function(object, newdata, newclass=NULL){

	U = object$U
	V = object$V
	M = object$Mean
	pi.list = object$pi.list
	out = MN_CLASSIFICATION(newdata, class=newclass, M = M, U=U, V=V, pi.list=pi.list)

	return(out)

}
