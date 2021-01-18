# ********************************************************************************************************************************************************************
# locStra: R and C++ implementations of local population stratification methods.
# 
# This file provides functions to compute:
# 
# 1) the covariance matrix.
# 2) the Jaccard similarity matrix.
# 3) the s-matrix as in the 'Stego' package, available on 'https://github.com/dschlauch/stego'.
# 4) the classic and robust versions of the genomic relationship matrix (grm).
# 5) the largest eigenvector via the power method (von Mises iteration), see 'https://en.wikipedia.org/wiki/Power_iteration'.
# 6) an automated full run of global and local correlations in population stratification data with the function 'fullscan'.
# 7) the windows to be used in the function 'fullscan'.
# 8) the k leading eigenvectors of the covariance matrix, Jaccard matrix, s-matrix, and (classic or robust) genomic relationship matrices without actually computing the similarity matrices. These are the functions "fastCovEVs", "fastJaccardEVs", "fastGrmEVs", and "fastSMatrixEVs".
# 
# All functions work with both standard (dense) R matrices and sparse matrix objects of the 'Matrix' class while never unpacking the sparse matrix during computation.
# All functions have flags to switch between R and C++ implementations, as well as separate implementations for dense or sparse matrices.
# Important note: For all functions the input matrix is always assumed to be oriented to contain the genomic data for one individual per column.
# ********************************************************************************************************************************************************************
#' @useDynLib locStra



# ******************************************************************************************
# auxiliary functions
# ******************************************************************************************

# NOT EXPORTED
# Auxiliary method to convert a sparse matrix of class 'Matrix' into a dense matrix (list) containing rows (i,j,x) to encode an entry x at position (i,j), which are then used as input to the C++ code.
sparseToList <- function(m) {
	L <- summary(Matrix(m,sparse=TRUE))
	# c++ code counts from 0, usually summary() orders according to columns, thus transpose to obtain row sorting
	cbind(L$i-1,L$j-1,L$x)
}



#' Auxiliary function to generate a two-column matrix of windows to be used in the function 'fullscan'.
#' 
#' @param len The overall length of the data which is to be scanned in windows.
#' @param size The window size.
#' @param offset The offset of the generated windows (e.g., if \code{offset=1} then sliding window, if \code{offset=size} then blocks).
#' 
#' @return A two-column matrix of sliding windows, with one window per row defined through start and end value.
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' require(locStra)
#' print(makeWindows(100,10,5))
#' 
#' @export
makeWindows <- function(len,size,offset) {
	if(len<size) return(matrix(c(1,len),nrow=1))
	w <- seq(1,len-(size-1),offset)
	wmatrix <- cbind(w,w+(size-1))
	lastw <- wmatrix[nrow(wmatrix),ncol(wmatrix)]
	if(lastw<len) {
		if(lastw+1<len) wmatrix <- rbind(wmatrix,c(lastw+1,len))
		else wmatrix[nrow(wmatrix),ncol(wmatrix)] <- wmatrix[nrow(wmatrix),ncol(wmatrix)]+1
	}
	unname(wmatrix)
}



#' C++ implementation of the power method (von Mises iteration) to compute the largest eigenvector of a dense input matrix.
#' 
#' @param m Symmetric matrix for which the largest eigenvector is sought.
#' @param initvector Optional vector compatible with the input matrix which serves as a starting value for the iteration. Default is zero.
#' 
#' @return The largest eigenvector of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Richard von Mises and Hilda Pollaczek-Geiringer (1929). Praktische Verfahren der Gleichungsaufloesung. ZAMM Zeitschrift fuer Angewandte Mathematik und Mechanik, 9:152-164.
#' 
#' @examples
#' require(locStra)
#' m <- matrix(1:9,3)
#' print(powerMethod(m))
#' 
#' @export
powerMethod <- function(m,initvector=0) {
	powerMethodCpp(m,initvector)
}



# ******************************************************************************************
# NOT EXPORTED: full R functions for dense and sparse matrices
# ******************************************************************************************

covMatrixR_dense <- function(X) {
	X = t(t(X) - colMeans(X))
	return( (t(X) %*% X)/(nrow(X)-1) )
}

covMatrixR_sparse <- function(X) {
	w <- colSums(X)
	return( (t(X) %*% X - outer(w,w)/nrow(X))/(nrow(X)-1) )
}

jaccardMatrixR <- function(X) {
	w <- colSums(X)
	matrix_and <- t(X) %*% X
	matrix_or <- abs(t(t(matrix_and - w) - w))
	matrix_and[matrix_or==0] <- 1
	matrix_or[matrix_or==0] <- 1
	res <- matrix_and/matrix_or
	return(res)
}

sMatrixR <- function(X, Djac=FALSE, phased=FALSE, minVariants=0) {
	numAlleles <- ifelse(phased, ncol(X), 2*ncol(X))
	sumVariants <- rowSums(X)
	invertMinorAllele <- sumVariants>(numAlleles/2)
	if(phased) {
		X[invertMinorAllele,] <- 1 - X[invertMinorAllele,]
	} else {
		X[invertMinorAllele,] <- 2 - X[invertMinorAllele,]
	}
	
	sumVariants <- rowSums(X)
	X <- X[sumVariants>=minVariants,]
	sumFilteredVariants <- rowSums(X)
	totalPossiblePairs <- numAlleles*(numAlleles-1)/2
	totalPairs <- sumFilteredVariants*(sumFilteredVariants-1)/2
	weights <- ifelse(totalPairs>0,totalPossiblePairs/totalPairs,0)
	if(Djac) {
		s_matrix_numerator <- t(X) %*% X
	} else {
		s_matrix_numerator <- t(X * weights) %*% X
	}
	s_matrix_denominator <- nrow(X)
	s_matrix_hap <- s_matrix_numerator/s_matrix_denominator
	
	if(phased) {
		s_matrix_dip <- (s_matrix_hap[c(T, F), c(T, F)] + s_matrix_hap[c(F,T), c(T, F)] + s_matrix_hap[c(T, F), c(F, T)] + s_matrix_hap[c(F, T), c(F, T)])/4
	} else {
		s_matrix_dip <- s_matrix_hap/4
	}
	return(s_matrix_dip)
}

grMatrixR_dense <- function(X,robust=TRUE) {
	p <- rowMeans(X)/2
	q <- 2*p*(1-p)
	X <- X - 2*p
	if(robust) {
		return( (t(X) %*% X)/sum(q) )
	}
	else {
		return( (t(X) %*% (X/q))/nrow(X) )
	}
}

grMatrixR_sparse <- function(X,robust=TRUE) {
	p <- rowMeans(X)/2
	q <- 2*p*(1-p)
	twop <- 2*p
	if(robust) {
		temp <- as.numeric(t(X) %*% twop)
		return( (t(t(t(X) %*% X - temp) - temp) + sum(twop*twop))/sum(q) )
	}
	else {
		temp <- as.numeric(t(X) %*% (twop/q))
		return( (t(t(t(X) %*% (X/q) - temp) - temp) + sum(twop*twop/q))/nrow(X) )
	}
}



# ******************************************************************************************
# NOT EXPORTED: fast eigenvector computation in R for both dense and sparse matrices
# ******************************************************************************************

# code to compute a randomized SVD according to
# N. Halko, P.G. Martinsson, and J.A. Tropp (2011). Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review: 53(2), pp. 217--288.
# The functions compute the k leading eigenvectors of X^t*X and X*X^t.
# The input is a matrix X=a*v*(A-w), where a is a scalar and v,w are vectors of length nrow(A). Implementations work for both dense and sparse matrices.
randomizedSVD_XtX <- function(a,v,A,w,k,q=2) {
	k0 <- k
	k <- max(k,2)
	Y <- matrix(rnorm(nrow(A)*2*k),nrow=nrow(A))
	Y <- a*( t(A) %*% (v*Y) - outer(rep(1,ncol(A)),colSums(v*w*Y)) )
	for(i in 1:q) {
		Y <- a*v*(A %*% Y - outer(w,colSums(Y)))
		Y <- a*( t(A) %*% (v*Y) - outer(rep(1,ncol(A)),colSums(v*w*Y)) )
	}
	temp <- qr(Y)
	Q <- qr.Q(temp)
	B <- t( a*v*(A %*% Q - outer(w,colSums(Q))) )
	temp <- svds(B,k=k,nu=k,nv=0)
	U <- Q %*% temp$u
	return(U[,1:k0])
}

randomizedSVD_XXt <- function(a,v,A,w,k,q=2) {
	k0 <- k
	k <- max(k,2)
	Y <- matrix(rnorm(ncol(A)*2*k),nrow=ncol(A))
	Y <- a*v*(A %*% Y - outer(w,colSums(Y)))
	for(i in 1:q) {
		Y <- a*( t(A) %*% (v*Y) - outer(rep(1,ncol(A)),colSums(v*w*Y)) )
		Y <- a*v*(A %*% Y - outer(w,colSums(Y)))
	}
	temp <- qr(Y)
	Q <- qr.Q(temp)
	B <- a*( t(Q*v) %*% A - colSums(Q*v*w) )
	temp <- svds(B,k=k,nu=k,nv=0)
	U <- Q %*% temp$u
	return(U[,1:k0])
}

fastCovEVsR <- function(X, k, q=2) {
	a <- 1/sqrt(nrow(X)-1)
	v <- 1
	w <- colMeans(X)
	randomizedSVD_XXt(a=a,v=v,A=t(X),w=w,k=k,q=q)
}

fastJaccardEVsR <- function(X, k, q=2) {
	a <- 1/sqrt(2*max(colSums(X)))
	v <- 1
	w <- rep(0,nrow(X))
	randomizedSVD_XtX(a=a,v=v,A=X,w=w,k=k,q=q)
}

# always phased=FALSE
fastSMatrixEVsR <- function(X, k, Djac=FALSE, minVariants=0, q=2) {
	numAlleles <- 2*ncol(X)
	sumVariants <- rowSums(X)
	invertMinorAllele <- sumVariants>(numAlleles/2)
	X[invertMinorAllele,] <- 2 - X[invertMinorAllele,]
	
	sumVariants <- rowSums(X)
	X <- X[sumVariants>=minVariants,]
	sumFilteredVariants <- rowSums(X)
	totalPossiblePairs <- numAlleles*(numAlleles-1)/2
	totalPairs <- sumFilteredVariants*(sumFilteredVariants-1)/2
	weights <- ifelse(totalPairs>0,totalPossiblePairs/totalPairs,0)
	s_matrix_denominator <- nrow(X)
	
	if(Djac) {
		a <- 1/sqrt(4*s_matrix_denominator)
		v <- 1
		w <- rep(0,nrow(X))
		randomizedSVD_XtX(a=a,v=v,A=X,w=w,k=k,q=q)
	} else {
		a <- 1/sqrt(4*s_matrix_denominator)
		v <- sqrt(weights)
		w <- rep(0,nrow(X))
		randomizedSVD_XtX(a=a,v=v,A=X,w=w,k=k,q=q)
	}
}

fastGrmEVsR <- function(X, k, robust=TRUE, q=2) {
	# compute population frequencies across rows
	p <- rowMeans(X)/2
	qv <- 2*p*(1-p)
	# compute grm
	if(robust) {
		a <- 1/sqrt(sum(qv))
		v <- 1
		w <- 2*p
		randomizedSVD_XtX(a=a,v=v,A=X,w=w,k=k,q=q)
	}
	else {
		a <- 1/sqrt(nrow(X))
		v <- 1/sqrt(qv)
		w <- 2*p
		randomizedSVD_XtX(a=a,v=v,A=X,w=w,k=k,q=q)
	}
}



# ******************************************************************************************
# main functions: similarity matrices and 'fullscan'
# ******************************************************************************************

#' C++ implementation to compute the covariance matrix for a (sparse) input matrix. The function is equivalent to the R command 'cov' applied to matrices.
#' 
#' @param m A (sparse) matrix for which the covariance matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' 
#' @return The covariance matrix of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references R Core Team (2014). R: A Language and Environment for Statistical Computing. R Foundation for Stat Comp, Vienna, Austria.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(covMatrix(sparseM))
#' 
#' @export
covMatrix <- function(m,useCpp=TRUE,sparse=TRUE) {
	if(useCpp) {
		if(sparse) covMatrixCpp_sparse(sparseToList(m),nrow(m),ncol(m))
		else covMatrixCpp_dense(as.matrix(m))
	}
	else {
		if(sparse) covMatrixR_sparse(m)
		else covMatrixR_dense(m)
	}
}



#' C++ implementation to compute the Jaccard similarity matrix for a (sparse) input matrix.
#' 
#' @param m A (sparse) matrix for which the Jaccard similarity matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' 
#' @return The Jaccard matrix of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Dmitry Prokopenko, Julian Hecker, Edwin Silverman, Marcello Pagano, Markus Noethen, Christian Dina, Christoph Lange and Heide Fier (2016). Utilizing the Jaccard index to reveal population stratification in sequencing data: a simulation study and an application to the 1000 Genomes Project. Bioinformatics, 32(9):1366-1372.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(jaccardMatrix(sparseM))
#' 
#' @export
jaccardMatrix <- function(m,useCpp=TRUE,sparse=TRUE) {
	if(useCpp) {
		if(sparse) jaccardMatrixCpp_sparse(sparseToList(m),nrow(m),ncol(m))
		else jaccardMatrixCpp_dense(as.matrix(m))
	}
	else {
		jaccardMatrixR(m)
	}
}



#' C++ implementation to compute the s-matrix (the weighted Jaccard similarity matrix) for a (sparse) input matrix as in the 'Stego' package: https://github.com/dschlauch/stego
#' 
#' @param m A (sparse) matrix for which the s-matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' @param Djac Flag to switch between the unweighted (\code{Djac=TRUE}) or weighted (\code{Djac=FALSE}) version. Default is \code{Djac=FALSE}.
#' @param phased Boolean flag to indicate if the input matrix is phased. Default is \code{phased=FALSE}.
#' @param minVariants Integer cutoff value for minimal number of variants. Default is \code{minVariants=0}.
#' 
#' @return The s-matrix (the weighted Jaccard matrix) of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Daniel Schlauch (2016). Implementation of the stego algorithm - Similarity Test for Estimating Genetic Outliers. https://github.com/dschlauch/stego
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(sMatrix(sparseM))
#' 
#' @export
sMatrix <- function(m,useCpp=TRUE,sparse=TRUE,Djac=FALSE,phased=FALSE,minVariants=0) {
	if(useCpp) {
		if(sparse) sMatrixCpp_sparse(sparseToList(m),nrow(m),ncol(m),Djac,phased,minVariants)
		else sMatrixCpp_dense(as.matrix(m),Djac,phased,minVariants)
	}
	else {
		sMatrixR(m,Djac,phased,minVariants)
	}
}



#' C++ implementation to compute the genomic relationship matrix (grm) for a (sparse) input matrix as defined in Yang et al. (2011).
#' 
#' @param m A (sparse) matrix for which the genomic relationship matrix is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' @param robust Flag to indicate if the classic (\code{robust=FALSE}) or robust (\code{robust=TRUE}) version of the genomic relationship matrix is desired. Default is \code{robust=TRUE}.
#' 
#' @return The genomic relationship matrix of \code{m}.
#' 
#' @importFrom Rdpack reprompt
#' @references Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet, 88(1):76-82.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,15,replace=TRUE),ncol=3)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(grMatrix(sparseM))
#' 
#' @export
grMatrix <- function(m,useCpp=TRUE,sparse=TRUE,robust=TRUE) {
	if(useCpp) {
		if(sparse) grmCpp_sparse(sparseToList(m),nrow(m),ncol(m),robust)
		else grmCpp_dense(as.matrix(m),robust)
	}
	else {
		if(sparse) grMatrixR_sparse(m,robust)
		else grMatrixR_dense(m,robust)
	}
}



#' A full scan of the input data \code{m} using a collection of windows given by the two-column matrix \code{windows}. For each window, the data is processed using the function \code{matrixFunction} (this could be, e.g., the \code{covMatrix} function), then the processed data is summarized using the function \code{summaryFunction} (e.g., the largest eigenvector computed with the function \code{powerMethod}), and finally the global and local summaries are compared using the function \code{comparisonFunction} (e.g., the vector correlation with R's function \code{cor}). The function returns a two-column matrix which contains per row the global summary statistics (e.g., the correlation between the global and local eigenvectors) and the local summary statistics (e.g., the correlation between the local eigenvectors of the previous and current windows) for each window.
#' 
#' @param m A (sparse) matrix for which the full scan is sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param windows A two-column matrix containing per column the windows on which the data is scanned. The windows can be overlapping. The windows can be computed using the function \code{makeWindows}.
#' @param matrixFunction Function on one matrix argument to process the data for each window (e.g., the covariance matrix).
#' @param summaryFunction Function on one argument to summarize the output of the function \code{matrixFunction} (e.g., the largest eigenvector).
#' @param comparisonFunction Function on two inputs to compute a comparison measure for the output of the function \code{summaryFunction} (e.g., vector correlation, or matrix norm).
#' 
#' @return A two-column matrix containing per row the global and local summary statistics for each window. Plotting the correlation data of the returned matrix gives a figure analogously to the figure shown here, which was generated with the example code below.
#' 
#' \figure{fig.pdf}{options: width=5in}
#' 
#' @importFrom Rdpack reprompt
#' @references Dmitry Prokopenko, Julian Hecker, Edwin Silverman, Marcello Pagano, Markus Noethen, Christian Dina, Christoph Lange and Heide Fier (2016). Utilizing the Jaccard index to reveal population stratification in sequencing data: a simulation study and an application to the 1000 Genomes Project. Bioinformatics, 32(9):1366-1372.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' data(testdata)
#' cor2 <- function(x,y) ifelse(sum(x)==0 | sum(y)==0, 0, cor(x,y))
#' windowSize <- 10000
#' w <- makeWindows(nrow(testdata),windowSize,windowSize)
#' resCov <- fullscan(testdata,w,covMatrix,powerMethod,cor2)
#' resJac <- fullscan(testdata,w,jaccardMatrix,powerMethod,cor2)
#' resSMx <- fullscan(testdata,w,sMatrix,powerMethod,cor2)
#' resGRM <- fullscan(testdata,w,grMatrix,powerMethod,cor2)
#' resAll <- cbind(resCov[,1], resJac[,1], resSMx[,1], resGRM[,1])
#' xlabel <- "SNP position"
#' ylabel <- "correlation between global and local eigenvectors"
#' mainlabel <- paste("window size",windowSize)
#' matplot(w[,1],abs(resAll),type="b",xlab=xlabel,ylab=ylabel,ylim=c(0,1),main=mainlabel)
#' legend("topright",legend=c("Cov","Jaccard","s-Matrix","GRM"),pch=paste(1:ncol(resAll)))
#' 
#' @export
fullscan <- function(m,windows,matrixFunction,summaryFunction,comparisonFunction) {
	matrix_global <- matrixFunction(m)
	summary_global <- summaryFunction(matrix_global)
	last_summary <- rep(0,length(summary_global))
	
	# go through sliding window
	res <- matrix(0, nrow=nrow(windows), ncol=2)
	for(i in 1:nrow(windows)) {
		matrix_local <- matrixFunction(m[windows[i,1]:windows[i,2],])
		summary_local <- summaryFunction(matrix_local)
		res[i,] <- c(comparisonFunction(summary_global,summary_local), comparisonFunction(last_summary,summary_local))
		last_summary <- summary_local
	}
	return(res)
}



# ******************************************************************************************
# main functions: fast eigenvector computation of similarity matrices
# ******************************************************************************************

#' Computation of the k leading eigenvectors of the covariance matrix for a (sparse) input matrix.
#' 
#' @param m A (sparse) matrix for which the eigenvectors of its covariance matrix are sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param k The number of leading eigenvectors.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' @param q The number of power iteration steps (default is \code{q=2}).
#' 
#' @return The k leading eigenvectors of the covariance matrix of \code{m} as a column matrix.
#' 
#' @importFrom Rdpack reprompt
#' @references R Core Team (2014). R: A Language and Environment for Statistical Computing. R Foundation for Stat Comp, Vienna, Austria.
#' @references N. Halko, P.G. Martinsson, and J.A. Tropp (2011). Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review: 53(2), pp. 217--288.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,100,replace=TRUE),ncol=5)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(fastCovEVs(sparseM,k=2,useCpp=FALSE))
#' 
#' @export
fastCovEVs <- function(m,k,useCpp=TRUE,sparse=TRUE,q=2) {
	if(useCpp) {
		if(sparse) fastCovEVsCpp_sparse(sparseToList(m),nrow(m),ncol(m),k,q)
		else fastCovEVsCpp_dense(as.matrix(m),k,q)
	}
	else {
		fastCovEVsR(m,k,q)
	}
}



#' Computation of the k leading eigenvectors of the Jaccard similarity matrix for a (sparse) input matrix. Note that this computation is only approximate and does not necessarily coincide with the result obtained by extracting the k leading eigenvectors of the Jaccard matrix computed with the function \code{jaccardMatrix}.
#' 
#' @param m A (sparse) matrix for which the eigenvectors of its Jaccard matrix are sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param k The number of leading eigenvectors.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' @param q The number of power iteration steps (default is \code{q=2}).
#' 
#' @return The k leading eigenvectors of the Jaccard matrix of \code{m} as a column matrix.
#' 
#' @importFrom Rdpack reprompt
#' @references Dmitry Prokopenko, Julian Hecker, Edwin Silverman, Marcello Pagano, Markus Noethen, Christian Dina, Christoph Lange and Heide Fier (2016). Utilizing the Jaccard index to reveal population stratification in sequencing data: a simulation study and an application to the 1000 Genomes Project. Bioinformatics, 32(9):1366-1372.
#' @references N. Halko, P.G. Martinsson, and J.A. Tropp (2011). Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review: 53(2), pp. 217--288.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,100,replace=TRUE),ncol=5)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(fastJaccardEVs(sparseM,k=2,useCpp=FALSE))
#' 
#' @export
fastJaccardEVs <- function(m,k,useCpp=TRUE,sparse=TRUE,q=2) {
	if(useCpp) {
		if(sparse) fastJaccardEVsCpp_sparse(sparseToList(m),nrow(m),ncol(m),k,q)
		else fastJaccardEVsCpp_dense(as.matrix(m),k,q)
	}
	else {
		fastJaccardEVsR(m,k,q)
	}
}



#' Computation of the k leading eigenvectors of the s-matrix (the weighted Jaccard similarity matrix) for a (sparse) input matrix. Note that in contrast to the parameters of the function \code{sMatrix}, the choice \code{phased=FALSE} cannot be modified for the fast eigenvector computation.
#' 
#' @param m A (sparse) matrix for which the eigenvectors of its s-matrix are sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param k The number of leading eigenvectors.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' @param Djac Flag to switch between the unweighted (\code{Djac=TRUE}) or weighted (\code{Djac=FALSE}) version. Default is \code{Djac=FALSE}.
#' @param minVariants Integer cutoff value for minimal number of variants. Default is \code{minVariants=0}.
#' @param q The number of power iteration steps (default is \code{q=2}).
#' 
#' @return The k leading eigenvectors of the s-matrix of \code{m} as a column matrix.
#' 
#' @importFrom Rdpack reprompt
#' @references Daniel Schlauch (2016). Implementation of the stego algorithm - Similarity Test for Estimating Genetic Outliers. https://github.com/dschlauch/stego
#' @references N. Halko, P.G. Martinsson, and J.A. Tropp (2011). Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review: 53(2), pp. 217--288.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,100,replace=TRUE),ncol=5)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(fastSMatrixEVs(sparseM,k=2,useCpp=FALSE))
#' 
#' @export
fastSMatrixEVs <- function(m,k,useCpp=TRUE,sparse=TRUE,Djac=FALSE,minVariants=0,q=2) {
	if(useCpp) {
		if(sparse) fastSMatrixEVsCpp_sparse(sparseToList(m),nrow(m),ncol(m),k,Djac,minVariants,q)
		else fastSMatrixEVsCpp_dense(as.matrix(m),k,Djac,minVariants,q)
	}
	else {
		fastSMatrixEVsR(m,k,Djac,minVariants,q)
	}
}



#' Computation of the k leading eigenvectors of the genomic relationship matrix, defined in Yang et al. (2011), for a (sparse) input matrix.
#' 
#' @param m A (sparse) matrix for which the eigenvectors of its genomic relationship matrix are sought. The input matrix is assumed to be oriented to contain the data for one individual per column.
#' @param k The number of leading eigenvectors.
#' @param useCpp Flag to switch between R or C++ implementations. Default is \code{useCpp=TRUE}.
#' @param sparse Flag to switch between purpose-built dense or sparse implementations. Default is \code{sparse=TRUE}.
#' @param robust Flag to indicate if the classic (\code{robust=FALSE}) or robust (\code{robust=TRUE}) version of the genomic relationship matrix is desired. Default is \code{robust=TRUE}.
#' @param q The number of power iteration steps (default is \code{q=2}).
#' 
#' @return The k leading eigenvectors of the genomic relationship matrix of \code{m} as a column matrix.
#' 
#' @importFrom Rdpack reprompt
#' @references Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet, 88(1):76-82.
#' @references N. Halko, P.G. Martinsson, and J.A. Tropp (2011). Finding Structure with Randomness: Probabilistic Algorithms for Constructing Approximate Matrix Decompositions. SIAM Review: 53(2), pp. 217--288.
#' 
#' @examples
#' require(locStra)
#' require(Matrix)
#' m <- matrix(sample(0:1,100,replace=TRUE),ncol=5)
#' sparseM <- Matrix(m,sparse=TRUE)
#' print(fastGrmEVs(sparseM,k=2,useCpp=FALSE))
#' 
#' @export
fastGrmEVs <- function(m,k,useCpp=TRUE,sparse=TRUE,robust=TRUE,q=2) {
	if(useCpp) {
		if(sparse) fastGrmEVsCpp_sparse(sparseToList(m),nrow(m),ncol(m),k,robust,q)
		else fastGrmEVsCpp_dense(m,k,robust,q)
	}
	else {
		fastGrmEVsR(m,k,robust,q)
	}
}
