#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;

/*
using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXi;
using Eigen::MatrixXd;
using Eigen::SparseMatrix;
*/

#include <random>
std::normal_distribution<double> normalDistribution(0.0,1.0);
std::random_device rd{};
std::mt19937 gen{rd()};



// ************************************************************************************
// power method
// ************************************************************************************

// power method to compute the largest eigenvalue, set v=0 for automatic rowStarting value
// [[Rcpp::export]]
VectorXd powerMethodCpp(MatrixXd& X, VectorXd& v, double eps=1e-6, int maxiter=100) {
	if(X.rows()!=X.cols()) {
		Rcpp::stop("powerMethod requires square numeric matrix.");
	}
	
	VectorXd v_old(v);
	if(v_old.size()<X.rows()) {
		v_old = VectorXd::Zero(X.rows()).array()+1;
	}
	
	if(v_old.size()!=X.rows()) {
		Rcpp::stop("powerMethod requires X and v to be compatible.");
	}
	
	if(eps<=0) {
		eps = 1e-6;
	}
	
	VectorXd v_new;
	VectorXd absdiff;
	for(int steps=0; steps<maxiter; steps++) {
		v_new = X * v_old;
		v_new.normalize();
		absdiff = v_new.array().abs() - v_old.array().abs();
		if(absdiff.norm() <= eps) {
			break;
		}
		v_old = v_new;
		steps = steps + 1;
	}
	return v_new;
}



// ************************************************************************************
// auxiliary functions
// ************************************************************************************

// cpp function to write a sparse matrix with integer entries
SparseMatrix<int> triplesToSparseIntMatrix(MatrixXi& T, int nrows, int ncols) {
	SparseMatrix<int> X(nrows,ncols);
	for(int i=0; i<T.rows(); i++) {
		X.insert(T(i,0),T(i,1)) = T(i,2);
	}
	return X;
}

// cpp function to write a sparse matrix with double entries
SparseMatrix<double> triplesToSparseDoubleMatrix(MatrixXd& T, int nrows, int ncols) {
	SparseMatrix<double> X(nrows,ncols);
	for(int i=0; i<T.rows(); i++) {
		X.insert(int(T(i,0)),int(T(i,1))) = T(i,2);
	}
	return X;
}

// multiply columns of sparse matrix with "v"
SparseMatrix<double> colMultiplySparseDoubleMatrix(SparseMatrix<double>& X, VectorXd& v) {
	SparseMatrix<double> Y(X.rows(),X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<double>::InnerIterator it(X,k); it; ++it) {
			Y.insert(it.row(),it.col()) = it.value()*v(it.row());
		}
	}
	return Y;
}

// subset of rows according to indicator vector "v" for dense matrix
MatrixXd subsetRowsDenseMatrix(MatrixXd& X, VectorXd& v) {
	MatrixXd res = MatrixXd::Zero(int(v.sum()),X.cols());
	int rowcounter = 0;
	for(int i=0; i<v.size(); i++) {
		if(int(v(i))==1) {
			res.row(rowcounter) = X.row(i);
			rowcounter++;
		}
	}
	return res;
}

// subset of rows according to indicator vector "v" for sparse matrix
SparseMatrix<double> subsetRowsSparseIntMatrix(SparseMatrix<int>& X, VectorXd& v) {
	// assign new row numbers
	VectorXi w = VectorXi::Zero(v.size());
	int counter = 0;
	for(int i=0; i<v.size(); i++) {
		if(int(v(i))==1) {
			w(i) = counter;
			counter++;
		}
	}
	// now copy only the rows in w
	SparseMatrix<double> Y(counter,X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			if(int(v(it.row()))==1) Y.insert(w(it.row()),it.col()) = it.value();
		}
	}
	return Y;
}

// submatrix of dense matrix according to indicator vectors "rows" and "cols"
MatrixXd subDenseMatrix(MatrixXd& X, VectorXd& rows, VectorXd& cols) {
	MatrixXd res = MatrixXd::Zero(int(rows.sum()), int(cols.sum()));
	int rowcounter = 0;
	int colcounter;
	for(int i=0; i<X.rows(); i++) {
		if(int(rows(i))==1) {
			colcounter = 0;
			for(int j=0; j<X.cols(); j++) {
				if(int(cols(j))==1) {
					res(rowcounter,colcounter) = X(i,j);
					colcounter++;
				}
			}
			rowcounter++;
		}
	}
	return res;
}

// row sums for sparse integer matrix
VectorXd rowSumsSparseIntMatrix(SparseMatrix<int>& X) {
	VectorXd v = VectorXd::Zero(X.rows());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			v[it.row()] += it.value();
		}
	}
	return v;
}

// row sums for sparse double matrix
VectorXd rowSumsSparseDoubleMatrix(SparseMatrix<double>& X) {
	VectorXd v = VectorXd::Zero(X.rows());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<double>::InnerIterator it(X,k); it; ++it) {
			v[it.row()] += it.value();
		}
	}
	return v;
}

// column sums for sparse integer matrix
VectorXd colSumsSparseIntMatrix(SparseMatrix<int>& X) {
	VectorXd v = VectorXd::Zero(X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			v[it.col()] += it.value();
		}
	}
	return v;
}

// column sums for sparse double matrix
VectorXd colSumsSparseDoubleMatrix(SparseMatrix<double>& X) {
	VectorXd v = VectorXd::Zero(X.cols());
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<double>::InnerIterator it(X,k); it; ++it) {
			v[it.col()] += it.value();
		}
	}
	return v;
}

MatrixXd normalRandomMatrix(int nrows, int ncols) {
	MatrixXd temp = MatrixXd::Zero(nrows,ncols);
	for(int i=0;i<nrows;i++) {
		for(int j=0;j<ncols;j++) {
			temp(i,j) = normalDistribution(gen);
		}
	}
	return temp;
}



// ************************************************************************************
// dense matrix functions
// ************************************************************************************

// cpp function to compute the covariance matrix "cov" in R
// [[Rcpp::export]]
MatrixXd covMatrixCpp_dense(MatrixXd& X) {
	// scale: subtract column means from rows
	MatrixXd s = X.rowwise() - X.colwise().mean();
	// return scaled cross product
	return s.transpose() * s * 1.0/(X.rows()-1.0);
}

// fast vectorized version of the Jaccard matrix computation
// [[Rcpp::export]]
MatrixXd jaccardMatrixCpp_dense(MatrixXd& X) {
	VectorXd colsums = X.colwise().sum();
	MatrixXd matrix_and = X.transpose() * X;
	MatrixXd matrix_or = ((matrix_and.rowwise() - colsums.transpose()).colwise() - colsums).cwiseAbs();
	// go through all entries of "or" matrix and if zero then ensure quotient will be one
	for(int i=0; i<matrix_or.rows(); i++) {
		for(int j=0; j<matrix_or.cols(); j++) {
			if(int(matrix_or(i,j))==0) {
				matrix_and(i,j) = 1.0;
				matrix_or(i,j) = 1.0;
			}
		}
	}
	return matrix_and.cwiseQuotient(matrix_or);
}

// cpp function to compute the s-matrix
// [[Rcpp::export]]
MatrixXd sMatrixCpp_dense(MatrixXd X, bool Djac=false, bool phased=false, int minVariants=0) {
	double numAlleles = X.cols();
	if(!phased) numAlleles *= 2.0;
	VectorXd sumVariants = X.rowwise().sum();
	
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>0.5*numAlleles) {
			if(phased) X.row(i) = 1.0 - X.row(i).array();
			else X.row(i) = 2.0 - X.row(i).array();
		}
	}
	
	sumVariants = X.rowwise().sum();
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>=minVariants) sumVariants(i) = 1.0;
		else sumVariants(i) = 0.0;
	}
	MatrixXd Y = subsetRowsDenseMatrix(X,sumVariants);
	
	VectorXd sumFilteredVariants = Y.rowwise().sum();
	double totalPossiblePairs = numAlleles*(numAlleles-1.0)/2.0;
	VectorXd totalPairs = sumFilteredVariants.array()*(sumFilteredVariants.array()-1.0)/2.0;
	VectorXd weights = VectorXd::Zero(totalPairs.size());
	for(int i=0; i<totalPairs.size(); i++) {
		if(totalPairs(i)>0.0) weights(i) = totalPossiblePairs/totalPairs(i);
		else weights(i) = 0.0;
	}
	
	MatrixXd s_matrix_numerator;
	if(Djac) {
		s_matrix_numerator = Y.transpose() * Y;
	} else {
		MatrixXd Z = Y.array().colwise() * weights.array();
		s_matrix_numerator = Z.transpose() * Y;
	}
	MatrixXd s_matrix_hap = s_matrix_numerator * (1.0/Y.rows());
	
	MatrixXd s_matrix_dip;
	if(phased) {
		// alternative binary vectors for rows and columns
		VectorXd rowsTF = VectorXd::Zero(s_matrix_hap.rows());
		for(int i=0; i<rowsTF.size(); i++) rowsTF(i) = (i+1)%2;
		VectorXd colsTF = VectorXd::Zero(s_matrix_hap.cols());
		for(int i=0; i<colsTF.size(); i++) colsTF(i) = (i+1)%2;
		VectorXd rowsFT = 1.0-rowsTF.array();
		VectorXd colsFT = 1.0-colsTF.array();
		s_matrix_dip = ( subDenseMatrix(s_matrix_hap,rowsTF,colsTF) + subDenseMatrix(s_matrix_hap,rowsFT,colsTF) +
							subDenseMatrix(s_matrix_hap,rowsTF,colsFT) + subDenseMatrix(s_matrix_hap,rowsFT,colsFT) )/4.0;
	} else {
		s_matrix_dip = s_matrix_hap/4.0;
	}
	return s_matrix_dip;
}

// dense version of the genomic relationship matrix: Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet, 88(1):76-82.
// [[Rcpp::export]]
MatrixXd grmCpp_dense(MatrixXd X, bool robust=true) {
	// compute population frequencies across rows
	VectorXd p = 0.5 * X.rowwise().mean();
	VectorXd counterp = 1.0-p.array();
	VectorXd q = 2.0*p.cwiseProduct(counterp);
	// compute grm
	X = X.colwise() - 2.0*p;
	if(robust) return X.transpose() * X * 1.0/q.sum();
	else return X.transpose() * (X.array().colwise() * q.array().inverse()).matrix() * 1.0/X.rows();
}



// ************************************************************************************
// sparse matrix functions
// 
// the input to sparse functions is always a dense n*3 index matrix which contains all
// non-zero entries as (x,y,value) tripels per row
// ************************************************************************************

// cpp function to compute the covariance matrix "cov" in R
// [[Rcpp::export]]
MatrixXd covMatrixCpp_sparse(MatrixXd& T, int nrows, int ncols) {
	SparseMatrix<double> X = triplesToSparseDoubleMatrix(T,nrows,ncols);
	// return cov matrix
	VectorXd w = colSumsSparseDoubleMatrix(X);
	MatrixXd temp = X.transpose()*X;
	return 1.0/(X.rows()-1.0) * ( temp - w*w.transpose()*1.0/X.rows() );
}

// fast vectorized version of the Jaccard matrix computation
// [[Rcpp::export]]
MatrixXd jaccardMatrixCpp_sparse(MatrixXi& T, int nrows, int ncols) {
	SparseMatrix<int> X = triplesToSparseIntMatrix(T,nrows,ncols);
	VectorXd colsums = colSumsSparseIntMatrix(X);
	MatrixXd matrix_and = (X.transpose() * X).cast<double>();
	MatrixXd matrix_or = ((matrix_and.rowwise() - colsums.transpose()).colwise() - colsums).cwiseAbs();
	// go through all entries of "or" matrix and if zero then ensure quotient will be one
	for(int i=0; i<matrix_or.rows(); i++) {
		for(int j=0; j<matrix_or.cols(); j++) {
			if(int(matrix_or(i,j))==0) {
				matrix_and(i,j) = 1.0;
				matrix_or(i,j) = 1.0;
			}
		}
	}
	return matrix_and.cwiseQuotient(matrix_or);
}

// cpp function to compute the s-matrix
// [[Rcpp::export]]
MatrixXd sMatrixCpp_sparse(MatrixXi& T, int nrows, int ncols, bool Djac=false, bool phased=false, int minVariants=0) {
	SparseMatrix<int> X = triplesToSparseIntMatrix(T,nrows,ncols);
	
	double numAlleles = X.cols();
	if(!phased) numAlleles *= 2.0;
	VectorXd sumVariants = rowSumsSparseIntMatrix(X);
	
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			if(sumVariants(it.row())>0.5*numAlleles) {
				if(phased) X.coeffRef(it.row(),it.col()) = 1.0 - it.value();
				else X.coeffRef(it.row(),it.col()) = 2.0 - it.value();
			}
		}
	}
	
	sumVariants = rowSumsSparseIntMatrix(X);
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>=minVariants) sumVariants(i) = 1.0;
		else sumVariants(i) = 0.0;
	}
	SparseMatrix<double> Y = subsetRowsSparseIntMatrix(X,sumVariants);
	
	VectorXd sumFilteredVariants = rowSumsSparseDoubleMatrix(Y);
	double totalPossiblePairs = numAlleles*(numAlleles-1.0)/2.0;
	VectorXd totalPairs = sumFilteredVariants.array()*(sumFilteredVariants.array()-1.0)/2.0;
	VectorXd weights = VectorXd::Zero(totalPairs.size());
	for(int i=0; i<totalPairs.size(); i++) {
		if(totalPairs(i)>0.0) weights(i) = totalPossiblePairs/totalPairs(i);
		else weights(i) = 0.0;
	}
	
	MatrixXd s_matrix_numerator;
	if(Djac) {
		s_matrix_numerator = Y.transpose() * Y;
	} else {
		SparseMatrix<double> Z = colMultiplySparseDoubleMatrix(Y,weights);
		s_matrix_numerator = Z.transpose() * Y;
	}
	MatrixXd s_matrix_hap = s_matrix_numerator * (1.0/Y.rows());
	
	MatrixXd s_matrix_dip;
	if(phased) {
		// alternative binary vectors for rows and columns
		VectorXd rowsTF = VectorXd::Zero(s_matrix_hap.rows());
		for(int i=0; i<rowsTF.size(); i++) rowsTF(i) = (i+1)%2;
		VectorXd colsTF = VectorXd::Zero(s_matrix_hap.cols());
		for(int i=0; i<colsTF.size(); i++) colsTF(i) = (i+1)%2;
		VectorXd rowsFT = 1.0-rowsTF.array();
		VectorXd colsFT = 1.0-colsTF.array();
		s_matrix_dip = ( subDenseMatrix(s_matrix_hap,rowsTF,colsTF) + subDenseMatrix(s_matrix_hap,rowsFT,colsTF) +
							subDenseMatrix(s_matrix_hap,rowsTF,colsFT) + subDenseMatrix(s_matrix_hap,rowsFT,colsFT) )/4.0;
	} else {
		s_matrix_dip = s_matrix_hap/4.0;
	}
	return s_matrix_dip;
}

// sparse version of the genomic relationship matrix: Yang J, Lee SH, Goddard ME, Visscher PM (2011). GCTA: a tool for genome-wide complex trait analysis. Am J Hum Genet, 88(1):76-82.
// [[Rcpp::export]]
MatrixXd grmCpp_sparse(MatrixXd& T, int nrows, int ncols, bool robust=true) {
	SparseMatrix<double> X = triplesToSparseDoubleMatrix(T,nrows,ncols);
	// compute population frequencies across rows
	VectorXd p = 0.5 * rowSumsSparseDoubleMatrix(X) * 1.0/X.cols();
	VectorXd counterp = 1.0-p.array();
	VectorXd q = 2.0*p.cwiseProduct(counterp);
	// compute grm
	VectorXd twop = 2.0*p;
	if(robust) {
		MatrixXd temp = X.transpose()*X;
		return ( ( (temp.colwise() - X.transpose()*twop).rowwise() - twop.transpose()*X ).array() + twop.dot(twop) ) * 1.0/q.sum();
	}
	else {
		VectorXd invq = q.array().inverse();
		SparseMatrix<double> Y = colMultiplySparseDoubleMatrix(X,invq);
		MatrixXd temp = X.transpose() * Y;
		return ( ( (temp.colwise() - X.transpose()*(twop.cwiseProduct(invq))).rowwise() - twop.transpose()*Y ).array() + twop.dot(twop.cwiseProduct(invq)) ) * 1.0/X.rows();
	}
}



// ************************************************************************************
// fast eigenvector computation for dense matrices
// ************************************************************************************

MatrixXd randomizedSVD_XtX_Cpp_dense(double a, VectorXd& v, MatrixXd& A, VectorXd& w, int k, int q=2) {
	// matrix with normally distributed entries
	MatrixXd Y = normalRandomMatrix(A.rows(),2*k);
	
	MatrixXd temp1 = Y.array().colwise() * v.array();
	MatrixXd temp2 = Y.array().colwise() * v.cwiseProduct(w).array();
	Y = a*( (A.transpose()*temp1).rowwise() - temp2.colwise().sum() );
	
	for(int i=0; i<q; i++) {
		Y = (A*Y - w*Y.colwise().sum()).array().colwise() * (a*v).array();
		
		temp1 = Y.array().colwise() * v.array();
		temp2 = Y.array().colwise() * v.cwiseProduct(w).array();
		Y = a*( (A.transpose()*temp1).rowwise() - temp2.colwise().sum() );
	}
	
	HouseholderQR<MatrixXd> qr(Y);
	MatrixXd Q = qr.householderQ();
	MatrixXd B = (A*Q - w*Q.colwise().sum()).array().colwise() * (a*v).array();
	JacobiSVD<MatrixXd> svd(B.transpose(), ComputeFullU);
	MatrixXd U = svd.matrixU();
	MatrixXd U2 = MatrixXd::Zero(U.rows(),k);
	for(int i=0; i<k; i++) {
		U2.col(i) = U.col(i);
	}
	return Q*U2;
}

MatrixXd randomizedSVD_XXt_Cpp_dense(double a, VectorXd& v, MatrixXd& A, VectorXd& w, int k, int q=2) {
	// matrix with normally distributed entries
	MatrixXd Y = normalRandomMatrix(A.rows(),2*k);
	
	Y = (A*Y - w*Y.colwise().sum()).array().colwise() * (a*v).array();
	
	MatrixXd temp1;
	MatrixXd temp2;
	for(int i=0; i<q; i++) {
		temp1 = Y.array().colwise() * v.array();
		temp2 = Y.array().colwise() * v.cwiseProduct(w).array();
		Y = a*( (A.transpose()*temp1).rowwise() - temp2.colwise().sum() );
		
		Y = (A*Y - w*Y.colwise().sum()).array().colwise() * (a*v).array();
	}
	
	HouseholderQR<MatrixXd> qr(Y);
	MatrixXd Q = qr.householderQ();
	temp1 = Q.array().colwise() * v.array();
	temp2 = Q.array().colwise() * v.cwiseProduct(w).array();
	MatrixXd B = a*( (temp1.transpose()*A).colwise() - temp2.colwise().sum().transpose() );
	JacobiSVD<MatrixXd> svd(B, ComputeFullU);
	MatrixXd U = svd.matrixU();
	MatrixXd U2 = MatrixXd::Zero(U.rows(),k);
	for(int i=0; i<k; i++) {
		U2.col(i) = U.col(i);
	}
	return Q*U2;
}

// [[Rcpp::export]]
MatrixXd fastCovEVsCpp_dense(MatrixXd& X, int k, int q=2) {
	double a = 1.0/sqrt(X.rows()-1.0);
	VectorXd v = VectorXd::Zero(X.cols()).array() + 1.0;
	VectorXd w = X.colwise().mean();
	X.transposeInPlace();
	return randomizedSVD_XXt_Cpp_dense(a,v,X,w,k,q);
}

// [[Rcpp::export]]
MatrixXd fastJaccardEVsCpp_dense(MatrixXd X, int k, int q=2) {
	double a = 1.0/sqrt(2.0*X.colwise().sum().maxCoeff());
	VectorXd v = VectorXd::Zero(X.rows()).array() + 1.0;
	VectorXd w = VectorXd::Zero(X.rows());
	return randomizedSVD_XtX_Cpp_dense(a,v,X,w,k,q);
}

// [[Rcpp::export]]
MatrixXd fastSMatrixEVsCpp_dense(MatrixXd X, int k, bool Djac=false, int minVariants=0, int q=2) {
	int numAlleles = 2.0*X.cols();
	VectorXd sumVariants = X.rowwise().sum();
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>0.5*numAlleles) {
			X.row(i) = 2.0 - X.row(i).array();
		}
	}
	
	sumVariants = X.rowwise().sum();
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>=minVariants) sumVariants(i) = 1.0;
		else sumVariants(i) = 0.0;
	}
	MatrixXd Y = subsetRowsDenseMatrix(X,sumVariants);
	
	VectorXd sumFilteredVariants = Y.rowwise().sum();
	double totalPossiblePairs = numAlleles*(numAlleles-1.0)/2.0;
	VectorXd totalPairs = sumFilteredVariants.array()*(sumFilteredVariants.array()-1.0)/2.0;
	VectorXd weights = VectorXd::Zero(totalPairs.size());
	for(int i=0; i<totalPairs.size(); i++) {
		if(totalPairs(i)>0.0) weights(i) = totalPossiblePairs/totalPairs(i);
		else weights(i) = 0.0;
	}
	int s_matrix_denominator = Y.rows();
	
	if(Djac) {
		double a = 1.0/sqrt(4.0*s_matrix_denominator);
		VectorXd v = VectorXd::Zero(X.rows()).array() + 1.0;
		VectorXd w = VectorXd::Zero(X.rows());
		return randomizedSVD_XtX_Cpp_dense(a,v,X,w,k,q);
	} else {
		double a = 1.0/sqrt(4.0*s_matrix_denominator);
		VectorXd v = weights.array().cwiseSqrt();
		VectorXd w = VectorXd::Zero(X.rows());
		return randomizedSVD_XtX_Cpp_dense(a,v,X,w,k,q);
	}
}

// [[Rcpp::export]]
MatrixXd fastGrmEVsCpp_dense(MatrixXd& X, int k, bool robust=true, int q=2) {
	// compute population frequencies across rows
	VectorXd p = X.rowwise().mean()/2.0;
	VectorXd counterp = 1.0-p.array();
	VectorXd qv = 2.0*p.cwiseProduct(counterp);
	// compute grm
	if(robust) {
		double a = 1.0/sqrt(qv.sum());
		VectorXd v = VectorXd::Zero(p.size()).array() + 1.0;
		VectorXd w = 2.0*p;
		return randomizedSVD_XtX_Cpp_dense(a,v,X,w,k,q);
	}
	else {
		double a = 1.0/sqrt(double(X.rows()));
		VectorXd v = qv.array().cwiseSqrt().inverse();
		VectorXd w = 2.0*p;
		return randomizedSVD_XtX_Cpp_dense(a,v,X,w,k,q);
	}
}



// ************************************************************************************
// fast eigenvector computation for sparse matrices
// ************************************************************************************

MatrixXd randomizedSVD_XtX_Cpp_sparse(double a, VectorXd& v, SparseMatrix<double>& A, VectorXd& w, int k, int q=2) {
	// matrix with normally distributed entries
	MatrixXd Y = normalRandomMatrix(A.rows(),2*k);
	
	MatrixXd temp1 = Y.array().colwise() * v.array();
	MatrixXd temp2 = Y.array().colwise() * v.cwiseProduct(w).array();
	Y = a*( (A.transpose()*temp1).rowwise() - temp2.colwise().sum() );
	
	for(int i=0; i<q; i++) {
		Y = (A*Y - w*Y.colwise().sum()).array().colwise() * (a*v).array();
		
		temp1 = Y.array().colwise() * v.array();
		temp2 = Y.array().colwise() * v.cwiseProduct(w).array();
		Y = a*( (A.transpose()*temp1).rowwise() - temp2.colwise().sum() );
	}
	
	HouseholderQR<MatrixXd> qr(Y);
	MatrixXd Q = qr.householderQ();
	MatrixXd B = (A*Q - w*Q.colwise().sum()).array().colwise() * (a*v).array();
	JacobiSVD<MatrixXd> svd(B.transpose(), ComputeFullU);
	MatrixXd U = svd.matrixU();
	MatrixXd U2 = MatrixXd::Zero(U.rows(),k);
	for(int i=0; i<k; i++) {
		U2.col(i) = U.col(i);
	}
	return Q*U2;
}

MatrixXd randomizedSVD_XXt_Cpp_sparse(double a, VectorXd& v, SparseMatrix<double>& A, VectorXd& w, int k, int q=2) {
	// matrix with normally distributed entries
	MatrixXd Y = normalRandomMatrix(A.rows(),2*k);
	
	Y = (A*Y - w*Y.colwise().sum()).array().colwise() * (a*v).array();
	
	MatrixXd temp1;
	MatrixXd temp2;
	for(int i=0; i<q; i++) {
		temp1 = Y.array().colwise() * v.array();
		temp2 = Y.array().colwise() * v.cwiseProduct(w).array();
		Y = a*( (A.transpose()*temp1).rowwise() - temp2.colwise().sum() );
		
		Y = (A*Y - w*Y.colwise().sum()).array().colwise() * (a*v).array();
	}
	
	HouseholderQR<MatrixXd> qr(Y);
	MatrixXd Q = qr.householderQ();
	temp1 = Q.array().colwise() * v.array();
	temp2 = Q.array().colwise() * v.cwiseProduct(w).array();
	MatrixXd B = a*( (temp1.transpose()*A).colwise() - temp2.colwise().sum().transpose() );
	JacobiSVD<MatrixXd> svd(B, ComputeFullU);
	MatrixXd U = svd.matrixU();
	MatrixXd U2 = MatrixXd::Zero(U.rows(),k);
	for(int i=0; i<k; i++) {
		U2.col(i) = U.col(i);
	}
	return Q*U2;
}

// [[Rcpp::export]]
MatrixXd fastCovEVsCpp_sparse(MatrixXd& T, int nrows, int ncols, int k, int q=2) {
	SparseMatrix<double> X = triplesToSparseDoubleMatrix(T,nrows,ncols);
	double a = 1.0/sqrt(X.rows()-1.0);
	VectorXd v = VectorXd::Zero(X.cols()).array() + 1.0;
	VectorXd w = colSumsSparseDoubleMatrix(X) * 1.0/X.rows();
	X = X.transpose();
	return randomizedSVD_XXt_Cpp_sparse(a,v,X,w,k,q);
}

// [[Rcpp::export]]
MatrixXd fastJaccardEVsCpp_sparse(MatrixXd& T, int nrows, int ncols, int k, int q=2) {
	SparseMatrix<double> X = triplesToSparseDoubleMatrix(T,nrows,ncols);
	double a = 1.0/sqrt(2.0*colSumsSparseDoubleMatrix(X).maxCoeff());
	VectorXd v = VectorXd::Zero(X.rows()).array() + 1.0;
	VectorXd w = VectorXd::Zero(X.rows());
	return randomizedSVD_XtX_Cpp_sparse(a,v,X,w,k,q);
}

// [[Rcpp::export]]
MatrixXd fastSMatrixEVsCpp_sparse(MatrixXi& T, int nrows, int ncols, int k, bool Djac=false, int minVariants=0, int q=2) {
	SparseMatrix<int> X = triplesToSparseIntMatrix(T,nrows,ncols);
	
	double numAlleles = X.cols();
	VectorXd sumVariants = rowSumsSparseIntMatrix(X);
	for(int k=0; k<X.outerSize(); k++) {
		for(SparseMatrix<int>::InnerIterator it(X,k); it; ++it) {
			if(sumVariants(it.row())>0.5*numAlleles) {
				X.coeffRef(it.row(),it.col()) = 2.0 - it.value();
			}
		}
	}
	
	sumVariants = rowSumsSparseIntMatrix(X);
	for(int i=0; i<sumVariants.size(); i++) {
		if(sumVariants(i)>=minVariants) sumVariants(i) = 1.0;
		else sumVariants(i) = 0.0;
	}
	SparseMatrix<double> Y = subsetRowsSparseIntMatrix(X,sumVariants);
	
	VectorXd sumFilteredVariants = rowSumsSparseDoubleMatrix(Y);
	double totalPossiblePairs = numAlleles*(numAlleles-1.0)/2.0;
	VectorXd totalPairs = sumFilteredVariants.array()*(sumFilteredVariants.array()-1.0)/2.0;
	VectorXd weights = VectorXd::Zero(totalPairs.size());
	for(int i=0; i<totalPairs.size(); i++) {
		if(totalPairs(i)>0.0) weights(i) = totalPossiblePairs/totalPairs(i);
		else weights(i) = 0.0;
	}
	int s_matrix_denominator = Y.rows();
	
	if(Djac) {
		double a = 1.0/sqrt(4.0*s_matrix_denominator);
		VectorXd v = VectorXd::Zero(X.rows()).array() + 1.0;
		VectorXd w = VectorXd::Zero(X.rows());
		return randomizedSVD_XtX_Cpp_sparse(a,v,Y,w,k,q);
	} else {
		double a = 1.0/sqrt(4.0*s_matrix_denominator);
		VectorXd v = weights.array().cwiseSqrt();
		VectorXd w = VectorXd::Zero(X.rows());
		return randomizedSVD_XtX_Cpp_sparse(a,v,Y,w,k,q);
	}
}

// [[Rcpp::export]]
MatrixXd fastGrmEVsCpp_sparse(MatrixXd& T, int nrows, int ncols, int k, bool robust=true, int q=2) {
	SparseMatrix<double> X = triplesToSparseDoubleMatrix(T,nrows,ncols);
	// compute population frequencies across rows
	VectorXd p = 0.5 * rowSumsSparseDoubleMatrix(X) * 1.0/X.cols();
	VectorXd counterp = 1.0-p.array();
	VectorXd qv = 2.0*p.cwiseProduct(counterp);
	// compute grm
	if(robust) {
		double a = 1.0/sqrt(qv.sum());
		VectorXd v = VectorXd::Zero(p.size()).array() + 1.0;
		VectorXd w = 2.0*p;
		return randomizedSVD_XtX_Cpp_sparse(a,v,X,w,k,q);
	}
	else {
		double a = 1.0/sqrt(double(X.rows()));
		VectorXd v = qv.array().cwiseSqrt().inverse();
		VectorXd w = 2.0*p;
		return randomizedSVD_XtX_Cpp_sparse(a,v,X,w,k,q);
	}
}
