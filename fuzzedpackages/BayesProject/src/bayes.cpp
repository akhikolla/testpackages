#include <R.h>
#include <Rcpp.h>
using namespace Rcpp;
#include <RcppEigen.h>
using namespace Eigen;

/*
#include <Eigen/SVD>
#include <Eigen/Core>
using Eigen::Map;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

static Rcpp::Environment base("package:base");
static Rcpp::Function Rsample = base["sample"];
*/


// *********************************************************************************
// R functions

// row-wise cumsum
MatrixXd rowCumsum(MatrixXd X) {
  MatrixXd Y = X;
  for(int i=1; i<Y.cols(); i++) {
    Y.col(i) += Y.col(i-1);
  }
  return Y;
}

VectorXd seqXd(int from, int to) {
  VectorXd v(to-from+1);
  for(int i=0; i<v.size(); i++) {
    v[i] = from+i;
  }
  return v;
}

// R command colSums for matrix
VectorXd colSums(MatrixXd X) {
  VectorXd v(X.cols());
  for(int i=0; i<X.cols(); i++) {
    v[i] = X.col(i).sum();
  }
  return v;
}

// return subset of columns, v must be an integer vector (numbered from 0)
MatrixXd subsetCols(MatrixXd X, VectorXd v) {
  MatrixXd R = MatrixXd::Zero(X.rows(),v.size());
  for(int i=0; i<v.size(); i++) {
    R.col(i) = X.col(int(v[i]));
  }
  return R;
}

// concatenate two vectors
VectorXd concatenate(VectorXd x, VectorXd y) {
  VectorXd z(x.size()+y.size());
  for(int i=0; i<x.size(); i++) {
    z[i] = x[i];
  }
  for(int i=0; i<y.size(); i++) {
    z[x.size()+i] = y[i];
  }
  return z;
}

// argmax of vector, returns first occurrence
int argmax(VectorXd v) {
  double m = v.maxCoeff();
  for(int i=0; i<v.size(); i++) {
    if(!(v[i]<m)) return(i);
  }
  return(0);
}

// compute difference of row-wise means for 0:t and (t+1):n, where n is the number of columns of matrix X
VectorXd rowMeanDiff(MatrixXd X, int t) {
  // compute means for left part only and substract from overall means, since this saves one copy operation
  MatrixXd leftX(X);
  leftX.conservativeResize(X.rows(),t+1);
  return leftX.rowwise().mean().array() - (X.rowwise().sum().array() - leftX.rowwise().sum().array())/(X.cols()-leftX.cols());
}

// cusum for vector for 0:t and (t+1):n, where n=length(v), thus 0 <= t <= n-2
double vectorCusumPoint(VectorXd v, int t) {
  int n = v.size();
  double leftsum = 0;
  for(int i=0; i<=t; i++) {
    leftsum += v[i];
  }
  return (leftsum/(t+1) - (v.array().sum() - leftsum)/(n-(t+1)))/sqrt(1.0/(t+1)+1.0/(n-(t+1)));
}

// same as rowMeanDiff but now compute row-wise means with different timepoints for each row given in vector "t"
VectorXd matrixCusumPoints(MatrixXd X, VectorXd t) {
  VectorXd ret(X.rows());
  for(int i=0; i<X.rows(); i++) {
    ret[i] = vectorCusumPoint(X.row(i), int(t[i]));
  }
  return ret;
}

// Rcpp version of the cusum transform of the Inspect package
MatrixXd cusum_transform(MatrixXd x) {
  int p = x.rows();
  int n = x.cols();
  MatrixXd leftsums = rowCumsum(x);
  MatrixXd rightsums = (-leftsums).array().colwise() + leftsums.array().col(leftsums.cols()-1);
  leftsums.transposeInPlace();
  rightsums.transposeInPlace();
  
  // vector t=1:(n-1)
  VectorXd t = seqXd(1,n-1);
  // vector tt=t*(n-t)/n
  VectorXd tt = t.array()*(n-t.array())/n;
  // delete last row
  leftsums.conservativeResize(n-1,p);
  rightsums.conservativeResize(n-1,p);
  
  rightsums = rightsums.array().colwise() * (n-t.array()).inverse();
  leftsums = leftsums.array().colwise() * t.array().inverse();
  return ((rightsums - leftsums).array().colwise() * tt.array().sqrt()).transpose();
}
// *********************************************************************************



// *********************************************************************************
// projection direction
// *********************************************************************************

// Bayes projection directions: returns matrix of posterior cusum vectors as columns (dimensions p * (n-1) if input is p*n matrix)
// projections are only computed at a subset of all time points given in timePoints=1:(n-1)
// [[Rcpp::export]]
MatrixXd bayes_vhat(MatrixXd x, VectorXd timePoints, double K=1.0/sqrt(2.0)) {
  int p = x.rows();
  int n = x.cols();
  
  // basic cusum as in cusum_transform without scaling
  MatrixXd cusumLeft = rowCumsum(x);
  MatrixXd cusumRight = (-cusumLeft).array().colwise() + cusumLeft.array().col(cusumLeft.cols()-1);
  
  // only consider subset of indices "timePoints" to save time
  VectorXd indL = timePoints.array().inverse();
  VectorXd indR = (n-timePoints.array()).inverse();
  cusumLeft = subsetCols(cusumLeft,timePoints.array()-1);
  cusumRight = subsetCols(cusumRight,timePoints.array()-1);
  
  MatrixXd Diff = cusumLeft.array().rowwise() * indL.array().transpose() - cusumRight.array().rowwise() * indR.array().transpose();
  VectorXd sigma2 = indL + indR;
  MatrixXd denominator = K + (Diff.array().pow(2).rowwise() * (-2*sigma2.array()).transpose().inverse()).exp();
  MatrixXd vhat = Diff.cwiseProduct(denominator.array().inverse().matrix());
  MatrixXd vnorm = vhat.array().rowwise() * colSums(vhat.array().pow(2)).array().transpose().sqrt().inverse();
  return vnorm;
}



// *********************************************************************************
// changepoint functions
// *********************************************************************************

// project at all time points, compute Bayes projection direction, project and return projection & cusum vector
// the returned projection direction is the one for the largest entry in the cusum vector
// the first p (=number of rows of x) entries of the returned vector is the projection direction, the rest is the cusum vector
// [[Rcpp::export]]
VectorXd bayes_cpt(MatrixXd x, VectorXd timePoints, double K=1.0/sqrt(2.0)) {
  MatrixXd vhats = bayes_vhat(x,timePoints,K);
  MatrixXd proj = vhats.transpose() * x;
  // return cusum vector and projection vector (column) for argmax(cusum) as *one* vector
  VectorXd cusum = matrixCusumPoints(proj,timePoints.array()-1);
  int amax = argmax(cusum);
  return concatenate(vhats.col(amax), cusum);
}

// sum-cusum and max-cusum
// the first p (=number of rows of x) entries of the returned vector is the projection direction, the rest is the cusum vector
// [[Rcpp::export]]
VectorXd sum_max_cusum(MatrixXd x, bool sum_cusum=true) {
  MatrixXd cusum_matrix = cusum_transform(x).array().abs();
  // compute column-wise sums (sum-cusum) or max (max-cusum)
  VectorXd cusum;
  if(sum_cusum) {
    cusum = cusum_matrix.colwise().sum();
  }
  else {
    cusum = cusum_matrix.colwise().maxCoeff();
  }
  // compute argmax of cusum vector and the projection direction as the change in means
  int amax = argmax(cusum);
  VectorXd vhat = rowMeanDiff(x,amax);
  return concatenate(vhat, cusum);
}
