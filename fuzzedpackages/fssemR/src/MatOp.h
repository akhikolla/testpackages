#ifndef MATOP_H_
#define MATOP_H_

#include <Eigen/Core>
#include <RcppEigen.h>

using Eigen::Map;
using Eigen::MatrixXf;
using Eigen::VectorXf;
using Eigen::ArrayXd;
using Eigen::ArrayXf;
using Eigen::Matrix;
using std::max;

// operator on matrix
template<class MatType, class RowIndexType>
inline MatType get_Rows(const MatType& X, const RowIndexType& ix) {
  MatType Rows(ix.size(), X.cols());
  for(unsigned int i = 0; i < ix.size(); ++i) {
    Rows.row(i) = X.row(ix(i) - 1);
  }
  return Rows;
}

template<class MatType, class ColIndexType>
inline MatType get_Cols(const MatType& X, const ColIndexType& ix) {
  MatType Cols(X.rows(), ix.size());
  for(unsigned int i = 0; i < ix.size(); ++i) {
    Cols.col(i) = X.col(ix(i) - 1);
  }
  return Cols;
}

template<class MatType>
inline MatType rm_Row(const MatType& X, const int ix) {
  MatType Mat(X.rows()-1, X.cols());
  int i = 0;
  int j = 0;
  while(j < X.rows()) {
    if(j != ix - 1) {
      Mat.row(i) = X.row(j);
      i++;
      j++;
    } else {
      j++;
    }
  }
  return Mat;
}

template<class MatType>
inline MatType rm_Col(const MatType& X, const int ix) {
  MatType Mat(X.rows(), X.cols()-1);
  int i = 0;
  int j = 0;
  while(j < X.cols()) {
    if (j != ix - 1) {
      Mat.col(i) = X.col(j);
      i++;
      j++;
    } else {
      j++;
    }
  }
  return Mat;
}

template<class VecType, class IndexType>
inline VecType get_Vec(const VecType& x, const IndexType& ix) {
  VecType vec(ix.size());
  for(unsigned int i = 0; i < ix.size(); ++i) {
    vec(i) = x(ix(i) - 1);
  }
  return vec;
}

template<class MatType>
inline void set_Row(MatType& X, const int ix, const MatType& src) {
  int i = 0;
  int j = 0;
  while(i < X.cols()) {
    if(i != ix - 1) {
      X(ix - 1, i) = src(j, 0);
      i++;
      j++;
    } else {
      X(ix - 1, i) = 0;
      i++;
    }
  }
}

template<class MatType>
inline void set_Col(MatType& X, const int ix, const MatType& src) {
  int i = 0;
  int j = 0;
  while(i < X.rows()) {
    if(i != ix - 1) {
      X(i, ix - 1) = src(j, 0);
      i++;
      j++;
    } else {
      X(i, ix - 1) = 0;
      i++;
    }
  }
}

template<class MatType>
inline void center(MatType& X, MatType& meanX, bool trans = false) {
  int p = X.cols();
  int n = X.rows();
  meanX.resize(p, 1);
  float *begin, *end;
#ifdef _OPENMP
  #pragma omp parallel for private(begin, end) shared(meanX)
#endif
  for(int i = 0; i < p; i++) {
    begin = &X(0, i);
    end   = begin + n;
    meanX(i, 0) = X.col(i).mean();
    std::transform(begin, end, begin, std::bind(std::minus<float>(), std::placeholders::_1, meanX(i, 0)));
  }
  if(trans) {
    X.transposeInPlace();
    meanX.transposeInPlace();
  }
}

template<class MatType, class RowIndexType>
inline void set_Row(MatType& F, const MatType& fi, const RowIndexType& ix, int i) {
  int s = fi.rows();
  for(int j = 0; j < s; j++) {
    F(ix(j) - 1, i) = fi(j, 0);
  }
}


template<class MatType, class ColIndexType>
inline void set_Col(MatType& F, const MatType& fi, const ColIndexType& ix, int i) {
  int s = fi.rows();
  for(int j = 0; j < s; j++) {
    F(i, ix(j) - 1) = fi(j, 0);
  }
}

template<class MatType, class VecType>
inline MatType get_Fs(std::vector<MatType>& fi, VecType S, int k) {
  int p = fi.size();
  MatType F(p, k);
  F.setZero();
#ifdef _OPENMP
  #pragma omp parallel for shared(F)
#endif
  for(int i = 0; i < p; i++) {
    set_Col(F, fi[i], S[i], i);
  }
  return F;
}

template<class MatType>
inline MatType inverse(MatType& X, double& determinant) {
  Eigen::FullPivLU<MatType> H(X);
  determinant = H.determinant();
  return H.inverse();
}

#endif // MATOP_H_
