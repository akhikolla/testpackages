#ifndef LR_H_
#define LR_H_

#include "MatOp.h"

template<class MatType>
inline MatType LR_Shrink(MatType& Xi, MatType& Yi, MatType& yi, double gamma, const int n, const int p, const int k) {
  Eigen::FullPivLU<MatType> CH(Xi.transpose() * Xi);
  MatType Pi  = MatType::Identity(n, n) - Xi * CH.solve(Xi.transpose());
  MatType I = Yi.transpose() * Pi * Yi;
  I.diagonal().array() += gamma;
  Eigen::FullPivLU<MatType> LT(I);
  MatType bi = LT.solve(Yi.transpose() * Pi * yi);
  return bi;
}

template<class MatType>
inline double LR_Trace(MatType& Xi, MatType& Yi, MatType& yi, const int n) {
  Eigen::FullPivLU<MatType> CH(Xi.transpose() * Xi);
  MatType Pi  = MatType::Identity(n, n) - Xi * CH.solve(Xi.transpose());
  MatType I = Yi.transpose() * Pi * Yi;
  double t = I.diagonal().sum();
  return t;
}


template<class MatType>
inline MatType LR_Fi(MatType& Xi, MatType& Yi, MatType& yi, MatType& bi) {
  Eigen::FullPivLU<MatType> XX(Xi.transpose() * Xi);
  MatType fi = XX.solve(Xi.transpose() * (yi - Yi * bi));
  return fi;
}

template<class MatType>
inline double LR_Objerr(MatType& X, MatType& Y, MatType& B, MatType& F, MatType& mu) {
  int p = B.cols();
  MatType err = Y * (MatType::Identity(p, p) - B).transpose() - X * F;
  for(int i = 0; i < p; i++) {
    err.col(i).array() -= mu(i,0);
  }
  return err.squaredNorm();
}

template<class MatType, class VecType, class IdxType>
inline double L2lr(MatType& X, MatType& Y, IdxType& S, MatType& B, VecType& fi, MatType& mu, double gamma, const int n, const int p, const int k) {
  MatrixXf xm, ym;
  center(X, xm);
  center(Y, ym);
  MatrixXf Xi, Yi, yi, bi, xi;
  fi.resize(p);
  double error = 0, sigma2;
#ifdef _OPENMP
  #pragma omp parallel for private(Xi, Yi, yi, bi) shared(B, fi) reduction(+:error)
#endif
  for(int i = 1; i <= p; i++)
  {
    Xi = get_Cols(X, S[i-1]);
    Yi = rm_Col(Y, i);
    yi = Y.col(i-1);
    bi = LR_Shrink(Xi, Yi, yi, gamma, n, p, k);
    set_Row(B, i, bi);
    fi[i-1] = LR_Fi(Xi, Yi, yi, bi);
    error += (yi - Yi * bi - Xi * fi[i-1]).squaredNorm();
  }
  mu  = (MatrixXf::Identity(p, p) - B) * ym;
#ifdef _OPENMP
  #pragma omp parallel for private(xi) shared(mu)
#endif
  for (int i = 1; i <= p; i++) {
    xi = get_Rows(xm, S[i-1]);
    mu(i-1, 0) -= (xi.transpose() * fi[i-1])(0,0);
  }
  return error;
}

template<class MatType, class IdxType>
inline double L2lamax(MatType& X, MatType& Y, IdxType& S, const int n, const int p, const int k) {
  MatrixXf xm, ym;
  center(X, xm);
  center(Y, ym);
  MatrixXf Xi, Yi, yi;
  double lambda = 0;
#ifdef _OPENMP
  #pragma omp parallel for private(Xi, Yi, yi) reduction(max:lambda)
#endif
  for(int i = 1; i <= p; i++)
  {
    Xi = get_Cols(X, S[i-1]);
    Yi = rm_Col(Y, i);
    yi = Y.col(i-1);
    lambda = max(lambda, LR_Trace(Xi, Yi, yi, n));
  }
  return lambda;
}


#endif // LR_H_
