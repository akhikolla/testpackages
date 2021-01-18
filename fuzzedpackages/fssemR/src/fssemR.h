#ifndef FSSEMR_H_
#define FSSEMR_H_

//#include "Spectra/SymEigsSolver.h"
#include "LR.h"
#ifdef _OPENMP
#include <omp.h>
#endif
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

#define epsilon 1e-12

inline double inert(double k)
{
  return (k - 1) / (k + 2);
}

inline double sign(double x)
{
  return (x > 0) ? 1 : ((x < 0) ? -1 : 0);
}

template<class MatType>
inline MatType ProjectMat(MatType& Xi)
{
  Eigen::FullPivLU<MatType> XX(Xi * Xi.transpose());
  return XX.solve(Xi);
}

template<class MatType>
inline void Lipschitz_Init(MatType& ImB, MatType& o, double& det, int i)
{
  MatType ImBi = rm_Col(ImB, i);
  Eigen::FullPivLU<MatType> t(ImBi.transpose() * ImBi);
  o   = t.inverse();
  det = t.determinant();
}

template<class MatType>
inline double Lipschitz(MatType& o, MatType& g, MatType& s, double c2, double det, double Y2norm, double sigma2, int n, int p)
{
  MatType og   = o * g.transpose();
  MatType io   = MatType::Identity(p-1, p-1) - g * og;
  MatType sg   = s * og;
  float lambda = epsilon;
  MatType c    = (1.0 - (s * o * s.transpose()).array()).matrix();
  Eigen::FullPivLU<MatType> H(p-1, p-1);
  H.compute(io + MatType::Identity(p-1, p-1) * lambda);
  MatType x    = -1.0 * H.solve(sg.transpose());
  MatType L    = n * c2 / (x.transpose() * io * x + 2 * sg * x + c).array() / (det + lambda) + Y2norm / sigma2;
  while(L(0, 0) < 0) {
    lambda *= 10;
    H.compute(io + MatType::Identity(p-1, p-1) * lambda);
    x  = -1.0 * H.solve(sg.transpose());
    L  = n * c2 / (x.transpose() * io * x + 2 * sg * x + c).array() / (det + lambda) + Y2norm / sigma2;
  }
  return L(0, 0);
}


// proximal operator
template<class MatType>
inline MatType prox_lasso(MatType& u, MatType& w, double c, double lambda)
{
  return ((u.array() - lambda * w.array() / c).cwiseMax(0) + (u.array() + lambda * w.array() / c).cwiseMin(0)).matrix();
}


// sigma2 update
template<class MatType>
inline double get_sigma2(MatType& X, MatType& Y, MatType& B, MatType& F) {
  int n = Y.cols();
  int p = Y.rows();
  MatType err = (MatType::Identity(p, p) - B) * Y - F.transpose() * X;
  return err.squaredNorm() / (n * p - 1);
}

// error
template<class MatType>
inline double get_err(MatType& Xn, MatType& Xp) {
  return (Xn - Xp).norm() / (1 + Xp.norm());
}

// object
template<class MatType>
inline double get_obj(MatType& B, MatType& W, double sigma2, double lambda, double detIB, int n, int p)
{
  double  L = -n / 2 * std::log(detIB * detIB) + n * p / 2 * std::log(sigma2);
  MatType R = lambda * (W.array() * B.array()).matrix();
  R.diagonal().array() = 0.0;
  return L + R.cwiseAbs().sum();
}


// loglikelihoodFSSEM
template<class MatType>
inline double loglikelihoodFSSEM(std::vector<MatType>& Bs, std::vector<MatType>& Wl, MatType& Wf,
                                 std::vector<double>& DetIBs, double lambda, double rho,
                                 const std::vector<int>& n, const int p, const int m, double sigma2)
{
  double loglik = 0;
  MatType l1(p, p);
  l1.setZero();
  MatType lf = rho * (Wf.array() * (Bs[1] - Bs[0]).cwiseAbs().array()).matrix();
  for(int i = 0; i < m; i++)
  {
    loglik = loglik - n[i] / 2 * log(DetIBs[i] * DetIBs[i]) + n[i] * p / 2 * std::log(sigma2);
    l1 = l1 + lambda * (Wl[i].array() * Bs[i].cwiseAbs().array()).matrix();
  }
  lf.diagonal().array() = 0;
  l1.diagonal().array() = 0;
  loglik = loglik + l1.sum() + lf.sum();
  return loglik;
}

template<class MatType>
inline void prox_flsa(std::vector<MatType>& u, std::vector<MatType>& x, double lambda, double rho, double c, std::vector<MatType>& w, MatType& r)
{
  MatType du = u[0] - u[1];
  int p = du.rows();
  for(int i = 0; i < 2; i++)
  {
    x[i].resize(p-1, 1);
  }
  for(int i = 0; i < p; i++)
  {
    if(std::abs(du(i, 0)) <= (2 * rho * r(i, 0) / c))
    {
      x[0](i, 0) = (u[0](i, 0) + u[1](i, 0)) / 2;
      x[1](i, 0) = (u[0](i, 0) + u[1](i, 0)) / 2;
      x[0](i, 0) = std::fmax(x[0](i, 0) - lambda * (w[0](i, 0) + w[1](i, 0)) / (2 * c), 0) +
        std::fmin(x[0](i, 0) + lambda * (w[0](i, 0) + w[1](i, 0)) / (2 * c), 0);
      x[1](i, 0) = std::fmax(x[1](i, 0) - lambda * (w[0](i, 0) + w[1](i, 0)) / (2 * c), 0) +
        std::fmin(x[1](i, 0) + lambda * (w[0](i, 0) + w[1](i, 0)) / (2 * c), 0);
    }
    else
    {
      x[0](i, 0) = u[0](i, 0) - sign(du(i, 0)) * rho * r(i, 0) / c;
      x[1](i, 0) = u[1](i, 0) + sign(du(i, 0)) * rho * r(i, 0) / c;
      x[0](i, 0) = std::fmax(x[0](i, 0) - lambda * w[0](i, 0) / c, 0) +
        std::fmin(x[0](i, 0) + lambda * w[0](i, 0) / c, 0);
      x[1](i, 0) = std::fmax(x[1](i, 0) - lambda * w[1](i, 0) / c, 0) +
        std::fmin(x[1](i, 0) + lambda * w[1](i, 0) / c, 0);
    }
  }
}

// sigma2 update
template<class MatType>
inline double get_sigma2(std::vector<MatType>& Xs, std::vector<MatType>& Ys, std::vector<MatType>& Bs, std::vector<MatType>& Fs, const int m) {
  int df = 0;
  double err2 = 0.0;
  for(int i = 0; i < m; i++)
  {
    err2 += (Ys[i] - Bs[i] * Ys[i] - Fs[i].transpose() * Xs[i]).squaredNorm();
    df += (Ys[i].cols() * Ys[i].rows());
  }
  return err2 / df;
}


#endif //FSSEMR_H_
