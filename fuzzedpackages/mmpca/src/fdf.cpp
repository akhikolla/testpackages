#include <RcppEigen.h>
#include <gsl/gsl_sf_trig.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <iostream>
using namespace std;

using Eigen::Array;
using Eigen::Dynamic;
using Eigen::Map;
using Eigen::Matrix;
using Eigen::MatrixBase;

typedef Matrix<double, Dynamic, Dynamic> matrix;
typedef Matrix<double, Dynamic, 1> vector;
typedef Array<Array<matrix, Dynamic, Dynamic>, Dynamic, 1> pvt;

void init_parallel() {
  Eigen::initParallel();
}

// small constant used in soft abs and soft sign
const double eps = 1e-6;
const double inveps = 1/eps;

template <typename Var, typename Const>
void partialvrs(pvt * pvs, MatrixBase<Var> const & V_,
    const MatrixBase<Const> & xi) {
  MatrixBase<Var> & V = const_cast<MatrixBase<Var>&>(V_);
  int k = xi.cols();
  int p = xi.rows();
  V.setIdentity();
  (*pvs)(0) = Array<matrix, Dynamic, Dynamic>(p, k);
  (*pvs)(1) = Array<matrix, Dynamic, Dynamic>(p, k);
  for (int j = k-1; j >= 0; j--) {
    for (int i = p-1; i >= j + 1; i--) {
      (*pvs)(0)(i, j) = matrix(V.row(i));
      (*pvs)(1)(i, j) = matrix(V.row(j));
      double cx = cos(xi(i, j)) - 1;
      double sx = sin(xi(i, j));
      V.row(i) += -sx * (*pvs)(1)(i, j) + cx * (*pvs)(0)(i, j);
      V.row(j) += cx * (*pvs)(1)(i, j) + sx * (*pvs)(0)(i, j);
    }
  }
}

template <typename Mat>
void partialvls(pvt * pvs, const MatrixBase<Mat> & xi) {
  int k = xi.cols();
  int p = xi.rows();
  matrix V(p, p);
  V.setIdentity();
  (*pvs)(2) = Array<matrix, Dynamic, Dynamic>(p, k);
  (*pvs)(3) = Array<matrix, Dynamic, Dynamic>(p, k);
  for (int j = 0; j < k; j++) {
    for (int i = j + 1; i < p; i++) {
      (*pvs)(2)(i, j) = matrix(V.block(0, i, i+1, 1));
      (*pvs)(3)(i, j) = matrix(V.block(0, j, i+1, 1));
      double cx = cos(xi(i, j)) - 1;
      double sx = sin(xi(i, j));
      V.block(0, i, i+1, 1) += cx * (*pvs)(2)(i, j) + sx * (*pvs)(3)(i, j);
      V.block(0, j, i+1, 1) += -sx * (*pvs)(2)(i, j) + cx * (*pvs)(3)(i, j);
    }
  }
}

template <typename Var, typename Const>
void Vxi(MatrixBase<Var> const & V_, const MatrixBase<Const> & xi) {
  MatrixBase<Var> & V = const_cast<MatrixBase<Var>&>(V_);
  V.setIdentity();
  int p = V.rows();
  for (int j = xi.cols()-1; j >= 0; j--) {
    for (int i = p-1; i >= j+1; i--) {
      matrix vi(V.row(i));
      double cx = cos(xi(i, j)) - 1;
      double sx = sin(xi(i, j));
      V.row(i) += -sx * V.row(j) + cx * vi;
      V.row(j) += cx * V.row(j) + sx * vi;
    }
  }
}

void f_vxi(double * v, const double * xi, int p, int k) {
  Map<matrix> eig_v(v, p, k);
  Map<const matrix> eig_x(xi, p, k);
  Vxi(eig_v, eig_x);
}

double f_obj(const double * theta, const double ** x, const double ** masks,
    const double * lambda, const int k, const int * inds, const int * p,
    const int m, const int n, const int len) {//, const int num_threads) {
  //omp_set_num_threads(num_threads);
  double loss = 0;
  int ix = 0;
  int * ixs = new int[n + 1];
  for (int i = 0; i < n; i++) {
    ixs[i] = ix;
    ix += k * p[i];
  }
  ixs[n] = ix;
  Map<const matrix> D(theta + ixs[n], k, n);
  Array<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> V(D.cols());
  //# pragma omp parallel for schedule(dynamic) reduction(+:loss)
  // undeterministic results when parallelizing this loop. dont know why.
  for (int i = 0; i < V.rows(); i++) {
    Map<const matrix> xi(theta + ixs[i], p[i], k);
    V(i) = Matrix<double, Dynamic, Dynamic>(p[i], k);
    Vxi(V(i), xi);
    matrix vd = V(i) * D.col(i).asDiagonal();
    loss += lambda[2] * vd.array().abs().sum();
    for (int j = 0; j < V(i).rows(); j++) {
      loss += lambda[3] * sqrt(vd.row(j).array().square().sum());
    }
  }
  //# pragma omp parallel for schedule(dynamic) reduction(+:loss)
  // undeterministic results when parallelizing this loop. dont know why.
  for (int i = 0; i < m; i++) {
    int row = inds[i];
    int col = inds[m + i];
    Map<const matrix> X(x[i], p[row], p[col]);
    Map<const matrix> mask(masks[i], p[row], p[col]);
    loss += (X - V(row) *
      (D.col(row).cwiseProduct(D.col(col))).asDiagonal() *
      V(col).transpose()).cwiseProduct(mask).array().square().sum();
  }
  loss += lambda[0] * (eps * (inveps * D.array())
    .unaryExpr([](double x) { return gsl_sf_lncosh(x); })).abs().sum();
  for (int i = 0; i < D.rows(); i++) {
    loss += lambda[1] * sqrt(D.row(i).array().square().sum());
  }
  return loss;
}

void d_obj(double * grad, const double * theta, const double ** x,
    const double ** masks, const double * lambda, const int k, const int * inds,
    const int * p, const int m, const int n, const int len,
    const double * indices, const int indices_len, const int num_threads) {
#ifdef _OPENMP
  omp_set_num_threads(num_threads);
#endif
  memset(grad, 0.0, sizeof(double) * len);
  int ix = 0;
  int * ixs = new int[n + 1];
  Array<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> V(n);
  Array<pvt, Dynamic, 1> pvss(n);
  for (int i = 0; i < n; i++) {
    ixs[i] = ix;
    ix += k * p[i];
    pvss(i) = pvt(4);
  }
  // precalculate rotations and partial rotations
  # pragma omp parallel for schedule(dynamic)
  for (int j = 0; j < 2*n; j++) {
    int i = j / 2;
    Map<const matrix> xi(theta + ixs[i], p[i], k);
    if (j % 2 == 0) {
      V(i) = Matrix<double, Dynamic, Dynamic>(p[i], k);
      partialvrs(&pvss(i), V(i), xi);
    } else {
      partialvls(&pvss(i), xi);
    }
  }
  ixs[n] = ix;
  Map<const matrix> D(theta + ixs[n], k, n);
  Map<matrix> dD(grad + ixs[n], k, n);
  dD.array() += lambda[0] * (inveps * D).array().tanh();
  // integration penalty gradient of D
  for (int i = 0; i < k; i++) {
    double a = sqrt(D.row(i).array().square().sum());
    if (a >= 1e-8) {
      dD.row(i).array() += lambda[1] * D.row(i).array() / a;
    }
  }
  Array<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> T1(n);
  Array<Matrix<double, Dynamic, Dynamic>, Dynamic, 1> T2(n);
  # pragma omp parallel for schedule(dynamic)
  for (int view = 0; view < n; view++) {
    matrix L1(k, p[view]);
    L1.setZero();
    // loss gradient of xi and loss gradient of D
    for (int block = 0; block < m; block++) {
      if (inds[block] == view) {
        int other = inds[block + m];
        Map<const matrix> X(x[block], p[view], p[other]);
        Map<const matrix> mask(masks[block], p[view], p[other]);
        matrix L2(V(other) * (D.col(other).cwiseProduct(D.col(view)))
          .asDiagonal());
        L1 += ((X - V(view) * L2.transpose()).cwiseProduct(mask) * L2)
          .transpose();
        dD.col(view) += (-2*V(view).transpose() * (X - V(view) *
          L2.transpose()).cwiseProduct(mask) * V(other) *
          D.col(other).asDiagonal()).diagonal();
      } else if (inds[block + m] == view) {
        int other = inds[block];
        Map<const matrix> X(x[block], p[other], p[view]);
        Map<const matrix> mask(masks[block], p[other], p[view]);
        matrix L2(V(other) * (D.col(other).cwiseProduct(D.col(view)))
          .asDiagonal());
        L1 += ((X.transpose() - V(view) * L2.transpose())
          .cwiseProduct(mask.transpose()) * L2).transpose();
        dD.col(view) += (-2*V(view).transpose() * (X.transpose() - V(view) *
          L2.transpose()).cwiseProduct(mask.transpose()) * V(other) *
          D.col(other).asDiagonal()).diagonal();
      }
    }
    matrix VD(V(view)*D.col(view).asDiagonal());
    matrix VD2((V(view)*(D.col(view).array().square().matrix()).asDiagonal())
      .transpose());
    if (lambda[3] > 0) {
      Array<double, 1, Dynamic> denoms(
        sqrt(VD.array().square().rowwise().sum()));
      for (int r = 0; r < p[view]; r++) {
        if (denoms(r) > 1e-8) {
          denoms(r) = 1/denoms(r);
        } else {
          denoms(r) = 0.0;
        }
      }
      VD2.array().rowwise() *= denoms;
      dD.col(view) += lambda[3] * ((V(view).array().colwise() *
        denoms.transpose()) *
        VD.array()).matrix().colwise().sum().transpose();
    }
    T1(view) = matrix((inveps * VD).array().tanh().matrix());
    // its much cheaper to make the following sum and then get just one big
    // matrix multiplication per view instead of one per term in objective
    T2(view) = matrix(lambda[2] * (T1(view) * D.col(view).asDiagonal())
      .transpose() + lambda[3] * VD2 - 2 * L1);
  }
  for (int view = 0; view < n; view++) {
    for (int i = 0; i < k; i++) {
      dD(i, view) += lambda[2] *
        (V(view).col(i).transpose() * T1(view).col(i))(0, 0);
    }
  }
  # pragma omp parallel for
  for (int c = 0; c < indices_len; c += 3) {
    int view = indices[c];
    int i = indices[c+1];
    int j = indices[c+2];
    Map<const matrix> xi(theta + ixs[view], p[view], k);
    Map<matrix> dxi(grad + ixs[view], p[view], k);
    double cx = cos(xi(i, j));
    double msx = -sin(xi(i, j));
    Matrix<double, Dynamic, 1> ti(T2(view).block(0, 0, k, i+1)
      * pvss(view)(2)(i, j));
    Matrix<double, Dynamic, 1> tj(T2(view).block(0, 0, k, i+1)
      * pvss(view)(3)(i, j));
    dxi(i, j) = (msx * (pvss(view)(0)(i, j) * ti) -
      cx * (pvss(view)(1)(i, j) * ti) +
      cx * (pvss(view)(0)(i, j) * tj) +
      msx * (pvss(view)(1)(i, j) * tj)).sum();
  }
}

void inv_v(double * xi, double * t, int n) {
  memset(xi, 0.0, n * n * sizeof(double));
  double * tnew = (double *) malloc(n * n * sizeof(double));
  double * f = (double *) malloc(n * n * sizeof(double));
  double * identity = (double *) calloc(n * n, sizeof(double));
  for (int i = 0; i < n; i++) {
    identity[i*n+i] = 1.0;
  }

  for (int v = n - 1; v >= 1; v--) {
    xi[v*n+(v-1)] = atan2(t[(v-1)*n+v], t[v*n+v]);
    if (v > 1) {
      for (int k = v - 2; k >= 0; k--) {
        double s = sin(xi[v*n+(k+1)]);
        if (abs(t[(k+1)*n+v] - s) < 1e-17) {
          if (v == n - 1 && abs(s) < 1e-17) {
            xi[v*n+k] = asin(t[k*n+v]);
          } else {
            xi[v*n+k] = 0.0;
          }
        } else {
          xi[v*n+k] = atan2(t[k*n+v], t[(k+1)*n+v]/s);
        }
        if (xi[v*n+k] > M_PI/2) xi[v*n+k] = xi[v*n+k] - M_PI;
        if (xi[v*n+k] < -M_PI/2) xi[v*n+k] = xi[v*n+k] + M_PI;
      }
      memcpy(f, t, sizeof(double) * n * n);
      memcpy(tnew, identity, sizeof(double) * n * n);
      for (int k = v - 1; k >= 0; k--) {
        double c = cos(xi[v*n+k]);
        double s = sin(xi[v*n+k]);
        for (int j = 0; j < n; j++) {
          tnew[k*n+j] = c * t[k*n+j] - s * f[(k+1)*n+j];
          f[k*n+j] = s * t[k*n+j] + c * f[(k+1)*n+j];
        }
      }
      memcpy(t, tnew, sizeof(double) * n * n);
    }
    R_CheckUserInterrupt();
  }

  free(identity);
  free(f);
  free(tnew);
}
