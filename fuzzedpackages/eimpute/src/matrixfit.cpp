// [[Rcpp::depends(RcppEigen)]]
#ifndef WIN_BUILD
// [[Rcpp::plugins(cpp14)]]
#else
// [[Rcpp::plugins(cpp11)]]
#endif
#include <iostream>
#include <cmath>
#include <algorithm>
#include <random>
#include <vector>
#include <RcppEigen.h>
#include <SymEigs.h>
#include <Rcpp.h>
#include "matops.h"
#include <R.h>
#include <Rinternals.h>

using namespace std;
using namespace Rcpp;
using namespace Eigen;
using namespace Spectra;
typedef Eigen::Map<const Eigen::VectorXd> MapConstVec;
#ifndef WIN_BUILD
#include <rsvd/Constants.hpp>
#include <rsvd/ErrorEstimators.hpp>
#include <rsvd/RandomizedSvd.hpp>
using namespace Rsvd;
#endif

struct value_index
{
  double value;
  int index;
};

bool smaller(const value_index &x, const value_index &y)
{
  return x.value < y.value;
}

bool bigger(const value_index &x, const value_index &y)
{
  return x.value > y.value;
}

#ifndef WIN_BUILD
Eigen::MatrixXd random_trun_svd(Eigen::MatrixXd X, int k)
{
  mt19937_64 randomEngine{};
  randomEngine.seed(1029);
  RandomizedSvd<MatrixXd, mt19937_64, SubspaceIterationConditioner::Lu> rsvd(randomEngine);
  rsvd.compute(X, k);
  return rsvd.matrixU() * rsvd.singularValues().asDiagonal() * rsvd.matrixV().adjoint();
}
#endif

// [[Rcpp::export]]
Eigen::MatrixXd trun_svd(Eigen::MatrixXd X, int k)
{
  int m = X.rows();
  int n = X.cols();
  MatrixXd Y(m, n);
  int K = k;

  Rcpp::NumericVector ctr_vec = 0;
  Rcpp::NumericVector scl_vec = 0;
  MapConstVec ctr_map(ctr_vec.begin(), m);
  MapConstVec scl_map(scl_vec.begin(), m);
  bool center = false;
  bool scale = false;
  SEXP A_mat = PROTECT(Rf_allocMatrix(REALSXP, m, n));
  double *rans = REAL(A_mat);

  for (int i = 0; i < m; i++)
  {
    for (int j = 0; j < n; j++)
      rans[i + m * j] = X(i, j);
  }
  // Operation for original matrix
  // MatProd* op_orig = get_mat_prod_op(A_mat, m, n, A_mat, 1);
  MatProd *op_orig = new MatProd_matrix(A_mat, m, n);
  // Operation for SVD
  MatProd *op;

  if (m > n)
  {
    op = new SVDTallOp(op_orig, center, scale, ctr_map, scl_map);
    SymEigsSolver<double, LARGEST_ALGE, MatProd> eig_r(op, K, 2 * K + 1 > n ? n : 2 * K + 1);
    // MatrixXd R = X.transpose() * X;
    // DenseSymMatProd<double> op_r(R);
    // SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eig_r(&op_r, K, 2 * K > m ? m : 2 * K);
    eig_r.init();
    UNPROTECT(1);
    int nconv = eig_r.compute();
    VectorXd evalues;
    if (eig_r.info() == SUCCESSFUL)
    {
      evalues = eig_r.eigenvalues();
      if (nconv < K)
      {
        Rcpp::warning("only %d singular values converged, less than K = %d", nconv, K);
        K = nconv;
      }
      VectorXd d = evalues.head(K);
      d = d.array().sqrt();

      MatrixXd v = eig_r.eigenvectors(K);
      MatrixXd u = X * v;
      MatrixXd D;
      D.setIdentity(K, K);

      for (int i = 0; i < K; i++)
      {
        u.col(i).array() /= d(i);
        D(i, i) = d(i);
      }
      Y = u * D * v.transpose();
    }
  }
  else
  {
    //MatProd* L;
    // L = new SVDWideOp(op_orig);
    op = new SVDWideOp(op_orig, center, scale, ctr_map, scl_map);
    SymEigsSolver<double, LARGEST_ALGE, MatProd> eig_l(op, K, 2 * K + 1 > n ? n : 2 * K + 1);
    // MatrixXd L = X * X.transpose();
    // DenseSymMatProd<double> op_l(L);
    // SymEigsSolver< double, LARGEST_ALGE, DenseSymMatProd<double> > eig_l(&op_l, K, 2 * K > m ? m : 2 * K);
    eig_l.init();
    UNPROTECT(1);
    int nconv = eig_l.compute();
    VectorXd evalues;
    if (eig_l.info() == SUCCESSFUL)
    {
      evalues = eig_l.eigenvalues();
      if (nconv < K)
      {
        Rcpp::warning("only %d singular values converged, less than K = %d", nconv, K);
        K = nconv;
      }

      VectorXd d = evalues.head(K);
      d = d.array().sqrt();

      MatrixXd u = eig_l.eigenvectors(K);
      MatrixXd v = X.transpose() * u;
      MatrixXd D;
      D.setIdentity(K, K);

      for (int i = 0; i < K; i++)
      {
        v.col(i).array() /= d(i);
        D(i, i) = d(i);
      }

      Y = u * D * v.transpose();
      // MatrixXd vec_l = eig_l.eigenvectors(K);
      // Y = vec_l * vec_l.transpose() * X;
    }
  }
  return Y;
}

MatrixXd DS(MatrixXd M, MatrixXd L, vector<value_index> imp, int s)
{
  int m = M.rows();
  int n = M.cols();
  MatrixXd S = MatrixXd::Zero(m, n);
  MatrixXd S_t = M - L;
  for (int k = 0; k < s; ++k)
  {
    int i = int(imp[k].index / n);
    int j = imp[k].index % n;
    S(i, j) = S_t(i, j);
  }
  return S;
}

//' @noRd
//' @param omega The matrix index of the observed value
//' @param X The obeserved value of the matrix
//' @param m, n The dimension of the matrix
//' @param rank The rank of matrix
//' @param max_it	 maximum number of iterations.
//' @param tol convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates.
//' @param type computing singular value decomposition, 1 is truncated singular value decomposition, 2 is randomized singular value decomposition
//' @description Use Rcpp to fit a low-rank matrix approximation to a matrix with two method computing singular value decomposition.
// [[Rcpp::export]]
List kkt_fix(Eigen::MatrixXi &omega, Eigen::VectorXd &X, int m, int n, int rank, int max_it, double tol, int type)
{
  // when rho = 1, it is equivalent to Hard Impute
  int l = omega.rows();
  double temp = X.mean();
  double train_error;
  Eigen::VectorXd X_train = X;
  Eigen::MatrixXd Z_old = MatrixXd::Constant(m, n, temp);
  double eps = 1;
  int count = 0;
  int r = 1;
  int c = 1;

  Eigen::MatrixXd (*svd_method)(Eigen::MatrixXd, int);
#ifndef WIN_BUILD
  if (type == 1)
  {
    svd_method = &trun_svd;
  }
  else
  {
    svd_method = &random_trun_svd;
  }
#else
  svd_method = &trun_svd;
#endif

  while (eps > tol && count < max_it)
  {
    for (int i = 0; i < l; i++)
    {
      r = omega(i, 0);
      c = omega(i, 1);
      X_train(i) = Z_old(r, c);
      Z_old(r, c) = X(i);
    }
    Z_old = svd_method(Z_old, rank);
    eps = 0;
    for (int i = 0; i < l; i++)
    {
      r = omega(i, 0);
      c = omega(i, 1);
      eps = (X_train(i) - Z_old(r, c)) * (X_train(i) - Z_old(r, c)) + eps;
    }
    eps = eps / X_train.squaredNorm();
    count++;
  }

  train_error = 0;
  for (int i = 0; i < l; i++)
  {
    r = omega(i, 0);
    c = omega(i, 1);
    train_error = (Z_old(r, c) - X(i)) * (Z_old(r, c) - X(i)) + train_error;
  }
  train_error = train_error / X.array().square().sum();
  // cout << count << ',' << eps << ',' << rank << endl;
  return Rcpp::List::create(Named("Z") = Z_old, Named("count") = count, Named("train_error") = train_error);
}

//' @noRd
//' @param omega The matrix index of the observed value
//' @param X The obeserved value of the matrix
//' @param m, n The dimension of the matrix
//' @param r_min The start rank for searching
//' @param r.max the max rank for searching.
//' @param n_fold number of folds in cross validation
//' @param max_it	maximum number of iterations.
//' @param tol convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates.
//' @param type computing singular value decomposition, 1 is truncated singular value decomposition, 2 is randomized singular value decomposition
//' @description Use Rcpp to search rank with cross validation.
// [[Rcpp::export]]
Eigen::VectorXd cv_rank(Eigen::MatrixXi &omega, Eigen::VectorXd &X, int m, int n, int r_min, int r_max, int n_fold, int max_it, double tol, int type)
{
  Rcpp::List res_list;
  int num = omega.rows();
  int r_num;
  r_num = r_max - r_min + 1;
  int num_fold = ceil((double)num / (double)n_fold);
  int test_stard, test_end;
  Eigen::MatrixXi train_set_ind;
  Eigen::VectorXd train_set;
  Eigen::MatrixXi test_set_ind;
  Eigen::VectorXd test_set;
  Eigen::MatrixXd cv_res_matrix;
  double test_error_old = 0;

  Eigen::VectorXd rank_test_error = Eigen::VectorXd::Zero(r_num);
  for (int i = 0; i < n_fold; i++) {
    test_stard = num_fold * i;
    test_end = num_fold * (i + 1);
    test_end = test_end > num ? num : test_end;
    int test_size = test_end - test_stard;
    int train_size = num - test_size;

    train_set = Eigen::VectorXd::Zero(train_size);
    train_set_ind = Eigen::MatrixXi::Zero(train_size, 2);
    test_set = Eigen::VectorXd::Zero(test_size);
    test_set_ind = Eigen::MatrixXi::Zero(test_size, 2);

    int train_temp = 0;
    int test_temp = 0;
    for (int s = 0; s < num; s++) {
      if (s >= test_stard && s < test_end) {
        test_set(test_temp) = X(s);
        test_set_ind.row(test_temp) = omega.row(s);
        test_temp++;
      } else {
        train_set(train_temp) = X(s);
        train_set_ind.row(train_temp) = omega.row(s);
        train_temp++;
      }
    }

    for (int k = 0; k < r_num; k++) {
      int temp_rank = r_min + k;
      res_list = kkt_fix(train_set_ind, train_set, m, n, temp_rank, max_it, tol, type);
      cv_res_matrix = res_list["Z"];
      test_error_old = 0;
      for (int j = 0; j < test_size; j++) {
        int row = test_set_ind(j, 0);
        int col = test_set_ind(j, 1);
        test_error_old += (cv_res_matrix(row, col) - test_set(j)) * (cv_res_matrix(row, col) - test_set(j));
      }
      test_error_old = test_error_old / test_set.squaredNorm();
      rank_test_error(k) = rank_test_error(k) + test_error_old;
    }
  }

  return rank_test_error / n_fold;
}

//' @noRd
//' @param omega The matrix index of the observed value
//' @param X The obeserved value of the matrix
//' @param m, n The dimension of the matrix
//' @param r_min The start rank for searching
//' @param r.max the max rank for searching.
//' @param max_it	 maximum number of iterations.
//' @param tol convergence threshold, measured as the relative change in the Frobenius norm between two successive estimates.
//' @param type computing singular value decomposition, 1 is truncated singular value decomposition, 2 is randomized singular value decomposition
//' @description Use Rcpp to search rank with information criterion rule.
// [[Rcpp::export]]
Eigen::VectorXd ic_rank(Eigen::MatrixXi &omega, Eigen::VectorXd &X, int m, int n, int r_min, int r_max, int max_it, double tol, int type)
{
  Rcpp::List res_list;
  int r_temp;
  int r_len;
  Eigen::VectorXd res_error;

  r_len = r_max - r_min + 1;
  res_error.setZero(r_len);

  for(int j = 0; j < r_len; j++){
    r_temp = r_min + j;
    res_list = kkt_fix(omega, X, m, n, r_temp, max_it, tol, type);
    res_error(j) = res_list["train_error"];
  }
  return res_error;
}
