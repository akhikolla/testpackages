/*
 *     file fRLR/src/frlr.cpp  Fit Repeated Linear Regressions
 *     Copyright (C) 2017 Lijun Wang
 *
 *     This program is free software: you can redistribute it and/or modify
 *     it under the terms of the GNU General Public License as published by
 *     the Free Software Foundation, either version 3 of the License, or
 *     (at your option) any later version.
 *
 *     This program is distributed in the hope that it will be useful,
 *     but WITHOUT ANY WARRANTY; without even the implied warranty of
 *     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *     GNU General Public License for more details.
 *
 *     You should have received a copy of the GNU General Public License
 *     along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
  */

#include<Rcpp.h>
#include<iostream>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<vector>
#include<string>
#include<stdlib.h>
#include<stdio.h>

#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_combination.h>
#include<gsl/gsl_statistics.h>
#include<gsl/gsl_fit.h>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_multifit.h>
#include<gsl/gsl_cdf.h>
#include<gsl/gsl_permutation.h>
#include<gsl/gsl_blas.h>
#ifdef _OPENMP
#include<omp.h>
#endif
using namespace std;

int get_col_from_r_matrix(Rcpp::NumericVector m, size_t nrow, size_t ncol, size_t i, gsl_vector* v)
{
  if (ncol <= i)
  {
    return 0;
  }
  for(int j = 0; j < nrow; j++)
  {
    gsl_vector_set(v, j, m[i*nrow+j]);
  }
  return 1;
}

//[[Rcpp::export]]
Rcpp::List frlr1(SEXP R_X, SEXP R_Y, SEXP R_COV)
{
  // convert data type
  Rcpp::NumericVector X(R_X);
  Rcpp::NumericVector Y(R_Y);
  Rcpp::NumericVector COV(R_COV);

  int nrow = Y.size();
  int nX = X.size();
  int ncol = nX/nrow;

  int nCOV = COV.size();
  int COV_COL = nCOV/nrow;

  int col_b = COV_COL+1;
  gsl_matrix *b = gsl_matrix_alloc(nrow, col_b);
  gsl_matrix *B = gsl_matrix_alloc(COV_COL+1, col_b);
  gsl_matrix *invB = gsl_matrix_alloc(COV_COL+1, col_b);

  gsl_vector *covariate_vec = gsl_vector_alloc(nrow);
  for (int ip = 0; ip < COV_COL; ip++)
  {
    get_col_from_r_matrix(COV, nrow, COV_COL, ip, covariate_vec);
    gsl_matrix_set_col(b, 1+ip, covariate_vec);
  }
  gsl_vector_free(covariate_vec);

  gsl_vector *Yv = gsl_vector_alloc(nrow);
  get_col_from_r_matrix(Y, nrow, 1, 0, Yv);

  // intercept
  gsl_vector *x0 = gsl_vector_alloc(nrow);
  gsl_vector_set_all(x0, 1);
  gsl_matrix_set_col(b, 0, x0);

  // degree of freedom
  int df = nrow - COV_COL - 1 - 1;

  vector<int> r;
  vector<double> r_p;

  gsl_permutation *permutation_B = gsl_permutation_alloc(B->size1);
  int status;

  // B = b'b
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, b, b, 0.0, B);

  // inv(B) by LU decomposition
  gsl_linalg_LU_decomp(B, permutation_B, &status);
  gsl_linalg_LU_invert(B, permutation_B, invB);

  # pragma omp parallel for schedule(dynamic) // faster!!!
  for (int j = 0; j < ncol; j++)
  {
    // the identical terms
    gsl_vector *x1 = gsl_vector_alloc(nrow);
    get_col_from_r_matrix(X, nrow, ncol, j, x1);

    gsl_vector *V_1i = gsl_vector_alloc(col_b);
    gsl_vector *invB_mul_V_1i = gsl_vector_alloc(col_b);
    double B_1, D;
    double invXX_11;
    gsl_matrix *invXX_22 = gsl_matrix_calloc(col_b, col_b); // all elements set to 0
    gsl_vector *invXX_21 = gsl_vector_alloc(col_b);
    double XY_1;
    gsl_vector *XY_2 = gsl_vector_alloc(col_b);
    double beta_1;
    gsl_vector *beta_2 = gsl_vector_alloc(col_b);
    gsl_vector *Yhat = gsl_vector_alloc(nrow);
    double rss, zscore1, pvalue1;

    // a
    double a;
    gsl_blas_ddot(x1, x1, &a);

    // V_1i = b'a_1i
    gsl_blas_dgemv(CblasTrans, 1.0, b, x1, 0.0, V_1i);

    // invB_mul_V_1i
    gsl_blas_dgemv(CblasNoTrans, 1.0, invB, V_1i, 0.0, invB_mul_V_1i);

    // B_1 = V_1i' mul invB mul V_1i
    gsl_blas_ddot(invB_mul_V_1i, V_1i, &B_1);

    // D = I - B_1 * invA_1i
    D =  a - B_1;

    // invXX_11
    invXX_11 = 1/a + 1/a*B_1/D;


    // invXX_22
    gsl_matrix_memcpy(invXX_22, invB);
    //gsl_blas_dsyr(CblasUpper, -1.0/D, invB_mul_V_1i, invXX_22);
    gsl_blas_dger(1.0/D, invB_mul_V_1i, invB_mul_V_1i, invXX_22);

    // invXX_21
    gsl_vector_memcpy(invXX_21, invB_mul_V_1i);
    gsl_vector_scale(invXX_21, -1.0/D);

    // X'Y
    gsl_blas_ddot(x1, Yv, &XY_1);
    //gsl_blas_dgemv(CblasTrans, 1.0, x1, Yv, 0.0, XY_1);
    gsl_blas_dgemv(CblasTrans, 1.0, b, Yv, 0.0, XY_2);

    // beta
    double tmp;
    gsl_blas_ddot(invXX_21, XY_2, &tmp);
    beta_1 = invXX_11*XY_1 + tmp;

    //gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_21, XY_1, 0.0, beta_2);
    gsl_vector_memcpy(beta_2, invXX_21);
    gsl_vector_scale(beta_2, XY_1);

    gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_22, XY_2, 1.0, beta_2);

    // RSS
    gsl_vector_memcpy(Yhat, x1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, b, beta_2, beta_1, Yhat);

    gsl_vector_sub(Yhat, Yv);

    gsl_blas_ddot(Yhat, Yhat, &rss);

    // zscore
    zscore1 = beta_1/(sqrt(rss/df*invXX_11));
    pvalue1 = 2*(zscore1 < 0 ? (1 - gsl_cdf_tdist_P(-zscore1, df)) : (1 - gsl_cdf_tdist_P(zscore1, df)));

    gsl_vector_free(x1);
    gsl_vector_free(XY_2);
    gsl_vector_free(beta_2);
    gsl_vector_free(Yhat);
    gsl_matrix_free(invXX_22);
    gsl_vector_free(invXX_21);
    gsl_vector_free(V_1i);
    gsl_vector_free(invB_mul_V_1i);
    #pragma omp critical
    {
      r.push_back(j);
      r_p.push_back(pvalue1);
    }
  }
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("r") = r,
                                                    Rcpp::Named("r.p.value") = r_p);

  gsl_vector_free(x0);
  return output;
}

//[[Rcpp::export]]
Rcpp::List frlr(SEXP R_X, SEXP R_Y, SEXP R_COV)
{
  // convert data type
  Rcpp::NumericVector X(R_X);
  Rcpp::NumericVector Y(R_Y);
  Rcpp::NumericVector COV(R_COV);
  
  int nrow = Y.size();
  int nX = X.size();
  int ncol = nX/nrow;
  
  int nCOV = COV.size();
  int COV_COL = nCOV/nrow;
  
  int col_b = COV_COL+1;
  gsl_matrix *b = gsl_matrix_alloc(nrow, col_b);
  gsl_matrix *B = gsl_matrix_alloc(COV_COL+1, col_b);
  gsl_matrix *invB = gsl_matrix_alloc(COV_COL+1, col_b);
  
  gsl_vector *covariate_vec = gsl_vector_alloc(nrow);
  for (int ip = 0; ip < COV_COL; ip++)
  {
    get_col_from_r_matrix(COV, nrow, COV_COL, ip, covariate_vec);
    gsl_matrix_set_col(b, 1+ip, covariate_vec);
  }
  gsl_vector_free(covariate_vec);
  
  gsl_vector *Yv = gsl_vector_alloc(nrow);
  get_col_from_r_matrix(Y, nrow, 1, 0, Yv);
  
  // intercept
  gsl_vector *x0 = gsl_vector_alloc(nrow);
  gsl_vector_set_all(x0, 1);
  gsl_matrix_set_col(b, 0, x0);
  
  // degree of freedom
  int df = nrow - COV_COL - 1 - 1;
  
  vector<int> r1;
  vector<int> r2;
  vector<double> r_p;
  
  gsl_permutation *permutation_B = gsl_permutation_alloc(B->size1);
  int status;
  
  // B = b'b
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, b, b, 0.0, B);
  
  // inv(B) by LU decomposition
  gsl_linalg_LU_decomp(B, permutation_B, &status);
  gsl_linalg_LU_invert(B, permutation_B, invB);
  for (int i = 0; i < ncol; i++)
  {
    # pragma omp parallel for schedule(dynamic) // faster!!!
    for (int j = i+1; j < ncol; j++)
    {
      // the identical terms
      gsl_vector *x1 = gsl_vector_alloc(nrow);
      get_col_from_r_matrix(X, nrow, ncol, i, x1);
      gsl_vector *x1tmp = gsl_vector_alloc(nrow);
      get_col_from_r_matrix(X, nrow, ncol, j, x1tmp);
      gsl_vector_mul(x1, x1tmp);
      gsl_vector_free(x1tmp);
      
      gsl_vector *V_1i = gsl_vector_alloc(col_b);
      gsl_vector *invB_mul_V_1i = gsl_vector_alloc(col_b);
      double B_1, D;
      double invXX_11;
      gsl_matrix *invXX_22 = gsl_matrix_calloc(col_b, col_b); // all elements set to 0
      gsl_vector *invXX_21 = gsl_vector_alloc(col_b);
      double XY_1;
      gsl_vector *XY_2 = gsl_vector_alloc(col_b);
      double beta_1;
      gsl_vector *beta_2 = gsl_vector_alloc(col_b);
      gsl_vector *Yhat = gsl_vector_alloc(nrow);
      double rss, zscore1, pvalue1;
      
      // a
      double a;
      gsl_blas_ddot(x1, x1, &a);
      
      // V_1i = b'a_1i
      gsl_blas_dgemv(CblasTrans, 1.0, b, x1, 0.0, V_1i);
      
      // invB_mul_V_1i
      gsl_blas_dgemv(CblasNoTrans, 1.0, invB, V_1i, 0.0, invB_mul_V_1i);
      
      // B_1 = V_1i' mul invB mul V_1i
      gsl_blas_ddot(invB_mul_V_1i, V_1i, &B_1);
      
      // D = I - B_1 * invA_1i
      D =  a - B_1;
      
      // invXX_11
      invXX_11 = 1/a + 1/a*B_1/D;
      
      
      // invXX_22
      gsl_matrix_memcpy(invXX_22, invB);
      //gsl_blas_dsyr(CblasUpper, -1.0/D, invB_mul_V_1i, invXX_22);
      gsl_blas_dger(1.0/D, invB_mul_V_1i, invB_mul_V_1i, invXX_22);
      
      // invXX_21
      gsl_vector_memcpy(invXX_21, invB_mul_V_1i);
      gsl_vector_scale(invXX_21, -1.0/D);
      
      // X'Y
      gsl_blas_ddot(x1, Yv, &XY_1);
      //gsl_blas_dgemv(CblasTrans, 1.0, x1, Yv, 0.0, XY_1);
      gsl_blas_dgemv(CblasTrans, 1.0, b, Yv, 0.0, XY_2);
      
      // beta
      double tmp;
      gsl_blas_ddot(invXX_21, XY_2, &tmp);
      beta_1 = invXX_11*XY_1 + tmp;
      
      //gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_21, XY_1, 0.0, beta_2);
      gsl_vector_memcpy(beta_2, invXX_21);
      gsl_vector_scale(beta_2, XY_1);
      
      gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_22, XY_2, 1.0, beta_2);
      
      // RSS
      gsl_vector_memcpy(Yhat, x1);
      gsl_blas_dgemv(CblasNoTrans, 1.0, b, beta_2, beta_1, Yhat);
      
      gsl_vector_sub(Yhat, Yv);
      
      gsl_blas_ddot(Yhat, Yhat, &rss);
      
      // zscore
      zscore1 = beta_1/(sqrt(rss/df*invXX_11));
      pvalue1 = 2*(zscore1 < 0 ? (1 - gsl_cdf_tdist_P(-zscore1, df)) : (1 - gsl_cdf_tdist_P(zscore1, df)));
      
      gsl_vector_free(x1);
      gsl_vector_free(XY_2);
      gsl_vector_free(beta_2);
      gsl_vector_free(Yhat);
      gsl_matrix_free(invXX_22);
      gsl_vector_free(invXX_21);
      gsl_vector_free(V_1i);
      gsl_vector_free(invB_mul_V_1i);
      #pragma omp critical
      {
        r1.push_back(i);
        r2.push_back(j);
        r_p.push_back(pvalue1);
      }
    }
  }
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("r1") = r1,
                                                   Rcpp::Named("r2") = r2,
                                                   Rcpp::Named("r.p.value") = r_p);
  
  gsl_vector_free(x0);
  return output;
}


//[[Rcpp::export]]
Rcpp::List frlr2(SEXP R_X, SEXP R_idx1, SEXP R_idx2, SEXP R_Y, SEXP R_COV)
{
  // gsl_matrix is row-major order

  // convert data type
  Rcpp::NumericVector X(R_X);
  Rcpp::NumericVector Y(R_Y);
  Rcpp::NumericVector COV(R_COV);
  Rcpp::IntegerVector idx1(R_idx1);
  Rcpp::IntegerVector idx2(R_idx2);

  int nrow = Y.size();
  int nX = X.size();
  int ncol = nX/nrow;

  int n = idx1.size();
  int nCOV = COV.size();
  int COV_COL = nCOV/nrow;

  int col_b = COV_COL+1;
  gsl_matrix *b = gsl_matrix_alloc(nrow, col_b);
  gsl_matrix *B = gsl_matrix_alloc(COV_COL+1, col_b);
  gsl_matrix *invB = gsl_matrix_alloc(COV_COL+1, col_b);

  gsl_vector *covariate_vec = gsl_vector_alloc(nrow);
  for (int ip = 0; ip < COV_COL; ip++)
  {
    get_col_from_r_matrix(COV, nrow, COV_COL, ip, covariate_vec);
    gsl_matrix_set_col(b, 1+ip, covariate_vec);
  }
  gsl_vector_free(covariate_vec);

  gsl_vector *Yv = gsl_vector_alloc(nrow);
  get_col_from_r_matrix(Y, nrow, 1, 0, Yv);

  // intercept
  gsl_vector *x0 = gsl_vector_alloc(nrow);
  gsl_vector_set_all(x0, 1);
  gsl_matrix_set_col(b, 0, x0);

  // degree of freedom
  int df = nrow - COV_COL - 2 - 1;


  vector<int> r1, r2;
  vector<double> r1_p, r2_p;

  gsl_permutation *permutation_B = gsl_permutation_alloc(B->size1);
  int status;

  // B = b'b
  gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, b, b, 0.0, B);

  // inv(B) by LU decomposition
  gsl_linalg_LU_decomp(B, permutation_B, &status);
  gsl_linalg_LU_invert(B, permutation_B, invB);

  # pragma omp parallel for schedule(dynamic) // faster!!!
  for (int j = 0; j < n; j++)
  {
    // the identical terms
    gsl_vector *x1 = gsl_vector_alloc(nrow);
    gsl_vector *x2 = gsl_vector_alloc(nrow);
    get_col_from_r_matrix(X, nrow, ncol, idx1[j]-1, x1);
    get_col_from_r_matrix(X, nrow, ncol, idx2[j]-1, x2);

    double A_1i_11, A_1i_12, A_1i_22;
    double a11, a12, a21, a22, a_det;
    gsl_matrix *invA_1i = gsl_matrix_alloc(2, 2);
    gsl_matrix *a_1i = gsl_matrix_alloc(nrow, 2);
    gsl_matrix *V_1i = gsl_matrix_alloc(col_b, 2);
    gsl_matrix *invB_mul_V_1i = gsl_matrix_alloc(col_b, 2);
    gsl_matrix *m_tmp = gsl_matrix_alloc(2, 2);
    gsl_matrix *B_1 = gsl_matrix_alloc(2, 2);
    gsl_matrix *invD = gsl_matrix_alloc(2, 2);
    gsl_matrix *m_tmp2 = gsl_matrix_alloc(col_b, 2);
    gsl_matrix *m_tmp3 = gsl_matrix_alloc(2, 2);
    gsl_matrix *m_tmp4 = gsl_matrix_alloc(col_b, 2);
    gsl_matrix *invXX_11 = gsl_matrix_alloc(2, 2);
    gsl_matrix *invXX_22 = gsl_matrix_alloc(col_b, col_b);
    gsl_matrix *invXX_21 = gsl_matrix_alloc(col_b, 2);
    gsl_vector *XY_1 = gsl_vector_alloc(2);
    gsl_vector *XY_2 = gsl_vector_alloc(col_b);
    gsl_vector *beta_1 = gsl_vector_alloc(2);
    gsl_vector *beta_2 = gsl_vector_alloc(col_b);
    gsl_vector *Yhat = gsl_vector_alloc(nrow);
    double rss, zscore1, zscore2, pvalue1, pvalue2;

    // A_1i
    gsl_blas_ddot(x1, x1, &A_1i_11);
    gsl_blas_ddot(x1, x2, &A_1i_12);
    gsl_blas_ddot(x2, x2, &A_1i_22);

    // invA_1i
    a_det = A_1i_11*A_1i_22-A_1i_12*A_1i_12;
    a11 = A_1i_22/a_det;
    a12 = -A_1i_12/a_det;
    a22 = A_1i_11/a_det;
    gsl_matrix_set(invA_1i, 0, 0, a11);
    gsl_matrix_set(invA_1i, 1, 1, a22);
    gsl_matrix_set(invA_1i, 0, 1, a12);
    gsl_matrix_set(invA_1i, 1, 0, a12);

    // construct a_1i
    gsl_matrix_set_col(a_1i, 0, x1);
    gsl_matrix_set_col(a_1i, 1, x2);

    // V_1i = b'a_1i
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, b, a_1i, 0.0, V_1i);

    // invB_mul_V_1i
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invB, V_1i, 0.0, invB_mul_V_1i);

    // B_1 = V_1i' mul invB mul V_1i
    gsl_blas_dgemm(CblasTrans, CblasNoTrans, 1.0, invB_mul_V_1i, V_1i, 0.0, B_1);

    // D = I - B_1 * invA_1i
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, B_1, invA_1i, 0.0, m_tmp);

    a11 = gsl_matrix_get(m_tmp, 0, 0);
    a22 = gsl_matrix_get(m_tmp, 1, 1);
    a12 = gsl_matrix_get(m_tmp, 0, 1);

    // D is noy symmetric !!
    a21 = gsl_matrix_get(m_tmp, 1, 0);
    a11 += 1;
    a22 += 1;
    a_det = a11*a22-a12*a21;

    gsl_matrix_set(invD, 0, 0, a22/a_det);
    gsl_matrix_set(invD, 1, 1, a11/a_det);
    gsl_matrix_set(invD, 1, 0, -a21/a_det);
    gsl_matrix_set(invD, 0, 1, -a12/a_det);

    // invXX_11
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invA_1i, B_1, 0.0, m_tmp3);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, m_tmp3, invD, 0.0, m_tmp); // Do not replace m_tmp with m_tmp3 !!
    gsl_matrix_memcpy(invXX_11, invA_1i);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, m_tmp, invA_1i, 1.0, invXX_11);

    // invXX_22
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invB_mul_V_1i, invA_1i, 0.0, m_tmp2);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, m_tmp2, invD, 0.0, m_tmp4);
    gsl_matrix_memcpy(invXX_22, invB);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, -1.0, m_tmp4, invB_mul_V_1i, 1.0, invXX_22);

    // invXX_21
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, invB, V_1i, 0.0, m_tmp2);
    gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, m_tmp2, invD, 0.0, m_tmp4);
    gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, -1.0, m_tmp4, invA_1i, 0.0, invXX_21);

    // X'Y
    gsl_blas_dgemv(CblasTrans, 1.0, a_1i, Yv, 0.0, XY_1);
    gsl_blas_dgemv(CblasTrans, 1.0, b, Yv, 0.0, XY_2);

    // beta
    gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_11, XY_1, 0.0, beta_1);
    gsl_blas_dgemv(CblasTrans, 1.0, invXX_21, XY_2, 1.0, beta_1);
    gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_21, XY_1, 0.0, beta_2);
    gsl_blas_dgemv(CblasNoTrans, 1.0, invXX_22, XY_2, 1.0, beta_2);

    // RSS
    gsl_blas_dgemv(CblasNoTrans, 1.0, a_1i, beta_1, 0.0, Yhat);
    gsl_blas_dgemv(CblasNoTrans, 1.0, b, beta_2, 1.0, Yhat);

    gsl_vector_sub(Yhat, Yv);

    gsl_blas_ddot(Yhat, Yhat, &rss);

    // zscore
    zscore1 = gsl_vector_get(beta_1, 0)/(sqrt(rss/df*gsl_matrix_get(invXX_11, 0, 0)));
    pvalue1 = 2*(zscore1 < 0 ? (1 - gsl_cdf_tdist_P(-zscore1, df)) : (1 - gsl_cdf_tdist_P(zscore1, df)));

    zscore2 = gsl_vector_get(beta_1, 1)/(sqrt(rss/df*gsl_matrix_get(invXX_11, 1, 1)));
    pvalue2 = 2*(zscore2 < 0 ? (1 - gsl_cdf_tdist_P(-zscore2, df)) : (1 - gsl_cdf_tdist_P(zscore2, df)));

    gsl_vector_free(x1);
    gsl_vector_free(x2);
    gsl_vector_free(XY_1);
    gsl_vector_free(XY_2);
    gsl_vector_free(beta_1);
    gsl_vector_free(beta_2);
    gsl_vector_free(Yhat);
    gsl_matrix_free(invXX_11);
    gsl_matrix_free(invXX_22);
    gsl_matrix_free(invXX_21);
    gsl_matrix_free(invA_1i);
    gsl_matrix_free(a_1i);
    gsl_matrix_free(V_1i);
    gsl_matrix_free(invB_mul_V_1i);
    gsl_matrix_free(m_tmp);
    gsl_matrix_free(B_1);
    gsl_matrix_free(invD);
    gsl_matrix_free(m_tmp2);
    gsl_matrix_free(m_tmp3);
    gsl_matrix_free(m_tmp4);
    #pragma omp critical
    {
      r1.push_back(idx1[j]);
      r2.push_back(idx2[j]);
      r1_p.push_back(pvalue1);
      r2_p.push_back(pvalue2);
    }
  }
  Rcpp::DataFrame output = Rcpp::DataFrame::create(Rcpp::Named("r1") = r1,
                                                    Rcpp::Named("r2") = r2,
                                                    Rcpp::Named("r1.p.value") = r1_p,
                                                    Rcpp::Named("r2.p.value") = r2_p);

  gsl_vector_free(x0);
  return output;
}
