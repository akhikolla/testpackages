/****************************************************************************
 * Copyright (C) 2014 by Taylor Arnold, Veeranjaneyulu Sadhanala,           *
 *                       Ryan Tibshirani                                    *
 *                                                                          *
 * This file is part of the glmgen library / package.                       *
 *                                                                          *
 *   glmgen is free software: you can redistribute it and/or modify it      *
 *   under the terms of the GNU Lesser General Public License as published  *
 *   by the Free Software Foundation, either version 2 of the License, or   *
 *   (at your option) any later version.                                    *
 *                                                                          *
 *   glmgen is distributed in the hope that it will be useful,              *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of         *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          *
 *   GNU Lesser General Public License for more details.                    *
 *                                                                          *
 *   You should have received a copy of the GNU Lesser General Public       *
 *   License along with glmgen. If not, see <http://www.gnu.org/licenses/>. *
 ****************************************************************************/

#include <string.h>
#include <stdio.h>
#include <math.h>

#include "lattice.h"
#include "tf.h"
#include "utils.h"
#include "cs.h"
#include "proxTV.h"

#include <R_ext/Print.h>

/**
 * @brief Calculates the left hand side of the beta0 step
 *        in the (optionally) constrained lattice. The
 *        weights should not have any missing values.
 *
 * @param N          number of observations
 * @param rho        tuning parameter for the ADMM algorithm
 * @param w          a vector of sample weights
 * @param E          a sparse matrix for supplying a set of linear constraints
 *                   on the fused lasso. Input a zero row matrix for the
 *                   unconstrained version
 * @return gqr
 */
gqr * calc_lattice_lhs(int N, double rho, double *w, cs *E)
{
  int k;
  cs *eye;
  cs *T;
  cs *Et;
  cs *EtE;
  cs *lhs;
  gqr * lhs_qr;

  T = cs_spalloc (N, N, N, 1, 1) ;
  for (k = 0 ; k < N; k++)
  {
    cs_entry (T, k, k, w[k] + rho);
  }
  eye = cs_compress (T);
  Et = cs_transpose(E, 1);
  EtE = cs_multiply(Et, E);

  lhs = cs_add(eye, EtE, 1, 1);
  lhs_qr = glmgen_qr(lhs);

  cs_spfree(T);
  cs_spfree(Et);
  cs_spfree(EtE);
  cs_spfree(lhs);

  return(lhs_qr);
}

/**
 * @brief Main wrapper for fitting fused lasso over a lattice.
 *   Calculates one value of lambda per cal.
 *
 * @param y             a vector of responses
 * @param w             a vector of sample weights
 * @param n             number of rows in the lattice
 * @param m             number of columns in the lattice
 * @param p             depth of the lattice; 1 for 2D grids
 * @param max_iter      maximum number of ADMM interations
 * @param lambda        lambda value in the fused lasso
 * @param rho           tuning parameter for the ADMM algorithm
 * @param eps           stopping criterion
 * @param verbose       0/1 flag for printing progress
 * @param na_flag       0/1 flag for whether any missing values are present
 * @param beta0         allocated space of size n*m*p for the output solution
 * @param beta1         allocated space of size n*m*p
 * @param beta2         allocated space of size n*m*p
 * @param beta3         allocated space of size n*m*p
 * @param thisy1        allocated space of size n*m*p
 * @param thisy2        allocated space of size n*m*p
 * @param thisy3        allocated space of size n*m*p
 * @param thisy4        allocated space of size n*m*p
 * @param u1            allocated space of size n*m*p
 * @param u2            allocated space of size n*m*p
 * @param u3            allocated space of size n*m*p
 * @param u4            allocated space of size d
 * @param E             a sparse matrix for supplying a set of linear constraints
 *                        on the fused lasso. Input a zero row matrix for the
 *                        unconstrained version
 * @param c             gives the right hand side of the constrain
 * @param d             number of constraints (number of rows of E)
 * @param buff          allocated space of size n*m*p
 * @param abuff         allocated space of size n*m*p
 * @param wbuff         allocated space of size n*m*p
 * @param lattice_type  integer code for the type of lattice
 * @param method_type   integer code; method for solving the 1D subproblem
 * @return void
 */
void do_lattice (double *y, double *w, int n, int m, int p,
                  int max_iter, double *lambda, double rho,
                  double eps, int verbose, int naflag,
                  double *beta0,
                  double *beta1, double *beta2, double *beta3,
                  double *thisy1, double *thisy2,  double *thisy3,
                  double *thisy4,
                  double *u1, double *u2, double *u3, double *u4,
                  cs *E, double *c, int d,
                  double *buff, double *abuff, double *wbuff,
                  int lattice_type, int method_type)
{
/*  int i;
  int j;*/
  int l;
  int k;
  int N;
  int iter;
  int rho_factor;
  double err0;
  double err1;
  double err2;
  double err3;
  double err4;
  cs *Et;
  gqr * lhs_qr;

  iter = 0;
  err0 = eps*2;
  err1 = eps*2;
  err2 = eps*2;
  err3 = eps*2;
  err4 = eps*2;
  N = n*m*p;
  rho_factor = (lattice_type == LATTICE_2D_GRID) ? 2 : 3;

  Et = cs_transpose(E, 1);
  lhs_qr = calc_lattice_lhs(N, rho_factor*rho, w, E);

  if (verbose)
  {
    Rprintf("lambda = %04.3f\n"
           "==================================================\n",
            lambda[0]);
  }

  while ((err0 >= eps || err1 >= eps || err2 >= eps || err3 >= eps || err4 >= eps ||
          iter < MIN_ITER) && iter < max_iter)
  {
    /* Update beta0 */
    for (k = 0; k < N; k++)
    {

      beta0[k] = (rho * (beta1[k] + beta2[k]) + y[k]*w[k] -
                   u1[k] - u2[k] - u3[k]);
    }

    for (k = 0; k < d; k++)
    {
      buff[k] = c[k] - u4[k];
    }

    cs_gaxpy(Et, buff, beta0);
    glmgen_qrsol(lhs_qr, beta0);

    /* Update beta1, beta2, beta3 */
    err0 = 0;
    for (k = 0; k < N; k++)
    {
      thisy1[k] = beta0[k] + u1[k] / rho;
      thisy2[k] = beta0[k] + u2[k] / rho;
      thisy3[k] = beta0[k] + u3[k] / rho;

      if (!isnan(y[k]))
      {
        err0 = MAX(err0, fabs(thisy4[k] - beta0[k]));
      }
      thisy4[k] = beta0[k];
    }

    switch (lattice_type) {
      case LATTICE_2D_GRID:
        if (!naflag) {
          do_d1_cols(thisy1, buff, abuff, beta1, n, m,
            lambda, method_type);
          do_d1_rows(thisy2, buff, abuff, beta2, n, m,
            lambda+(n-1)*m, method_type);
        } else {
          do_d1_cols_na(thisy1, buff, abuff, wbuff, beta1, n, m,
            lambda, method_type);
          do_d1_rows_na(thisy2, buff, abuff, wbuff, beta2, n, m,
            lambda+(n-1)*m, method_type);
        }
        break;

      case LATTICE_HEX_GRID:
        if (!naflag) {
          do_d1_cols(thisy1, buff, abuff, beta1, n, m,
            lambda, method_type);
          do_d1_rows(thisy2, buff, abuff, beta2, n, m,
            lambda+(n-1)*m, method_type);
          do_d1_hexs(thisy3, buff, abuff, beta3, n, m,
            lambda+(n-1)*m+n*(m-1), method_type);
        } else {
          do_d1_cols_na(thisy1, buff, abuff, wbuff, beta1, n, m,
            lambda, method_type);
          do_d1_rows_na(thisy2, buff, abuff, wbuff, beta2, n, m,
            lambda+(n-1)*m, method_type);
          do_d1_hexs_na(thisy3, buff, abuff, wbuff, beta3, n, m,
            lambda+(n-1)*m+n*(m-1), method_type);
        }
        break;

      case LATTICE_3D_GRID:
        for (l = 0; l < p; l++)
        {
          if (!naflag) {
            do_d1_cols(thisy1 + l*n*m, buff, abuff, beta1 + l*n*m,
                        n, m, lambda, method_type);
            do_d1_rows(thisy2 + l*n*m, buff, abuff, beta2 + l*n*m,
                        n, m, lambda+(n-1)*m, method_type);
          } else {
            do_d1_cols_na(thisy1 + l*n*m, buff, abuff, wbuff,
                        beta1 + l*n*m, n, m, lambda, method_type);
            do_d1_rows_na(thisy2 + l*n*m, buff, abuff, wbuff,
                        beta2 + l*n*m, n, m, lambda+(n-1)*m,
                        method_type);
          }
        }

        if (!naflag) {
          do_d1_rows(thisy3, buff, abuff, beta3, n*m, p,
            lambda+(n-1)*m+n*(m-1), method_type);
        } else {
          do_d1_rows_na(thisy3, buff, abuff, wbuff, beta3, n*m, p,
            lambda+(n-1)*m+n*(m-1), method_type);
        }
        break;
    }

    /* Update u1, u2, u3 */
    err1 = err2 = err3 = 0;
    for (k = 0; k < N; k++)
    {
      u1[k] += rho * (beta0[k] - beta1[k]);
      u2[k] += rho * (beta0[k] - beta2[k]);
      if (lattice_type != LATTICE_2D_GRID)
      {
        u3[k] += rho * (beta0[k] - beta3[k]);
      }
      if (!(isnan(*(y + k))))
      {
        err1 = MAX(err1, fabs(beta0[k] - beta1[k]));
        err2 = MAX(err2, fabs(beta0[k] - beta2[k]));
        if (lattice_type != LATTICE_2D_GRID)
        {
          err3 = MAX(err3, fabs(beta0[k] - beta3[k]));
        }
      }
    }

    /* Update u4 */
    for (k = 0; k < d; k++)
    {
      buff[k] = -1 * c[k];
    }
    cs_gaxpy(E, beta0, buff);

    err4 = 0;
    for (k = 0; k < d; k++)
    {
      u4[k] += rho * buff[k];
      err4 = MAX(err4, fabs(buff[k]));
    }

    if (verbose)
    {
      Rprintf("iter #%03d =>\n  ||bold - bnew||: %02.4f  "
             "||b0 - b1||: %02.4f  "
             "||b0 - b2||: %02.4f  "
             "||E b0 - c||: %02.4f\n",
              iter, err0, err1, err2, err4);
    }
    iter++;
  }

  if (verbose)
  {
    Rprintf("\n");
  }

  glmgen_gqr_free(lhs_qr);
  cs_spfree(Et);
}

/**
 * @brief Low-level function for calculating fused lasso
 *   across the columns of a grid using DP algorithm
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n*m
 * @param abuff        allocated space of size n*m
 * @param ans          allocated space of size n*m;
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_cols (double *y, double *buff,
                  double *abuff, double *ans, int n, int m,
                  double *lambda, int method_type)
{
  int i;
  /*int j;*/

  for (i = 0; i < m; i++)
  {
    memcpy(buff,  y + i*n, sizeof(double) * n);
    switch(method_type)
    {
      case LATTICE_DP:
      tf_dp(n, buff, lambda[0], abuff);
      break;

      case LATTICE_PROX:
      tautString_TV1(buff, lambda[0], abuff, n);
      break;

      case LATTICE_PROX_W:
      tautString_TV1_Weighted(buff, lambda+i*(n-1), abuff, n);
      break;
    }
    memcpy(ans + i*n, abuff, sizeof(double) * n);
  }

}


/**
 * @brief Low-level function for calculating fused lasso
 *   across the rows of a grid using DP algorithm
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n*m
 * @param abuff        allocated space of size n*m
 * @param ans          allocated space of size n*m;
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_rows (double *y, double *buff,
                  double *abuff, double *ans, int n, int m,
                  double *lambda, int method_type)
{
  int i;
  int j;

  for (j = 0; j < n; j++)
  {
    for (i = 0; i < m; i++)
    {
      memcpy(buff + i,  y + j + i*n, sizeof(double));
    }
    switch(method_type)
    {
      case LATTICE_DP:
      tf_dp(m, buff, lambda[0], abuff);
      break;

      case LATTICE_PROX:
      tautString_TV1(buff, lambda[0], abuff, m);
      break;

      case LATTICE_PROX_W:
      tautString_TV1_Weighted(buff, lambda+j*(m-1), abuff, m);
      break;
    }
    for (i = 0; i < m; i++)
    {
      memcpy(ans + j + i*n, abuff + i, sizeof(double));
    }
  }

}

/**
 * @brief Low-level function for calculating fused lasso
 *   across the diagonals of a hexagonal grid using
 *   DP algorithm
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n*m
 * @param abuff        allocated space of size n*m
 * @param ans          allocated space of size n*m;
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_hexs (double *y, double *buff,
                  double *abuff, double *ans, int n, int m,
                  double *lambda, int method_type)
{
  int i;
  int j;

  for (i = 1; i < m; i++)
  {
    for (j = 0; j < n; j++)
    {
      if ((j % 2 == 0))
      {
        memcpy(buff + j,  y + j + i*n,     sizeof(double));
      } else {
        memcpy(buff + j,  y + j + (i-1)*n, sizeof(double));
      }
    }
    switch(method_type)
    {
      case LATTICE_DP:
      tf_dp(n, buff, lambda[0], abuff);
      break;

      case LATTICE_PROX:
      tautString_TV1(buff, lambda[0], abuff, n);
      break;

      case LATTICE_PROX_W:
      tautString_TV1_Weighted(buff, lambda+i*(n-1), abuff, n);
      break;
    }
    for (j = 0; j < n; j++)
    {
      if ((j % 2 == 0))
      {
        memcpy(ans + j + i*n,     abuff + j, sizeof(double));
      } else {
        memcpy(ans + j + (i-1)*n, abuff + j, sizeof(double));
      }
    }
  }

  /* Need to fill in the start of the odd rows and ends
     of the even rows */
  for (j = 0; j < n; j++)
  {
    if ((j % 2 == 0))
    {
      memcpy(ans + j, y + j, sizeof(double));
    } else {
      memcpy(ans + j + (m-1)*n, y + j + (m-1)*n, sizeof(double));
    }
  }
}

/**
 * @brief Low-level function for calculating fused lasso
 *   across the columns of a square grid using DP algorithm
 *   where some values of y are NaN
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n*m
 * @param abuff        allocated space of size n*m
 * @param wbuff        allocated space of size n*m
 * @param ans          allocated space of size n*m;
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_cols_na (double *y, double *buff, double *abuff,
                    double *wbuff, double *ans, int n, int m,
                    double *lambda, int method_type)
{
  int i;
  int j;
  int k;
  int len;

  len = 0;
  for (i = 0; i < m; i++)
  {
    for (j = 0; j <= n; j++) {
      if ( j == n || isnan(*(y + i*n + j)) ) {
        if (len != 0)
        {
          switch(method_type)
          {
            case LATTICE_DP:
            tf_dp(len, buff, lambda[0], abuff);
            break;

            case LATTICE_PROX:
            tautString_TV1(buff, lambda[0], abuff, len);
            break;

            case LATTICE_PROX_W:
            tautString_TV1_Weighted(buff, wbuff, abuff, len);
            break;
          }
          for (k = j - len; k < j; k++)
          {
            memcpy(ans + i*n + k, abuff + k - (j - len),
                    sizeof(double));
          }
          len = 0;
        }
      } else {
        memcpy(buff + len,  y + i*n + j, sizeof(double));
        memcpy(wbuff + len, lambda + i*(n-1) + j, sizeof(double));
        len++;
      }
    }
  }

}

/**
 * @brief Low-level function for calculating fused lasso
 *   across the rows of a grid using DP algorithm
 *   where some values of y are NaN
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n*m
 * @param abuff        allocated space of size n*m
 * @param wbuff        allocated space of size n*m
 * @param ans          allocated space of size n*m;
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_rows_na (double *y, double *buff, double *abuff,
                    double *wbuff, double *ans, int n, int m,
                    double *lambda, int method_type)
{
  int i;
  int j;
  int k;
  int len;

  len = 0;
  for (j = 0; j < n; j++)
  {
    for (i = 0; i <= m; i++) {
      if ( i == m || isnan(*(y + i*n + j))) {
        if (len != 0)
        {
          switch(method_type)
          {
            case LATTICE_DP:
            tf_dp(len, buff, lambda[0], abuff);
            break;

            case LATTICE_PROX:
            tautString_TV1(buff, lambda[0], abuff, len);
            break;

            case LATTICE_PROX_W:
            tautString_TV1_Weighted(buff, wbuff, abuff, len);
            break;
          }
          for (k = i - len; k < i; k++)
          {
            memcpy(ans + k*n + j, abuff + k - (i - len),
                    sizeof(double));
          }
          len = 0;
        }
      } else {
        memcpy(buff + len,  y + i*n + j, sizeof(double));
        memcpy(wbuff + len, lambda + i + j*(m-1), sizeof(double));
        len++;
      }
    }
  }

}

/**
 * @brief Low-level function for calculating fused lasso
 *   across the diagonals of a hexagonal grid using
 *   DP algorithm where some values of y are NaN
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n*m
 * @param abuff        allocated space of size n*m
 * @param wbuff        allocated space of size n*m
 * @param ans          allocated space of size n*m
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_hexs_na (double *y, double *buff, double *abuff,
                    double *wbuff, double *ans, int n, int m,
                    double *lambda, int method_type)
{
  int i;
  int j;
  int k;
  int len;
  double this_val;

  len = 0;
  for (i = 1; i < m; i++)
  {
    for (j = 0; j <= n; j++)
    {
      if ((j % 2 == 0))
      {
        this_val = *(y + j + i*n);
      } else {
        this_val = *(y + j + (i-1)*n);
      }
      if ( j == n || isnan(this_val) ) {
        if (len != 0)
        {
          switch(method_type)
          {
            case LATTICE_DP:
            tf_dp(len, buff, lambda[0], abuff);
            break;

            case LATTICE_PROX:
            tautString_TV1(buff, lambda[0], abuff, len);
            break;

            case LATTICE_PROX_W:
            tautString_TV1_Weighted(buff, wbuff, abuff, len);
            break;
          }
          for (k = j - len; k < j; k++)
          {
            if ((k % 2 == 0))
            {
              memcpy(ans + k + i*n, abuff + k - (j - len),
                      sizeof(double));
            } else {
              memcpy(ans + k + (i-1)*n, abuff + k - (j - len),
                      sizeof(double));
            }
          }
          len = 0;
        }
      } else {
        *(buff + len) = this_val;
        memcpy(wbuff + len, lambda + (i-1)*(n-1) + j, sizeof(double));
        len++;
      }
    }
  }

  /* Need to fill in the start of the odd rows and ends
     of the even rows */
  for (j = 0; j < n; j++)
  {
    if ((j % 2 == 0))
    {
      memcpy(ans + j,           y + j,           sizeof(double));
    } else {
      memcpy(ans + j + (m-1)*n, y + j + (m-1)*n, sizeof(double));
    }
  }
}

/**
 * @brief Low-level function for calculating fused lasso
 *   across a chain of a graph
 *
 * @param y            a vector of responses
 * @param buff         allocated space of size n
 * @param abuff        allocated space of size n
 * @param wbuff        allocated space of size n
 * @param ebuff        allocated space of size n
 * @param ans          allocated space of size n
 *                       contains the solution on return
 * @param n            number of rows in the lattice
 * @param m            number of columns in the lattice
 * @param lambda       lambda value for fused lasso
 * @param method_type  integer code for the sub-method
 * @return void
 */
void do_d1_chain (double *y, double *buff, double *abuff,
                  double *wbuff, int *ebuff, double *ans, int n,
                  int *e, int edge_len,
                  double *lambda, int method_type)
{
  int j;
  int k;
  int len;

  len = 0;
  for (j = 0; j <= edge_len; j++)
  {
    if (j == edge_len || e[j] < 0 || isnan(y[e[j]])) {
      if (len != 0)
      {
        switch(method_type)
        {
          case LATTICE_DP:
          tf_dp(len, buff, lambda[0], abuff);
          break;

          case LATTICE_PROX:
          tautString_TV1(buff, lambda[0], abuff, len);
          break;

          case LATTICE_PROX_W:
          tautString_TV1_Weighted(buff, wbuff, abuff, len);
          break;
        }
        for (k = 0; k < len; k++)
        {
          ans[ebuff[k]] = abuff[k];
        }
        len = 0;
      }
    } else {
      buff[len] = y[e[j]];
      wbuff[len] = lambda[j];
      ebuff[len] = e[j];
      len++;
    }
  }

}

/**
 * @brief Main wrapper for fitting fused lasso over a graph
 *   that has been decomposed into chains
 *
 * @param y             a vector of responses
 * @param w             a vector of sample weights
 * @param e             a collapsed array of the edges
 * @param elen          an array of length num_chains giving the length
 *                       of the edge lists.
 * @param n             number of nodes in the graph
 * @param num_chains    number of chains in the decomposition
 * @param max_iter      maximum number of ADMM interations
 * @param lambda        lambda value in the fused lasso
 * @param rho           tuning parameter for the ADMM algorithm
 * @param eps           stopping criterion
 * @param verbose       0/1 flag for printing progress
 * @param beta0         allocated space of size n for the output solution
 * @param B             allocated space of size n*num_chains
 * @param thisY         allocated space of size n*(num_chains+1)
 * @param U             allocated space of size n*num_chains
 * @param u4            allocated space of size d
 * @param E             a sparse matrix for supplying a set of linear constraints
 *                        on the fused lasso. Input a zero row matrix for the
 *                        unconstrained version
 * @param c             gives the right hand side of the constrain
 * @param d             number of constraints (number of rows of E)
 * @param buff          allocated space of size n
 * @param abuff         allocated space of size n
 * @param wbuff         allocated space of size n
 * @param ebuff         allocated space of size n
 * @param method_type   integer code; method for solving the 1D subproblem
 * @return void
 */
void do_fused_graph (double *y, double *w, int *e, int *elen,
                      int n, int num_chains,
                      int max_iter, double *lambda, double rho,
                      double eps, int verbose,
                      double *beta0, double *B, double *thisY,
                      double *U, double *u4,
                      cs *E, double *c, int d,
                      double *buff, double *abuff, double *wbuff,
                      int *ebuff, int method_type)
{
  /*int i;*/
  /*int j;*/
  /*int l;*/
  int k;
  int m;
  int iter;
  int e_offset;
  double err0;
  double err1;
  cs *Et;
  gqr * lhs_qr;

  iter = 0;
  err0 = eps*2;
  err1 = eps*2;

  Et = cs_transpose(E, 1);
  lhs_qr = calc_lattice_lhs(n, num_chains*rho, w, E);

  if (verbose)
  {
    Rprintf("lambda = %04.3f\n"
           "==================================================\n",
            lambda[0]);
  }

  while ((err0 >= eps || err1 >= eps || iter < MIN_ITER) && iter < max_iter)
  {
    /* Update beta0 */
    for (k = 0; k < n; k++)
    {
      beta0[k] = y[k]*w[k];
      for (m = 0; m < num_chains; m++)
      {
        beta0[k] += rho * B[k+m*n] - U[k+m*n];
      }
    }

    for (k = 0; k < d; k++)
    {
      buff[k] = c[k] - u4[k];
    }

    cs_gaxpy(Et, buff, beta0);
    glmgen_qrsol(lhs_qr, beta0);

    /* Update beta's */
    err0 = 0;
    for (k = 0; k < n; k++)
    {
      for (m = 0; m < num_chains; m++)
      {
        thisY[k+m*n] = beta0[k] + U[k+m*n] / rho;
      }
      if (!isnan(y[k]))
      {
        err0 = MAX(err0, fabs(beta0[k] - thisY[k+num_chains*n]));
      }
      thisY[k+num_chains*n] = beta0[k];
    }

    e_offset = 0;
    for (m = 0; m < num_chains; m++)
    {
      do_d1_chain (thisY+m*n, buff, abuff, wbuff, ebuff, B + m*n, n, e + e_offset,
                   elen[m], lambda + e_offset, method_type);
      e_offset += elen[m];
    }

    /* Update U */
    err1 = 0;
    for (k = 0; k < n; k++)
    {
      for (m = 0; m < num_chains; m++)
      {
        U[k+m*n] += rho * (beta0[k] - B[k+m*n]);
        if (!(isnan(*(y + k))))
        {
          err1 = MAX(err1, fabs(beta0[k] - B[k+m*n]));
        }
      }
    }

    /* Update u4 */
    for (k = 0; k < d; k++)
    {
      buff[k] = -1 * c[k];
    }
    cs_gaxpy(E, beta0, buff);

    for (k = 0; k < d; k++)
    {
      u4[k] += rho * buff[k];
      err1 = MAX(err1, fabs(buff[k]));
    }

    if (verbose)
    {
      Rprintf("iter #%03d =>\n  max delta: %02.6f  max error: %02.6f\n", iter, err0, err1);
    }
    iter++;
  }

  if (verbose)
  {
    Rprintf("\n");
  }

  glmgen_gqr_free(lhs_qr);
  cs_spfree(Et);
}
