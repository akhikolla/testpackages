#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_blas.h>
#include <time.h>
#include "fdf.h"

double dx_nrm2(gsl_multimin_fdfminimizer * s) {
  return gsl_blas_dnrm2(gsl_multimin_fdfminimizer_dx(s));
}

//int test_dx_nrm2(gsl_multimin_fdfminimizer * s, double epsabs) {
//  if (dx_nrm2(s) < epsabs) return GSL_SUCCESS;
//  else return GSL_CONTINUE;
//}

void minimize(double * theta, int * iter, int * status, char * msg,
    double * stepsize, void * params, int len,
    gsl_multimin_function_fdf * gsl_multimin_params, const double * start,
    const bool trace) {
  gsl_vector_view gsl_opt = gsl_vector_view_array(theta, len);
  gsl_vector_const_view gsl_start = gsl_vector_const_view_array(start, len);

  gsl_multimin_fdfminimizer * s = gsl_multimin_fdfminimizer_alloc(
    gsl_multimin_fdfminimizer_vector_bfgs2, gsl_multimin_params->n);
  gsl_multimin_fdfminimizer_set(s, gsl_multimin_params, &gsl_start.vector, 0.01,
    1);
  
  int status_iter;//, status_test;
  //int successes = 0;
  for (*iter = 0; *iter < 20001; (*iter)++) {
    status_iter = gsl_multimin_fdfminimizer_iterate(s);
    if (trace && *iter % 10 == 0) {
      time_t t;
      time(&t);
      Rprintf("%.19s\t%d\t%f\t%f\t%f\r", ctime(&t), *iter, s->f, dx_nrm2(s),
        gsl_blas_dnrm2(s->gradient));
      R_FlushConsole();
      if (*iter % 1000 == 0) Rprintf("\n");
    }
    if (status_iter) break;
    //status_test = test_dx_nrm2(s, 1e-6);
    //if (status_test == GSL_SUCCESS) {
    //  successes++;
    //  if (successes == 10) break;
    //} else {
    //  successes = 0;
    //}
    R_CheckUserInterrupt();
  }
  *stepsize = dx_nrm2(s);
  /*if (status_test == GSL_SUCCESS) {
    strcpy(msg, "success, norm of dx < 1e-6");
    *status = 0;
  } else*/ if (status_iter) {
    strcpy(msg, gsl_strerror(status_iter));
    *status = 1;
  } else {
    strcpy(msg, "maximum iterations reached");
    *status = 2;
  }

  if (trace && *iter % 10 != 0) {
    time_t t;
    time(&t);
    Rprintf("%.19s\t%d\t%f\t%f\t%f\n", ctime(&t), *iter, s->f, dx_nrm2(s),
      gsl_blas_dnrm2(s->gradient));
    R_FlushConsole();
  }
  if (trace && *iter % 10 == 0 && *iter % 1000 >= 10) Rprintf("\n");

  gsl_vector_memcpy(&gsl_opt.vector, s->x);

  gsl_multimin_fdfminimizer_free(s);
}

struct parameters {
  const double ** x, ** masks, * lambda;
  const int k, * inds, * p;
  const size_t m, n, len;
  const double * indices;
  const int indices_len;
  const int num_threads;
};

double gsl_obj(const gsl_vector * theta, void * params) {
  parameters * p = (parameters *) params;
  return f_obj(theta->data, p->x, p->masks, p->lambda, p->k, p->inds, p->p,
    p->m, p->n, p->len);
}

void gsl_d_obj(const gsl_vector * theta, void * params, gsl_vector * grad) {
  parameters * p = (parameters *) params;
  d_obj(grad->data, theta->data, p->x, p->masks, p->lambda, p->k, p->inds, p->p,
    p->m, p->n, p->len, p->indices, p->indices_len, p->num_threads);
}

void gsl_fd_obj(const gsl_vector * theta, void * params, double * f,
    gsl_vector * grad) {
  parameters * p = (parameters *) params;
  d_obj(grad->data, theta->data, p->x, p->masks, p->lambda, p->k, p->inds, p->p,
    p->m, p->n, p->len, p->indices, p->indices_len, p->num_threads);
  *f = f_obj(theta->data, p->x, p->masks, p->lambda, p->k, p->inds, p->p,
    p->m, p->n, p->len);
}

int prep_indices_len(const int k, size_t n, const int * p) {
  int indices_len = 0;
  for (size_t view = 0; view < n; view++) {
    for (int j = 0; j < k && j < p[view]-1; j++) {
      for (int i = j+1; i < p[view]; i++) {
        indices_len += 3;
      }
    }
  }
  return indices_len;
}

void prep_indices(double * indices, const int k, size_t n, const int * p) {
  int indices_len = 0;
  for (size_t view = 0; view < n; view++) {
    for (int j = 0; j < k && j < p[view]-1; j++) {
      for (int i = j+1; i < p[view]; i++) {
        indices[indices_len++] = (double) view;
        indices[indices_len++] = i;
        indices[indices_len++] = j;
      }
    }
  }
}

void optim(double * theta, const double * start, size_t len, const double ** x,
    const double ** masks, const int * inds, const int k, size_t m, size_t n,
    const int * p, const double * lambda, int * iter, int * status, char * msg,
    double * upval, double * stepsize, const bool trace,
    const int num_threads) {
  int indices_len = prep_indices_len(k, n, p);
  double * indices = (double *) malloc(indices_len * sizeof(double));
  prep_indices(indices, k, n, p);
  parameters params {x, masks, lambda, k, inds, p, m, n, len, indices,
    indices_len, num_threads};
  gsl_multimin_function_fdf gsl_multimin_params {gsl_obj, gsl_d_obj, gsl_fd_obj,
    params.len, (void *)&params};
  minimize(theta, iter, status, msg, stepsize, &params, params.len,
    &gsl_multimin_params, start, trace);
  double lambda_zero[4] = {0.0, 0.0, 0.0, 0.0};
  parameters up_params {x, masks, lambda_zero, k, inds, p, m, n, len, indices,
    indices_len, num_threads};
  gsl_vector_const_view gsl_opt = gsl_vector_const_view_array(theta,
    params.len);
  *upval = gsl_obj(&gsl_opt.vector, &up_params);
  if (trace) {
    time_t t;
    time(&t);
    Rprintf("%.19s\tUnpenalized loss:\t%f\n", ctime(&t), *upval);
    R_FlushConsole();
  }
  free(indices);
}

size_t nrow(SEXP x) {
  return INTEGER(getAttrib(x, R_DimSymbol))[0];
}

size_t ncol(SEXP x) {
  return INTEGER(getAttrib(x, R_DimSymbol))[1];
}

extern "C" {

  SEXP r_obj(SEXP theta, SEXP x, SEXP masks, SEXP inds, SEXP k, SEXP p,
      SEXP lambda) {
    const double ** c_x = (const double **) malloc(length(x) * sizeof(double*));
    const double ** c_masks = (const double **) malloc(length(x) * sizeof(double*));
    for (int i = 0; i < length(x); i++) {
      c_x[i] = REAL(VECTOR_ELT(x, i));
      c_masks[i] = REAL(VECTOR_ELT(masks, i));
    }
    return ScalarReal(f_obj(REAL(theta), c_x, c_masks, REAL(lambda),
      *INTEGER(k), INTEGER(inds), INTEGER(p), length(x), length(p),
      length(theta)));
  }

  SEXP r_grad(SEXP theta, SEXP x, SEXP masks, SEXP inds, SEXP k, SEXP p,
      SEXP lambda, SEXP num_threads) {
    const double ** c_x = (const double **) malloc(length(x) * sizeof(double*));
    const double ** c_masks = (const double **) malloc(length(x) * sizeof(double*));
    for (int i = 0; i < length(x); i++) {
      c_x[i] = REAL(VECTOR_ELT(x, i));
      c_masks[i] = REAL(VECTOR_ELT(masks, i));
    }
    int indices_len = prep_indices_len(*INTEGER(k), length(p), INTEGER(p));
    double * indices = (double *) malloc(indices_len * sizeof(double));
    prep_indices(indices, *INTEGER(k), length(p), INTEGER(p));
    SEXP grad = PROTECT(allocVector(REALSXP, length(theta)));
    d_obj(REAL(grad), REAL(theta), c_x, c_masks, REAL(lambda), *INTEGER(k),
      INTEGER(inds), INTEGER(p), length(x), length(p), length(theta), indices,
      indices_len, *INTEGER(num_threads));
    free(indices);
    UNPROTECT(1);
    return grad;
  }

  SEXP r_vxi(SEXP xi) {
    SEXP u = PROTECT(allocMatrix(REALSXP, nrow(xi), ncol(xi)));
    f_vxi(REAL(u), REAL(xi), nrow(xi), ncol(xi));
    UNPROTECT(1);
    return u;
  }

  SEXP r_optim(SEXP start, SEXP x, SEXP masks, SEXP inds, SEXP k, SEXP p,
      SEXP lambda, SEXP trace, SEXP num_threads) {
    const double ** c_x = (const double **) malloc(length(x) * sizeof(double*));
    const double ** c_masks = (const double **) malloc(length(x) * sizeof(double*));
    for (int i = 0; i < length(x); i++) {
      c_x[i] = REAL(VECTOR_ELT(x, i));
      c_masks[i] = REAL(VECTOR_ELT(masks, i));
    }
    SEXP theta = PROTECT(allocVector(REALSXP, length(start)));
    int iter, status;
    double upval, stepsize;
    char msg[100];
    optim(REAL(theta), REAL(start), length(theta), c_x, c_masks, INTEGER(inds),
      *INTEGER(k), length(x), length(p), INTEGER(p), REAL(lambda), &iter,
      &status, msg, &upval, &stepsize, *LOGICAL(trace), *INTEGER(num_threads));
    SEXP res = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(res, 0, theta);
    SET_VECTOR_ELT(res, 1, ScalarInteger(iter));
    SET_VECTOR_ELT(res, 2, ScalarInteger(status));
    SET_VECTOR_ELT(res, 3, ScalarReal(upval));
    SET_VECTOR_ELT(res, 4, ScalarReal(stepsize));
    SET_VECTOR_ELT(res, 5, mkString(msg));
    UNPROTECT(2);
    return res;
  }
 
  SEXP r_inv_v(SEXP t) {
    const int n = nrow(t);
    SEXP xi = PROTECT(allocMatrix(REALSXP, n, n));
    inv_v(REAL(xi), REAL(t), n);
    UNPROTECT(1);
    return xi;
  }

  SEXP r_init_parallel() {
    init_parallel();
    return R_NilValue;
  }

  static const R_CallMethodDef callMethods[]  = {
    {"r_obj", (DL_FUNC) &r_obj, 7},
    {"r_grad", (DL_FUNC) &r_grad, 8},
    {"r_vxi", (DL_FUNC) &r_vxi, 1},
    {"r_optim", (DL_FUNC) &r_optim, 9},
    {"r_inv_v", (DL_FUNC) &r_inv_v, 1},
    {"r_init_parallel", (DL_FUNC) &r_init_parallel, 0},
    {NULL, NULL, 0}
  };

  void R_init_mmpca(DllInfo *info) {
    R_registerRoutines(info, NULL, callMethods, NULL, NULL);
    R_useDynamicSymbols(info, FALSE);
    R_forceSymbols(info, TRUE);
  }
}
