#ifndef R2SUNDIALS_H
#define R2SUNDIALS_H

#include <cvodes/cvodes.h> // prototypes for CVODE fcts., consts.
#include <nvector/nvector_serial.h> // access to serial N_Vector
#include <sunmatrix/sunmatrix_dense.h> // access to dense SUNMatrix
#include <sunlinsol/sunlinsol_dense.h> // access to dense SUNLinearSolver
#include <sunmatrix/sunmatrix_sparse.h> // access to sparse SUNMatrix
#include <sundials/sundials_types.h> // defs. of realtype, sunindextype

#define R2SUNDIALS_EVENT_IGNORE 0
#define R2SUNDIALS_EVENT_HOLD 1
#define R2SUNDIALS_EVENT_STOP -1

#include <sunlinsol_rmumps.h>

#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

#define getmem(y, f) {(y)=(f); if ((y) == NULL) stop("no memory for " #y);}
#define check_retval(expr) {int retval=(expr); if (retval != CV_SUCCESS) stop("r2sundials: call: %s\nraised flag: %s", #expr, CVodeGetReturnFlagName(retval));}

typedef struct {
  NumericVector psens;
  List lp;
} UserData;

typedef void (*funfree)(void *);
typedef void (*funfreep)(void **);
template <typename T>
using  funfree1=void (*)(void *, T);

template <typename T>
class Sunmem {
public:
  Sunmem() {}; // constructor
  ~Sunmem(); // destructor
  void add(void **pptr, funfree f);
  void add(void **pptr, funfreep f);
  void add(void **pptr, funfree1<T> f, T arg);
  void freeall();
private:
  std::vector<void**> vecptr;
  std::vector<void**> vecptrp;
  std::vector<void**> vecptr1;
  std::vector<funfree> vecf;
  std::vector<funfreep> vecfp;
  std::vector<funfree1<T>> vecf1;
  std::vector<T> vecarg;
};
// helper class for memory freeing
// serial destructor
template<typename T>
Sunmem<T>::~Sunmem<T>() {
  freeall();
}
template<typename T>
void Sunmem<T>::freeall() {
/*
Rcout << "call ~Sunmem\n";
Rcout << "freeing\t" << vecptr.size() << " simple pointers\n";
Rcout << "freeing\t" << vecptr1.size() << " simple pointers with one argument\n";
Rcout << "freeing\t" << vecptrp.size() << " ref pointers\n";
*/
  // free simple pointers
  for (int i=vecptr.size()-1; i >= 0; i--) {
    (vecf[i])(*(vecptr[i]));
  }
  vecptr.resize(0);
  vecf.resize(0);
  // free simple pointers with an argument
  for (int i=vecptr1.size()-1; i >= 0; i--)
    (vecf1[i])(*(vecptr1[i]), vecarg[i]);
  vecptr1.resize(0);
  vecarg.resize(0);
  vecf1.resize(0);
  // free pointers by ref
  for (int i=vecptrp.size()-1; i >= 0; i--)
    (vecfp[i])(vecptrp[i]);
  vecptrp.resize(0);
  vecfp.resize(0);
}
template<typename T>
void Sunmem<T>::add(void **pptr, funfree f) {
  vecptr.push_back(pptr);
  vecf.push_back(f);
}
template<typename T>
void Sunmem<T>::add(void **pptr, funfreep f) {
  vecptrp.push_back(pptr);
  vecfp.push_back(f);
}
template <typename T>
void Sunmem<T>::add(void **pptr, funfree1<T> f, T arg) {
  vecptr1.push_back(pptr);
  vecf1.push_back(f);
  vecarg.push_back(arg);
}
//template class Sunmem<int>;

// define a type for user supplied function rhs
typedef int (*rsunRhsFn)(double t, const vec &y, vec &ydot, RObject &param, NumericVector &psens);
typedef int (*rsunJacFn)(double t, const vec &y, const vec &ydot, mat &J, RObject &param, NumericVector &psens, vec &tmp1, vec &tmp2, vec &tmp3);
typedef int (*rsunSpJacFn)(double t, const vec &y, const vec &ydot, uvec &i, uvec &p, vec &v, int n, int nz, RObject &param, NumericVector &psens, vec &tmp1, vec &tmp2, vec &tmp3);
typedef int (*rsunRootFn)(double t, const vec &y, vec &vroot, RObject &param, NumericVector &psens);
typedef int (*rsunEventFn)(double t, const vec &y, vec &ynew, int Ns, std::vector<vec> &ySv, const ivec &rootsfound, RObject &param, NumericVector &psens);
typedef int (*rsunSensFn)(int Ns, double t, const vec &yv, const vec &ydotv, const std::vector<vec> &ySv, std::vector<vec> &ySdotv, RObject &param, NumericVector &psens, const vec &tmp1v, const vec &tmp2v);
typedef int (*rsunSens1Fn)(int Ns, double t, const vec &yv, const vec &ydotv, int iS, const vec &ySv, vec &ySdotv, RObject &param, NumericVector &psens, vec &tmp1v, vec &tmp2v);

int rhswrap(realtype t, N_Vector y, N_Vector ydot, void *user_data);
int jacwrap(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int spjacwrap(realtype t, N_Vector y, N_Vector fy, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);
int rootwrap(realtype t, N_Vector y, realtype *rootout, void *user_data);
int sensrhswrap(int Ns, realtype t, N_Vector y, N_Vector ydot, N_Vector *yS, N_Vector *ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2);
int sensrhs1wrap(int Ns, realtype t, N_Vector y, N_Vector ydot, int iS, N_Vector yS, N_Vector ySdot, void *user_data, N_Vector tmp1, N_Vector tmp2);

// error handler
void rsunerr(int error_code, const char *module, const char *function, char *msg, void *eh_data);

// [[Rcpp::plugins(cpp11)]]
template <typename... Args>
inline void warningNoCall(const char* fmt, Args&&... args ) {
    Rf_warningcall(R_NilValue, tfm::format(fmt, std::forward<Args>(args)... ).c_str());
}
#endif
