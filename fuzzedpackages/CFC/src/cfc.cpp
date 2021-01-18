// [[Rcpp::plugins(openmp)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <Rcpp.h>
using namespace Rcpp;
//using namespace arma;
#ifdef _OPENMP
#include <omp.h>
#endif

typedef arma::vec (*func)(arma::vec x, void* arg, int n);
typedef void* (*initfunc)(List arg);
typedef void (*freefunc)(void *arg);

double trapezoidal_step(arma::vec fvec, arma::vec gvec) {
  double h = gvec[2] - gvec[0];
  double ret = 0.5 * h * (fvec[0] + fvec[2]);
  return ret;
}

double simpson_step(arma::vec fvec, arma::vec gvec, double epsilon = 1e-16) {
  arma::vec gdiff_abs = abs(diff(gvec));
  if (gdiff_abs[0] < epsilon) {
    arma::vec ftmp(3); ftmp[0] = fvec[1]; ftmp[1] = 0.0; ftmp[2] = fvec[2];
    arma::vec gtmp(3); gtmp[0] = gvec[1]; gtmp[1] = 0.0; gtmp[2] = gvec[2];
    return trapezoidal_step(ftmp, gtmp);
  } else if (gdiff_abs[1] < epsilon) {
    arma::vec ftmp(3); ftmp[0] = fvec[0]; ftmp[1] = 0.0; ftmp[2] = fvec[1];
    arma::vec gtmp(3); gtmp[0] = gvec[0]; gtmp[1] = 0.0; gtmp[2] = gvec[1];
    return trapezoidal_step(ftmp, gtmp);
  }
  
  double h = gvec[2] - gvec[0];
  double r = (2*gvec[1] - gvec[0] - gvec[2]) / h;
  double ret = (h / 6) * ((fvec[0] + 4 * fvec[1] + fvec[2]) + 2 * r * (fvec[0] - fvec[2]) - 3 * r*r * (fvec[0] + fvec[2])) / (1 - r*r);
  return ret;
}

double simpson_step_bak(arma::vec fvec, arma::vec gvec) {
  double h = gvec[2] - gvec[0];
  double r = (2* gvec[1] - gvec[0] - gvec[2]) / h;
  double ret = (h / 6) * ((fvec[0] + 4 * fvec[1] + fvec[2]) + 2 * r * (fvec[0] - fvec[2]) - 3 * r*r * (fvec[0] + fvec[2])) / (1 - r*r);
  return ret;
}

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
// [[Rcpp::export]]
List cscr_samples_Cpp(List func_list, List init_list, List free_list, List arg_list, arma::vec tout
  , int Nmax, double rel_tol, int nsmp, int ncore) {
  if (Nmax < 2) stop("Maximum number of subdivisions cannot be less than 2");
  int K = func_list.length();
  int nout = tout.size();
  
  Progress p(0, false);
  
  //initialize
  //func fsurv_list[K];
  //initfunc finit_list[K];
  //freefunc ffree_list[K];
  
  func *fsurv_list = new func[K];
  initfunc *finit_list = new initfunc[K];
  freefunc *ffree_list = new freefunc[K];

  for (int k = 0; k < K; k++) {
    fsurv_list[k] = *(as< XPtr<func> >(func_list[k]));
    finit_list[k] = *(as< XPtr<initfunc> >(init_list[k]));
    ffree_list[k] = *(as< XPtr<freefunc> >(free_list[k]));
  }
  
  //void* args_list[K];
  void **args_list = new void*[K];

  for (int k = 0; k < K; k++) args_list[k] = finit_list[k](arg_list[k]);
  
  // body of function
  double tmax = max(tout);
  arma::cube I_out_value(nout, K, nsmp);
  arma::uvec is_maxiter(nsmp);
  arma::cube scube(nout, K, nsmp);
  
  //std::cout << "number of cores: " << ncore << std::endl;
#ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic) num_threads(ncore) 
#endif
  for (int j = 0; j < nsmp; j++) {
    if (!Progress::check_abort()) {
    int N;
    arma::vec xvec_new(3);
    arma::mat fmat_new(3, K); // fix inefficiency of recalculating midpoints
    arma::mat fprodmat_new(3, K);
    arma::vec prodtmp_new(3);
    arma::rowvec I_trap_int_new(K), I_simp_int_new(K);
    arma::rowvec I_trap_int_1(K), I_trap_int_2(K);
    arma::rowvec I_simp_int_1(K), I_simp_int_2(K);
    arma::vec xvec(2*Nmax + 1);
    arma::mat f_mat(2*Nmax + 1, K);
    arma::mat fprod_mat(2*Nmax + 1, K);
    arma::mat I_trap_int(Nmax, K);
    arma::mat I_simp_int(Nmax, K);
    arma::mat I_trap_cum(Nmax + 1, K);
    arma::mat I_simp_cum(Nmax + 1, K);
    arma::mat err_abs_int(Nmax, K);
    arma::mat err_abs_cum(Nmax, K);
    arma::mat err_rel_cum(Nmax, K);
    arma::rowvec err_abs(K);
    arma::rowvec err_rel(K);
    double err_rel_max;
    arma::uword idx, idx2;
    arma::uvec idx_nodes(Nmax + 1);
    for (int n = 0; n < Nmax + 1; n++) idx_nodes(n) = 2*n;
    arma::vec I_interp(Nmax + 1);
    arma::vec I_gen_out(nout - 2);
    arma::mat I_trap_out_all(nout, K);
    arma::mat I_simp_out_all(tout.n_elem, K);
    arma::mat smat(nout, K);
    
    N = 1;
    xvec_new[0] = 0.0; xvec_new[1] = tmax/2; xvec_new[2] = tmax;
    prodtmp_new.fill(1.0);
    
    for (int k = 0; k < K; k++) fmat_new.col(k) = fsurv_list[k](xvec_new, args_list[k], j);
    for (int k = 0; k < K; k++) prodtmp_new = prodtmp_new % fmat_new.col(k);
    for (int k = 0; k < K; k++) fprodmat_new.col(k) = prodtmp_new / fmat_new.col(k);
    
    for (int k = 0; k < K; k++) {
      I_trap_int_new(k) = -trapezoidal_step(fprodmat_new.col(k), fmat_new.col(k));
      I_simp_int_new(k) = -simpson_step(fprodmat_new.col(k), fmat_new.col(k)); 
    }
    
    xvec.fill(0.0);
    xvec(arma::span(0, 2*N)) = xvec_new;
    
    f_mat.fill(0.0);
    f_mat(arma::span(0, 2), arma::span::all) = fmat_new;
    fprod_mat.fill(0.0);
    fprod_mat(arma::span(0, 2), arma::span::all) = fprodmat_new;
    
    I_trap_int.fill(0.0);
    I_simp_int.fill(0.0);
    I_trap_int.row(0) = I_trap_int_new;
    I_simp_int.row(0) = I_simp_int_new;
    
    I_trap_cum.fill(0.0);
    I_trap_cum(1, arma::span::all) = I_trap_int_new;
    I_simp_cum.fill(0.0);
    I_simp_cum(1, arma::span::all) = I_simp_int_new;
    
    err_abs_int = abs(I_simp_int - I_trap_int);
    err_abs_cum = abs(I_simp_cum - I_trap_cum);
    err_rel_cum = err_abs_cum / abs(I_simp_cum);
    err_abs = err_abs_cum.row(N);
    err_rel = err_rel_cum.row(N);
    err_rel_max = max(err_rel);
    
    while (err_rel_max > rel_tol && N < Nmax) {
      err_abs_int.max(idx, idx2);
      
      xvec_new[0] = mean(xvec(arma::span(2*idx, 2*idx + 1)));
      xvec_new[1] = xvec(2*idx + 1);
      xvec_new[2] = mean(xvec(arma::span(2*idx + 1, 2*idx + 2)));
      xvec(arma::span(2*idx + 4, 2*N + 2)) = xvec(arma::span(2*idx + 2, 2*N));
      xvec(arma::span(2*idx + 1, 2*idx + 3)) = xvec_new;
      
      for (int k = 0; k < K; k++) fmat_new.col(k) = fsurv_list[k](xvec_new, args_list[k], j);
      prodtmp_new.ones();
      for (int k = 0; k < K; k++) prodtmp_new = prodtmp_new % fmat_new.col(k);
      for (int k = 0; k < K; k++) fprodmat_new.col(k) = prodtmp_new / fmat_new.col(k);
      
      f_mat(arma::span(2*idx + 4, 2*N + 2), arma::span::all) = f_mat(arma::span(2*idx + 2, 2*N), arma::span::all);
      f_mat(arma::span(2*idx + 1, 2*idx + 3), arma::span::all) = fmat_new;
      fprod_mat(arma::span(2*idx + 4, 2*N + 2), arma::span::all) = fprod_mat(arma::span(2*idx + 2, 2*N), arma::span::all);
      fprod_mat(arma::span(2*idx + 1, 2*idx + 3), arma::span::all) = fprodmat_new;
      
      for (int k = 0; k < K; k++) {
        I_trap_int_1(k) = -trapezoidal_step(fprod_mat(arma::span(2*idx, 2*idx + 2), k), f_mat(arma::span(2*idx, 2*idx + 2), k));
        I_simp_int_1(k) = -simpson_step(fprod_mat(arma::span(2*idx, 2*idx + 2), k), f_mat(arma::span(2*idx, 2*idx + 2), k));
        I_trap_int_2(k) = -trapezoidal_step(fprod_mat(arma::span(2*idx + 2, 2*idx + 4), k), f_mat(arma::span(2*idx + 2, 2*idx + 4), k));
        I_simp_int_2(k) = -simpson_step(fprod_mat(arma::span(2*idx + 2, 2*idx + 4), k), f_mat(arma::span(2*idx + 2, 2*idx + 4), k));
      }
      
      if (idx < N - 1) { // right segment has to move
        I_trap_int(arma::span(idx + 2, N), arma::span::all) = I_trap_int(arma::span(idx + 1, N - 1), arma::span::all);
        I_simp_int(arma::span(idx + 2, N), arma::span::all) = I_simp_int(arma::span(idx + 1, N - 1), arma::span::all);
      }
      I_trap_int(idx, arma::span::all) = I_trap_int_1;
      I_trap_int(idx + 1, arma::span::all) = I_trap_int_2;
      I_simp_int(idx, arma::span::all) = I_simp_int_1;
      I_simp_int(idx + 1, arma::span::all) = I_simp_int_2;
      
      I_trap_cum(arma::span(1, N + 1), arma::span::all) = cumsum(I_trap_int(arma::span(0, N), arma::span::all));
      I_simp_cum(arma::span(1, N + 1), arma::span::all) = cumsum(I_simp_int(arma::span(0, N), arma::span::all));
      
      err_abs_int = abs(I_simp_int - I_trap_int);
      err_abs_cum = abs(I_simp_cum - I_trap_cum);
      err_rel_cum = err_abs_cum / abs(I_simp_cum);
      err_abs = err_abs_cum.row(N + 1);
      err_rel = err_rel_cum.row(N + 1);
      err_rel_max = max(err_rel);
      
      ++N;
    }
    
    I_trap_out_all.fill(0.0);
    I_trap_out_all.row(nout - 1) = I_trap_cum.row(N);
    I_simp_out_all.fill(0.0);
    I_simp_out_all.row(nout - 1) = I_simp_cum.row(N);

    //std::cout << "N=" << N << std::endl;
    //std::cout << "Nmax=" << Nmax << std::endl;
    is_maxiter[j] = (N == Nmax ? 1 : 0);
    
    for (int k = 0; k < K; k++) {
      I_interp = I_trap_cum(arma::span(0, N), k);
      //std::cout << I_interp.size() << std::endl;
      //I_interp.print("I_interp:");
      
      interp1(xvec(idx_nodes(arma::span(0, N))), I_interp, tout(arma::span(1, nout - 2)), I_gen_out, "linear");
      //I_gen_out.print("I_gen_out:");
      I_trap_out_all(arma::span(1, nout - 2), k) = I_gen_out;

      I_interp = I_simp_cum(arma::span(0, N), k);
      interp1(xvec(idx_nodes(arma::span(0, N))), I_interp, tout(arma::span(1, nout - 2)), I_gen_out, "linear");
      I_simp_out_all(arma::span(1, nout - 2), k) = I_gen_out;
      //tout.print();
      //xvec(idx_nodes).print();
      smat.col(k) = fsurv_list[k](tout, args_list[k], j);
    }
    
    //mat I_out_buff(nout, K, fill::ones);
    //I_out_value.slice(j) = I_out_buff;
    I_out_value.slice(j) = I_simp_out_all;
    scube.slice(j) = smat;
    }
  }
  
  // free
  for (int k = 0; k < K; k++) ffree_list[k](args_list[k]);

  delete [] fsurv_list;
  delete [] finit_list;
  delete [] ffree_list;
  delete [] args_list;
  
  return List::create(Named("ci") = I_out_value, Named("s") = scube, Named("is.maxiter") = is_maxiter, Named("n.maxiter") = sum(is_maxiter));
}
