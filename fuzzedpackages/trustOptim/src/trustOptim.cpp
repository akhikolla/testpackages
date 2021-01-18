// trustOptim.cpp.  Part of the trustOptim package for the R programming language.
// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2014 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.
// See the trustOptim LICENSE file for more information.


#ifndef __TRUST_OPTIM_RINTERFACE
#define __TRUST_OPTIM_RINTERFACE


#include <common_R.hpp>

#include <CG-sparse.h>
#include <CG-quasi.h>
#include <Rfunc.cpp>
#include <RfuncHess.cpp>

using Rcpp::Function;
using Rcpp::List;
using Rcpp::as;
using Rcpp::NumericVector;
using Rcpp::IntegerVector;

//[[Rcpp::export]]
List sparseTR(NumericVector start,
	      Function fn,
	      Function gr,
	      Function hs,
	      const List control) {
  

    
  // use this version when the user supplies his own Hessian function.
  //  Hessian function must return a Matrix of class dgCMatrix
  
  using Eigen::VectorXi;
  using Eigen::Map;
    
  typedef SparseMatrix<double> optHessType;
  typedef SimplicialLLT<optHessType> optPrecondType; 

  int nvars = start.size();
  if (nvars<=0) throw MyException("Number of variables (starting values) must be positive\n",__FILE__, __LINE__);
  
    // Control parameters for optimizer

    double rad = as<double>(control["start.trust.radius"]);
    const double min_rad = as<double>(control["stop.trust.radius"]);
    const double tol = as<double>(control["cg.tol"]);
    const double prec = as<double>(control["prec"]);
    const int report_freq = as<int>(control["report.freq"]);
    const int report_level = as<int>(control["report.level"]);
    const int report_precision = as<int>(control["report.precision"]);
    const int maxit = as<int>(control["maxit"]);
    const double contract_factor = as<double>(control["contract.factor"]);
    const double expand_factor = as<double>(control["expand.factor"]);
    const double contract_threshold = as<double>(control["contract.threshold"]);
    const double expand_threshold_rad = as<double>(control["expand.threshold.radius"]);
    const double expand_threshold_ap = as<double>(control["expand.threshold.ap"]);
    const double function_scale_factor = as<double>(control["function.scale.factor"]);
    const int precond_refresh_freq = as<int>(control["precond.refresh.freq"]);
    const int precond_ID = as<int>(control["preconditioner"]);
    const int trust_iter = as<int>(control["trust.iter"]);

    Map<VectorXd> startX(start.begin(),nvars); 
     
    Rcpp::S4 sh_  = hs(startX);
    Map<SparseMatrix<double> > sh = Rcpp::as<Map<SparseMatrix<double>>>(sh_);
    int nnz = (sh.nonZeros() + nvars)/2;

    RfuncHess func(nvars, nnz, fn, gr, hs);
  
    Trust_CG_Sparse<Map<VectorXd>, RfuncHess,optHessType, optPrecondType> 
	opt(func, startX, rad, min_rad, tol, prec,
	    report_freq, report_level, report_precision,
	    maxit, contract_factor, expand_factor,
	    contract_threshold, expand_threshold_rad,
	    expand_threshold_ap, function_scale_factor,
	    precond_refresh_freq, precond_ID, trust_iter);

  
    opt.run();

    // collect results and return

    VectorXd P(nvars);
    VectorXd grad(nvars);
    optHessType hess(nvars,nvars);
    hess.reserve(nnz);

    double fval, radius;
    int iterations;
    MB_Status status;

    status = opt.get_current_state(P, fval, grad, hess, iterations, radius);

    List res;
    res = List::create(Rcpp::Named("fval") = Rcpp::wrap(fval),
		       Rcpp::Named("solution") = Rcpp::wrap(P),
		       Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		       Rcpp::Named("hessian") = Rcpp::wrap(hess),
		       Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		       Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		       Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		       Rcpp::Named("nnz") = Rcpp::wrap(nnz),
		       Rcpp::Named("method") = Rcpp::wrap("Sparse")
		       );
   
    return(res);
  
	}




//[[Rcpp::export]]
List  quasiTR(NumericVector start,
	      Function fn,
	      Function gr,
	      const List control) {
  
  
  using Eigen::VectorXi;
  using Eigen::Map;

  typedef MatrixXd optHessType;
  typedef LLT<optHessType> optPrecondType;
  
  int nvars = start.size();
  
  double rad = as<double>(control["start.trust.radius"]);
  const double min_rad = as<double>(control["stop.trust.radius"]);
  const double tol = as<double>(control["cg.tol"]);
  const double prec = as<double>(control["prec"]);
  const int report_freq = as<int>(control["report.freq"]);
  const int report_level = as<int>(control["report.level"]);
  const int report_precision = as<int>(control["report.precision"]);
  const int maxit = as<int>(control["maxit"]);
  const double contract_factor = as<double>(control["contract.factor"]);
  const double expand_factor = as<double>(control["expand.factor"]);
  const double contract_threshold = as<double>(control["contract.threshold"]);
  const double expand_threshold_rad = as<double>(control["expand.threshold.radius"]);
  const double expand_threshold_ap = as<double>(control["expand.threshold.ap"]);
  const double function_scale_factor = as<double>(control["function.scale.factor"]);
  const int precond_refresh_freq = as<int>(control["precond.refresh.freq"]);
  const int precond_ID = as<int>(control["preconditioner"]);
  const int quasi_newton_method = as<int>(control["quasi.newton.method"]);
  const int trust_iter = as<int>(control["trust.iter"]);
  
  std::string method_string;
  
  if (quasi_newton_method==1) {
    method_string = "SR1";
  } else {
    method_string = "BFGS";
  }
  
  Rfunc func(nvars, fn, gr);
  
  Map<VectorXd> startX(start.begin(),nvars); 
  
  // Control parameters for optimizer
  
  Trust_CG_Optimizer<Map<VectorXd>, Rfunc,
		     optHessType, optPrecondType> opt(func,
						      startX, rad, min_rad, tol,
						      prec, report_freq,
						      report_level,
						      report_precision,
						      maxit, contract_factor,
						      expand_factor,
						      contract_threshold,
						      expand_threshold_rad,
						      expand_threshold_ap,
						      function_scale_factor,
						      precond_refresh_freq,
						      precond_ID,
						      quasi_newton_method,
						      trust_iter);
  
  opt.run();
  
  // collect results and return
  
  VectorXd P(nvars);
  VectorXd grad(nvars);
  
  // return sparse hessian information in CSC format
  
  double fval, radius;
  int iterations;
  MB_Status status;
  
  status = opt.get_current_state(P, fval, grad,
				 iterations, radius);
  
  List res;
  res = List::create(Rcpp::Named("fval") = Rcpp::wrap(fval),
		     Rcpp::Named("solution") = Rcpp::wrap(P),
		     Rcpp::Named("gradient") = Rcpp::wrap(grad),	
		     Rcpp::Named("iterations") = Rcpp::wrap(iterations),
		     Rcpp::Named("status") = Rcpp::wrap((std::string) MB_strerror(status)),
		     Rcpp::Named("trust.radius") = Rcpp::wrap(radius),
		     Rcpp::Named("method") = Rcpp::wrap(method_string),
		     Rcpp::Named("hessian.update.method") = Rcpp::wrap(quasi_newton_method)
		     );
  
  return(res);
  
}




#endif
