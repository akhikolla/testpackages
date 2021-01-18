// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2013-2015 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, You can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef __TRUST_OPTIM_RFUNC_HESS__
#define __TRUST_OPTIM_RFUNC_HESS__

// uses a user-supplied function for the Hessian

#include <common_R.hpp>

using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;

class RfuncHess {

  int nvars; 
  const Rcpp::Function & fn;
  const Rcpp::Function & gr;
  const Rcpp::Function & hs;
  
protected:
    
  int nnz;

public:

  RfuncHess(const int nvars_,
	    const int nnz_,
	    const Rcpp::Function& fn_,
	    const Rcpp::Function& gr_,
	    const Rcpp::Function& hs_) :
    nvars(nvars_), fn(fn_), gr(gr_), hs(hs_), nnz(nnz_){}


  ~RfuncHess(){};

  template <typename Tpars>
  void get_f(const Eigen::MatrixBase<Tpars>& P_,
	     double& f) {
    
    Eigen::MatrixBase<Tpars>& P = const_cast<Eigen::MatrixBase<Tpars>& >(P_);
    if (P.size()!=nvars) {
      throw MyException("Incorrect number of parameters\n",
			__FILE__, __LINE__);
    }
    
    Rcpp::NumericVector pars(P.derived().data(),
			     P.derived().data() + P.derived().size());  
    f = Rcpp::as<double>(fn(pars));
    return;
  }
  
  template <typename Tpars, typename Tgrad>
  void get_df(const MatrixBase<Tpars>& P_,
	      const MatrixBase<Tgrad>& df_) {
    
    using Rcpp::NumericVector;
    
    Eigen::MatrixBase<Tpars>&
      P = const_cast<Eigen::MatrixBase<Tpars>& >(P_);
    Eigen::MatrixBase<Tgrad> &
      df = const_cast<Eigen::MatrixBase<Tgrad>& >(df_);
    
    if (P.size()!=nvars){
      throw MyException("Incorrect number of parameters\n",
			__FILE__, __LINE__);
    }

    if (df.size()!=nvars) {
      throw MyException("Incorrect gradient length\n",
			__FILE__, __LINE__);
    }
    
    NumericVector pars(P.derived().data(),
		       P.derived().data() + P.size());
  
    NumericVector grad_  = gr(pars);
    
    VectorXd grad = VectorXd::Map(grad_.begin(),
				  nvars);
    
    df = grad;
  
    return;  
  }
  
  template <typename Tpars, typename Tgrad>
  void get_fdf(const Eigen::MatrixBase<Tpars>& P_,
	       double& f_,
	       Eigen::MatrixBase<Tgrad>& df_)
  {
    get_f(P_, f_);
    get_df(P_, df_);
    return;     
  }
  

  template<typename Tpars, typename Tout>
  void get_hessian(const Eigen::MatrixBase<Tpars>& P,
		   const Eigen::SparseMatrixBase<Tout>& out_){

    SparseMatrixBase<Tout>& out = const_cast<SparseMatrixBase<Tout>&>(out_);
    
    using Rcpp::NumericVector;
    using Eigen::VectorXd;
    
    if (P.size()!=nvars) {
      throw MyException("Incorrect number of parameters\n",
			__FILE__, __LINE__);
    }
    
    NumericVector pars(P.derived().data(),
		       P.derived().data() + P.size());
    
    Rcpp::S4 sh_ = hs(pars);
    Map<SparseMatrix<double> > sh(Rcpp::as<Map<SparseMatrix<double> > >(sh_));
    out = sh.selfadjointView<Lower>();
    
  }

  inline int get_nnz() {return(nnz);}  
};

#endif




