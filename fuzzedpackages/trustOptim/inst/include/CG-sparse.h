// CG-sparse.cpp.  Part of the trustOptim package for the R programming language.
//
// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2013 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.


#ifndef __TRUST_OPTIM_CGSPARSE
#define __TRUST_OPTIM_CGSPARSE

#include <CG-base.h>

using Eigen::Dynamic;

template<typename TP, typename TFunc, typename THess, typename TPreLLt>  // TP is type for the parameter vector
class Trust_CG_Sparse : public Trust_CG_Base<TP, TFunc, THess, TPreLLt> {

    typedef typename TP::Index Index;
    typedef Trust_CG_Base<TP, TFunc, THess, TPreLLt> BaseType;

    typedef SparseMatrix<double> SparseMatrixXd;

    using BaseType::Bk;
    using BaseType::nvars;
    using BaseType::func;
    using BaseType::function_scale_factor;
    using BaseType::xk;
    using BaseType::precond_ID;
    using BaseType::PrecondLLt;
    using BaseType::status;
    using BaseType::rad;
    using BaseType::iter;

private:

    Index nnz;

    void update_precond() {

	switch (precond_ID) {
	case 0:
	    update_precond_identity();
	    break;
	case 1:
	    update_precond_modified_Cholesky();
	    break;
	default:
	    update_precond_identity();   
	}	
	return;
    }

    void update_precond_identity() {
	// Do nothing.  Preconditioner will always be identity matrix
	return;
    }

    void update_precond_modified_Cholesky() {

	bool success = FALSE;
	double alpha, beta, bmin;
	
	VectorXd T(nvars);
	SparseMatrixXd BB(nvars, nvars);
	BB = Bk.template selfadjointView<Lower>();
	
	for (int j=0; j<BB.outerSize(); j++) {
	    T(j) = sqrt(BB.innerVector(j).dot(BB.innerVector(j)));
	}
	
	for (int j=0; j<BB.outerSize(); j++) {
	    for (typename SparseMatrixXd::InnerIterator it(BB, j); it; ++it) {
		BB.coeffRef(it.row(),j) *= 1/sqrt(T(it.row()) * T(j));
	    }
	}
	
	beta = sqrt(BB.cwiseAbs2().sum());
	bmin = BB.coeff(0,0);
	for (int j=0; j<nvars; j++) {
	    bmin = std::min(bmin,BB.coeff(j,j));
	}
	
	if (bmin > 0) alpha = 0; else alpha = beta/2;
	
	int ii = 0;
	do {
	    ii++;	    
	    PrecondLLt.factorize(BB); // try decomp of modified Hessian
	    if (PrecondLLt.info() == Eigen::Success) {
		success = TRUE;
	    } else {
		alpha = std::max(2*alpha, beta/2.) - alpha;
		for (int j = 0; j<nvars; j++) {
		    BB.coeffRef(j,j) += alpha;
		}
	    }	    
	}
	while (!success);
    }

    void init_precond() {
	switch (precond_ID) {
	case 0:
	    init_precond_identity();
	    break;
	    break;
	case 1:
	    init_precond_modified_Cholesky();
	    break;
	default:
	    init_precond_identity();
	}	
	return;
    }

    void init_precond_identity() {
	SparseMatrixXd tmp(nvars, nvars);
	tmp.reserve(nvars);
	
	for (int j=0; j<nvars; j++) {
	    tmp.startVec(j);
	    tmp.insertBack(j,j) = 1.;
	}
	tmp.makeCompressed();
	
	PrecondLLt.compute(tmp);
	
	return;
    }

    inline void init_precond_modified_Cholesky() {
	// this is a special case of a hessian preconditioner
	// in general, preconditioner must be lower triangular and positive definite
	PrecondLLt.analyzePattern(Bk); // set sparsity pattern same as Hessian    
	update_precond_modified_Cholesky();
    }

    void init_sparse_hessian() {
	Bk.resize(nvars, nvars);
	Bk.reserve(nnz);
	func.get_hessian(xk, Bk);
    }

    inline void update_sparse_hessian() {
	func.get_hessian(xk, Bk);
    }

    inline void update_hessian() {
	update_sparse_hessian(); // get a brand new Hessian
	Bk *= function_scale_factor;
    }

public:
    // Constructor
    Trust_CG_Sparse(TFunc & func_,
		    const MatrixBase<TP>& startX_,
		    const double & rad_,
		    const double & min_rad_,
		    const double & tol_,
		    const double & prec_,
		    const int & report_freq_,
		    const int  report_level_,
		    const int & report_precision_,
		    const int & maxit_,
		    const double & contract_factor_,
		    const double & expand_factor_,
		    const double & contract_threshold_,
		    const double & expand_threshold_rad_,
		    const double & expand_threshold_ap_,
		    const double & function_scale_factor_,
		    const int & precond_refresh_freq_,
		    const int & precond_ID_,
		    const int & trust_iter_) :
	BaseType(func_, startX_, rad_, min_rad_, tol_, prec_, report_freq_,
		 report_level_, report_precision_, maxit_, contract_factor_,
		 expand_factor_, contract_threshold_, expand_threshold_rad_,
		 expand_threshold_ap_, function_scale_factor_,
		 precond_refresh_freq_, precond_ID_, trust_iter_),
	nnz(func_.get_nnz())
    {
	init_sparse_hessian();
	Bk *= function_scale_factor; 
	init_precond();     
    }

    template<typename Tvec, typename Tout>
    MB_Status get_current_state(const MatrixBase<Tvec> & pars_,
				double & fval, 
				const MatrixBase<Tvec> & grad_,
				const SparseMatrixBase<Tout> & hess_,
				int & iterations,
				double & radius) {
	
	MatrixBase<Tvec> & pars = const_cast<MatrixBase<Tvec>& >(pars_);
	MatrixBase<Tvec> & grad = const_cast<MatrixBase<Tvec>& >(grad_);
	SparseMatrixBase<Tout> & hess = const_cast<SparseMatrixBase<Tout>&>(hess_);
	
	pars = xk;
	func.get_fdf(xk, fval, grad);
	func.get_hessian(xk, hess);
	
	iterations = iter;
	radius = rad;
	
	return status;
    }

}; // end class Trust_CG_Sparse

#endif


