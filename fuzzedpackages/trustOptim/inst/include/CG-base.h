// CG-base.cpp.  Part of the trustOptim package for the R programming language.

// This file is part of trustOptim, a nonlinear optimization package
// for the R statistical programming platform.
//
// Copyright (C) 2013 Michael Braun
//
// This Source Code Form is subject to the license terms of the Mozilla
// Public License v. 2.0. If a copy of the MPL was not distributed
// with this file, you can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef __TRUST_OPTIM_CGBASE
#define __TRUST_OPTIM_CGBASE

#include <math.h>
#include <iostream>
#include <iomanip>
#include <exception>
#include <stdexcept>
#include <string>

#include "common_R.hpp"

#ifndef REPORT_HEADER_FREQ
#define REPORT_HEADER_FREQ 25
#endif

using Eigen::Matrix;
using Eigen::MatrixBase;
using Eigen::Dynamic;
using Eigen::VectorXi;
using Eigen::VectorXd;
using Eigen::MatrixXd;
using Eigen::LLT;
using Eigen::SimplicialLLT;
using Eigen::Upper;
using Eigen::Lower;
using Eigen::SparseMatrix;
using Eigen::SparseMatrixBase;

using std::endl;

template<typename TP, typename TFunc, typename THess, typename TPreLLt>  // TP is type for the parameter vector
class Trust_CG_Base {

public:

    // Constructor
    Trust_CG_Base(TFunc&, const MatrixBase<TP>&,	
		  const double &,
		  const double &,
		  const double &,
		  const double &,
		  const int &,
		  int,
		  const int &,
		  const int &,
		  const double &,
		  const double &,
		  const double &,
		  const double &,
		  const double &,
		  const double &,
		  const int &,
		  const int &,
		  const int &);


    int run();

protected:

    TFunc & func; // functor that contains function to be optimized
    const MatrixBase<TP>& startX;
    double rad; // trust region radius ("step.size" in R List)
    const double & min_rad;
    const double & tol;
    const double & prec;
    const int & report_freq;
    int  report_level;
    const int & report_precision;
    const int & maxit;
    const double & contract_factor;
    const double & expand_factor;
    const double & contract_threshold;
    const double & expand_threshold_rad;
    const double & expand_threshold_ap;
    const double & function_scale_factor; 
    const int & precond_refresh_freq;
    const int & precond_ID;
    int nvars;
    const int & trust_iter;


    THess Bk;
    TPreLLt PrecondLLt;

    int iter;

    double f; // current function value
    VectorXd xk; // current value of x
    VectorXd gk; // current gradient
    VectorXd sk; // x_{k} - x_{k-1}
    VectorXd yk; // g_{k} - g_{k-1}
    MB_Status status;

    VectorXd try_x; // x for proposed step
    VectorXd try_g; // gradient at proposed step

    double try_f, gs;
    double nrm_gk, ared, pred, sBs, ap, nrm_sk_scaled;
    int i, j;
    int header_freq, page_count, f_width, g_width, r_width;

    int num_CG_iters;
    std::string CG_stop_reason;


    // internal workspace for trust region computation

    VectorXd zj;
    VectorXd rj;
    VectorXd dj;
    VectorXd zj_old;
    VectorXd yj;
    VectorXd wd;
    VectorXd wz;

    virtual void update_hessian() =0;
    virtual void update_precond() =0;

    template<typename Tvec>
    void solve_trust_CG(const MatrixBase<Tvec>&);

    MB_Status update_one_step();

    double find_tau(const VectorXd& z,
		    const VectorXd& d)
    {
	UPz(PrecondLLt, d, wd);
	UPz(PrecondLLt, z, wz);
    
	double d2 = wd.squaredNorm();
	double z2 = wz.squaredNorm();
	double zd = wd.dot(wz);
    
	double root = zd*zd - d2*(z2 - (rad*rad));  
	double tau = (sqrt(root) - zd) / d2;
    
	return(tau);
    }

    inline void UPz(const LLT<THess> & X,
		    const VectorXd & v,
		    const VectorXd & out_) {

	VectorXd & out = const_cast<VectorXd&>(out_);
	out = X.matrixU() *  v;
    }
  
  
    inline void UPz(const SimplicialLLT<THess> & X,
		    const VectorXd& v,
		    const VectorXd& out_) {

	VectorXd & out = const_cast<VectorXd&>(out_);
	out = X.permutationP() * v;
	out =  X.matrixU().template triangularView<Upper>() * out;
	
    }
  
    inline double get_nrm_sk(const SimplicialLLT<THess>& X) {
	double res = (X.matrixU() * (X.permutationPinv() * sk).eval()).norm();
	return(res);
    }


    inline double get_nrm_sk(const LLT<THess>& X) {
	double res = (X.matrixU() * sk).norm();
	return(res);
    }

    void report_state(const int&);
    void report_header();

}; // end Trust_CG_Base


// Constructor
template<typename TP, typename TFunc, typename THess, typename TPreLLt>  
Trust_CG_Base<TP, TFunc, THess, TPreLLt>
::Trust_CG_Base(TFunc & func_,
		const MatrixBase<TP>& startX_,
		const double & rad_,
		const double & min_rad_,
		const double & tol_,
		const double & prec_,
		const int & report_freq_,
		int report_level_,
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
		const int & trust_iter_
		) :
    func(func_), startX(startX_),
    rad(rad_), min_rad(min_rad_),
    tol(tol_), prec(prec_),
    report_freq(report_freq_), report_level(report_level_),
    report_precision(report_precision_), maxit(maxit_),
    contract_factor(contract_factor_), expand_factor(expand_factor_),
    contract_threshold(contract_threshold_),
    expand_threshold_rad(expand_threshold_rad_),
    expand_threshold_ap(expand_threshold_ap_),
    function_scale_factor(function_scale_factor_), 
    precond_refresh_freq(precond_refresh_freq_),
    precond_ID(precond_ID_), nvars(startX_.size()),
    trust_iter(trust_iter_)
{

    if ( (function_scale_factor==0) || !my_finite(function_scale_factor) ) {
	throw MyException("Invalid function.scale.factor",
			  __FILE__,__LINE__);
    }
    
    xk = startX; // current x = initial x
    gk.resize(nvars);
    sk.resize(nvars);
    yk.resize(nvars);
    try_x.resize(nvars);
    try_g.resize(nvars);

    // Get initial values, gradients and Hessians
    func.get_fdf(xk, f, gk); // get initial function value and gradient

    if (!my_finite(f)) {
	throw MyException("Function value at starting point is not finite.",
			  __FILE__, __LINE__);
    }

    f *= function_scale_factor;
    gk *= function_scale_factor;
    nrm_gk = gk.norm();

    f_width = std::max(log10(std::abs(f)),1. ) + report_precision + 5;
    g_width = std::max(log10(std::abs(nrm_gk)),1. ) + report_precision + 5;
    // modified_gr
    r_width = std::max(log10(std::abs(rad)),1. ) + report_precision + 5;

    if (!my_finite(nrm_gk)) {
	throw MyException("Norm of gradient at starting point is not finite",
			  __FILE__,__LINE__); 
    }
    
    // allocate workspace for trust region iteration
    
    zj.setZero(nvars);
    rj.setZero(nvars);
    dj.setZero(nvars);
    zj_old.setZero(nvars);
    yj.setZero(nvars);
    wd.resize(nvars);
    wz.resize(nvars);

    header_freq = REPORT_HEADER_FREQ;
    page_count = header_freq;

} // end constructor



template<typename TP, typename TFunc, typename THess, typename TPreLLt>
template<typename Tvec>
void Trust_CG_Base<TP, TFunc, THess, TPreLLt>
::solve_trust_CG(const MatrixBase<Tvec>& pk_){

    // solves trust region subproblem using Steinhaug 1983, Section 2 algorithm
    // outputs:  proposed step direction p

 
    MatrixBase<Tvec>& pk = const_cast<MatrixBase<Tvec>&>(pk_);

    double norm_rj, dot_ry, dot_ry_old, norm_zj, aj, bj, tau, dBd, norm_gk;
    int j;
    double crit;

    zj.setZero();
    rj = -gk;

    UPz(PrecondLLt, rj, wd); // wd is workspace
    norm_rj = wd.norm();

    UPz(PrecondLLt, gk, wd); // wd is workspace
    norm_gk = wd.norm();

    // solving LL'y = r 
    yj = PrecondLLt.solve(rj);
    dj = yj;

    std::stringstream reason;  
    for (j=0; j<trust_iter; j++) {
	dBd = dj.dot(Bk.template selfadjointView<Lower>() * dj);
	if (dBd <= 0) {

	    // find tau that minimizes quadratic model, and return p
	    tau = find_tau(zj, dj);
	    pk = zj + tau*dj;
	    num_CG_iters = j+1;
	    reason << "Negative curvature";   
	    break;
	}
    
	aj = rj.dot(yj) / dBd;
	zj_old = zj;
	zj.noalias() += aj*dj; // now z_{j+1}

	UPz(PrecondLLt, zj, wd); // wd is workspace
	norm_zj = wd.norm();

	if (norm_zj >= rad) {
 
	    //find tau>=0 s.t. p intersects trust region
	    tau = find_tau(zj_old, dj);  
	    pk = zj_old + tau*dj;
	    num_CG_iters = j+1;
	    reason << "Intersect TR bound";   
	    break;
	}

	dot_ry = rj.dot(yj);

	rj.noalias() -= aj * (Bk.template selfadjointView<Lower>() * dj).eval();    

	UPz(PrecondLLt, rj, wd);
	norm_rj = wd.norm();
    
	crit = norm_rj / norm_gk;

	if (crit < tol) {
	    pk = zj;
	    num_CG_iters = j+1;
	    reason << "Reached tolerance";
	    break;
	}

	dot_ry_old = dot_ry;

	// updating yj
	yj = PrecondLLt.solve(rj);
	dot_ry = rj.dot(yj);
    	bj = dot_ry/dot_ry_old;
	dj *= bj;
	dj.noalias() += yj;
    } // end loop

    if (j>=trust_iter) {
	pk=zj;
	num_CG_iters = j;
	reason << "Exceeded max CG iters";
    }

    CG_stop_reason = reason.str();

    return;
  
} // end solve_trust_CG

template<typename TP, typename TFunc, typename THess, typename TPreLLt>
MB_Status Trust_CG_Base<TP, TFunc, THess, TPreLLt>::update_one_step() {

    using std::endl;
    // Solve Trust Region Subproblem. Output direction sk.
    // Keep even if not accepted (for SR1 updates).

    MB_Status step_status = UNKNOWN;
    solve_trust_CG(sk);

    nrm_sk_scaled = get_nrm_sk(PrecondLLt);

    if (!my_finite(nrm_sk_scaled)) {
	step_status = FAILEDCG;
    } else {
	
	// check proposed direction.  Flag if no movement.  
	
	try_x = xk+sk;
  	func.get_f(try_x, try_f); //  get f(x+s)

	if (my_finite(try_f)) {
	    try_f *= function_scale_factor;
	    ared = f - try_f;
	    gs = gk.dot(sk);
	    sBs = sk.dot(Bk.template selfadjointView<Lower>() * sk);    
	    pred = -(gs + sBs/2);    
	    if (pred < 0) step_status = ENEGMOVE;
	    ap = ared/pred;
	  	  
	} else {
	    step_status = FAILEDCG; // new function value is not finite
	}
    }

    if ( (step_status != FAILEDCG) && (step_status != ENEGMOVE) ) {
    	if (ap > contract_threshold) {
	    func.get_df(try_x, try_g); //  get grad(x+s)	  
	    if (my_finite(try_g.norm())){		
		try_g *= function_scale_factor;
		// used for quasi-newton.  Should eventually move to CG-quasi.
		yk = try_g - gk;		  
		f = try_f;
		xk += sk;
		gk = try_g;
		nrm_gk = gk.norm();
		
		if ((ap > expand_threshold_ap)
		     && (nrm_sk_scaled >= expand_threshold_rad*rad))
		    {  
			step_status = EXPAND;
		    } else {
		    step_status = MOVED;
		}
	    } else {
		step_status = FAILEDCG; // new gradient is not finite
	    }	
	} else { //Contracting trust region
	    step_status = CONTRACT;
	}
    }


    if ( (step_status == CONTRACT)
	 || (step_status == FAILEDCG)
	 || (step_status == ENEGMOVE) ) {
	rad *= contract_factor;
    }
  
    if (step_status == EXPAND) {
	rad *= expand_factor;
    }
  
    return step_status;

} // end update_one_step



template<typename TP, typename TFunc, typename THess, typename TPreLLt>
void Trust_CG_Base<TP, TFunc, THess, TPreLLt>::report_header() {

    using std::setw;
    using std::left;
    using std::right;
    using std::ios;
    using std::setiosflags;
    using std::setprecision;

    if (report_level >= 1) {  
	TRUST_COUT  <<  endl << setw(floor(log10(double(maxit)))+1) << right << "iter";
	TRUST_COUT  << setw(f_width) << right << "f  ";
    }
    if (report_level >=2) {
	TRUST_COUT  << right << setw(g_width) << right << "nrm_gr";
	TRUST_COUT  << setw(27) << right << "status";  
    }
    if (report_level >= 3) {
	TRUST_COUT  << setw(r_width) << right << "rad";
    }
    if (report_level >=4) {
	TRUST_COUT << setw(floor(log10(double(trust_iter)))+6) << right << "CG iter";
	TRUST_COUT << setw(27) << "CG result";
    }
    if (report_level >= 1) {
	TRUST_COUT << endl;
    }
}

template<typename TP, typename TFunc, typename THess, typename TPreLLt>
void Trust_CG_Base<TP, TFunc, THess, TPreLLt>::report_state(const int& iter) {

    using std::setw;
    using std::left;
    using std::right;
    using std::ios;
    using std::setiosflags;
    using std::setprecision;

    if ( page_count == header_freq ) {
	report_header();
	page_count = 0;
    }

    page_count++;

    if (report_level >= 1) {
	TRUST_COUT << setiosflags(ios::fixed) << setprecision(report_precision);
	TRUST_COUT  << setw(floor(log10(double(maxit)))+1) << right << iter;
	TRUST_COUT  << setw(f_width) << right << f;
    }
    if (report_level >=2) {
	TRUST_COUT  << setw(g_width) << right << nrm_gk;
	TRUST_COUT  << setw(27) << right << MB_strerror(status);
    }
    if (report_level >= 3) {
	TRUST_COUT  << setprecision(report_precision) << setw(r_width) << right << rad;
    }
    if (report_level >=4) {
	TRUST_COUT << setw(floor(log10(double(trust_iter)))+6) << right << num_CG_iters;
	TRUST_COUT << setw(27) << right << CG_stop_reason;
    }
    if (report_level >= 1) {
	TRUST_COUT << endl;
    }

}


template<typename TP, typename TFunc, typename THess, typename TPreLLt> 
int Trust_CG_Base<TP, TFunc, THess, TPreLLt>::run() {

    using std::endl;
    iter = 0;
    status = CONTINUE;
    TRUST_COUT << "Beginning optimization\n";
 
    do
	{
	    iter++;

	    if (check_interrupt()) {
		throw MyException("Interrupt detected",__FILE__,__LINE__);
	    }

	    status = update_one_step();

	    // report state
	    if ( (report_freq>0) && (iter % report_freq == 0)) {
		report_state(iter);
	    }

	    if ( (status == FAILEDCG) || (status == ENEGMOVE) ) {
		status = CONTINUE;
	    }

	    // check converge of algorithm
	    if (nrm_gk/sqrt(double(nvars)) <= prec) {
		status = SUCCESS; // successful convergence of algorithm
	    }

	    if (iter >= maxit) {
		status = EMAXITER;
	    }

	    if ( rad < min_rad) { // trust radius has collapsed
		status = ETOLG; // can't move anymore.  we're done
	    }

	    // update Hessian

	    if ( (status == MOVED) || (status == EXPAND) ) {
		update_hessian();
		if ( iter % precond_refresh_freq == 0 ) {  
		    update_precond();
		}	
		status = CONTINUE;
	    }      
	    if (status == CONTRACT) status = CONTINUE;      
	}
    while (status==CONTINUE);
  
    TRUST_COUT << "\nIteration has terminated\n";
    report_level = 2;
    report_state(iter);
    TRUST_COUT << endl;
    
    return status;
}


#endif
