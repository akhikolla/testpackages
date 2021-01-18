#include "gse-rocke.h"
using namespace Rcpp;
using namespace arma;
using namespace std;


/***************************************************/
/*               Function prototypes               */
/***************************************************/
vec rho2(vec x, vec gam);
double rho2p(double x, double gam);
double rho2pp(double x, double gam);
double solve_scales( vec maj, vec cc1, vec gam1, double tol, int miter, double bdp );
mat pmd_adj(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group);
double scales(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group, double tol, int miter, double bdp);
mat iterS( mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, double sk, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group, double* wgts_mem, double* wgtsp_mem, double* ximp_mem, int* error_code_mem);
mat gradient( mat x, mat sigma1, mat sigma2, mat mu1, mat mu2, mat omega, double s1, 
	umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group,
	double tol_scale, int miter_scale, double bdp);
int update( mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, double sk, umat miss_group_unique, 
	uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group, 
	double* wgts_mem, double* wgtsp_mem, double* ximp_mem);
SEXP GSE_Rocke(SEXP X, SEXP N, SEXP P, SEXP Mu0, SEXP S0, SEXP Tol, SEXP Maxiter, SEXP Tol_scale, SEXP Miter_scale, SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Tuning_const_group, SEXP Gamma_tune_group, SEXP Print_step, SEXP Bdp);


/***************************************************/
/*                Main GSE function                */
/***************************************************/
SEXP GSE_Rocke(SEXP X, SEXP N, SEXP P, SEXP Mu0, SEXP S0, SEXP Tol, SEXP Maxiter, SEXP Tol_scale, SEXP Miter_scale, SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Tuning_const_group, SEXP Gamma_tune_group, SEXP Print_step, SEXP Bdp)
{
	try{	
	// Input argument
	mat x = as<mat>(X);
	int n = as<int>(N);
	int p = as<int>(P);
	mat mu0 = as<mat>(Mu0);
	mat Sigma0 = as<mat>(S0);
	double tol = as<double>(Tol);
	int maxiter = as<int>(Maxiter);
	double tol_scale = as<double>(Tol_scale);
	int miter_scale = as<int>(Miter_scale);
	umat miss_group_unique = as<umat>(Miss_group_unique);
	uvec miss_group_counts = as<uvec>(Miss_group_counts);
	vec tuning_const_group = as<vec>(Tuning_const_group);
	vec gamma_tune_group = as<vec>(Gamma_tune_group);
	int print_step = as<int>(Print_step);
	double bdp = as<double>(Bdp);
	
	// individual weights 
	vec wgts(n);
	double* wgts_mem = wgts.memptr();
	vec wgtsp(n);
	double* wgtsp_mem = wgtsp.memptr();
	
	// predicted data 
	mat ximp(n,p);
	double* ximp_mem = ximp.memptr();
	
	// Error code
	int error_code = 0; // 0 = no error, 1 = non-positive definite est scatter
	int* error_code_mem;
	error_code_mem = &error_code;
	
	// Basic setup
	double stilde0 = scales(x, Sigma0, Sigma0, true, mu0, miss_group_unique, miss_group_counts, tuning_const_group, gamma_tune_group, tol_scale, miter_scale, bdp);
	double ep = 1.0;
	int iter = 0;

	// Proceed only if no error
	Sigma0 = stilde0 * Sigma0; 
	mat Omega = Sigma0;
	stilde0 = 1; 
	
	// Start iteration
	mat mu1(1,p);
	mat Sigma1(p,p);
	do
	{
		iter++;
		mat iter_res = iterS(x, Omega, Sigma0, false, mu0, stilde0, miss_group_unique, miss_group_counts, tuning_const_group, gamma_tune_group, wgts_mem, wgtsp_mem, ximp_mem, error_code_mem);
		mu1 = iter_res.row(0);
		iter_res.shed_row(0);
		mat Stilde1 = iter_res;

		// Ensure symmetric
		Stilde1 = symmatl(Stilde1);

		// Check singularity
		mat Stilde1_chol(p,p);
		bool error_code_chol = chol(Stilde1_chol, Stilde1);
		if( !error_code_chol ) error_code = 1; 
		
		if( error_code == 0){
			double s1 = scales(x,Stilde1,Stilde1,true,mu1,miss_group_unique,miss_group_counts,tuning_const_group, gamma_tune_group,tol_scale, miter_scale, bdp);
			Sigma1 = s1 * Stilde1;
			double stilde1 = scales(x,Omega,Stilde1,false,mu1,miss_group_unique,miss_group_counts,tuning_const_group, gamma_tune_group,tol_scale, miter_scale, bdp);
			ep = 1 - stilde1/stilde0;

			if( ep <= 0){
				// gradient step 
				mat grad_res = gradient(x, Sigma0, Sigma1, mu0, mu1, Omega, stilde0, 
					miss_group_unique, miss_group_counts, tuning_const_group, gamma_tune_group,
					tol_scale, miter_scale, bdp);
				double stilde2 = grad_res(0,0);
				grad_res.shed_row(0);
				mat mu2 = grad_res.row(0);
				grad_res.shed_row(0);
				mat Sigma2 = grad_res;
				double ep2 = 1 - stilde2/stilde0;
				if( stilde2 >= stilde1 ){
					Sigma0 = Sigma1;
					mu0 = mu1;
					stilde0 = stilde1;
				} else{
					Sigma0 = Sigma2;
					mu0 = mu2;
					stilde0 = stilde2;
					ep = ep2;
				}
			} else{
				// Updates
				Sigma0 = Sigma1;
				mu0 = mu1;
				stilde0 = stilde1;
			}
		}
		//if( print_step == 2 ) Rcout << "iter=" << iter << "; tol=" << ep << "; scale=" << stilde0 << std::endl;
	}
	while( ( ep > tol) && (error_code == 0) && (iter <= maxiter) );
	
	int ximp_update = update(x, Omega, Sigma0, false, mu0, stilde0, miss_group_unique, miss_group_counts, tuning_const_group, gamma_tune_group, wgts_mem, wgtsp_mem, ximp_mem);

	return List::create( Named("S")=Sigma0, Named("mu")=mu0, Named("stilde0")=stilde0, 
		Named("weights")=wgts, Named("weightsprm")=wgtsp,Named("ximp")=ximp,
		Named("iter")=iter, Named("ep")=ep, Named("error_code")=error_code);

	} catch( std::exception& ex ){
		forward_exception_to_r( ex );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	return wrap(NA_REAL);
}






/***************************************************/
/*               Iterative step                    */
/***************************************************/
mat iterS( mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, double sk, umat miss_group_unique, 
	uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group, 
	double* wgts_mem, double* wgtsp_mem, double* ximp_mem, int* error_code_mem)
{ 
	/*********************************/	
	/* Some declaration of variables */
	int n_counts = miss_group_unique.n_rows;
	int n = x.n_rows;
	unsigned int p = x.n_cols;
	int error_code = *error_code_mem; // 0 = no error, 1 = non-positive scale, 2 = non-positive definite est scatter
	
	try{
		// individual weights and weights prime (new July 28, 2014)
		vec wgts(wgts_mem, n, false, true);
		vec wgtsp(wgtsp_mem, n, false, true);
		mat ximp(ximp_mem, n, p, false, true);
		
		uvec pp_grp = sum(miss_group_unique, 1);  // vector of observed variables for each missingness group
		int rowid_start = 0;

		double sweight1 = 0.0;				// sum of weight 1, w
		double sweight2 = 0.0;				// sum of weight 2, w X w*

		mat mu1 = zeros<mat>(1,p);			// mu_(k+1)
		mat Sigma1 = zeros<mat>(p,p);			// Sigma^tilde_(k+1)

		// things help to cross check (can be deleted)
		//mat xpredlist(n, p);				// n X p matrix of best predicted values 

		// Start computing
		for(int i = 0; i < n_counts; i++){
			mat mu_nonmiss(1, pp_grp(i));
			mat sigma0_nonmiss( pp_grp(i), pp_grp(i) );
			mat sigmak_nonmiss( pp_grp(i), pp_grp(i) );
			mat sigmaku_nonmiss( p, pp_grp(i) );
			mat xi( miss_group_counts(i) , pp_grp(i) );
			mat xpred(1, p);
			mat Ck = sigmak;
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp_grp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma0_nonmiss(mm, nn) = sigma0(j,k);
								sigma0_nonmiss(nn, mm) = sigma0(k,j);
								sigmak_nonmiss(mm, nn) = sigmak(j,k);
								sigmak_nonmiss(nn, mm) = sigmak(k,j);
								nn++;
							}
						}
						xi.col(mm) = x( span(rowid_start, rowid_end ), j );
						sigmaku_nonmiss.col(mm) = sigmak.col(j);	
						mu_nonmiss(0, mm) = mu(0, j);
						mm++;
					}
				}
			} else{
				sigmak_nonmiss = sigmak;
				sigma0_nonmiss = sigma0;
				xi = x.rows(rowid_start, rowid_end);
				mu_nonmiss = mu;
				//xpredlist.rows(rowid_start, rowid_end) = x.rows(rowid_start, rowid_end);
			}

			// Check singularity
			mat sigmak_nonmiss_chol(p,p);
			bool error_code_chol = chol(sigmak_nonmiss_chol, sigmak_nonmiss);
			if( !error_code_chol ) error_code = 2; 

			mat A = ones<mat>( pp_grp(i), pp_grp(i));
			mat diagA = diagmat(A);
			mat sigmak_nonmiss_inv = solve( sigmak_nonmiss, diagA );
			mat betas = sigmaku_nonmiss * sigmak_nonmiss_inv;
			if( pp_grp(i) < p){
				Ck = sigmaku_nonmiss * sigmak_nonmiss_inv * trans(sigmaku_nonmiss);
			}
			Ck = sigmak - Ck;
			
			double ppi = (double) pp_grp(i);
			
			// Added May 10, 2012, if fixed pt, no need to calculate dee or de0, as they are the same and ratio is always 1
			double dee = 1, de0 = 1, dee_ratio = 1;
			if( !equalsig ){
				dee = det( sigmak_nonmiss );
				dee = pow(dee, 1/ppi );
				de0 = det( sigma0_nonmiss );
				de0 = pow(de0, 1/ppi );
				dee_ratio = dee / de0;
			}

			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m) - mu_nonmiss;
				double md = as_scalar(xii * sigmak_nonmiss_inv * trans(xii));
				double md2 = md*dee_ratio / (tuning_const_group(i) * sk );
				double weight1 = rho2p( md2, gamma_tune_group(i) ) * dee_ratio;
				double weight2 = weight1 * (md / ppi);

				wgts(rowid_start + m) = weight1;
				wgtsp(rowid_start + m) = rho2pp( md2, gamma_tune_group(i))*dee_ratio/tuning_const_group(i);
				
				sweight1 += weight1;
				sweight2 += weight2;
				if( pp_grp(i) < p ){
					xpred = trans(betas * trans(xii)) + mu;
				} else{
					xpred = x.row(rowid_start + m);	
				}
				for(int jj=0; jj < p; jj++) ximp(rowid_start + m, jj) = xpred(0, jj); // new July 28, 2014 for outputing the imputed data matrix
				mu1 += weight1 * xpred;
				xpred -= mu;
				Sigma1 += weight1 * (trans(xpred) * xpred) + weight2 * Ck;
			}
			rowid_start = rowid_start + miss_group_counts(i);
		}

		mu1 = mu1/sweight1;
		Sigma1 = Sigma1/sweight2;

		mat res = join_cols(mu1, Sigma1);
		
		return res; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	mat res(p, p+1);
	res.fill(NA_REAL);
	return res;
}


/*****************************************************************/
/* Gradient step                                                 */
/*****************************************************************/
mat gradient( mat x, mat sigma1, mat sigma2, mat mu1, mat mu2, mat omega, double s1, 
	umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group,
	double tol_scale, int miter_scale, double bdp)
{ 
	int n = x.n_rows;
	unsigned int p = x.n_cols;
	try{
	/*********************************/	
	/* Some declaration of variables */
	mat mumin = mu1;
	mat sigmamin = sigma1;
	double smin = s1;
	for(int i=1; i<=9; i++){
		double xi = (double) 0.1*i;
		mat mu3 = (1.0 - xi) * mu1 + xi * mu2;
		mat sigma3 = (1.0 - xi) * sigma1 + xi * sigma2;
		sigma3 = symmatl(sigma3);
		double stilde3 = scales(x,omega,sigma3,false,mu3,miss_group_unique,miss_group_counts,tuning_const_group, gamma_tune_group,tol_scale, miter_scale, bdp);
		double s3 = scales(x,sigma3,sigma3,true,mu3,miss_group_unique,miss_group_counts,tuning_const_group, gamma_tune_group,tol_scale, miter_scale, bdp);
		sigma3 = s3 * sigma3;
		if( stilde3 < smin){
			mumin = mu3;
			sigmamin = sigma3;
			smin = stilde3;
		}
	}
	rowvec smin_tmp(p);
	smin_tmp(0) = smin;
	mat res = join_cols(mumin, sigmamin);
	res = join_cols(smin_tmp, res);
	return res; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	mat res(p, p+2);
	res.fill(NA_REAL);
	return res;
}


int update( mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, double sk, umat miss_group_unique, 
	uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group, 
	double* wgts_mem, double* wgtsp_mem, double* ximp_mem)
{ 
	/*********************************/	
	/* Some declaration of variables */
	int n_counts = miss_group_unique.n_rows;
	int n = x.n_rows;
	unsigned int p = x.n_cols;
	
	try{
		// individual weights and weights prime (new July 28, 2014)
		vec wgts(wgts_mem, n, false, true);
		vec wgtsp(wgtsp_mem, n, false, true);
		mat ximp(ximp_mem, n, p, false, true);
		
		uvec pp_grp = sum(miss_group_unique, 1);  // vector of observed variables for each missingness group
		int rowid_start = 0;

		// Start computing
		for(int i = 0; i < n_counts; i++){
			mat mu_nonmiss(1, pp_grp(i));
			mat sigma0_nonmiss( pp_grp(i), pp_grp(i) );
			mat sigmak_nonmiss( pp_grp(i), pp_grp(i) );
			mat sigmaku_nonmiss( p, pp_grp(i) );
			mat xi( miss_group_counts(i) , pp_grp(i) );
			mat xpred(1, p);
			mat Ck = sigmak;
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp_grp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma0_nonmiss(mm, nn) = sigma0(j,k);
								sigma0_nonmiss(nn, mm) = sigma0(k,j);
								sigmak_nonmiss(mm, nn) = sigmak(j,k);
								sigmak_nonmiss(nn, mm) = sigmak(k,j);
								nn++;
							}
						}
						xi.col(mm) = x( span(rowid_start, rowid_end ), j );
						sigmaku_nonmiss.col(mm) = sigmak.col(j);	
						mu_nonmiss(0, mm) = mu(0, j);
						mm++;
					}
				}
			} else{
				sigmak_nonmiss = sigmak;
				sigma0_nonmiss = sigma0;
				xi = x.rows(rowid_start, rowid_end);
				mu_nonmiss = mu;
				//xpredlist.rows(rowid_start, rowid_end) = x.rows(rowid_start, rowid_end);
			}

			mat A = ones<mat>( pp_grp(i), pp_grp(i));
			mat diagA = diagmat(A);
			mat sigmak_nonmiss_inv = solve( sigmak_nonmiss, diagA );
			mat betas = sigmaku_nonmiss * sigmak_nonmiss_inv;
			if( pp_grp(i) < p){
				Ck = sigmaku_nonmiss * sigmak_nonmiss_inv * trans(sigmaku_nonmiss);
			}
			Ck = sigmak - Ck;

			double ppi = (double) pp_grp(i);
			
			// Added May 10, 2012, if fixed pt, no need to calculate dee or de0, as they are the same and ratio is always 1
			double dee = 1, de0 = 1, dee_ratio = 1;
			if( !equalsig ){
				dee = det( sigmak_nonmiss );
				dee = pow(dee, 1/ppi );
				de0 = det( sigma0_nonmiss );
				de0 = pow(de0, 1/ppi );
				dee_ratio = dee / de0;
			}
			
			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m) - mu_nonmiss;
				double md = as_scalar(xii * sigmak_nonmiss_inv * trans(xii));
				double md2 = md*dee_ratio / (tuning_const_group(i) * sk );
				double weight1 = rho2p( md2, gamma_tune_group(i) ) * dee_ratio;
				double weight2 = weight1 * (md / ppi);

				wgts(rowid_start + m) = weight1;
				wgtsp(rowid_start + m) = rho2pp( md2, gamma_tune_group(i))*dee_ratio/tuning_const_group(i);

				if( pp_grp(i) < p ){
					xpred = trans(betas * trans(xii)) + mu;
				} else{
					xpred = x.row(rowid_start + m);
				}
				for(int jj=0; jj < p; jj++) ximp(rowid_start + m, jj) = xpred(0, jj); // new July 28, 2014 for outputing the imputed data matrix
				xpred -= mu;
			}
			rowid_start = rowid_start + miss_group_counts(i);
		}
		return 0; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	return 0;
}

/*****************************************************************/
/* Compute scales or S* for given scaled md and tuning constants */
/*****************************************************************/
double scales(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group, double tol, int miter, double bdp)
{
	try{
		mat res = pmd_adj(x, sigma0, sigmak, equalsig, mu, miss_group_unique, miss_group_counts, tuning_const_group, gamma_tune_group);
		vec partialVec = res.col(0);
		vec cc1 = res.col(1);
		vec gam1 = res.col(2);
		double s0 = solve_scales( partialVec, cc1, gam1, tol, miter, bdp );
		return s0; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	return NA_REAL;
}



mat pmd_adj(mat x, mat sigma0, mat sigmak, bool equalsig, mat mu, umat miss_group_unique, uvec miss_group_counts, vec tuning_const_group, vec gamma_tune_group)
{
	int n_counts = miss_group_unique.n_rows;
	int n = x.n_rows;
	unsigned int p = x.n_cols;
	try{
		vec partialVec(n);
		vec cc1(n);
		vec gam1(n);
		uvec pp = sum(miss_group_unique, 1);
		int rowid_start = 0;

		for(int i = 0; i < n_counts; i++){
			mat mu_nonmiss(1, pp(i));
			mat sigma0_nonmiss( pp(i), pp(i) );
			mat sigmak_nonmiss( pp(i), pp(i) );
			mat xi( miss_group_counts(i) , pp(i) );
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma0_nonmiss(mm, nn) = sigma0(j,k);
								sigma0_nonmiss(nn, mm) = sigma0(k,j);
								sigmak_nonmiss(mm, nn) = sigmak(j,k);
								sigmak_nonmiss(nn, mm) = sigmak(k,j);
								nn++;
							}
						}
						xi.col(mm) = x( span(rowid_start, rowid_end ), j );
						mu_nonmiss(0, mm) = mu(0, j);
						mm++;
					}
				}
			} else{
				sigmak_nonmiss = sigmak;
				sigma0_nonmiss = sigma0;
				xi = x.rows(rowid_start, rowid_end);
				mu_nonmiss = mu;
			}

			mat A = ones<mat>( pp(i), pp(i));
			mat diagA = diagmat(A);
			mat sigmak_nonmiss_inv = solve( sigmak_nonmiss, diagA );

			double ppi = (double) pp(i);
			double aabb_ratio = 1;
			if( !equalsig ){
				double aa = det( sigmak_nonmiss );
				aa = pow( aa, 1/ppi );
				double bb = det( sigma0_nonmiss );
				bb = pow( bb, 1/ppi );
				aabb_ratio = aa / bb;
			}
			
			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m) - mu_nonmiss;
				double md = as_scalar(xii * sigmak_nonmiss_inv * trans(xii));
				partialVec(rowid_start + m) = md * aabb_ratio / tuning_const_group(i);
				cc1(rowid_start + m) = tuning_const_group(i);
				gam1(rowid_start + m) = gamma_tune_group(i);
				}
			rowid_start = rowid_start + miss_group_counts(i);
		}	
		
		mat res( n, 3);
		res.col(0) = partialVec;
		res.col(1) = cc1;
		res.col(2) = gam1;
		return res; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	mat res( n, 2);
	res.fill(NA_REAL);
	return res;
}


// May 17, 2013: Add stopping criteria for solving the scale
double solve_scales( vec maj, vec cc1, vec gam1, double tol, int miter, double bdp )
{
	int n = maj.n_elem;
	uvec maj_ind = sort_index( maj );
	vec maj_sort(n);
	vec cc1_sort(n);
	vec cc1_sort_cum(n);
	double scc1_half = sum(cc1) * bdp;
	double scc1_tmp = 0;
	int mid_ind = 0;
	bool found = false;
	for(int i = 0; i < n; i++){
		maj_sort(i) = maj( maj_ind(i) );
		cc1_sort(i) = cc1( maj_ind(i) );
		scc1_tmp += cc1_sort(i);
		cc1_sort_cum(i) = scc1_tmp;
		if( cc1_sort_cum(i) >= scc1_half  && !found){
			mid_ind = i;
			found = true;
		}
	}
	
	double s0 = maj_sort(mid_ind);
	double eps = 1.0;
	int iter= 0;
	while( eps > tol){
		vec zz;
		zz = rho2(maj_sort/s0, gam1);
		vec zzcc1 = cc1_sort % zz;
		double s1 = s0 * sum(zzcc1) / scc1_half;
		eps = fabs(s1 - s0)/s0;
		s0 = s1;
		iter++;
		if( iter > miter ) break;
	}
	return s0;
}

/***************************************************/
/*                        Rho                      */
/***************************************************/
// Rho 2 = modified Rocke's biflat function (Maronna et al, 2006)
vec rho2(vec x, vec gam){
	x = min(x, 1 + gam);
	x = max(x, 1 - gam);
	vec y = (x - 1)/gam;
	vec z = 0.25 * (y % (3 - y % y)) + 0.5;
	return z; 
}

double rho2p(double x, double gam){
	double y = 0.0;
	if( x >= 1 - gam && x <= 1 + gam){
		double z = (x - 1)/gam;
		y = 0.75/gam * (1 - z*z);
	}
	return y; 
}

double rho2pp(double x, double gam){
	double y = 0.0;
	if( x >= 1 - gam && x <= 1 + gam){
		y = -1.5 * (x - 1)/ pow(gam, 3);
	}
	return y; 
}
