#include "cov-em.h"
using namespace Rcpp;
using namespace arma;
using namespace std;


/***************************************************/
/*               Function prototypes               */
/***************************************************/
// Rcpp export functions
SEXP CovEM_Rcpp(SEXP X, SEXP N, SEXP P, SEXP Theta0, SEXP G, SEXP D, 
	SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Miss_group_obs_col, SEXP Miss_group_mis_col,
	SEXP Miss_group_p, SEXP Miss_group_n, SEXP Tol, SEXP Maxiter);

mat CovEM(mat x, int n, int p, vec theta0, mat G, int d, umat miss_group_unique, uvec miss_group_counts, 
	umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, int miss_group_n,
	double tol, int maxiter);
void sweep(double* theta_mem, int d, double* G_mem, int G_ncol, int k, int rev); // OKAY, NO BUG
void sweepobs(double* theta_mem, int d, double* G_mem, int G_ncol, int p, 
	umat miss_group_unique, int miss_group_i); // OKAY, NO BUG
void preEM( double* theta_mem, int d, double* G_mem, int G_ncol, 
	mat x, int n, int p, umat miss_group_unique, uvec miss_group_counts, 
	umat miss_group_obs_col, uvec miss_group_p, int miss_group_n); // OKAY, NO BUG
vec iterEM( double* theta_mem, double* tobs_mem, int d, double* G_mem, int G_ncol, 
	mat x, int n, int p, umat miss_grp_unique, uvec miss_grp_counts, 
	umat miss_grp_obs_col, umat miss_grp_mis_col, uvec miss_grp_p, int miss_grp_n);


/***************************************************/
/*                Main GSE function                */
/***************************************************/
SEXP CovEM_Rcpp(SEXP X, SEXP N, SEXP P, SEXP Theta0, SEXP GG, SEXP D, 
	SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Miss_group_obs_col, SEXP Miss_group_mis_col,
	SEXP Miss_group_p, SEXP Miss_group_n, SEXP Tol, SEXP Maxiter)
{
	try{	
	using namespace Rcpp;
	using namespace arma;
	mat x = as<mat>(X);
	int n = as<int>(N);
	int p = as<int>(P);
	vec theta0 = as<vec>(Theta0);
	mat G = as<mat>(GG);
	int d = as<int>(D);
	umat miss_group_unique = as<umat>(Miss_group_unique);
	uvec miss_group_counts = as<uvec>(Miss_group_counts);
	umat miss_group_obs_col = as<umat>(Miss_group_obs_col);
	umat miss_group_mis_col = as<umat>(Miss_group_mis_col);
	uvec miss_group_p = as<uvec>(Miss_group_p);
	int miss_group_n = as<int>(Miss_group_n);
	double tol = as<double>(Tol);
	int maxiter = as<int>(Maxiter);

	mat res = CovEM(x, n, p, theta0, G, d, miss_group_unique, miss_group_counts, 
		miss_group_obs_col, miss_group_mis_col, miss_group_p, miss_group_n,tol, maxiter);
	return wrap(res);
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	return wrap(NA_REAL);
}






/***************************************************/
/*               Iterative step                    */
/***************************************************/
mat CovEM(mat x, int n, int p, vec theta0, mat G, int d, umat miss_group_unique, uvec miss_group_counts, 
	umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, 
	int miss_group_n, double tol, int maxiter)
{
	// initial setup
	double* theta0_mem = theta0.memptr();
	double* G_mem = G.memptr();
	vec tobs( theta0_mem, d);
	double* tobs_mem = tobs.memptr();

	// preliminary calculation before EM iteration
	preEM( tobs_mem, d, G_mem, p+1, 
		x, n, p, miss_group_unique, miss_group_counts, 
		miss_group_obs_col, miss_group_p, miss_group_n); 

	// Iterative step EM
	int iter = 0; bool converged = false;
	double eps = 1.0;
	vec told( theta0_mem, d);
	vec told_tmp(theta0_mem, d);
	while( !converged){
		double* told_mem = told.memptr();
		vec tnew = iterEM( told_mem, tobs_mem, d, G_mem, p+1, x, n, p, 
			miss_group_unique, miss_group_counts, miss_group_obs_col, 
			miss_group_mis_col,miss_group_p, miss_group_n);
		vec tdiff = abs(tnew - told_tmp);
		eps = max(tdiff);
		if( eps <= tol) converged  = true;
		told = tnew;
		told_tmp = tnew;
		if( iter > maxiter) break;
		iter++;
	}

	mat res(p+3, p); res.zeros();
	res.at(0,0) = iter; res.at(1,0) = eps; 
	for(int i = 0; i < p; i++){
		res.at(2,i) = told.at(G.at(0,i+1)); // extract the est loc
		for(int j = i; j < p; j++){
			// extract the est cov
			res.at(3+i,j) = told.at(G.at(i+1,j+1));
			res.at(3+j,i) = res.at(3+i,j);
		}
	}
	return res;
}


// sweeping operator
// Input:
// d: dimension of theta (the parameters in G)
// theta_mem: pointer to theta, parameters
// G_mem: pointer to G matrix to be swept (G is only an index to theta)
// G_ncol: number of column in G
// k: int, the row and column to be swept
// rev: int, 1 means ordinary and -1 means reverse
void sweep(double* theta_mem, int d, double* G_mem, int G_ncol, int k, int rev){
	// reconstruct the matrix without copying
	vec theta(theta_mem, d, false, true);
	mat G(G_mem, G_ncol, G_ncol, false, true);

	// Step 1: replace the pivot cell
	double akk = theta(G(k,k));
	theta(G(k,k)) = -1/akk;
	// Step 2: replace the corresponding row and col
	for(int j=0; j < G_ncol; j++)
		if(j != k) theta(G(j,k)) = rev * theta(G(j,k))/akk;
			
	// Step 3: replace the remaining ones
	for(int j=0; j < G_ncol; j++){
		for(int l=j; l < G_ncol; l++){
			if( (j != k) & (l != k) ){
				double bb = theta(G(j,k));
				double cc = theta(G(l,k));
				theta(G(j,l))=theta(G(j,l)) - akk * bb * cc;
			}
		}
	}
}


// sweeping iteration based on observed entries
// Input 
// theta_mem, d, G_mem, G_ncol: (see above)
// miss_group_unique: matrix of missing pattern with 1=obs, 0=miss
// miss_group_i: index of current missing pattern = 0,...,nrow(miss_group_unique)-1
void sweepobs(double* theta_mem, int d, double* G_mem, int G_ncol, int p,
	umat miss_group_unique, int miss_group_i)
{
	// reconstruct the matrix without copying
	vec theta(theta_mem, d, false, true);
	mat G(G_mem, G_ncol, G_ncol, false, true);
	for(int i=0; i < p; i++){
		// if the variable is observed and the corresponding column
		// has not been swept (as indicated by the corresponding sign 
		// of diagonal entry), sweep it
		if( (miss_group_unique(miss_group_i, i) == 1) & (theta(G(i+1,i+1)) > 0.0)){
			sweep(theta_mem, d, G_mem, G_ncol, i+1, 1);
		} 
		// if the variable is missing and the corresponding column
		// has been swept, reverse it
		else if( (miss_group_unique(miss_group_i, i) == 0) & (theta(G(i+1,i+1)) < 0.0)){
			sweep(theta_mem, d, G_mem, G_ncol, i+1, -1);
		}
	}
}

// Compute sufficient statistics on only the observed entries
void preEM( double* theta_mem, int d, double* G_mem, int G_ncol, 
	mat x, int n, int p, umat miss_group_unique, uvec miss_group_counts, 
	umat miss_group_obs_col, uvec miss_group_p, int miss_group_n)
{
	// Setup
	vec theta(theta_mem, d, false, true);
	mat G(G_mem, G_ncol, G_ncol, false, true);
	// Initializing theta
	theta(0) = 1.0;
	for(int i = 1; i < d; i++) theta(i) = 0.0;
	// Start computing
	int rowid_st = 0;
	for(int i = 0; i < miss_group_n; i++){
		// Loop over each to sufficient statistics
		for(unsigned int m = 0; m < miss_group_counts(i); m++){
			for(unsigned int j = 0; j < miss_group_p(i); j++){
				theta(G(0,miss_group_obs_col(i,j))) += x(rowid_st + m, miss_group_obs_col(i,j)-1);
				for(unsigned int k= j; k < miss_group_p(i); k++)
					theta(G(miss_group_obs_col(i,j),miss_group_obs_col(i,k))) += x(rowid_st + m, miss_group_obs_col(i,j)-1)*x(rowid_st + m, miss_group_obs_col(i,k)-1);
			}
		}
		rowid_st = rowid_st + miss_group_counts(i);
	}
}




vec iterEM( double* theta_mem, double* tobs_mem, int d, double* G_mem, int G_ncol, 
	mat x, int n, int p, umat miss_grp_unique, uvec miss_grp_counts, 
	umat miss_grp_obs_col, umat miss_grp_mis_col, uvec miss_grp_p, int miss_grp_n)
{ 	
	// Input
	vec theta(theta_mem, d, false, true);
	vec tobs(tobs_mem, d, false, true);
	mat G(G_mem, G_ncol, G_ncol, false, true);
	vec cmean(p);

	// output
	vec tnew(tobs_mem, d);
	
	try{		
		// Start computing
		int rowid_st = 0;
		for(int i = 0; i < miss_grp_n; i++){
			sweepobs(theta_mem, d, G_mem, G_ncol, p, miss_grp_unique, i);

			// Loop over each to compute conditional mean and covariances
			int noc = miss_grp_p(i); 	// number of obs
			int nmc = p - noc;		// number of miss

			for(unsigned int m = 0; m < miss_grp_counts(i); m++){
				// compute cmean = E(xmis_i|xobs_i,theta)
				for(int j=0; j < nmc; j++){
					cmean( miss_grp_mis_col(i,j)-1) = theta(G(0,miss_grp_mis_col(i,j)));
					for(int k=0; k < noc; k++){
						cmean( miss_grp_mis_col(i,j)-1) += theta(G(miss_grp_obs_col(i,k),miss_grp_mis_col(i,j)))*x(rowid_st + m,miss_grp_obs_col(i,k)-1);
					}
				}
				//update location vector (corresponding to missing part)
				for(int j=0; j < nmc; j++){
					tnew(G(0,miss_grp_mis_col(i,j))) += cmean( miss_grp_mis_col(i,j)-1);
					// covariance of missing vs observed: cross-product
					for(int k=0; k < noc; k++){
						tnew(G(miss_grp_obs_col(i,k),miss_grp_mis_col(i,j))) += x(rowid_st+m, miss_grp_obs_col(i,k)-1)*cmean( miss_grp_mis_col(i,j)-1);
					}
					// covariance of missing: cross-product + Ci
					for(int k=j; k < nmc; k++){
						tnew(G(miss_grp_mis_col(i,k),miss_grp_mis_col(i,j))) += cmean( miss_grp_mis_col(i,k)-1)*cmean( miss_grp_mis_col(i,j)-1) + theta(G(miss_grp_mis_col(i,k),miss_grp_mis_col(i,j)));
					}
				}
			}
			rowid_st = rowid_st + miss_grp_counts(i);
		}
		for(int i = 1; i < d; i++) tnew(i) = tnew(i)/n;
		double* tnew_mem = tnew.memptr();
		sweep(tnew_mem, d, G_mem, G_ncol, 0, 1);
		return tnew;

	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	tnew.fill(NA_REAL);
	return tnew; 
}




