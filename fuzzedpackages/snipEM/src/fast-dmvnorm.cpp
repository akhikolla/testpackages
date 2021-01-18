#include "fast-dmvnorm.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

/***************************************************/
/*                  Some constants                 */
/***************************************************/
#define log2PI 1.837877

/***************************************************/
/*               Function prototypes               */
/***************************************************/
// Export to R
SEXP fast_mvnorm_density(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts);

// Function to compute partial mahalanobis distances
// Input:
// x_mu_diff:	centralized data matrix with missing values filled with coordinate median
// sigma: 	est scatter
// miss_group_unique, miss_group_counts: matrix and vec indicating the missing patterns
vec fast_dmvnorm(mat x_mu_diff, mat sigma, umat miss_group_unique, uvec miss_group_counts)
{
	int n_counts = miss_group_unique.n_rows;
	int n = x_mu_diff.n_rows;
	unsigned int p = x_mu_diff.n_cols;

	vec partialVec(n);
	uvec pp = sum(miss_group_unique, 1);
	try{
		int rowid_start = 0;
		for(int i = 0; i < n_counts; i++){
			mat sigma_nonmiss( pp(i), pp(i) );
			mat xi( miss_group_counts(i) , pp(i) );
			int rowid_end = rowid_start + miss_group_counts(i) - 1;
			if( pp(i) < p ){
				int mm = 0;
				for(unsigned int j=0; j<p; j++){
					int nn=mm;
					if(miss_group_unique(i,j) == 1){
						for(unsigned int k=j; k<p; k++){
							if( miss_group_unique(i,k) == 1 ){
								sigma_nonmiss(mm, nn) = sigma(j,k);
								sigma_nonmiss(nn, mm) = sigma(k,j);
								nn++;
							}
						}
						xi.col(mm) = x_mu_diff( span(rowid_start, rowid_end ), j );
						mm++;
					}
				}
			} else{
				sigma_nonmiss = sigma;
				xi = x_mu_diff.rows(rowid_start, rowid_end);
			}

			mat A = ones<mat>( pp(i), pp(i));
			mat diagA = diagmat(A);
			mat sigma_nonmiss_inv = solve( sigma_nonmiss, diagA );
			
			double sigma_nonmiss_logdet, sigma_nonmiss_logdet_sign;
			log_det(sigma_nonmiss_logdet, sigma_nonmiss_logdet_sign, sigma_nonmiss);
			
			for(unsigned int m = 0; m < miss_group_counts(i); m++){
				mat xii = xi.row(m);
				partialVec(rowid_start + m) = -0.5*pp(i)*log2PI - 0.5*sigma_nonmiss_logdet - 0.5*as_scalar(xii * sigma_nonmiss_inv * trans(xii));
			}
			rowid_start = rowid_start + miss_group_counts(i);
		}	
		return partialVec; 
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	partialVec.fill(NA_REAL);
	return partialVec;
}

	
SEXP fast_mvnorm_density(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts)
{
	try{
		mat x_mu_diff = as<mat>(X_mu_diff);
		mat sigma = as<mat>(Sigma);
		umat miss_group_unique = as<umat>(Miss_group_unique);
		uvec miss_group_counts = as<uvec>(Miss_group_counts);
		vec partialVec = fast_dmvnorm(x_mu_diff, sigma, miss_group_unique, miss_group_counts);
		return wrap(partialVec);
	} catch( std::exception& __ex__ ){
		forward_exception_to_r( __ex__ );
	} catch(...){
		::Rf_error( "c++ exception " "(unknown reason)" );
	}
	return wrap(NA_REAL);
}
