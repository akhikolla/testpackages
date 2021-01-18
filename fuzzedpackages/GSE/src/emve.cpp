#include <math.h>

#include "R.h"
#include "emve.h"
#include "cov-em.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

/***************************************************/
/*               Function prototypes               */
/***************************************************/
// Export to R
SEXP emve_Rcpp(SEXP X, SEXP X_nonmiss, SEXP Pu, SEXP N, SEXP P, SEXP Theta0, SEXP GG, SEXP D, SEXP X_miss_group_match,
    SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Miss_group_obs_col, SEXP Miss_group_mis_col,
    SEXP Miss_group_p, SEXP Miss_group_n, SEXP NResample, SEXP NSubsampleSize, SEXP Subsample_ID, SEXP CC, SEXP CK, SEXP Maxits);
SEXP fast_partial_mahalanobis(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts);

// Internal C++ functions
cube emve_resamp(mat x, umat x_nonmiss, vec pu, int n, int p, vec theta0, mat G, int d,
    uvec x_miss_group_match, umat miss_group_unique, uvec miss_group_counts,
    umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, int miss_group_n,
    int nResample, int nSubsampleSize, uvec subsample_id, vec cc, vec ck, int EM_maxits);
double emve_scale_constraint( mat sigma, umat miss_group_unique, uvec miss_group_counts );//ok
double emve_scale(vec d, vec cc, vec ck); //ok
vec fast_pmd(mat x_mu_diff, mat sigma, umat miss_group_unique, uvec miss_group_counts); //ok
void subsampling(double* subsample_mem, unsigned int* subsamp_nonmis_mem,
    mat x, umat x_nonmiss, int nSubsampleSize, int p, uvec subsample_id);
mat concentrate_step(mat x, umat x_nonmiss, vec pu, vec pmd, uvec x_miss_group_match,
    umat miss_group_unique, uvec miss_group_counts,
    umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, int miss_group_n,
    int n, int n_half, int p, vec theta0, mat G, int d, int EM_maxits,
    double* xh_mem, unsigned  int* misgrpuqh_mem, unsigned  int* msgrpcth_mem,
    unsigned int* msgrpoch_mem, unsigned int* msgrpmch_mem, unsigned int* msgrpph_mem); //ok


// Compute EM MLE
mat CovEM(mat x, int n, int p, vec theta0, mat G, int d, umat miss_group_unique, uvec miss_group_counts,
    umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, int miss_group_n,
    double tol, int maxiter); // OKAY
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
/*                Called functions                 */
/***************************************************/
SEXP emve_Rcpp(SEXP X, SEXP X_nonmiss, SEXP Pu, SEXP N, SEXP P, SEXP Theta0, SEXP GG, SEXP D, SEXP X_miss_group_match,
    SEXP Miss_group_unique, SEXP Miss_group_counts, SEXP Miss_group_obs_col, SEXP Miss_group_mis_col,
    SEXP Miss_group_p, SEXP Miss_group_n, SEXP NResample, SEXP NSubsampleSize, SEXP Subsample_ID, SEXP CC, SEXP CK, SEXP Maxits)
{
    try{
        mat x = as<mat>(X);
        umat x_nonmiss = as<umat>(X_nonmiss);
        vec pu = as<vec>(Pu);
        int n = as<int>(N);
        int p = as<int>(P);

        vec theta0 = as<vec>(Theta0);
        mat G = as<mat>(GG);
        int d = as<int>(D);
        uvec x_miss_group_match = as<uvec>(X_miss_group_match);
        umat miss_group_unique = as<umat>(Miss_group_unique);
        uvec miss_group_counts = as<uvec>(Miss_group_counts);
        umat miss_group_obs_col = as<umat>(Miss_group_obs_col);
        umat miss_group_mis_col = as<umat>(Miss_group_mis_col);
        uvec miss_group_p = as<uvec>(Miss_group_p);
        int miss_group_n = as<int>(Miss_group_n);

        int nResample = as<int>(NResample);
        int nSubsampleSize = as<int>(NSubsampleSize);
        uvec subsample_id = as<uvec>(Subsample_ID);

        vec cc = as<vec>(CC);
        vec ck = as<vec>(CK);
        int EM_maxits = as<int>(Maxits);

        // resampling an obtain a list of candidates to be further condensed
        cube cand_list = emve_resamp(x, x_nonmiss, pu, n, p, theta0, G, d,
                    x_miss_group_match, miss_group_unique, miss_group_counts,
                    miss_group_obs_col, miss_group_mis_col, miss_group_p, miss_group_n,
                    nResample, nSubsampleSize, subsample_id, cc, ck, EM_maxits);

        return wrap(cand_list);
    } catch( std::exception& __ex__ ){
        forward_exception_to_r( __ex__ );
    } catch(...){
        ::Rf_error( "c++ exception " "(unknown reason)" );
    }
    return wrap(NA_REAL);
}

SEXP fast_partial_mahalanobis(SEXP X_mu_diff, SEXP Sigma, SEXP Miss_group_unique, SEXP Miss_group_counts)
{
    try{
        mat x_mu_diff = as<mat>(X_mu_diff);
        mat sigma = as<mat>(Sigma);
        umat miss_group_unique = as<umat>(Miss_group_unique);
        uvec miss_group_counts = as<uvec>(Miss_group_counts);
        vec partialVec = fast_pmd(x_mu_diff, sigma, miss_group_unique, miss_group_counts);
        return wrap(partialVec);
    } catch( std::exception& __ex__ ){
        forward_exception_to_r( __ex__ );
    } catch(...){
        ::Rf_error( "c++ exception " "(unknown reason)" );
    }
    return wrap(NA_REAL);
}




// Function to perform initial resampling and condensing
// Input:
// x:			data matrix with missing values filled with coordinate median
// x_nonmiss: 	a matrix of non-missing indicator with 1 means non-missing
// n, p:			dimension of the data
// theta0, G, d:	parameters for CovEM calculations
// miss_group_...: information about the missing patterns
// nResample:	number of subsamples
// nSubsampleSize:	size for each subsample
// cc, kk: 	tuning constants when computing the MVE scale
// EM_maxits: max num iteration for EM
cube emve_resamp(mat x, umat x_nonmiss, vec pu, int n, int p, vec theta0, mat G, int d,
    uvec x_miss_group_match, umat miss_group_unique, uvec miss_group_counts,
    umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, int miss_group_n,
    int nResample, int nSubsampleSize, uvec subsample_id, vec cc, vec ck, int EM_maxits)
{
    int nCand = 1;
    try{
        // an array of intermediate results for each good conditioned subsample
        // for each nCand, it contains a matrix of p+2 X p
        // first row: (mve scale, 0, ..., 0) of p X 1
        // second row: (cand location) of p X 1
        // remaining: (cand scatter) of p X p
        cube cand_list(p + 2, p, nCand);
        mat cand_list_pmd(n, nCand);
        rowvec cand_mu(p);
        mat cand_S(p,p);
        double cand_sc;
        double sc_constr;
        vec cand_scale(nCand); // vector of scale for each candidate location and scatter
        cand_scale.fill(1e+20);
        double max_cand_scale;
        double min_cand_scale;
        uword curMaxScaleInd = 0; // index of the candidate with largest scale
        uword curMinScaleInd = 0;
        mat x_mu_diff(n,p);
        vec pmd(n);

        // for subsampling
        mat subsample(nSubsampleSize, p); double* subsample_mem = subsample.memptr();
        umat subsample_nonmiss(nSubsampleSize, p); unsigned int* subsamp_nonmis_mem = subsample_nonmiss.memptr();
        int nsubsamp = 0;
        uword rk = 0;

        // for concentration steps
        int n_half = (int) n/2;
        mat x_half( n_half, p);
        double* xh_mem = x_half.memptr();
        umat miss_grp_uq_half( miss_group_n, p); unsigned int* misgrpuqh_mem = miss_grp_uq_half.memptr();
        uvec miss_grp_cts_half( miss_group_n); unsigned int* msgrpcth_mem = miss_grp_cts_half.memptr();
        umat miss_grp_oc_half(miss_group_n, p); unsigned int* msgrpoch_mem = miss_grp_oc_half.memptr();
        umat miss_grp_mc_half(miss_group_n, p); unsigned int* msgrpmch_mem = miss_grp_mc_half.memptr();
        uvec miss_grp_p(n_half); unsigned int* msgrpph_mem = miss_grp_p.memptr();

        // Start resampling
        GetRNGstate(); /* set the seed from R? */
        for(int i=0; i < nResample; i++)
        {
            bool keep_subsamp = true;
            rk = 0;
            // keep resample until obtain a good subsample
            while( keep_subsamp)
            {
                subsampling(subsample_mem, subsamp_nonmis_mem, x, x_nonmiss, nSubsampleSize, p, subsample_id);
                urowvec b = sum( subsample_nonmiss );
                if( b.min() > 0 ){
                    cand_S = cov(subsample);
                    keep_subsamp = false;
                }
            }
            rk = arma::rank(cand_S);
            if( rk == p){
                // subsamples location (coordinate median)
                cand_mu = mean(subsample);
                // rescale the subsamples scatter matrix based on the subsample (not the whole sample)
                sc_constr = emve_scale_constraint(cand_S, miss_group_unique, miss_group_counts );
                cand_S = cand_S * sc_constr;
                // compute EMVE scale
                x_mu_diff = x;
                x_mu_diff.each_row() -= cand_mu;
                pmd = fast_pmd(x_mu_diff, cand_S, miss_group_unique, miss_group_counts);
                cand_sc = emve_scale(pmd, cc, ck);

                // concentration step: select 50% of data points with smallest adj pmd
                pmd = pmd/cand_sc;
                mat cand_res_concentrate = concentrate_step( x, x_nonmiss, pu, pmd,
                    x_miss_group_match, miss_group_unique, miss_group_counts,
                    miss_group_obs_col, miss_group_mis_col, miss_group_p, miss_group_n,
                    n, n_half, p, theta0, G, d, EM_maxits,
                    xh_mem, misgrpuqh_mem, msgrpcth_mem,
                    msgrpoch_mem, msgrpmch_mem, msgrpph_mem);
                cand_mu = cand_res_concentrate.row(0);
                cand_S =  cand_res_concentrate.rows(1,p);
                rk = arma::rank(cand_S);
                if( rk == p){
                // rescale the subsamples scatter matrix based on the subsample (not the whole sample)
                sc_constr = emve_scale_constraint(cand_S, miss_group_unique, miss_group_counts );
                cand_S = cand_S * sc_constr;
                // compute EMVE scale
                x_mu_diff = x; x_mu_diff.each_row() -= cand_mu;
                pmd = fast_pmd(x_mu_diff, cand_S, miss_group_unique, miss_group_counts);
                cand_sc = emve_scale(pmd, cc, ck);

                // store the candidates if the scale is larger than the current maximum scale
                max_cand_scale = cand_scale.max(curMaxScaleInd);
                if(cand_sc < cand_scale(curMaxScaleInd) ) {
                    cand_list.slice(curMaxScaleInd)(0,0) = cand_sc;
                    cand_list.slice(curMaxScaleInd).row(1) = cand_mu;
                    cand_list.slice(curMaxScaleInd).rows(2,p+1) = cand_S;
                    cand_scale(curMaxScaleInd) = cand_sc;
                    cand_list_pmd.col(curMaxScaleInd) = pmd;
                    if(nsubsamp < nCand) nsubsamp++;
                }
                }
            }
        }

        /*****************************************************************/
        // candidate Gaussian MLE based on the entire sample
        // mat cand_res_all = CovEM(x, n, p, theta0, G, d, miss_group_unique, miss_group_counts,
            // miss_group_obs_col, miss_group_mis_col, miss_group_p, miss_group_n, 0.0001, EM_maxits);
        // cand_res_all.shed_rows(0,1);
        // cand_mu = cand_res_all.row(0);
        // cand_S =  cand_res_all.rows(1,p);

        // rk = arma::rank(cand_S);
        // if( rk == p){
        // // rescale the subsamples scatter matrix based on the subsample (not the whole sample)
        // sc_constr = emve_scale_constraint(cand_S, miss_group_unique, miss_group_counts );
        // cand_S = cand_S * sc_constr;

        // // compute EMVE scale
        // x_mu_diff = x; x_mu_diff.each_row() -= cand_mu;
        // pmd = fast_pmd(x_mu_diff, cand_S, miss_group_unique, miss_group_counts);
        // cand_sc = emve_scale(pmd, cc, ck);
        // max_cand_scale = cand_scale.max(curMaxScaleInd);
        // if(cand_sc < cand_scale(curMaxScaleInd) ) {
            // cand_list.slice(curMaxScaleInd)(0,0) = cand_sc;
            // cand_list.slice(curMaxScaleInd).row(1) = cand_mu;
            // cand_list.slice(curMaxScaleInd).rows(2,p+1) = cand_S;
            // cand_scale(curMaxScaleInd) = cand_sc;
            // cand_list_pmd.col(curMaxScaleInd) = pmd;
            // if(nsubsamp < nCand) nsubsamp++;
        // }
        // }

        // final concentration step
        // for(int i = 0; i < nsubsamp; i++){
            // cand_sc = cand_list.slice(i)(0,0);
            // cand_mu = cand_list.slice(i).row(1);
            // cand_S = cand_list.slice(i).rows(2,p+1);
            // pmd = cand_list_pmd.col(i);

            // // select 50% of data points with smallest adj pmd
            // pmd = pmd/cand_sc;
            // mat cand_res_concentrate = concentrate_step( x, x_nonmiss, pu, pmd_mem,
                // x_miss_group_match, miss_group_unique, miss_group_counts,
                // miss_group_obs_col, miss_group_mis_col, miss_group_p, miss_group_n,
                // n, n_half, p, theta0, G, d, EM_maxits,
                // xh_mem, misgrpuqh_mem, msgrpcth_mem,
                // msgrpoch_mem, msgrpmch_mem, msgrpph_mem);
            // cand_mu = cand_res_concentrate.row(0);
            // cand_S =  cand_res_concentrate.rows(1,p);

            // rk = arma::rank(cand_S);
            // if( rk == p){
            // // rescale the subsamples scatter matrix based on the subsample
            // sc_constr = emve_scale_constraint(cand_S, miss_group_unique, miss_group_counts );
            // cand_S = cand_S * sc_constr;

            // // compute EMVE scale
            // x_mu_diff = x; x_mu_diff.each_row() -= cand_mu;
            // pmd = fast_pmd(x_mu_diff, cand_S, miss_group_unique, miss_group_counts);
            // cand_sc = emve_scale(pmd, cc, ck);
            // // store the candidates if the scale is larger than the current maximum scale
            // if( cand_sc < cand_list.slice(i)(0,0) ){
                // cand_list.slice(i)(0,0) = cand_sc;
                // cand_list.slice(i).row(1) = cand_mu;
                // cand_list.slice(i).rows(2,p+1) = cand_S;
                // cand_scale(i) = cand_sc;
            // }
            // }
        // }

        // Export only the min scale candidate
        cube cand_list_final(p+2,p,1);
        min_cand_scale = cand_scale.min(curMinScaleInd);
        cand_sc = cand_list.slice(curMinScaleInd)(0,0);
        cand_S = cand_list.slice(curMinScaleInd).rows(2,p+1);
        cand_list.slice(curMinScaleInd).rows(2,p+1) = cand_S * cand_sc;
        cand_list_final.slice(0) = cand_list.slice(curMinScaleInd);

        return cand_list_final;
    } catch( std::exception& __ex__ ){
        forward_exception_to_r( __ex__ );
    } catch(...){
        ::Rf_error( "c++ exception " "(unknown reason)" );
    }
    cube cand_list_final(p+2,p,nCand);
    cand_list_final.fill(NA_REAL);
    return cand_list_final;
}


mat concentrate_step(mat x, umat x_nonmiss, vec pu, vec pmd, uvec x_miss_group_match,
    umat miss_group_unique, uvec miss_group_counts,
    umat miss_group_obs_col, umat miss_group_mis_col, uvec miss_group_p, int miss_group_n,
    int n, int n_half, int p, vec theta0, mat G, int d, int EM_maxits,
    double* xh_mem, unsigned  int* misgrpuqh_mem, unsigned  int* msgrpcth_mem,
    unsigned int* msgrpoch_mem, unsigned int* msgrpmch_mem, unsigned int* msgrpph_mem)
{
    // information required for subsampling
    mat x_half(xh_mem, n_half, p, false, true);
    umat miss_grp_uq_half( misgrpuqh_mem, miss_group_n, p, false, true); miss_grp_uq_half.zeros();
    uvec miss_grp_cts_half( msgrpcth_mem, miss_group_n, false, true); miss_grp_cts_half.zeros();
    umat miss_grp_oc_half( msgrpoch_mem, miss_group_n, p, false, true); miss_grp_oc_half.zeros();
    umat miss_grp_mc_half( msgrpmch_mem, miss_group_n, p, false, true); miss_grp_mc_half.zeros();
    uvec miss_grp_p_half( msgrpph_mem, n_half, false, true);miss_grp_p_half.zeros();
    int miss_grp_n_half;

    // compute adjusted pmd
    vec pmd_adj(n);
  	for(int i = 0; i < n; i++) pmd_adj(i) = R::pchisq(pmd(i), (double) pu(i), 1, 0);
	  double pmd_adj_med = median(pmd_adj);
	  uvec x_half_ind = find( pmd_adj <= pmd_adj_med, n_half);

    // find new missing pattern based on the halved samples
    int miss_new_add = -1;
    for(int i = 0; i < n_half; i++){
        int index = (int) x_half_ind(i);
        x_half.row(i) = x.row( index );

        // check whether we found a new missing pattern
        urowvec check_new_miss(p);
        if( i == 0 ){ check_new_miss.zeros();
        } else{  check_new_miss = (x_nonmiss.row( index ) == miss_grp_uq_half.row( miss_new_add )); }

        if( sum(check_new_miss)== (unsigned int) p){
            // if found same missing pattern
            miss_grp_cts_half( miss_new_add )++;
        }else{
            // if found a new missing pattern
            miss_grp_uq_half.row( miss_new_add + 1) = x_nonmiss.row( index );
            miss_grp_cts_half( miss_new_add + 1)++;
            miss_grp_oc_half.row( miss_new_add + 1) = miss_group_obs_col.row( x_miss_group_match( index)-1 );
            miss_grp_mc_half.row( miss_new_add + 1) = miss_group_mis_col.row( x_miss_group_match( index)-1 );
            miss_grp_p_half( miss_new_add + 1) = miss_group_p( x_miss_group_match( index)-1);
            miss_new_add++;
        }

    }
    miss_grp_n_half = miss_new_add + 1;

    mat res = CovEM(x_half, n_half, p, theta0, G, d, miss_grp_uq_half, miss_grp_cts_half,
        miss_grp_oc_half, miss_grp_mc_half, miss_grp_p_half, miss_grp_n_half,
        0.0001, EM_maxits);
    res.shed_rows(0,1);

    return res;
}







// Function to compute partial mahalanobis distances
// Input:
// x_mu_diff:	centralized data matrix with missing values filled with coordinate median
// sigma: 	est scatter
// miss_group_unique, miss_group_counts: matrix and vec indicating the missing patterns
vec fast_pmd(mat x_mu_diff, mat sigma, umat miss_group_unique, uvec miss_group_counts)
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

            for(unsigned int m = 0; m < miss_group_counts(i); m++){
                mat xii = xi.row(m);
                partialVec(rowid_start + m) = as_scalar(xii * sigma_nonmiss_inv * trans(xii));
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


double emve_scale_constraint( mat sigma, umat miss_group_unique, uvec miss_group_counts )
{
    try{
        int n_counts = miss_group_unique.n_rows;
        mat sigma_nonmiss;
        double a = 0;
        double val, sign;
        unsigned int p = sigma.n_cols;
        uvec pp = sum( miss_group_unique, 1);
        for(int i = 0; i < n_counts; i++){
            mat sigma_nonmiss( pp(i), pp(i) );
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
                        mm++;
                    }
                }
            } else{
                sigma_nonmiss = sigma;
            }
            log_det(val, sign, sigma_nonmiss);
            a += val * miss_group_counts(i);
        }
        double ppp = as_scalar(sum(pp % miss_group_counts));
        a = exp(-1 * a / ppp);
        return a;
    } catch( std::exception& __ex__ ){
        forward_exception_to_r( __ex__ );
    } catch(...){
        ::Rf_error( "c++ exception " "(unknown reason)" );
    }
    return NA_REAL;
}

double emve_scale(vec d, vec cc, vec ck)
{
    int n = d.n_elem;
    vec ee = d/cc;
    uvec ff = sort_index(ee);
    double cct = sum(ck);
    ck = ck/cct;
    vec ee_sort(n);
    vec ckcumsum(n);
    ee_sort(0) = ee(ff(0));
    ckcumsum(0) = ck(ff(0));
    int mm = 0;
    for(int i = 1; i < n; i++){
        ee_sort(i) = ee(ff(i));
        ckcumsum(i) = ckcumsum(i-1) + ck(ff(i));
        if(ckcumsum(i) >= 0.5 && mm == 0) mm = i;
    }
    double me = 0.0;
    if(mm == 0){
        me = ee_sort(mm);
    } else if(ckcumsum(mm)==0.5){
        me = ee_sort(mm);
    } else{
        double aa = 1.0/(0.5-ckcumsum(mm-1));
        double bb = 1.0/(ckcumsum(mm)-0.5);
        double dd = aa+bb;
        aa = aa/dd;
        bb = bb/dd;
        me = aa*ee_sort(mm-1)+bb*ee_sort(mm);
    }
    return me;
}


/********************************************************************/
// Misc functions
/********************************************************************/
void subsampling(double* subsample_mem, unsigned int* subsamp_nonmis_mem,
    mat x, umat x_nonmiss, int nSubsampleSize, int p, uvec subsample_id)
{
    mat subsamp( subsample_mem, nSubsampleSize, p, false, true);
    umat subsamp_nonmiss( subsamp_nonmis_mem, nSubsampleSize, p, false, true);
    uvec subsample_id_shuff = shuffle(subsample_id);
    for(int i = 0; i < nSubsampleSize; i++)
    {
        subsamp.row(i) = x.row( subsample_id_shuff(i) );
        subsamp_nonmiss.row(i) = x_nonmiss.row( subsample_id_shuff(i) );
    }
}
