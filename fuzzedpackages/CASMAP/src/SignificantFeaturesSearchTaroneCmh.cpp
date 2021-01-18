/*
 * SignificantFeaturesSearchTaroneCmh.cpp
 *
 *  Created on: 8 May 2017
 *      Author: mikolajr
 */

#include "SignificantFeaturesSearchTaroneCmh.h"

#include<math.h> //pow
#include<algorithm> //sort

#include "chi2.h"
#include "pval.h"



namespace SignificantPattern
{

SignificantFeaturesSearchTaroneCmh::SignificantFeaturesSearchTaroneCmh()
    : SignificantFeaturesSearchWithCovariates()
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh()\n");
    #endif
    execute_constructor_taronecmh();
}
SignificantFeaturesSearchTaroneCmh::~SignificantFeaturesSearchTaroneCmh() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantFeaturesSearchTaroneCmh()\n");
    #endif
    execute_destructor_taronecmh();
}

void SignificantFeaturesSearchTaroneCmh::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::execute_constructor()\n");
    #endif
    execute_constructor_taronecmh();
}
void SignificantFeaturesSearchTaroneCmh::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::execute_destructor()\n");
    #endif
    execute_destructor_taronecmh();
    super::execute_destructor();
}

void SignificantFeaturesSearchTaroneCmh::execute_constructor_taronecmh() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::execute_constructor_taronecmh()\n");
    #endif
    // Initialise grid of candidate corrected significance thresholds
    double log10_p; unsigned short j;
    for(log10_p=0, j=0; j<=NGRID; log10_p-=log10_p_step, j++)
        pgrid[j] = pow( ((double)10.0),log10_p );

    //NOTE: fix above 10->10.0, to support
    K = 0;

    freq_constructor();
}
void SignificantFeaturesSearchTaroneCmh::execute_destructor_taronecmh(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::execute_destructor_taronecmh()\n");
    #endif
    freq_destructor();
}

void SignificantFeaturesSearchTaroneCmh::freq_init(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::freq_init()\n");
    #endif
    K = getNumCovariates();
    freq_par_cov = new longint*[L];
    freq_par_cov[0] = new longint[L*K];
    for(longint j=1; j<L; j++) freq_par_cov[j] = freq_par_cov[0] + j*K;
    freq_cnt_cmh = new longint[NGRID+1];
    freq_clear();
    // Do not call super to make sure the ununsed freq_{par,cnt} are not alloc'ed
}
void SignificantFeaturesSearchTaroneCmh::freq_clear(){
    std::fill_n(freq_par_cov[0], L*K, 0); //memset(freq_par_cov[0],0,L*K*sizeof(longint));
    std::fill_n(freq_cnt_cmh, NGRID+1, 0);
    // Do not call super to make sure the alloc'ed freq_{par,cnt} are not accessed
}
void SignificantFeaturesSearchTaroneCmh::freq_constructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::freq_constructor()\n");
    #endif
    freq_par_cov = 0; freq_cnt_cmh = 0;
}
void SignificantFeaturesSearchTaroneCmh::freq_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::freq_destructor()\n");
    fprintf(stderr, "\tfreq_par_cov=%p\n", (void *) freq_par_cov);
    #endif
    if (freq_par_cov) {
        #ifdef DEBUG
        fprintf(stderr, "\tfreq_par_cov[0]=%p\n", (void *) freq_par_cov[0]);
        #endif
        if (freq_par_cov[0]) delete [] freq_par_cov[0]; freq_par_cov[0] = 0;
        delete [] freq_par_cov;
    };
    #ifdef DEBUG
    fprintf(stderr, "\tfreq_cnt_cmh=%p\n", (void *) freq_cnt_cmh);
    #endif
    if (freq_cnt_cmh) delete [] freq_cnt_cmh;
    freq_constructor();
}


void SignificantFeaturesSearchTaroneCmh::algorithm_init(){
    super::algorithm_init();

    // Initialise threshold value
    idx_pgrid = 1; delta = pgrid[idx_pgrid];

    freq_init();

    K = getNumCovariates(); // Note: currently redundant - K is also set in freq_init()
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::algorithm_init(), K=%d, haveCovariatesChanged? %d \n", K, haveCovariatesChanged());
    #endif
    // Initialise tables of number of observations per class, cumulative number
    // of observations, and number of observations in positive class
    if (haveCovariatesChanged()) {
        Nt = getNumObservationsInClasses();
        nt = getNumPositiveObservationsInClasses();
        #ifdef DEBUG
        fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::algorithm_init(), covariates changed: Nt[0]=%ld, nt[0]=%ld\n", Nt[0], nt[0]);
        #endif

        // Note: resize(d,0) will set 0 value only for new elements
        cum_Nt.resize(K+1); std::fill(cum_Nt.begin(),cum_Nt.end(),0);

        Nt_nt.resize(K); std::fill(Nt_nt.begin(),Nt_nt.end(),0);
        hypercorner_bnd.resize(K);
        gammat.resize(K);
        gammabint.resize(K);
        for(unsigned short j=0; j<K; j++){
            cum_Nt[j+1] = cum_Nt[j] + Nt[j];
            Nt_nt[j] = Nt[j]-nt[j];
            hypercorner_bnd[j] = hypercorner_bnd_k(j);
            gammat[j] = ((double)nt[j])/Nt[j];
            gammabint[j] = gammat[j]*(1-gammat[j]);
        }

        //isprunable() vectors
        f_vals.resize(K);
        g_vals.resize(K);
        betas.resize(K);
        idx_betas_sorted.resize(K);

        markCovariatesSeen();
    }
}
void SignificantFeaturesSearchTaroneCmh::algorithm_end(){
    freq_destructor();
    super::algorithm_end();
}





/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

double SignificantFeaturesSearchTaroneCmh::compute_pval(longint a, longint* x){
    return cmh_pval(a, x, K, Nt.data(), gammat.data(), gammabint.data());
}

double SignificantFeaturesSearchTaroneCmh::compute_score(longint a, longint* x){
    return cmh_score(a, x, K, Nt.data(), gammat.data(), gammabint.data());
}

double SignificantFeaturesSearchTaroneCmh::compute_odds_ratio(longint* at, longint* x){
    return cmh_odds_ratio(at, x, K, Nt.data(), nt.data());
}

double SignificantFeaturesSearchTaroneCmh::compute_minpval(longint *x){
    return cmh_minpval(x, K, Nt.data(), nt.data(), Nt_nt.data(), gammat.data(),
            gammabint.data());
}

unsigned short SignificantFeaturesSearchTaroneCmh::bucket_idx(double pval){
    double idx;
    idx = floor(-log10(pval)/log10_p_step);
    if(idx<0) idx = 0;
    else if(idx>NGRID) idx = NGRID;
    // cast only after the target type range has been assured
    return (unsigned short) idx;
}

//int SignificantFeaturesSearchTaroneCmh::qsort_cmp_betas(const void *x, const void *y){
//    if ( betas[*((longint*)x)] < betas[*((longint*)y)] ) return (-1);
//    else return 1;
//}
void SignificantFeaturesSearchTaroneCmh::idx_betas_sort(unsigned short j)
{
    //qsort(idx_betas_sorted.data(),j,sizeof(unsigned short),&qsort_cmp_betas); //Equivalent to argsort(betas[0:j])
    std::sort(idx_betas_sorted.begin(), idx_betas_sorted.begin()+j,
        [&](unsigned short x, unsigned short y) { return  betas[x] < betas[y]; }
    );
}
double SignificantFeaturesSearchTaroneCmh::compute_envelope_minpval(longint* x) {
    double f_sum, g_sum, Tcmh_max_corner_l, Tcmh_max_corner_r, Tcmh_aux_corner;
    unsigned short j,k;
    longint dim_margin;

    // If for any of the K tables, its margin x is not prunable, then we cannot
    // prune the interval
    for(k=0; k<K; k++) if(notprunable_k(x[k], k)) return 0;

    // Compute the maximum value of left handside function
    for(j=0,k=0; k<K; k++){
        // Discard all dimensions for which we're at the border, as they don't
        // contribute to the function neither in the numerator nor the
        // denominator
        dim_margin = dim_margin_k(x[k], k);
        if(dim_margin > 0) {
            f_vals[j] = gammat[k]*dim_margin;
            g_vals[j] = gammabint[k]*x[k]*(1-((double)x[k])/Nt[k]);
            betas[j] = g_vals[j]/f_vals[j];
            idx_betas_sorted[j] = j;
            j++;
        }
    } // betas[0...j-1] > 0 and betas[j...K] == 0, with #{k: dim_margin <= 0} == K-j


    idx_betas_sort(j);
    f_sum = 0; g_sum = 0; Tcmh_max_corner_l = 0;
    for(k=0; k<j; k++){
        f_sum += f_vals[idx_betas_sorted[k]];
        g_sum += g_vals[idx_betas_sorted[k]];
        Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
        Tcmh_max_corner_l = (Tcmh_max_corner_l >= Tcmh_aux_corner) ? Tcmh_max_corner_l : Tcmh_aux_corner; //Tcmh_max_corner_l=max(Tcmh_max_corner_l,Tcmh_aux_corner)
    }
    // Compute the maximum value of right handside function
    for(j=0,k=0; k<K; k++){
        // Discard all dimensions for which we're at the border, as they don't
        // contribute to the function neither in the numerator nor the
        // denominator
        dim_margin = dim_margin_k(x[k], k);
        if(dim_margin > 0) {
            f_vals[j] = (1-gammat[k])*dim_margin;
            //g_vals doesn't change, hence it does not need to be recomputed
            betas[j] = g_vals[j]/f_vals[j];
            idx_betas_sorted[j] = j;
            j++;
        }
    } // betas[0...j-1] > 0 and betas[j...K] == 0, with #{k: dim_margin <= 0} == K-j

    idx_betas_sort(j);
    f_sum = 0; g_sum = 0; Tcmh_max_corner_r = 0;
    for(k=0; k<j; k++){
        f_sum += f_vals[idx_betas_sorted[k]];
        g_sum += g_vals[idx_betas_sorted[k]];
        Tcmh_aux_corner = (f_sum*f_sum)/g_sum;
        Tcmh_max_corner_r = (Tcmh_max_corner_r >= Tcmh_aux_corner) ? Tcmh_max_corner_r : Tcmh_aux_corner; //Tcmh_max_corner_r=max(Tcmh_max_corner_r,Tcmh_aux_corner)
    }
    return Chi2_sf(((Tcmh_max_corner_r >= Tcmh_max_corner_l) ? Tcmh_max_corner_r : Tcmh_max_corner_l),1);
}

void SignificantFeaturesSearchTaroneCmh::decrease_threshold(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantFeaturesSearchTaroneCmh::decrease_threshold(), pre m=%ld, delta=%g, idx_pgrid=%d, freq_cnt[idx_pgrid]=%ld\n", m, delta, idx_pgrid, freq_cnt_cmh[idx_pgrid]);
    #endif
    // Remove the intervals which become untestable after the change
    m -= freq_cnt_cmh[idx_pgrid];
    // Change threshold
    idx_pgrid++; delta = pgrid[idx_pgrid];
}

} /* namespace SignificantPattern */
