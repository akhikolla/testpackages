/*
 * SignificantIntervalSearchFais.cpp
 *
 *  Created on: 31 Mar 2017
 *      Author: mikolajr
 */

#include "SignificantIntervalSearchFais.h"


/* CONSTANT DEFINES */
#define NO_VERBOSE 1



namespace SignificantPattern {

SignificantIntervalSearchFais::SignificantIntervalSearchFais() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais()\n");
    #endif
    execute_constructor_fais();
}
SignificantIntervalSearchFais::~SignificantIntervalSearchFais() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchFais()\n");
    #endif
    execute_destructor_fais();
}

void SignificantIntervalSearchFais::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::execute_constructor()\n");
    #endif
    execute_constructor_fais();
}
void SignificantIntervalSearchFais::execute_destructor() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::execute_destructor()\n");
    #endif
    execute_destructor_fais();
    super::execute_destructor();
}
void SignificantIntervalSearchFais::execute_constructor_fais() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::execute_constructor_fais()\n");
    #endif
    pValsTestableInts = IntervalSetWithOddsRatio();
    pValsSigInts = IntervalSetWithOddsRatio();

    sl1=0; sl2=0; su1=0; su2=0;
    flag=0;

    psi_constructor();

    freq_constructor();
}
void SignificantIntervalSearchFais::execute_destructor_fais() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::execute_destructor_fais()\n");
    #endif
    psi_destructor();
    freq_destructor();
}

void SignificantIntervalSearchFais::algorithm_init(){
    super::algorithm_init();

    sl1 = 1; sl2 = N_over_2; su1 = N-sl2; su2 = N-sl1;
    flag = 1;

    psi_init();
    freq_init();

}
void SignificantIntervalSearchFais::algorithm_end(){
    freq_destructor();
    super::algorithm_end();

    SummaryFais& summary = getSummaryRef();
    summary.setTestabilityRegion(sl1, sl2, su1, su2);
}

void SignificantIntervalSearchFais::psi_init(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::psi_init()\n");
    #endif
    psi = new double[N+1];
    psi_clear();
}
void SignificantIntervalSearchFais::psi_constructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::psi_constructor()\n");
    #endif
    psi = 0;
}
void SignificantIntervalSearchFais::psi_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::psi_destructor()\n");
    fprintf(stderr, "\tpsi=%p\n", (void *) psi);
    #endif
    if (psi) delete [] psi;
    psi_constructor();
}

void SignificantIntervalSearchFais::freq_init(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::freq_init()\n");
    #endif
    freq_par = new longint[L];
    freq_cnt = new longint[N+1];
    freq_clear();
}
void SignificantIntervalSearchFais::freq_clear(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::freq_clear()\n");
    #endif
    std::fill_n(freq_par, L, 0); //memset(freq_par,0,L*sizeof(longint));
    std::fill_n(freq_cnt, N+1, 0);
}
void SignificantIntervalSearchFais::freq_constructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::freq_constructor()\n");
    #endif
    freq_par = 0;
    freq_cnt = 0;
}
void SignificantIntervalSearchFais::freq_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFais::freq_destructor()\n");
    fprintf(stderr, "\tfreq_par=%p\n", (void *) freq_par);
    #endif
    if (freq_par) delete [] freq_par;
    #ifdef DEBUG
    fprintf(stderr, "\tfreq_cnt=%p\n", (void *) freq_cnt);
    #endif
    if (freq_cnt) delete [] freq_cnt;
    freq_constructor();
}

void SignificantIntervalSearchFais::decrease_threshold(){
    if(flag){ // Flag==1 means the last call to decrease_threshold() shrunk "the W" on the left side
        // Update number of testable intervals
        m -= freq_cnt[sl1]; m -= freq_cnt[su2];
        sl1++; su2--; // Shrink Sigma_k on extremes of the W
        // Check what the new case will be
        if (psi[sl1] >= psi[sl2]) delta = psi[sl1];
        else{ delta = psi[sl2]; flag = 0; }
    }else{ // Flag==0 means the last call to decrease_threshold() shrunk "the W" on the right side
        // Update number of testable intervals
        if(sl2==su1) m -= freq_cnt[sl2];//(beware of case sl2==su1 since it could lead to discounting the same thing twice!)
        else {m -= freq_cnt[sl2]; m -= freq_cnt[su1];}
        sl2--; su1++; // Shrink Sigma_k on center of the W
        // Check what the new case will be
        if (psi[sl1] >= psi[sl2]){ delta = psi[sl1]; flag = 1; }
        else delta = psi[sl2];
        //No need to update LCM minimum support in this case, since sl1 remains the same
    }
}

void SignificantIntervalSearchFais::process_first_layer_pvalues()
{
    longint tau, j, queue_idx, a;
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char *Y_tr = phenotype.getVectorPtr();
    unsigned char *X_tr_aux;
    double score, odds_ratio, pval;
    // Clear the current layer frequency counters
    freq_clear();
    // Process each length 1 interval
    for(tau=0; tau<L; tau++)
    {
        // Compute number of 1s in the interval
        X_tr_aux = X_tr[tau];
        // Compute number of 1s in the interval
        for(j=0; j<N; j++) freq_par[tau] += X_tr_aux[j];
        // If the interval is testable...
        // Update frequency-buckets and number of testable intervals
        #ifndef NO_SINGLE_FEATURES
        if(istestable_int(tau))
        {
            // Compute cell count
            a = 0;
            for(j=0; j<N; j++) if(X_tr_aux[j])
                a += Y_tr[j];
            // Compute the p-value
            //pval = compute_interval_pval(a, tau); n_pvalues_computed++;
            score = compute_interval_score(a, tau);
            pval = score_to_pval(score);  // Convert score into P-value
            odds_ratio = compute_interval_odds_ratio(a, tau);  // Compute odds ratio
            n_pvalues_computed++;
            //precompute_pvals(freq_par[tau]);
            // Check if the P-value is significant
            testAndSaveInterval(delta_opt, score, odds_ratio, pval, tau, l, a);
            //testAndSaveInterval(delta_opt, pval, tau, l, a);
        }
        #endif
        // If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
        // then do NOT append the left-child to the testable queue (i.e. prune it)
        if((tau==0) || isprunable_int(tau) || isprunable_int(tau-1)) continue;
        // Compute index of array position to be used, wrapping around if necessary
        queue_idx = testable_queue_front + testable_queue_length;
        queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
        // Actually append children to testable queue
        testable_queue[queue_idx] = tau-1;
        // Update queue length
        testable_queue_length++;
    }
}

void SignificantIntervalSearchFais::process_intervals_pvalues(){
    longint tau, j, queue_idx, a;
    unsigned char *X_tr_aux, *X_par_aux;
    unsigned char *Y_tr = phenotype.getVectorPtr();
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char **X_par = genotype_par.getMatrixPtr();
    double score, odds_ratio, pval;
    // While testable-interval queue is not empty, continue to process intervals
    while(testable_queue_length){
        // Pop a testable interval from the queue
        tau = testable_queue[testable_queue_front];
        testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
        testable_queue_length--;
        // Check if we have started processing a new layer by detecting non-monotonicity in tau
        if(tau < last_tau) {
            l++;
            #ifndef NO_VERBOSE
            printf("\tProcessing layer %lld...\n", (long long) l+1);
            #endif
        }
        if((L_max>0) && ((l+1) > L_max)) {
            #ifndef NO_VERBOSE
            printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n", (long long) l+1);
            #endif
            break;
        }
        last_tau = tau;
        // In this case, the testable region does not change, so we don't need to check if the interval
        // has to be pruned now. If it was appended to the queue, then it has to be processed for sure
        // Compute OR and frequency of the interval
        X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
        for(j=0; j<N; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; freq_par[tau]++;}
        // If the interval is testable, increase counter of testable items and frequency-buckets and
        // check if the corrected significance threshold must be reduced
        if(istestable_int(tau)){
            // Compute cell count
            a = 0;
            for(j=0; j<N; j++) if(X_par_aux[j]) a += Y_tr[j];
            // Compute the P-value
            //pval = compute_interval_pval(a, tau); n_pvalues_computed++;
            score = compute_interval_score(a, tau);
            pval = score_to_pval(score);
            odds_ratio = compute_interval_odds_ratio(a, tau);  // Compute odds ratio
            n_pvalues_computed++;
            //precompute_pvals(freq_par[tau]);
            // Check if the P-value is significant
            testAndSaveInterval(delta_opt, score, odds_ratio, pval, tau, l, a);
            //testAndSaveInterval(delta_opt, pval, tau, l, a);
        }
        // If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
        // then do NOT append the left-child to the testable queue (i.e. prune it)
        if((tau==0) || isprunable_int(tau) || isprunable_int(tau-1)) continue;
        // Compute index of array position to be used, wrapping around if necessary
        queue_idx = testable_queue_front + testable_queue_length;
        queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
        // Actually append children to testable queue
        testable_queue[queue_idx] = tau-1;
        // Update queue length
        testable_queue_length++;
    }
}

void SignificantIntervalSearchFais::process_first_layer_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    longint tau, j, queue_idx;
    // Process each length 1 interval
    for(tau=0; tau<L; tau++){
        n_featuresets_processed++;
        // Compute number of 1s in the interval
        for(j=0; j<N; j++) freq_par[tau] += X_tr[tau][j];
        // If the interval is testable...
        // Update frequency-buckets and number of testable intervals
        #ifndef NO_SINGLE_FEATURES
        if(istestable_int(tau)){
            freq_cnt[freq_par[tau]]++; m++;
            update_threshold();
        }
        #endif
        // If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
        // then do NOT append the left-child to the testable queue (i.e. prune it)
        if((tau==0) || isprunable_int(tau) || isprunable_int(tau-1)) continue;
        // Compute index of array position to be used, wrapping around if necessary
        queue_idx = testable_queue_front + testable_queue_length;
        queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
        // Actually append children to testable queue
        testable_queue[queue_idx] = tau-1;
        // Update queue length
        testable_queue_length++;
    }
}

void SignificantIntervalSearchFais::process_intervals_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char **X_par = genotype_par.getMatrixPtr();
    longint tau, j, queue_idx;
    // While testable-interval queue is not empty, continue to process intervals
    while(testable_queue_length){
        // Pop a testable interval from the queue
        tau = testable_queue[testable_queue_front];
        testable_queue_front = (testable_queue_front<(L-1)) ? testable_queue_front + 1: 0;
        testable_queue_length--;
        // Check if we have started processing a new layer by detecting non-monotonicity in tau
        if(tau < last_tau) {
            l++;
            #ifndef NO_VERBOSE
            printf("\tProcessing layer %lld...\n", (long long) l+1);
            #endif
        }
        if((L_max>0) && ((l+1) > L_max)) {
            #ifndef NO_VERBOSE
            printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n", (long long) l+1);
            #endif
            break;
        }
        last_tau = tau;
        // Check any of the two parents is prunable, stop processing. Notice that this check is necessary
        // even if the current interval was appended to the testable queue, because the threshold and
        // testability regions might have been modified between the time in which the current interval
        // was appended to the queue and the time in which it is being processed
        if(isprunable_int(tau) || isprunable_int(tau+1)) continue;
        n_featuresets_processed++;
        // Compute OR and frequency of the interval
        for(j=0; j<N; j++) if((!X_par[tau][j]) && X_tr[tau+l][j]){ X_par[tau][j] = 1; freq_par[tau]++;}
        // If the interval is testable, increase counter of testable items and frequency-buckets and
        // check if the corrected significance threshold must be reduced
        if(istestable_int(tau)){
            // Update frequency-buckets and number of testable items
            freq_cnt[freq_par[tau]]++; m++;
            update_threshold();
        }
        // If either the current interval or the previous one are prunable (i.e. have more than su2 ones)
        // then do NOT append the left-child to the testable queue (i.e. prune it)
        if((tau==0) || isprunable_int(tau) || isprunable_int(tau-1)) continue;
        // Compute index of array position to be used, wrapping around if necessary
        queue_idx = testable_queue_front + testable_queue_length;
        queue_idx = (queue_idx < L) ? queue_idx : queue_idx - L;
        // Actually append children to testable queue
        testable_queue[queue_idx] = tau-1;
        // Update queue length
        testable_queue_length++;
    }
}

} /* namespace SignificantPattern */
