/*
 * SignificantIntervalSearchFastCmh.cpp
 *
 *  Created on: 2 Mar 2017
 *      Author: mikolajr
 */

#include "SignificantIntervalSearchFastCmh.h"
#include "chi2.h" // Chi2_sf

#include <vector>

/* CONSTANT DEFINES */
#define NO_VERBOSE 1



using namespace std;


namespace SignificantPattern
{

SignificantIntervalSearchFastCmh::SignificantIntervalSearchFastCmh()
    : SignificantFeaturesSearch(), // explicitly construct a virtual base class
      SignificantIntervalSearch(), SignificantFeaturesSearchTaroneCmh()

{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFastCmh()\n");
    #endif
    execute_constructor_fastcmh();
}
SignificantIntervalSearchFastCmh::~SignificantIntervalSearchFastCmh() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchFastCmh()\n");
    #endif
    execute_destructor_fastcmh();
}

void SignificantIntervalSearchFastCmh::execute_constructor() {
    super_cov::execute_constructor();
    super_feat::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFastCmh::execute_constructor()\n");
    #endif
    execute_constructor_fastcmh();
}
void SignificantIntervalSearchFastCmh::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFastCmh::execute_destructor()\n");
    #endif
    execute_destructor_fastcmh();
    super_feat::execute_destructor();
    super_cov::execute_destructor();
}

void SignificantIntervalSearchFastCmh::execute_constructor_fastcmh() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFastCmh::execute_constructor_fastcmh()\n");
    #endif
    pValsTestableInts = IntervalSetWithOddsRatio();
    pValsSigInts = IntervalSetWithOddsRatio();
}
void SignificantIntervalSearchFastCmh::execute_destructor_fastcmh(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchFastCmh::execute_destructor_fastcmh()\n");
    #endif
}




/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */

/// @todo
/// factor out diff's wrt super class implementations to enable override only on
/// "for each of the K tables" parts of the code; this requires definition of
/// accessors methods to freq_par{,_cov} and freq_cnt_cmh
void SignificantIntervalSearchFastCmh::process_first_layer_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    longint tau, j, k, queue_idx;
    unsigned char *X_tr_aux;
    longint *aux_ptr;
    double pmh_min_val;
    // Process each length 1 interval
    for(tau=0; tau<L; tau++){
        n_featuresets_processed++;
        // Compute number of 1s in the interval for each of the K tables
        X_tr_aux = X_tr[tau];
        for(k=0; k<K; k++){
            aux_ptr = &freq_par_cov[tau][k];
            for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) *aux_ptr += X_tr_aux[j];
        }
        #ifndef NO_SINGLE_FEATURES
            // If the interval is testable...
            // Update frequency-buckets and number of testable intervals
            pmh_min_val = compute_minpval(freq_par_cov[tau]);
            if(istestable(pmh_min_val)){//Definition of testability in terms of critical values
                // Increase counter of appropriate bucket and overall number of testable intervals
                freq_cnt_cmh[bucket_idx(pmh_min_val)]++; m++;
                update_threshold();
            }
        #endif
        // If either the current interval or the previous one are prunable
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

void SignificantIntervalSearchFastCmh::process_intervals_threshold(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char **X_par = genotype_par.getMatrixPtr();
    unsigned short k;
    longint tau, j, queue_idx;
    unsigned char *X_tr_aux, *X_par_aux;
    longint *aux_ptr;
    double pmh_min_val;
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
            printf("\tProcessing layer %lld...\n",l+1);
            #endif
        }
        if((L_max>0) && ((l+1) > L_max)) {
            #ifndef NO_VERBOSE
            printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
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
        // Compute OR and frequency of the interval for each of the K tables
        X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
        for(k=0; k<K; k++){
            aux_ptr = &freq_par_cov[tau][k];
            for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; (*aux_ptr)++;}
        }
        // If the interval is testable, increase counter of testable items and frequency-buckets and
        // check if the corrected significance threshold must be reduced
        pmh_min_val = compute_minpval(freq_par_cov[tau]);
        if(istestable(pmh_min_val)){//Definition of testability in terms of critical values
            // Increase counter of appropriate bucket and overall number of testable intervals
            freq_cnt_cmh[bucket_idx(pmh_min_val)]++; m++;
            update_threshold();
        }
        // If either the current interval or the previous one are prunable
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

void SignificantIntervalSearchFastCmh::process_first_layer_pvalues(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char *Y_tr = phenotype.getVectorPtr();
    unsigned short k;
    longint tau, j, queue_idx, a;
    std::vector<longint> at(K);
    unsigned char *X_tr_aux;
    longint *aux_ptr;
    double score, odds_ratio, pval;
    // Clear the current layer frequency counters
    freq_clear();
    // Process each length 1 interval
    for(tau=0; tau<L; tau++){
        // Direct pointer to the relevant feature vector
        X_tr_aux = X_tr[tau];
        // Compute number of 1s in the interval for each of the K tables
        for(k=0; k<K; k++){
            aux_ptr = &freq_par_cov[tau][k];
            for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) *aux_ptr += X_tr_aux[j];
        }
        // If the interval is testable...
        // Update frequency-buckets and number of testable intervals
        #ifndef NO_SINGLE_FEATURES
            //here we check if pval_min of test is smaller than threshold
            if(istestable_int(tau)){//Definition of testability in terms of critical values
                // Compute global cell count considering all tables together
                //a = 0;
                //for(j=0; j<N; j++) if(X_tr_aux[j]) a += Y_tr[j];
                a = 0;
                std::fill(at.begin(), at.end(), 0);
                for(k=0; k<K; k++){
                    for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) if(X_tr_aux[j]) at[k] += Y_tr[j];
                    a += at[k];
                }

                // Compute the p-value
                //pval = compute_interval_pval(a, tau); n_pvalues_computed++;
                score = compute_interval_score(a, tau);
                pval = score_to_pval(score);
                odds_ratio = compute_interval_odds_ratio(at, tau);  // Compute odds ratio
                n_pvalues_computed++;
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

//Instead of writing to file, saved to vector
void SignificantIntervalSearchFastCmh::process_intervals_pvalues(){
    unsigned char **X_tr = genotype.getMatrixPtr();
    unsigned char **X_par = genotype_par.getMatrixPtr();
    unsigned char *Y_tr = phenotype.getVectorPtr();
    unsigned short k;
    longint tau, j, queue_idx, a;
    std::vector<longint> at(K);
    unsigned char *X_tr_aux, *X_par_aux;
    longint *aux_ptr;
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
            printf("\tProcessing layer %lld...\n",l+1);
            #endif
        }
        if((L_max>0) && ((l+1) > L_max)) {
            #ifndef NO_VERBOSE
            printf("\tMaximum interval length achieved at l=%lld. Stopping enumeration...\n",l+1);
            #endif
            break;
        }
        last_tau = tau;
        // In this case, the testable region does not change, so we don't need to check if the interval
        // has to be pruned now. If it was appended to the queue, then it has to be processed for sure
        // Compute OR and frequency of the interval for each of the K tables
        X_tr_aux = X_tr[tau+l]; X_par_aux = X_par[tau];
        for(k=0; k<K; k++){
            aux_ptr = &freq_par_cov[tau][k];
            for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) if((!X_par_aux[j]) && X_tr_aux[j]){ X_par_aux[j] = 1; (*aux_ptr)++;}
        }
        // If the interval is testable, increase counter of testable items and frequency-buckets and
        // check if the corrected significance threshold must be reduced
        if(istestable_int(tau)){//Definition of testability in terms of critical values
            // Compute global cell count considering all tables together
            //a = 0;
            //for(j=0; j<N; j++) if(X_par_aux[j]) a += Y_tr[j];
            a = 0;
            std::fill(at.begin(), at.end(), 0);
            for(k=0; k<K; k++){
                for(j=cum_Nt[k]; j<cum_Nt[k+1]; j++) if(X_par_aux[j]) at[k] += Y_tr[j];
                a += at[k];
            }

            // Compute the P-value
            //pval = compute_interval_pval(a, tau); n_pvalues_computed++;
            score = compute_interval_score(a, tau);
            pval = score_to_pval(score);
            odds_ratio = compute_interval_odds_ratio(at, tau);  // Compute odds ratio
            n_pvalues_computed++;
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


} /* namespace SignificantPattern */
