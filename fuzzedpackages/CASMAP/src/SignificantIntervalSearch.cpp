/**
 * @file    SignificantIntervalSearch.cpp
 * @author  mabaker
 * @date    Aug 18, 2016
 *
 * Single class source file.
 *
 */
#include "SignificantIntervalSearch.h"



/* CONSTANT DEFINES */
#define NO_VERBOSE 1



using namespace std;

namespace SignificantPattern
{

SignificantIntervalSearch::SignificantIntervalSearch()
    : SignificantFeaturesSearch() // explicitly construct a virtual base class
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch()\n");
    #endif
    execute_constructor_int();
}
SignificantIntervalSearch::~SignificantIntervalSearch() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearch()\n");
    #endif
    execute_destructor_int();
}

void SignificantIntervalSearch::execute_constructor(){
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::execute_constructor()\n");
    #endif
    execute_constructor_int();
}
void SignificantIntervalSearch::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::execute_destructor()\n");
    #endif
    execute_destructor_int();
    super::execute_destructor();
}
void SignificantIntervalSearch::execute_constructor_int() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::execute_constructor_int()\n");
    #endif
    genotype_par = Genotype();
    last_tau=0;
    testable_queue_constructor();
}
void SignificantIntervalSearch::execute_destructor_int() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::execute_destructor_int()\n");
    #endif
    testable_queue_destructor();
}

void SignificantIntervalSearch::algorithm_init(){
    super::algorithm_init();
    // Initialise parent matrix
    genotype_par.initialiseMatrix(L, N);
    testable_queue_init();
}
void SignificantIntervalSearch::algorithm_end(){
    super::algorithm_end();

    SummaryInt& summary = getSummaryRef();
    summary.setMaxTestableIntervalLength(l+1);
}



void SignificantIntervalSearch::testable_queue_init(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::testable_queue_init()\n");
    #endif
    testable_queue = new longint[L];
    testable_queue_clear();
}
void SignificantIntervalSearch::testable_queue_clear(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::testable_queue_clear()\n");
    #endif
    std::fill_n(testable_queue, L, 0);
    testable_queue_front = 0; testable_queue_length = 0;
}
void SignificantIntervalSearch::testable_queue_constructor() {
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::testable_queue_constructor()\n");
    #endif
    testable_queue = 0;
    testable_queue_front = 0; testable_queue_length = 0;
}
void SignificantIntervalSearch::testable_queue_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearch::testable_queue_destructor()\n");
    fprintf(stderr, "\ttestable_queue=%p\n", (void *) testable_queue);
    #endif
    if (testable_queue) delete [] testable_queue;
    testable_queue_constructor();
}



void SignificantIntervalSearch::process_significant_features() {
    std::vector<longint> lVector;
    std::vector<longint> tauVector;
    getPValsSigInts().getLAndTauVectors(lVector, tauVector);
    std::vector<double> pValVector = getPValsSigInts().getPValueVector();
    std::vector<double> oddsRatioVector = getPValsSigInts().getOddsRatioVector();
    std::vector<double> scoreVector = getPValsSigInts().getScoreVector();

    intervals.cpp_intervalsFromMemory(tauVector, lVector, scoreVector, oddsRatioVector, pValVector);
    // Filter significant intervals
    filter.cpp_filterIntervalsFromMemory(tauVector, lVector, scoreVector, oddsRatioVector, pValVector);
}

    bool SignificantIntervalSearch::testAndSaveInterval(double threshold, double score, double odds_ratio, double pval, longint tau, longint l, longint a)
    {
        if(outputForTestableInts)
        {
            saveTestableInterval(pval, score, odds_ratio, tau, l, a);
        }
        bool isSignificant = (pval <= threshold);
        if(isSignificant)
        {
            saveSignificantInterval(pval, score, odds_ratio, tau, l, a);
            n_significant_featuresets++;
        }
        return isSignificant;
    }

void SignificantIntervalSearch::compute_corrected_significance_threshold(){
    // Give feedback to user
    #ifndef NO_VERBOSE
    printf("COMPUTING CORRECTED SIGNIFICANCE THRESHOLD...\n");
    #endif
    // Initialise the queue as empty
    testable_queue_clear();
    // Initialise current layer index, current number of testable intervals and current number of intervals processed to 0
    l = 0; m = 0; n_featuresets_processed = 0;
    // Initialise the value of the OR vectors of current layer to original dataset
    genotype_par = genotype;
    // Process the upper-most layer (i.e. the layer composed of length 1 intervals)
    #ifndef NO_VERBOSE
    printf("\tProcessing layer %lld...\n",(long long) l+1);
    #endif
    process_first_layer_threshold();
    // Artificially initialise last_tau to L-1 to ensure the first iteration of process_intervals()
    // increases the number of layers processed if the testable queue is non-empty
    last_tau = L-1;
    // Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
    process_intervals_threshold();
    // Set final corrected significance threshold
    compute_final_corrected_significance_threshold();
}

void SignificantIntervalSearch::find_significant_features(){
    // Give feedback to user
    #ifndef NO_VERBOSE
    printf("\n\nSCANNING DATASET FOR SIGNIFICANT INTERVALS...\n\n");
    #endif
    // Initialise the queue as empty
    testable_queue_clear();
    // Initialise current layer index and current number of computed p-values to 0
    l = 0; n_pvalues_computed = 0; n_significant_featuresets = 0;
    // Initialise the value of the OR vectors of current layer to original dataset
    genotype_par = genotype;
    // Process the upper-most layer (i.e. the layer composed of length 1 intervals)
    #ifndef NO_VERBOSE
    printf("\tProcessing layer %lld...\n", (long long) l+1);
    #endif
    process_first_layer_pvalues();
    // Artificially initialise last_tau to L-1 to ensure the first iteration of process_intervals()
    // increases the number of layers processed if the testable queue is non-empty
    last_tau = L-1;
    // Process the rest of layers (i.e. intervals of length > 1) until the pruning naturally stops the execution
    process_intervals_pvalues();
}

} /* namespace SignificantPattern */
