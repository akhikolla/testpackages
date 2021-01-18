/*
 * SignificantIntervalSearchChi.cpp
 */

#include "SignificantIntervalSearchChi.h"

#include <stdlib.h>

#include "chi2.h"
#include "pval.h"



/* CONSTANT DEFINES */
#define NO_VERBOSE 1



using namespace std;

extern double Chi2_sf(double,double);

namespace SignificantPattern
{

SignificantIntervalSearchChi::SignificantIntervalSearchChi()
    : SignificantIntervalSearchFais()
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchChi()\n");
    #endif
    execute_constructor_chi();
}
SignificantIntervalSearchChi::~SignificantIntervalSearchChi() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchChi()\n");
    #endif
    execute_destructor_chi();
}

void SignificantIntervalSearchChi::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchChi::execute_constructor()\n");
    #endif
    execute_constructor_chi();
}
void SignificantIntervalSearchChi::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchChi::execute_destructor()\n");
    #endif
    execute_destructor_chi();
    super::execute_destructor();
}

void SignificantIntervalSearchChi::execute_constructor_chi() {
    class_ratio = 0; class_ratio_bin = 0;
}
void SignificantIntervalSearchChi::execute_destructor_chi(){
}


void SignificantIntervalSearchChi::algorithm_init(){
    // Precompute constants for psi_init (called from algorithm_init)
    class_ratio = ((double)n)/N; class_ratio_bin = class_ratio*(1-class_ratio);
    super::algorithm_init();
    // Initialise threshold value
    delta = psi[sl1];
}
void SignificantIntervalSearchChi::algorithm_end(){
    super::algorithm_end();
}


void SignificantIntervalSearchChi::psi_clear(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchChi::psi_clear()\n");
    #endif
    chi2_minpvals(N, n, N_over_2, class_ratio, class_ratio_bin, psi);
}


/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
double SignificantIntervalSearchChi::compute_pval(longint a, longint x){
    return chi2_pval(a, x, N, n, class_ratio_bin);
}

double SignificantIntervalSearchChi::compute_score(longint a, longint x){
    return chi2_score(a, x, N, n, class_ratio_bin);
}

double SignificantIntervalSearchChi::compute_odds_ratio(longint a, longint x){
    return odds_ratio(a, x, N, n);
}

double SignificantIntervalSearchChi::score_to_pval(double score){
    return Chi2_sf(score, 1);
}

} /* namespace SignificantPattern */
