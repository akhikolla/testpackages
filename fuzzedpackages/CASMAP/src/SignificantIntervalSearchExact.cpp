/*
 * SignificantIntervalSearchExact.cpp
 *
 *  Created on: Aug 18, 2016
 *      Author: mabaker
 */

#include "SignificantIntervalSearchExact.h"

#include <stdlib.h>

#include "pval.h"



/* CONSTANT DEFINES */
#define NO_VERBOSE 1
//#define RATIO_TH -23
#define LOGPVAL_EQ_RELTOL 1E-12 // exp(-27.6) \approx 1E-12


using namespace std;

namespace SignificantPattern
{

SignificantIntervalSearchExact::SignificantIntervalSearchExact()
    : SignificantIntervalSearchFais()
{
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchExact()\n");
    #endif
    execute_constructor_exact();
}
SignificantIntervalSearchExact::~SignificantIntervalSearchExact() {
    #ifdef DEBUG
    fprintf(stderr, "~SignificantIntervalSearchExact()\n");
    #endif
    execute_destructor_exact();
}

void SignificantIntervalSearchExact::execute_constructor() {
    super::execute_constructor();
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchExact::execute_constructor()\n");
    #endif
    execute_constructor_exact();
}
void SignificantIntervalSearchExact::execute_destructor(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchExact::execute_destructor()\n");
    #endif
    execute_destructor_exact();
    super::execute_destructor();
}

void SignificantIntervalSearchExact::execute_constructor_exact() {
    loggamma_constructor();
}
void SignificantIntervalSearchExact::execute_destructor_exact() {
    loggamma_destructor();
}

void SignificantIntervalSearchExact::algorithm_init(){
    super::algorithm_init();

    // Initialise threshold value
    delta = ((double) n)/N; //$\psi(1)=\frac{n}{N}$
    loggamma_init();
}



void SignificantIntervalSearchExact::loggamma_init(){
    if (!loggamma) {
        loggamma = new double[N+1];
        loggamma_clear();
    }
}
void SignificantIntervalSearchExact::loggamma_clear(){
    longint x;
    // Initialise cache with appropriate values
    for(x=0;x<=N;x++) loggamma[x] = lgamma(x+1);//Gamma(x) = (x-1)!
    // Initialise log_inv_binom_N_n
    log_inv_binom_N_n = loggamma[n] + loggamma[N-n] - loggamma[N];
}
void SignificantIntervalSearchExact::loggamma_constructor(){
     loggamma = 0;
     log_inv_binom_N_n = 0;
}
void SignificantIntervalSearchExact::loggamma_destructor() {
    if (loggamma) delete [] loggamma;
    loggamma_constructor();
}

void SignificantIntervalSearchExact::psi_clear(){
    #ifdef DEBUG
    fprintf(stderr, "SignificantIntervalSearchExact::psi_clear()\n");
    #endif
    fisher_minpvals(N, n, N_over_2, psi);
}


/* ---------------------------------------FUNCTIONS TO FIND THE SIGNIFICANT INTERVALS-------------------------------- */
double SignificantIntervalSearchExact::compute_pval(longint a, longint x){
    return fisher_pval(a, x, N, n, loggamma, log_inv_binom_N_n);
}

double SignificantIntervalSearchExact::compute_score(longint a, longint x){
    return fisher_pval(a, x, N, n, loggamma, log_inv_binom_N_n);  // TODO: Change dummy function
}

double SignificantIntervalSearchExact::compute_odds_ratio(longint a, longint x){
    return odds_ratio(a, x, N, n);
}

double SignificantIntervalSearchExact::score_to_pval(double score){
    return score;  // TODO: Change dummy function
}


} /* namespace SignificantPattern */
