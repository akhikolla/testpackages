#include "functions.h"

// bootstrap particle filter
double bootstrapPartFilter (int N, arma::vec pars, IntegerMatrix state, IntegerMatrix stateNew, 
                            NumericVector weights, NumericVector weightsNew, NumericMatrix dataset, SEXP func_)
{
    // N is number of particles
    // pars is vector of parameters
    // state and stateNew are integer vectors
    // weights* are vector of particle weights
    // dataset is matrix of time and then time series counts
    // func_ is simulation function
    
    // initialise variables
    int j, k, l, m, t, r;
    double totWeight = 0.0, u = 0.0, maxWeight = 0.0, cumWeight = 0.0;
    
    // initialise log-likelihood
    double LL = 0.0;
    
    // setup output
    NumericVector out(state.ncol() + 1);
    
    // extract function pointer
    funcPtr func = *XPtr<funcPtr>(func_);
    
    // setup data
    IntegerVector counts(dataset.ncol() - 1);
    
    // set weights
    totWeight = 1.0;
    for(k = 0; k < N; k++){
        weights[k] = 1.0 / ((double) N);
    }
    
    // loop over time series
    t = 0;
    while(t < (dataset.nrow() - 1) && R_finite(LL) != 0){
        
        // set data
        for(k = 0; k < counts.size(); k++) {
            counts[k] = (int) dataset(t + 1, k + 1);
        }
        
        // loop over particles
        for(k = 0; k < N; k++) {
            // simulate forward
            if(t == 0) {
                r = k;
            } else {
                // resample a particle from the previous time step
                u = R::runif(0.0, 1.0);
                r = 0;
                cumWeight = weights[r];
                while(u > cumWeight){
                    r++;
                    cumWeight += weights[r];
                }
            }
            out = core_processing<funcPtr>(func, as<NumericVector>(wrap(pars)), 
                    dataset(t, 0), dataset(t + 1, 0), state(r, _), counts);
            for(j = 0; j < state.ncol(); j++){
                stateNew(k, j) = (int) out[j + 1];
            }
            // set new weight (on log-scale)
            weightsNew[k] = out[0];
            // Rprintf("count = %d stateNew = %d wn[%d] = %f\n", counts[0], stateNew(k, 1), k, weightsNew[k]);
            // Rprintf("w = %f\n", R::dpois(counts[0], 1e-5, 1));
            // for(j = 0; j < 4; j++) {
            //     Rprintf("u[%d] = %d ", j, stateNew(k, j));
            // }
            // Rprintf("\n");
        }
        
        // update states and weights (deep copy)
        state = clone(stateNew);
        weights = clone(weightsNew);
        
        if(any(is_finite(weights))){
            l = 0;
            while(R_finite(weights[l]) == 0){
                l++;
            }
            if(l >= weights.size()){
                stop("Error in bootstrap particle filter - too many non-finite weights");
            }
            // normalise weights
            maxWeight = weights[l];
            for(m = l; m < weights.size(); m++){
                if(R_finite(weights[m]) != 0){
                    maxWeight = (maxWeight > weights[m] ? maxWeight:weights[m]);
                }
            }
            totWeight = 0.0;
            for(m = 0; m < weights.size(); m++){
                if(R_finite(weights[m]) != 0){
                    weightsNew[m] = exp(weightsNew[m] - maxWeight);
                    totWeight += weightsNew[m];
                }
            }
            totWeight = maxWeight + log(totWeight);
            for(m = 0; m < weights.size(); m++){
                if(R_finite(weights[m]) == 0){
                    weights[m] = 0.0;
                } else {
                    weights[m] = exp(weights[m] - totWeight);
                }
                // Rprintf("weights[%d] = %f\n", m, weights[m]);
            }
            // Rprintf("max = %f tot = %f totsum = %f\n", maxWeight, totWeight, sum(weights));
            // getchar();
            
            // update log-likelihood
            LL += totWeight - log((double) N);
        } else {
            LL = NA_REAL;
            // Rprintf("LL = %f\n", LL);
        }
        
        // increment time step
        t++;
    }
    // return log-likelihood
    return(LL);
}
