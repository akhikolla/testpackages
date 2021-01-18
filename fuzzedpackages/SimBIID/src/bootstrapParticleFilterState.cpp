#include "functions.h"

// bootstrap particle filter, returns states
List bootstrapPartFilterState (int N, NumericMatrix pars, NumericMatrix dataset, IntegerVector iniStates, SEXP func_)
{
    // N is number of particles
    // pars is vector of parameters
    // dataset is matrix of time and then time series counts
    // func_ is simulation function
    
    // initialise variables
    int i, j, k, l, m, t, r, valid;
    double totWeight = 0.0, u = 0.0, maxWeight = 0.0, cumWeight = 0.0;
    
    // extract function pointer
    funcPtr func = *XPtr<funcPtr>(func_);
    
    // setup data
    IntegerVector counts(dataset.ncol() - 1);
    
    // set up output matrix
    int nclass = iniStates.size();
    List outlist(pars.nrow());
    
    // set up intermediary matrices
    IntegerMatrix states(dataset.nrow() * N, nclass);
    IntegerMatrix state(N, nclass);
    IntegerMatrix stateNew(N, nclass);
    NumericVector weights(N);
    NumericVector weightsNew(N);
    IntegerMatrix trajectories(N, dataset.nrow());
    IntegerMatrix tempout(dataset.nrow(), nclass);
    NumericVector tempsim(state.ncol() + 1);
    
    // loop over parameter sets
    for(i = 0; i < pars.nrow(); i++) {
        
        // set weights
        totWeight = 1.0;
        for(k = 0; k < N; k++){
            weights[k] = 1.0 / ((double) N);
        }
        
        // set initial states
        for(j = 0; j < N; j++){
            for(k = 0; k < nclass; k++){
                state(j, k) = iniStates[k];
            }
        }
        for(j = 0; j < N; j++){
            for(k = 0; k < nclass; k++){
                states(j, k) = state(j, k);
            }
        }
        
        // set initial trajectories
        for(k = 0; k < N; k++){
            trajectories(k, 0) = k;
        }
        
        // loop over time series
        t = 0; valid = 0;
        while(t < (dataset.nrow() - 1) && valid == 0){
            
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
                    if(r >= N) {
                        stop("Incorrect sampling of particle in bootPartFilterState\n");
                    }
                }
                
                // update trajectories
                trajectories(k, t + 1) = r;
                
                // run simulation
                tempsim = core_processing<funcPtr>(func, pars(i, _), 
                        dataset(t, 0), dataset(t + 1, 0), state(r, _), counts);
                        
                // update states
                for(j = 0; j < state.ncol(); j++) {
                    stateNew(k, j) = (int) tempsim[j + 1];
                }
                // set new weight (on log-scale)
                weightsNew[k] = tempsim[0];
            }
            
            // update states and weights (deep copy)
            state = clone(stateNew);
            weights = clone(weightsNew);
            for(j = 0; j < N; j++){
                for(k = 0; k < nclass; k++){
                    states((t + 1) * N + j, k) = state(j, k);
                }
            }
            
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
            } else {
                // stop("No non-zero weights\n");
                valid = 1;
            }
            t++;
        }
        if(valid == 1) {
            outlist[0] = NA_REAL;
        } else {
        
            // resample trajectory based on final weights
            u = R::runif(0.0, 1.0);
            r = 0;
            cumWeight = weights[r];
            while(u > cumWeight){
                r++;
                cumWeight += weights[r];
            }
            if(r >= N) {
                stop("Incorrect sampling of particle in bootPartFilterState\n");
            }
            // Rprintf("t = %d k = %d N = %d r = %d\n", t, k, N, r);
            for(k = 0; k < nclass; k++){
                tempout(t, k) = states(t * N + r, k);
            }
            while(t > 0){
                r = trajectories(r, t);
                t--;
                for(k = 0; k < nclass; k++){
                    tempout(t, k) = states(t * N + r, k);
                }
            }
            // for(k = 0; k < nclass; k++){
            //     tempout(0, k) = iniStates[k];
            // }
            
            // save output
            outlist[i] = clone(tempout);
        }
    }
    
    // return sampled states
    return outlist;
}
