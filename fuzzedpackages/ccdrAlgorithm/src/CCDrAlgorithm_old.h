//
//  CCDrAlgorithm.h
//  ccdr_proj
//
//  Created by Bryon Aragam on 3/20/14.
//  Copyright (c) 2014-2015 Bryon Aragam. All rights reserved.
//

#ifndef CCDrAlgorithm_h
#define CCDrAlgorithm_h

#include <vector>
#include <math.h>

#include "defines.h"

//------------------------------------------------------------------------------/
//   CCDR ALGORITHM CLASS
//------------------------------------------------------------------------------/

//
// Keeps track of data that is relevant to the running of the CCDr algorithm, but not strictly necessary
//   for keeping track of the data associated with the betas matrix (this data is tracked by SparseBlockMatrix).
//   Examples of such data includes the maximum number of iterations per sweep, the error threshold, and whether
//   or not the algorithm should terminate based on these criteria (among others).
//
// The most important parameters are supplied by the user:
//   
//   maxIters = the maximum number of iterations 
//   eps = the error threshold used to determine convergence
//   alpha = the multiplier used to terminate the algorithm once the active set gets too large 
//           (i.e. activeSetLength > alpha * pp)
//   pp = the number of nodes in the model being estimated by the algorithm
//
// Other important parameters are tracked automatically:
//
//   numSweeps = the total number of sweeps run so far
//   maxAbsError = the total accumulated error from each single parameter update run so far
//
// There is also a vector called 'stopFlags' which is used to keep track of the various reasons for terminating the
//   algorithm. We made this a vector so that if new stopping conditions are added (or we want to just keep track
//   of something we currently aren't), it will be easy to simply push_back this value.
//
// The main purpose of this class is to keep track of the data that ultimately decides whether or not the algorithm
//   keeps running. We recall how this decision is made here:
//
//   1) First, a full sweep of every parameter in the model is made
//          a) If, at any point, activeSetLength > maxEdges = alpha*pp, terminate immediately and disregard this model
//          ***When run over a grid of lambdas, the estimated model for each lambda strictly BEFORE termination is kept
//   2) Based on the active set computed in (1), we iterate over this active set performing SPUs
//       for each active edge until:
//          a) The updates have converged (i.e. the maximum absolute error after updating all the active edges is < eps)
//          b) The maximum number of iterations (maxIters) has been reached 
//   3) Go back to (1), until (1a) or the total number of complete sweeps (over all parameters) is > maxIters
//
// There are thus two functions for checking a stopping condition:
// 
//   keepGoing() checks (3) => model for this lambda is completely finished
//   moar() checks (2a) and (2b) => model for this active set is finished, but another complete sweep will follow
//
class CCDrAlgorithm{
    
public:
    // user-defined input
    unsigned int maxIters;  // maximum number of iterations for the algorithm
    double eps;             // convergence threshold
    
    //
    // Constructors
    //
    CCDrAlgorithm(unsigned int m, double e, double a, unsigned int p);
    
    //
    // Member functions
    //
    bool keepGoing() const;         // determines whether or not the algorithm should continue or terminate
    bool moar(int iters) const;     // determines whether or not to continue iterating over the current active set
    int edgeThreshold() const;      // get the maximum number of edges allowed (maxEdges = alpha*pp)
    double getError() const;        // get the current accumulated error
    int getStopFlag(int f) const;   // return the value of one of the stop flags
    void activeSetChanged();        // set the flag that indicates the active set has changed (and hence need to continue)
    void belowThreshold();          // set the flag that indicates the size of the active set is below the threshold level alpha*pp (and hence need to continue)
    void resetFlags();              // reset all stop flags to zero
    void updateError(double e);     // add a value to the error term
    void resetError();              // reset the error term (maxAbsError) to zero
    void addSweep();                // increment numSweeps
    
private:
    // 
    // This vector keeps track of whether or not to continue iterating the coordinate descent
    //  procedure (0 = stop, 1 = continue). The various slots refer to different reasons for
    //  terminating the algorithm:
    // stopFlags[0] = active set has changed
    // stopFlags[1] = active set is smaller than threshold set by user (alpha * pp)
    // 
    std::vector<int> stopFlags; // to keep track of whether or not to continue iterating
    double alpha;
    unsigned int maxEdges;
    
    // thresholds
    unsigned int numSweeps; // to keep track of how many full sweeps we have performed, including each check of the active set
    double maxAbsError;     // to store the error from each iteration of the CCDr algorithm
    
};

// Explicit constructor
CCDrAlgorithm::CCDrAlgorithm(unsigned int m, double e, double a, unsigned int p){
    maxIters = m;
    eps = e;
    alpha = a;
    maxEdges = round(a * p);
    numSweeps = 0;
    maxAbsError = 0;
    stopFlags = std::vector<int>(2, 0);
}

//
// Checks stopFlags to determine whether or not to continue running more complete sweeps
//
bool CCDrAlgorithm::keepGoing() const{
#ifdef _DEBUG_ON_
    if(stopFlags[0] == 0){
        FILE_LOG(logDEBUG1) << "After running concaveCDInit, active set has not changed: numSweeps = " << numSweeps;
    } else if(stopFlags[0] == 1){
        FILE_LOG(logDEBUG1) << "Active set has changed: numSweeps = " << numSweeps;
    }
    if(stopFlags[1] == 0){
        FILE_LOG(logDEBUG1) << "After running concaveCDInit, active set has exceeded edge threshold: numSweeps = " << numSweeps;
    }
    if(numSweeps > maxIters){
        FILE_LOG(logDEBUG1) << "Maximum number of iterations of concaveCDInit reached with maxAbsError = " << maxAbsError << ": numSweeps = " << numSweeps << " > " << maxIters;
    }
#endif
    
    int prod = 1;
    
    // check the stop flags: if any element is zero, we terminate
    for(int i = 0; i < stopFlags.size(); ++i){
        prod *= stopFlags[i];
    }
    
    // check if maxIters has been exceeded
    if(numSweeps > maxIters) prod = 0;
    
    // if prod = 1, keep going, if prod = 0, stop
    return (prod > 0);
}

// 
// Checks to see if the iterations over a fixed active set should stop
//
bool CCDrAlgorithm::moar(int iters) const{
    #ifdef _DEBUG_ON_
    if(maxAbsError <= eps){
        FILE_LOG(logDEBUG1) << "Parameter values converged in concaveCD after " << iters << " iterations: maxAbsError = " << maxAbsError << " <= " << eps;
    } else if(iters > maxIters){
        FILE_LOG(logDEBUG1) << "Maximum number of iterations of concaveCD reached with maxAbsError = " << maxAbsError << ": iters = " << iters << " > " << maxIters;
    }
    #endif
    
    return (maxAbsError > eps && iters <= maxIters);
}

int CCDrAlgorithm::edgeThreshold() const{
    return maxEdges;
}

double CCDrAlgorithm::getError() const{
    return maxAbsError;
}

int CCDrAlgorithm::getStopFlag(int f) const{
    return stopFlags[f];
}

// See comments for stopFlags variable in class declaration
void CCDrAlgorithm::activeSetChanged(){
    stopFlags[0] = 1;
}

// See comments for stopFlags variable in class declaration
void CCDrAlgorithm::belowThreshold(){
    stopFlags[1] = 1;
}

// Set all flags to be zero
void CCDrAlgorithm::resetFlags(){
    // fast / clever way to reset everything to zero
    std::fill(stopFlags.begin(), stopFlags.end(), 0);
}

// Update the total error
//
// NOTE: The argument passed to this function is "raw", i.e. it can be nonzero and represents the raw difference
//       between the old and new values. This function enforces the implementation of the particular error function
//       used in the algorithm. By default, we use the L1 error, i.e. the sum of absolute differences of all the 
//       updates. 
//
//       Other candidates include maximum absolute error and the L2 error.
void CCDrAlgorithm::updateError(double e){
    maxAbsError += fabs(e); // L1 error
}

void CCDrAlgorithm::resetError(){
    maxAbsError = 0;
}

void CCDrAlgorithm::addSweep(){
    numSweeps++;
}

#endif
