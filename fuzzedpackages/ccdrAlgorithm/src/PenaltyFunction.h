//
//  PenaltyFunction.h
//  ccdr_proj
//
//  Created by Bryon Aragam on 3/19/14.
//  Copyright (c) 2014-2015 Bryon Aragam. All rights reserved.
//

#ifndef PenaltyFunction_h
#define PenaltyFunction_h

#include <math.h>

#include "defines.h"
#include "penalties.h"

//------------------------------------------------------------------------------/
//   PENALTY FUNCTION CLASS
//------------------------------------------------------------------------------/

//
// Implements various functions associated with the minimax concave penalty (MCP) 
//  and L1 penalty functions:
//
//   p = penalty
//   Dp = derivative
//   DDp = second derivative
//   threshold = associated threshold function (see SparseNet paper)
//
// Strictly speaking, other functions could be used here (SCAD, L1, capped L1, etc). 
//
class PenaltyFunction{
    
public:
    //
    // Constructor
    //
    PenaltyFunction(double g);
    
    double threshold(double z, double lambda) const{
        return thresholdPtr(z, lambda, gamma);
    }
    
    double p(double z, double lambda) const{
        return pPtr(z, lambda, gamma);
    }
    
private:
    // This is the shape parameter for the penalty function, denoted by gamma for the MCP.
    //  Note that for other penalty functions, this parameter is often denoted by a different
    //  letter (e.g. 'a' for SCAD).
    double gamma;
    double (*pPtr)(double, double, double);
    double (*thresholdPtr)(double, double, double);
};

// Explicit constructor
PenaltyFunction::PenaltyFunction(double g){
    if(g >= 0){
        gamma = g;
        pPtr = &MCPPenalty;
        thresholdPtr = &MCPThreshold;
    } else{
        gamma = g;
        pPtr = &LassoPenalty;
        thresholdPtr = &LassoThreshold;
    }
}

#endif
