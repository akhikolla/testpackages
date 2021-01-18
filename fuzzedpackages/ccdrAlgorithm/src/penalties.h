//
//  penalties.h
//  ccdr_proj
//
//  Created by Bryon Aragam on 6/18/14.
//  Copyright (c) 2014-2015 Bryon Aragam. All rights reserved.
//

#ifndef penalties_h
#define penalties_h

#include <math.h>

#include "defines.h"

//------------------------------------------------------------------------------/
//   PENALTY FUNCTION DEFINITIONS
//------------------------------------------------------------------------------/

//
// sign
//
//   Returns the usual sign function for a real number
//
double sign(double x){
	if(x > 0) return 1;
	else if(x < 0) return -1;
	else return 0;
}

//
// MCPPenalty
//
//   Minimax concave penalty function; equivalent to p_lambda(t; gamma) for the MCP.
//
//   Input:
//      t = value to penalize (e.g. coefficient)
//      lambda = regularization parameter
//      gamma = concavity paramater
//   Output: value of p_lambda(t; gamma)
//
double MCPPenalty(double b, double lambda, double gamma){
    if(b < gamma * lambda)
        return lambda * (b - 0.5 * b * b / (gamma * lambda));
    else
        return 0.5 * lambda *lambda * gamma;
}

//
// MCPThreshold
//
//   Thresholding function for MCP, equivalent to solving
//          argmin_t (t-z)^2 / 2 + p_lambda(t; gamma)
//
//   See Section 5.2.1 and Mazumder et al (2011) for details of this function and its derivation.
//
//   Update: we included the case where gamma <= 1.
//   This case is not discussed in Mazumder et al (2011), but we added it in case it happens.
//   In Zhang (2010), see the paragraph following formula (7.1) for a discussion on the minimizer by cases gamma > 1, and gamma <= 1.
double MCPThreshold(double z, double lambda, double gamma){
    if(gamma > 1) {
        if(fabs(z) <= lambda){
            return 0;
        } else if(lambda < fabs(z) && fabs(z) <= lambda * gamma){
            return sign(z) * gamma * (fabs(z) - lambda) / (gamma - 1.0);
        } else if(fabs(z) > lambda * gamma){
            return z;
        }
    } else {
        // this may happen when absorbing penalty weight alpha_j into penalty term
        // so gamma becomes alpha_j*gamma
        if(fabs(z) <= lambda * gamma) {
            return 0;
        } else return z;
    }

    #ifdef _DEBUG_ON_
        FILE_LOG(logERROR) << "There was a problem calculating the threshold function: z = " << z;
    #endif

    return 0;
}

//
// LassoPenalty
//
//   Lasso penalty function; equivalent to p_lambda(t) = lambda*t.
//
//   Input:
//      t = value to penalize (e.g. coefficient)
//      lambda = regularization parameter
//   Output: lambda * t
//
//   NOTE: The Lasso penalty is a convex function and hence has no concavity parameter
//
double LassoPenalty(double b, double lambda, double gamma = 0){
    return lambda * b;
}

//
// LassoThreshold
//
//   Thresholding function for Lasso, equivalent to solving
//          argmin_t (t-z)^2 / 2 + lambda * |t|
//
//   See Mazumder et al (2011) for the derivation of this function.
//
double LassoThreshold(double z, double lambda, double gamma = 0){
    if(fabs(z) <= lambda){
        return 0;
    } else {
        return sign(z) * (fabs(z) - lambda);
    }

    #ifdef _DEBUG_ON_
        FILE_LOG(logERROR) << "There was a problem calculating the threshold function: z = " << z;
    #endif

    return 0;
}
#endif
