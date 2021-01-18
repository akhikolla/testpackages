#include <hyperg.h>
#include <cmath>

// 22/02/2012: corrected Laplace approximation formula
// compute the log of the psi function,
// either by calling the hyp2f1 function or by a Laplace approximation
double logPsi(double b, double c, int n, int p, double R2){

	/* R2 = usual coefficient of determination
	 n = sample size
	 p = number of columns of X (including (!) intercept)
	 b,c = parameters of psi function
	 */

	double pmodel = p - 1; // now only the covariates without intercept
	
	double fterm = hyp2f1((n - 1) / 2.0, b, (pmodel + c) / 2.0, R2);
	double ret;
	
	if (R_FINITE(fterm)){
    	
    	// hyp2f1 return value was finite, so we return just:
    	ret = lbeta(b, (pmodel + c) / 2.0 - b) + log(fterm);    	
    
    } else {
    	
    	// there were numerical problems with hyp2f1, so we do a Laplace approximation:
    	
    	// components of quadratic equation in exp(tau) are
    	double alpha = (2 * b - pmodel - c) * (1 - R2);
    	
    	double beta = 4 * b - pmodel - c + R2 * (n - 1 - 2 * b);
    	
    	double gamma = 2 * b;
    	
    	// exp(hat(tau)) is the mode of g:
    	double gMode = (- beta - sqrt(pow(beta, 2.0) - 4 * alpha * gamma)) / (2 * alpha);

    	// hat(tau)
    	double tauMode = log(gMode);

    	// the log integrand evaluated at the mode is
    	double logIntegrandAtTauMode = b * tauMode 
    				       + (n - 1 - pmodel - c) / 2.0 * log1p(gMode)
    				       - (n - 1) / 2.0 * log1p((1 - R2) * gMode);
    	
    	// the log variance is
    	double logSigma2 = 
    		M_LN2
    		- tauMode
    		- logspace_sub(
    				log(static_cast<double>(n - 1)) + log1p(-R2) - 2 * log1p((1 - R2) * gMode),
    				log(n - 1 - pmodel - c) - 2 * log1p(gMode) 
    			);
//   		- log(
//    				(n - 1) * (1 - R2) / pow(1 + (1 - R2) * gMode, 2)
//    				- (n - 1 - pmodel - c)  / pow(1 + gMode, 2) 
//    			);
//    	
    	// so the Laplace approximation is
    	ret = M_LN_SQRT_2PI + 0.5 * logSigma2 + logIntegrandAtTauMode;    	
    }
	
    return(ret);
}


// compute the log of formula (17), p. 415 in Liang et al.:
double logBF_hyperg(double R2, int n, int p, double alpha)
{  
	/* log Marginal under reference prior for phi & intercept
     hyperg prior for regression coefficients; alpha > 2 
     log marginal for null model is 0 */

    double ret;
	
    if (p == 1) 
    	ret = 0.0;
    else 
    	ret = log(alpha/2.0 - 1.0) + logPsi(1.0, alpha, n, p, R2);
    
    return(ret);
}


// compute formula (18), the posterior expected g under hyper-g prior
double posteriorExpectedg_hyperg(double R2, int n, int p, double alpha, double logBF)
{ 
	double ret;
	
	if (p == 1) // null model: no g here...
		ret = 0.0;
	else  // non-null model: use the Bayes factor, because it is the normalizing constant
			// of the posterior density of g
		ret = exp(log(alpha/2.0 - 1.0) + logPsi(2.0, alpha, n, p, R2) - logBF);
	
	return(ret); 
}

// compute formula (19), the posterior expected shrinkage factor g/(1+g) under hyper-g prior
double posteriorExpectedShrinkage_hyperg(double R2, int n, int p, double alpha, double logBF)
{ 
	double ret;
	
	if (p == 1) // null model: no g here...
		ret = 0.0;
	else  // non-null model: use the Bayes factor, because it is the normalizing constant
			// of the posterior density of g
		ret = exp(log(alpha/2.0 - 1.0) + logPsi(2.0, alpha + 2, n, p, R2) - logBF);
	
	return(ret); 
}
