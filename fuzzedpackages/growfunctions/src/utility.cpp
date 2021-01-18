#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

SEXP wishrnd(mat& L, const mat& V, double nu)
{
        BEGIN_RCPP
        int p = L.n_rows;
        int i, j;
        mat R = chol(V);
        mat A(p,p); A.zeros();
        for(i = 0; i < p; i++)
        {
            A(i,i) = sqrt( rchisq( 1, (nu-double(i)) )[0] );
            for(j = 0; j < i; j++)
            {
                if(i  >  0)
                {
                    A(i,j) = rnorm(1,0,1)[0];
                }
            }
        }
        L = trans(R) * A * trans(A) * R;
       END_RCPP
} /* end function rmvnqr for sampling mvn using QR decomposition */

colvec rdirich(colvec& shape)
{
     // input a vector and add elements at the last row
     int M = shape.n_rows;
     colvec gam_var(M); gam_var.zeros();
     for(int m = 0; m < M; m++)
     {
          gam_var(m) = rgamma(1,shape(m),1)[0];
     }
     colvec dir_var = gam_var / sum(gam_var);
	return dir_var;
} /* end function rdirch() generating dirichlet random variables */

double logmvndens(const colvec& b_i, const colvec& m, const mat& Q)
{
        // Compute multivariate Gaussian density for matrix, B
        // with mean, m, and precision matrix, Q
        double val_q, sign_q;
        log_det(val_q,sign_q,symmatl(Q));
        int T = b_i.n_rows;
        double c = log(2*M_PI)*(-0.5*T);
        colvec bi_center = b_i - m;
	      double logdens = c + 0.5*val_q - 
                         0.5*as_scalar( bi_center.t()*symmatl(Q)*bi_center );
        return(logdens);
}

double loggmrfdens_full(const colvec& b_i, const colvec& m, const mat& Q, const vec& eigraw, double kappa)
{
        // Compute multivariate Gaussian density for matrix, B
        // with mean, m, and rank deficient precision matrix, Q, 
        // of rank = df
        int df                = eigraw.n_elem; /* Q is rank-deficient */
        double c              = log(kappa/(2*M_PI))*(0.5*double(df));
        colvec bi_center      = b_i - m;
        double logdens        = c + 0.5*sum(log(eigraw)) - 
                                   0.5*kappa*as_scalar( bi_center.t()*symmatl(Q)*bi_center );
        return(logdens);
}

double loggmrfdens(const colvec& b_i, const colvec& m, const mat& Q, int df, double kappa)
{
        // Compute multivariate Gaussian density for matrix, B
        // with mean, m, and rank deficient precision matrix, Q, 
        // of rank = df
        // don't need to compute the eigenvalues if used in normalized weights
        double c              = log(kappa/(2*M_PI))*(0.5*double(df));
        colvec bi_center      = b_i - m;
        double logdens        = c  - 
                                   0.5*kappa*as_scalar( bi_center.t()*symmatl(Q)*bi_center );
        return(logdens);
}

double log_dnorm_vec(const rowvec& y, const rowvec& b, double tau_e)
{
        // Compute sum of univariate Gaussian densities, logged.
        int T = y.n_elem;
        int i = 0;
        double log_dens = 0;
        for( i = 0; i < T; i++ )
        {
             log_dens    += R::dnorm(y(i),b(i),(1/sqrt(tau_e)),true); /* true for log_density */
        }
        return(log_dens);
}

double logmatrixdens(const mat& B_i, const mat& P, const mat& Q)
{
        // NumericVector r(resid);
        // devvec is an nc x 1 vector
        double val_p, val_q, sign_p, sign_q;
        log_det(val_p,sign_p,P); log_det(val_q,sign_q,Q);
	   int T = B_i.n_cols;
	   int S = B_i.n_rows;
        double c = log(2*M_PI)*(-0.5*T*S);
	   double logdens = c + 0.5*T*val_p + 0.5*S*val_q - 
                         0.5*as_scalar( trace(Q*B_i.t()*P*B_i) );
        return(logdens);
}

double dev(const colvec& resid, double taue)
{
        // NumericVector r(resid);
        // devvec is an nc x 1 vector
        NumericVector r = wrap(resid);
        NumericVector devvec  = dnorm( r, 0.0, sqrt(1/taue), true ); /* logD */
        double deviance = accumulate(devvec.begin(),devvec.end(),0.0);
        deviance *= -2;
        return(deviance);
}

SEXP dmarg(const colvec& resid, double taue, rowvec& devmarg)
{
        BEGIN_RCPP
        // NumericVector r(resid);
        // devvec is an nc x 1 vector
        int nc = resid.n_elem;
        NumericVector r = wrap(resid); 
        NumericVector devvec  = dnorm( r, 0.0, sqrt(1/taue), false ); /* f(y|theta)*/
        rowvec devcond(devvec.begin(), nc, false);  devmarg = devcond;
        END_RCPP
}

SEXP cpo(const mat& Devmarg, rowvec& logcpo, double& lpml)
{
        BEGIN_RCPP
        mat invDmarg = pow(Devmarg,-1); /* invert each member of Devmarg */
        logcpo = mean(invDmarg,0); /* 1 x nc */
        logcpo = pow(logcpo,-1); /* invert again to achieve f(y_i|y_-i) */
        logcpo = log(logcpo); /* return vector of log(cpo_i) */
        lpml = sum(logcpo); /* sum each elem of log(cpo_i) */
        END_RCPP
}

SEXP dic3comp(const colvec& Deviance, const mat& Devmarg, colvec& devres)
{
        BEGIN_RCPP
        double dbar = mean(Deviance);
        rowvec devmean = mean(Devmarg,0);
        rowvec dens = log(devmean); /* 1 x nc */
        double dhat = sum( dens ); dhat *= -2;
        double dic = 2*dbar - dhat;
        double pd = dic - dbar;
        devres(0) = dic; devres(1) = dbar; devres(2) = dhat; devres(3) = pd;
        END_RCPP
}

SEXP rmvncov(const mat& phi_inv, const colvec& h, colvec& b)
{    // use to sample from prior where phi_inv is the P x P covariance (not precision)
     BEGIN_RCPP
     // build posterior variance and mean
     int p      	= phi_inv.n_cols;
	colvec noise 	= randn<colvec>(p);
	b 		     = trimatu(chol(symmatl(phi_inv))).t() * noise + h;
     END_RCPP
} /* end function rmvnbasic for drawing a single mvn sample */

SEXP rmvnchol(const mat& U, const colvec& h, colvec& b)
{    // use to sample from prior where phi_inv is the P x P variance (not precision)
     BEGIN_RCPP
     // U is the cholesky decomposition of a covariance matrix, C = U'U
     // build posterior variance and mean
     int p          = U.n_cols;
	colvec noise 	= randn<colvec>(p);
	b 		     = trimatu(U).t() * noise + h;
     END_RCPP
} /* end function rmvnbasic for drawing a single mvn sample */

SEXP rmvnbasic(const mat& phi, const colvec& e, colvec& b)
{    // use to sample posterior (by solving for h = phi-1*e) 
     // phi is a precision matrix
     BEGIN_RCPP
     // phi is a precision matrix
     // build posterior variance and mean
     int p 		= phi.n_cols;
     colvec h   	= solve(phi,e);
	colvec noise 	= randn<colvec>(p);
	b 		     = solve(trimatu(chol(symmatl(phi))),noise) + h;
     END_RCPP
} /* end function rmvnbasic for drawing a single mvn sample */

SEXP rmvnsample(const mat& phi, const colvec& h, colvec& b)
{
        BEGIN_RCPP
        // phi is a precision matrix
        // build posterior variance and mean
        // S = P^-1 = (U_p' * U_p) ^-1 = U_p^-1 * (U_p^-1)' = U'U
        // note: U' does not equal U_p^-1
	      // b = U_p^-1 * z + h -> cov(b) = U_p^-1 * (U_p^-1)' = S
        int p      	      = phi.n_cols;
        colvec noise 	 = randn<colvec>(p);
	   b 		      = solve(trimatu(chol(symmatl(phi))),noise) + h;
        END_RCPP
} /* end function rmvnbasic for drawing a single mvn sample */


unsigned long rdrawone(const colvec& pr, unsigned long k)
{
        // extract sort index of pr in descending order
        uvec pOrderIndex = sort_index(pr,"descend");
        // draw uniform random variate we will compare
        // to cdf composed from categories in descending
        // order
        double drawP = runif(1,0,1)[0];
        double pSum = 0.0;
        unsigned long i, x;
        for(i = 0; i < k; i++)
        {
          x = pOrderIndex(i);
          pSum += pr(x);
          if(pSum > drawP)
          {
            return x;
          } /* end conditional test on category i */
        } /* end loop over k categories */
        return pOrderIndex(k-1);
} /* end function rdrawone */
