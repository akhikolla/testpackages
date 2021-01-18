#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

// update vector of cluster membership indicators, s(i),....,s(N)
SEXP clusterstep(const cube& B, mat& kappa_star, mat& B1, const uvec& o,
             const field<mat>& C, const mat& D, ucolvec& s, 
             //const field<sp_mat>& C,
             ucolvec& num, unsigned int& M, double& conc, int a, int b,
             const vec& ipr, colvec& Num)
    {
        BEGIN_RCPP
      
        // sample cluster assignments, s(1), ..., s(N)
        // B = (B_1,...,B_K), where B_k is N x T matrix for iGMRF term k
        // Q = (Q_1,...,Q_K), where Q_k is a T x T de-scaled iGMRF precision matrix
        // C = (C_1,...,C_K), where C_k = D_k^-1 * Omega_k, 
        // where Omega_k is the T x T adjacency matrix for iGMRF term, k
        // D is a K x T matrix where row k contains T diagonal elements of Q_k
        // K x M matrix, kappa_star records locations for each iGMRF term
        // o = (o_1,...,o_k) is a vector where each entry denotes the order of term K.
        // e.g. RW(1) -> o = 2, RW(2) -> o = 3, seas(3) -> o = 3
        int N = B.slice(0).n_rows;
        int T = B.slice(0).n_cols;
        int K = C.n_rows;
        double sweights = 0;
        // zro is the zeros.T vector 
        colvec zro(T); zro.zeros();
        uvec o_adjust = o; 
        //o_adjust.zeros();
        // capture quadratic product for rate kernel of posterior gamma
        // posterior for kappa_star(k,i).  
        // save B1 to latter (in another function) compute posterior for kappa_star
        // mat B1(K,N); 
        double a1k; /* posterior shape for kappa_star(k,i) under 1 obs */
        B1.zeros();
        int i, j, k;
	      unsigned int l;
           
        /* 
        mat D_k(T,T), Omega_k(T,T);
        cube Q(T,T,K);
        for(k = 0; k < k; k++)
        {
           D_k.zeros(); D_k.diag()        = D.row(k);
           Omega_k                      = D_k * C(k,0); 
           Q.slice(k)                   = D_k - Omega_k;
        } // end loop K over iGMRF terms 
        */
        
        for(i = 0; i < N; i++)
        {
            // check if _i assigned to singleton cluster
            // if so, remove the cluster associated to _i
            // and decrement the cluster labels for m > s(i)
            if(num(s(i)) == 1) /* remove singleton cluster */
            {
                kappa_star.shed_col(s(i));
                num.shed_row(s(i));
                Num.shed_row(s(i));
                M -= 1; /* decrement cluster count */

                //decrement cluster tracking values by 1 for tossed cluster
                s( find(s > s(i)) )          -= 1;
                
            } /* end cluster accounting adjustment for singleton cluster */
            else /* cluster contains more than one unit */
            {
                num(s(i))                    -= 1;
                /* scale up num to population totals, Num, based on H-T inverse probability estimator */
                Num(s(i))                    -= 1/ipr(i);
            } /* decrement non-singleton cluster count by one */

            // construct normalization constant, q0i, to sample s(i)
            // build loqq0 and exponentiate
            colvec bki(T), bbar_ki(T); /* T x 1, D_k^-1*Omega_k*b_ki = C(k,0)*b_ki */
            mat bbar_i(K,T); bbar_i.zeros();
            double logd_dk = 0; /* set of T 0 mean gaussian densities for term k */
            double logq0ki = 0, logq0i = 0, q0i = 0;
            // accumulate weight, q0i, for s(i) over K iGMRF terms  
            for( k = 0; k < K; k++)
            {
                 logq0ki       = 0; /* reset k-indexed log-like on each k */
                 //a1k           = 0.5*(double(T)) + a;
                 a1k           = 0.5*(double(T)-double(o_adjust(k))) + a;
                 bki           = B.slice(k).row(i).t();
                 bbar_ki       = C(k,0) * bki; /* T x 1 */
                 bbar_i.row(k) = bbar_ki.t();
                 B1(k,i)       = 0.5*dot( D.row(k), pow((bki-bbar_ki),2) ); /* no b */
                 logd_dk       = 0; /* set of T gaussian densities for term k */
                 /* dmvn(zro|m,Q.slice(k),true) */
                 for( j = 0; j < T; j++ )
                 {
                    logd_dk   += R::dnorm(0.0,0.0,double(1/sqrt(D(k,j))),true);
                 }
                 logq0ki      = logd_dk + lgamma(a1k) + a*log(b) -
                                   lgamma(a) - a1k*trunc_log(B1(k,i)+b);
                 logq0i       += logq0ki;
            } /* end loop k over iGMRF terms */
            q0i = trunc_exp(logq0i);

            // construct posterior sampling weights to sample s(i)
            colvec weights(M+1); weights.zeros();
            /* evaluate likelihood under kappa_star(k,i) */
            double lweights_l;
            for(l = 0; l < M; l++) /* cycle through all clusters for s(i) */
            {
                s(i)          = l; /* will compute likelihoods for every cluster */  
                lweights_l = 0; /* hold log densities for K computations */
                for(k = 0; k < K; k++)
                {
                    bki            = B.slice(k).row(i).t();
                    for( j = 0; j < T; j++ )
                    {
                      /* effectively making assignment, s(i) = l */
                      lweights_l   += trunc_log(R::dnorm(bki(j),bbar_i(k,j),
                                    double(1/sqrt(kappa_star(k,l)*D(k,j))),false));
                    } /* end loop j over time index */
                } /* end loop k over iGMRF terms */
                //if(lweights_l < -300){lweights_l = -300;}
                weights(l)          = trunc_exp(lweights_l);
                weights(l)          *= double(Num(s(i)))/(double(N) - 1/ipr(i) + conc);
            } /* end loop l over existing or populated clusters */
            /* M+1 or new component sampled from F_{0} */
            weights(M)              = conc/(double(N) - 1/ipr(i) + conc)*q0i;

            // normalize weights
            sweights = sum(weights);
            if(sweights == 0)
            {
                weights.ones(); weights *= 1/(double(M)+1);
            }
            else
            {
                weights /= sweights;
            }

            // conduct discrete posterior draw for s(j)
            unsigned long MplusOne = M + 1;
            s(i) = rdrawone(weights, MplusOne);

            // if new cluster chosen, generate new location
            if(s(i) == M)
            {
                /* sample posterior of ksi_star[k,m] for 1 (vs. n_m) observation */
                double a_star_k; /* shape for 1 obs */
                double bstar_ki;
                kappa_star.insert_cols(M,1); /* add K vector new location to kappa_star */
                num.insert_rows(M,1);
                Num.insert_rows(M,1);
                for(k = 0; k < K; k++)
                {
                     a_star_k         = 0.5*(double(T) - double(o_adjust(k))) + a; /* shape for 1 obs */
                     bstar_ki         = B1(k,i) + b; /* B1(k,i) is a scalar quadratic product */
                     /*
                     bki              = B.slice(k).row(i).t();
                     bstar_ki         = 0.5*( as_scalar(bki.t()*symmatl(Q.slice(k))*bki) ) + b;
                     */
                     kappa_star(k,M)  = rgamma(1, a_star_k, (1/bstar_ki))[0];
                }
                num(M)   = 1;
                Num(M)   = 1/ipr(i);
                M        = MplusOne;
            }
            else
            {
                num(s(i)) += 1;
                Num(s(i)) += 1/ipr(i);
            }
            
        } /* end loop i for cluster assignment to unit i = 1,...,N */
        END_RCPP
    } /* end function bstep for cluster assignments, s, and computing zb */



// update vector of cluster membership indicators, s(i),....,s(N)
SEXP clusterstep_alt(const cube& B, mat& kappa_star, mat& B1, const uvec& o,
            const cube& Q, ucolvec& s, 
            ucolvec& num, unsigned int& M, double& conc, int a, int b,
            const vec& ipr, colvec& Num)
    {
        BEGIN_RCPP
      
        // sample cluster assignments, s(1), ..., s(N)
        // B = (B_1,...,B_K), where B_k is N x T matrix for iGMRF term k
        // Q = (Q_1,...,Q_K), where Q_k is a T x T de-scaled iGMRF precision matrix
        // C = (C_1,...,C_K), where C_k = D_k^-1 * Omega_k, 
        // where Omega_k is the T x T adjacency matrix for iGMRF term, k
        // D is a K x T matrix where row k contains T diagonal elements of Q_k
        // K x M matrix, kappa_star records locations for each iGMRF term
        // o = (o_1,...,o_k) is a vector where each entry denotes the order of term K.
        // e.g. RW(1) -> o = 2, RW(2) -> o = 3, seas(3) -> o = 3
        int N = B.slice(0).n_rows;
        int T = B.slice(0).n_cols;
        uvec df = T - o; /* rank of each slice of Q */
        int K = Q.n_slices;
        double sweights = 0;
        // zro is the zeros.T vector 
        colvec zro(T); zro.zeros();
        // capture quadratic product for rate kernel of posterior gamma
        // posterior for kappa_star(k,i).  
        // save B1 to latter (in another function) compute posterior for kappa_star
        // mat B1(K,N); 
        double a1k; /* posterior shape for kappa_star(k,i) under 1 obs */
        B1.zeros();
        int i, k;
	      unsigned int l;
        
        for(i = 0; i < N; i++)
        {
            // check if _i assigned to singleton cluster
            // if so, remove the cluster associated to _i
            // and decrement the cluster labels for m > s(i)
            if(num(s(i)) == 1) /* remove singleton cluster */
            {
                kappa_star.shed_col(s(i));
                num.shed_row(s(i));
                Num.shed_row(s(i));
                M -= 1; /* decrement cluster count */

                //decrement cluster tracking values by 1 for tossed cluster
                s( find(s > s(i)) )          -= 1;
                
            } /* end cluster accounting adjustment for singleton cluster */
            else /* cluster contains more than one unit */
            {
                num(s(i))                    -= 1;
                /* scale up num to population totals, Num, based on H-T inverse probability estimator */
                Num(s(i))                    -= 1/ipr(i);
            } /* decrement non-singleton cluster count by one */
            
            // construct normalization constant, q0i, to sample s(i)
            // build loqq0i and exponentiate
            colvec b_ki(T), zro(T); zro.zeros();
            long double logq0ki, logq0i = 0, q0i = 0;
            // accumulate weight, q0i, for s(i) over K iGMRF terms  
            for( k = 0; k < K; k++)
            {
                 logq0ki       = 0; /* reset k-indexed log-like on each k */
                 a1k           = 0.5*(double(T)-double(o(k))) + a;
                 b_ki          = B.slice(k).row(i).t();
                 B1(k,i)       = 0.5*( as_scalar(b_ki.t()*symmatl(Q.slice(k))*b_ki) ); /* no b */
                 /* dmvn(zro|m,Q.slice(k),true) */
                 logq0ki      = loggmrfdens(zro,zro,Q.slice(k),df(k),1) + lgamma(a1k) +
                                   a*log(b) - lgamma(a) - a1k*log(B1(k,i)+b);
                 logq0i       += logq0ki;
            } /* end loop k over iGMRF terms */
            q0i = exp(logq0i);
            // construct posterior sampling weights to sample s(i)
            colvec weights(M+1); weights.zeros();
            /* evaluate likelihood under kappa_star(k,i) */
            long double lweights_l;
            for(l = 0; l < M; l++) /* cycle through all clusters for s(i) */
            {
                s(i)          = l; /* will compute likelihoods for every cluster */
                lweights_l    = 0; /* hold log densities for K computations */
                for(k = 0; k < K; k++)
                {
                    b_ki            = B.slice(k).row(i).t();
                    lweights_l     += loggmrfdens( b_ki, zro, Q.slice(k), df(k),
                                                  kappa_star(k,l) );
                } /* end loop k over iGMRF terms */
                weights(l)          = exp(lweights_l);
                weights(l)          *= double(Num(l))/(double(N) - (1/ipr(i)) + conc);
            } /* end loop l over existing or populated clusters */
            /* M+1 or new component sampled from F_{0} */
            weights(M)              = conc/(double(N) - (1/ipr(i)) + conc)*q0i;

            // normalize weights
            sweights = sum(weights);
            if(sweights == 0)
            {
                weights.ones(); weights *= 1/(double(M)+1);
                //weights += 1e-4;
                //weights /= sum(weights);
            }
            else
            {
                weights /=  sweights;
            }
            
            // conduct discrete posterior draw for s(j)
            unsigned long MplusOne = M + 1;
            s(i) = rdrawone(weights, MplusOne);
            // if new cluster chosen, generate new location
            if(s(i) == M)
            {
                /* sample posterior of ksi_star[k,m] for 1 (vs. n_m) observation */
                double a_star_k; /* shape for 1 obs */
                double bstar_ki;
                kappa_star.insert_cols(M,1); /* add K vector new location to kappa_star */
                num.insert_rows(M,1);
                Num.insert_rows(M,1);
                for(k = 0; k < K; k++)
                {
                     a_star_k         = 0.5*(double(T) - double(o(k))) + a; /* shape for 1 obs */
                     bstar_ki         = B1(k,i) + b; /* B1(k,i) is a scalar quadratic product */
                     kappa_star(k,M)  = rgamma(1, a_star_k, (1/bstar_ki))[0];
                }
                num(M)   = 1; 
                Num(M)   = 1/ipr(i);
                M        = MplusOne;
            }
            else
            {
                num(s(i)) += 1;
                Num(s(i)) += 1/ipr(i);
            }
            
        } /* end loop i for cluster assignment to unit i = 1,...,N */
        END_RCPP
    } /* end function bstep for cluster assignments, s, and computing zb */
    
SEXP move_kappastar(mat& kappa_star, const mat& B1, const ucolvec& s, uvec& o, 
                    int T, int a, int b, const vec& ipr)
{
     BEGIN_RCPP
     // K x M matrix, kappa_star, records location values 
     // K x N matrix, B1 contains quadratic products used for rate parameter
     // on posterior for kappa_star.  
     int K = kappa_star.n_rows; /* number of iGMRF terms */
     int M = kappa_star.n_cols; /* number of clusters */
     int N = B1.n_cols; /* number of units */
     int k, m; 
     double num_m = 0;
     uvec pos_m;
     uvec o_adjust = o; 
     //o_adjust.zeros();
     long double a1_mk; /* posterior shape parameter */
     long double b1_mk; /* posterior rate parameter */
     rowvec b1_k(N); b1_k.zeros();
     for(m = 0; m < M; m++)
     {
          pos_m          = find( s == m ); /* s is vector length N */
          //num_m          = pos_m.n_elem;
          num_m          = sum( 1/ipr(pos_m) );
          /* sample posterior for kappa_star(k,) for each iGMRF term, k = 1,...,K */
          for( k = 0; k < K; k++ )
          {
               b1_k                = B1.row(k); /* rowvec(N) */
               b1_mk               = sum(b1_k(pos_m)/ipr(pos_m)) + b; /* weighted totals for H-T adjustment */
               a1_mk               = 0.5*num_m*(double(T)-double(o_adjust(k))) + a;
               kappa_star(k,m)     = rgamma(1, a1_mk, (1/b1_mk))[0];
          } /* end loop over K iGMRF terms */ 
          
     } /* end loop m over clusters */
     END_RCPP
} /* end function move_kappastar() to sample cluster locations */

SEXP move_kappastar_alt(mat& kappa_star, const cube& B, const cube& Q, 
                    const ucolvec& s, uvec& o, 
                    int T, int a, int b, const vec& ipr)
{
     BEGIN_RCPP
     // K x M matrix, kappa_star, records location values 
     // N x T matrices, B_1,...,B_K contains the de-noised functions
     // on posterior for kappa_star.  
     int K = kappa_star.n_rows; /* number of iGMRF terms */
     int M = kappa_star.n_cols; /* number of clusters */
     int k = 0, m = 0, i = 0, count_m;
     double num_m;
     uvec pos_m;
     double a1_mk; /* posterior shape parameter */
     double b1_mk; /* posterior rate parameter */
     colvec b_ki(T);
     for(m = 0; m < M; m++)
     {
          pos_m          = find( s == m ); /* s is vector length N */
          count_m        = pos_m.n_elem;
          num_m          = sum( 1/ipr(pos_m) );
          /* sample posterior for kappa_star(k,) for each iGMRF term, k = 1,...,K */
          for( k = 0; k < K; k++ )
          {
               b1_mk           = 0;
               for( i = 0; i < count_m; i++ )
               {
                    b_ki                  = B.slice(k).row(pos_m(i)).t(); /* T x 1 */
                    b1_mk                 += 0.5*( as_scalar(b_ki.t()*symmatl(Q.slice(k))*b_ki)
                                                       / ipr(pos_m(i)) );       
               } /* end loop i over num_m weighted units in cluster m */
               b1_mk               += b;
               a1_mk               = 0.5*num_m*(double(T)-double(o(k))) + a;
               kappa_star(k,m)     = rgamma(1, a1_mk, (1/b1_mk))[0];
          } /* end loop over K iGMRF terms */ 
          
     } /* end loop m over clusters */
     /* add a bumper */
     END_RCPP
} /* end function move_kappastar() to sample cluster locations */

SEXP move_B(const mat& y, cube& B, const mat& kappa_star, const field<mat>& C,
            //const field<sp_mat>& C,
            mat& gamma, const mat& D, const ucolvec& s, double tau_e)
{
     BEGIN_RCPP
     // N x T matrix of standardized counts, y
     // B = (B_1,...,B_K), where B_k is N x T matrix for iGMRF term k
     // N x T matrix, gamma = sum_{k=1^K}(B_k)
     // K x T, D, where row k contains T x 1 diagonal elements of Q_k
     // sample cluster assignments, s(1), ..., s(N)  
     // Q = (Q_1,...,Q_K), where Q_k is a T x T de-scaled iGMRF precision matrix
     // C = (C_1,...,C_K), where C_k = D_k^-1 * Omega_k, 
     // where Omega_k is the T x T adjacency matrix for iGMRF term, k
     // D is a K x T matrix where row k contains T diagonal elements of Q_k
     // K x M matrix, kappa_star records locations for each iGMRF term 
     int K     = B.n_slices;
     int N     = y.n_rows;
     int T     = D.n_cols;
     colvec bbar_ki(T); bbar_ki.zeros();
     rowvec gammatilde_ki(T); gammatilde_ki.zeros();
     rowvec ytilde_ki(T); ytilde_ki.zeros();
     rowvec d_k(T);
     vec zro(T); zro.zeros();
     double e_ij, phi_ij, bkij, h_ij;
     int k = 0, i = 0, j = 0; /* loop variables */
     for( k = 0; k < K; k++ ) /* over iGMRF terms */
     {
          d_k                           = D.row(k);
          for( i = 0; i < N; i++ ) /* over units */
          {
               /* take out T x 1, b_ki, to be sampled from gamma */
               //gammatilde_ki        = gamma.row(i) - B.slice(k).row(i); /* 1 x T */
               //ytilde_ki            = y.row(i) - gammatilde_ki;
               gammatilde_ki            = gamma.row(i); /* 1 x T */
               // mean of univariate iGMRF, b_kij = 1/d_kj * (omega_kj(-j) * b_ki(-j))
               bbar_ki                  = C(k,0) * B.slice(k).row(i).t(); /* T x 1 */
               for( j = 0; j < T; j++ ) /* over time points */
               {
                    gammatilde_ki(j)    -= B.slice(k)(i,j);
                    ytilde_ki(j)        = y(i,j) - gammatilde_ki(j);
                    // mean of univariate iGMRF, b_kij = 1/d_kj * (omega_kj(-j) * b_ki(-j))
                    B.slice(k)(i,j)     = 0;
                    //bbar_ki(j)          = dot( C(k,0).row(j), B.slice(k).row(i) );
                    e_ij                = tau_e*ytilde_ki(j) 
                                             + d_k(j)*kappa_star(k,s(i)) * bbar_ki(j);
                    phi_ij              = tau_e + d_k(j)*kappa_star(k,s(i));
                    h_ij                = (e_ij / phi_ij);
                    bkij                = rnorm( 1, (h_ij), sqrt(1/phi_ij) )[0];
                    B.slice(k)(i,j)     = bkij;
                    /* put back new sampled values for b_kij */
                    gammatilde_ki(j)    += B.slice(k)(i,j);
               } /* end loop j over time points */
               /* put back new sampled values for b_ki */
               //gamma.row(i)         = gammatilde_ki + B.slice(k).row(i); 
               gamma.row(i)        = gammatilde_ki;
          } /* end loop i over units */
     } /* end loop k over iGMRF terms */
     END_RCPP
} /* end function ustep to sample B_1,...,B_K */


SEXP move_B_alt(const mat& y, cube& B, const mat& kappa_star, const field<mat>& C, 
            //const field<sp_mat>& C,
            mat& gamma, const mat& D, const ucolvec& s, double tau_e)
{
     BEGIN_RCPP
     // N x T matrix of standardized counts, y
     // B = (B_1,...,B_K), where B_k is N x T matrix for iGMRF term k
     // N x T matrix, gamma = sum_{k=1^K}(B_k)
     // K x T, D, where row k contains T x 1 diagonal elements of Q_k
     // sample cluster assignments, s(1), ..., s(N)  
     // Q = (Q_1,...,Q_K), where Q_k is a T x T de-scaled iGMRF precision matrix
     // C = (C_1,...,C_K), where C_k = D_k^-1 * Omega_k, 
     // where Omega_k is the T x T adjacency matrix for iGMRF term, k
     // D is a K x T matrix where row k contains T diagonal elements of Q_k
     // K x M matrix, kappa_star records locations for each iGMRF term 
     int K     = B.n_slices;
     int N     = y.n_rows;
     int T     = D.n_cols;
     double bbar_kij = 0; 
     rowvec gammatilde_ki(T); gammatilde_ki.zeros();
     rowvec ytilde_ki(T); ytilde_ki.zeros();
     rowvec d_k(T);
     double e_ij, phi_ij, bkij, h_ij;
     int k = 0, i = 0, j = 0; /* loop variables */
     for( k = 0; k < K; k++ ) /* over iGMRF terms */
     {
          d_k     = D.row(k);
          for( i = 0; i < N; i++ ) /* over units */
          {
               /* take out T x 1, b_ki, to be sampled from gamma */
               gammatilde_ki        = gamma.row(i); /* 1 x T */
               for( j = 0; j < T; j++ ) /* over time points */
               {
                    gammatilde_ki(j)    -= B.slice(k)(i,j);
                    ytilde_ki(j)        = y(i,j) - gammatilde_ki(j);
                    // mean of univariate iGMRF, b_kij = 1/d_kj * (omega_kj(-j) * b_ki(-j))
                    B.slice(k)(i,j)     = 0;
                    bbar_kij            = as_scalar(dot( C(k,0).row(j), B.slice(k).row(i) ));
                    e_ij                = tau_e*ytilde_ki(j) 
                                             + d_k(j)*kappa_star(k,s(i)) * bbar_kij;
                    phi_ij              = tau_e + d_k(j)*kappa_star(k,s(i));
                    h_ij                = (e_ij / phi_ij);
                    bkij                = rnorm( 1, (h_ij), sqrt(1/phi_ij) )[0];
                    B.slice(k)(i,j)     = bkij;
                    /* put back new sampled values for b_kij */
                    gammatilde_ki(j)    += B.slice(k)(i,j);
               } /* end loop j over time points */
               /* reset gamma.row(i) */
               gamma.row(i)         = gammatilde_ki;
          } /* end loop i over units */
     } /* end loop k over iGMRF terms */
     END_RCPP
} /* end function ustep to sample B_1,...,B_K */

SEXP move_taue(const mat& y, const mat& gamma, double& tau_e, double a, double b, const vec& ipr)
{
     BEGIN_RCPP
     // N x T matrix, gamma, includes de-noised functions
     // N x T matrix, y, holds the data for observation i at time t
     double N            = sum( 1/ipr );
     int T               = y.n_cols;
     //mat Ipr             = repmat(ipr,1,T);
     /* sum over columns (0) - time points - and then over rows (1) - units - weighting each unit */
     double b_1          = 0.5*sum( sum(pow( (y - gamma),2 ),1)/ipr ) + b; /* element-wise division */
     double a_1          = 0.5*(double(N)*double(T)) + a;
     tau_e               = rgamma(1, a_1, (1/b_1))[0];
     END_RCPP
} /* end function move_taue() to sample the global precision parameter */

SEXP move_taue_jitter(const mat& y, const mat& gamma, double& tau_e, double a, double b, double jitter,
                      const vec& ipr)
{
     BEGIN_RCPP
     // N x T matrix, gamma, includes de-noised functions
     // N x T matrix, y, holds the data for observation i at time t
     //int N               = y.n_rows;
     int N               = sum( 1/ipr );
     int T               = y.n_cols;
     //mat Ipr             = repmat(ipr,1,T);
     /* sum over columns (0) - time points - and then over rows (1) - units - weighting each unit */
     double b_1          = 0.5*sum( sum(pow( (y - gamma),2 ),1)/ipr ) + b - jitter;
     double a_1          = 0.5*(double(N)*double(T)) + a - jitter;
     //if( b_1 > 3.0*a_1 ){b_1 = 3.0*a_1;} /* add a bumper */
     tau_e               = rgamma(1, a_1, (1/b_1))[0];
     END_RCPP
} /* end function move_taue() to sample the global precision parameter */

SEXP miss_ystep(mat& y_rep, const mat& y, const mat& gamma, double tau_e)
{
     BEGIN_RCPP
     // N x T matrix, gamma, includes de-noised functions
     // N x T matrix, y, holds the data for observation i at time t
     int N               = y.n_rows;
     int i = 0, j = 0, n_posi = 0;
     uvec pos_i; pos_i.zeros();
     for( i = 0; i < N; i++ )
     {
          /* find missing elements in row(i) of y */
          pos_i               = find( y.row(i) == -9 ); /* y are the data, including missing vals */
          /* n_posi <= T, the number of time points */
          n_posi              = pos_i.n_elem;
          for( j = 0; j < n_posi; j++ )
          {
               /* y_rep is the completed data matrix, with missing vals of y filled in */
               y_rep(i,pos_i(j))  = rnorm( 1, gamma(i,pos_i(j)), sqrt(1/tau_e) )[0];
          }
     }
     END_RCPP
} /* end function miss_ystep() to sample missing data values */
