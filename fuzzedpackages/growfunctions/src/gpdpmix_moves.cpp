#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

SEXP gen_P(uvec& gp_mod, uvec& P_vec)
{
     BEGIN_RCPP
     int L     = gp_mod.n_elem;
     int l = 0;
     P_vec.zeros();
     for(l = 0; l < L; l++)
     {
          if( gp_mod(l) == 1 ) /* RQ */
          {
               P_vec(l)       = 3;
          }else{
               if( gp_mod(l) == 2 ) /* SE */
               {
                    P_vec(l)  = 2;
               }else{ /* Seasonal */
                         P_vec(l)  = 3;
               }/* end conditional statement on term is seasonal */
          }/*end conditional statement whether term is RQ */
     }/* end loop l over the terms */
     END_RCPP
} /* end function gen_P() */

mat gen_C(const colvec& thetastar_m, double tau_e, const mat& Omega_t, const cube& Omega_s,
          double jitter, uvec& gp_mod, uvec& n_parms, uvec& pos_s, int noise)
{
     // generates T x T GP covariance matrix, C, without global noise parameter, tau_e
     // use this where have to co-sample the GP function, z_{1},..,z_{N}
     // gp_mod indexes the GP covariance functional form for each of L terms
     // n_parms is an (L+1) vector indexes the starting position for each covariance term
     // pos_s is an L vector that indexes the seasonal covariance function in 1,..,L_s
     // - where pos_s(j) returns the correct term in 1,...,L_s if term j is seasonal
     // L = gp_mod.n_elem = n_parms.n_elem
     // Same Euclidean distance matrix, Omega_t, used for all trend terms
     int L     = gp_mod.n_elem;
     int T     = Omega_t.n_rows;
     mat C_m(T,T); C_m.zeros();
     mat eye_T = eye(T,T);
     int j = 0; /* loop variable */
     //int iter = 0;
     
     // trend covariance 
     for( j = 0; j < L; j++)
     {
          if( gp_mod(j) == 1 ) /* rational quadratic */
          {
               // theta_star[p,m] -> col m, p=0 := kappa; p=1 := rho; p=2 := alpha
               // T x T gp rational quadratic covariance matrix, C
               C_m  +=  (1/thetastar_m((0+n_parms(j)))) * 
                         pow( (1+(Omega_t/(thetastar_m((1+n_parms(j)))*thetastar_m((2+n_parms(j)))))), 
                         -thetastar_m((2+n_parms(j))) );
          }else{ /* squared exponential */
                    if( gp_mod(j) == 2 )
                    {
                         C_m  += (1/thetastar_m((0+n_parms(j)))) * exp(-(1/thetastar_m((1+n_parms(j))))
                                                                 *Omega_t);
                    }else{
                         if( gp_mod(j) == 3 ) /*seasonal*/
                         {
                              // p=3 := kappa_s; p=4 := rho_s; p=5 := rho_t
                              // memo: each seasonal term, j, has its own cov mat, */
                              /* Omega_s.slice(j) */
                              C_m  += (1/thetastar_m(0+n_parms(j))) *
                                        exp(-(1/thetastar_m(1+n_parms(j)))*Omega_s.slice(pos_s(j)) -
                                        (1/thetastar_m(2+n_parms(j)))*Omega_t);
                         } /*end conditional statement on whether j is seasonal covariance */
                    } /* end conditional on whether j is an SE covariance */
                    
          } /* end conditional statement on whether covariance trend term j is rational quadratic */            
     } /* end loop j over n_trend additive trend covariance terms */
     
     /* add noise */
     if( noise > 0 )
     {
          C_m  += (1/(tau_e))*eye_T;
     }
     
     /* add jitter */
     C_m            +=  jitter*eye_T;
     mat(C_ms)      = symmatl(C_m);
     //while( (cond(C_ms) > 1e6) & (iter < 5) )
     //{
     //     C_ms    += jitter*(iter+1)*eye_T;
     //     iter   += 1;
     //}
     
     return C_ms;
} /* end function gen_C() */

mat gen_Cterm(const colvec& thetastar_m, double tau_e, const mat& Omega_t, 
               const mat& Omega_s, double jitter, 
               int gp_mod_term, int noise)
{
        // generate T x T GP covariance matrix from a single formulation
        // gp_mod = 1 -> RQ, gp_mod = 2 -> SE, else gp_mod -> Seasonality
        // thetastar_m = theta_star.col(m) - P x 1
        int T            = Omega_t.n_rows;
        mat C_m(T,T); C_m.zeros();
        if( gp_mod_term == 1 ) /* rational quadratic */
        {
             // theta_star[3,m] -> row 1 := kappa; row 2 := rho; row 3 := alpha
             // T x T gp rational quadratic covariance matrix, C
                    C_m   = (1/thetastar_m(0)) * pow( (1 + (Omega_t/(thetastar_m(1)*thetastar_m(2)))), 
                                        -thetastar_m(2) );
        }else{
             if( gp_mod_term == 2 ) /* squared exponential */
             {
                  // theta_star[0,m] = kappa, theta_star[1,m] = rho
                  C_m    = (1/thetastar_m(0)) * exp(-(1/thetastar_m(1))*Omega_t);
             }else{ /* seasonal */
                  
                  C_m    = (1/thetastar_m(0)) *
                              exp(-(1/thetastar_m(1))*Omega_s -
                                   (1/thetastar_m(2))*Omega_t);
             } /* end computation of seasonal C_m */     
        } /* end conditional statement on whether C_m is 1. RQ, 2. SE, 2. seasonal */
        
        /* add noise */
        if( noise > 0 )
        {
               C_m  += (1/(tau_e))*eye(T,T);
        }
        
        /* add jitter */
        C_m         += jitter*eye(T,T);
        
        mat(C_ms)   = symmatl(C_m);
        
        return C_ms;
}

mat gen_Casym(const colvec& theta_i, const mat& Omega_t, const cube& Omega_s,
               uvec& gp_mod, uvec& n_parms, uvec& pos_s)
{
     // generates T_pred x T GP covariance matrix, C, to use for prediction
     // of T_pred x 1 function, f
     // Rest of inputs relate to P x 1 set of sampled GP covariance parameters, theta_i
     // gp_mod indexes the GP covariance functional form for each of L terms
     // n_parms is an (L+1) vector indexes the starting position for each covariance term
     // pos_s is an L vector that indexes the seasonal covariance function in 1,..,L_s
     // - where pos_s(j) returns the correct term in 1,...,L_s if term j is seasonal
     // L = gp_mod.n_elem = n_parms.n_elem
     // Same Euclidean distance matrix, Omega_t, used for all trend terms
     int L          = gp_mod.n_elem;
     int T_pred     = Omega_t.n_rows; /* Omega_t is T_pred x T */
     int T          = Omega_t.n_cols;
     mat C_i(T_pred,T); C_i.zeros();
     int j = 0; /* loop variable */
     
     // trend covariance 
     for( j = 0; j < L; j++)
     {
          if( gp_mod(j) == 1 ) /* rational quadratic */
          {
               // theta_i[p] ->  p=0 := kappa; p=1 := rho; p=2 := alpha
               // T_pred x T gp rational quadratic covariance matrix, C_i
               C_i  +=  (1/theta_i((0+n_parms(j)))) * 
                              pow( (1 + (Omega_t/(theta_i((1+n_parms(j)))*theta_i((2+n_parms(j)))))), 
                              -theta_i((2+n_parms(j))) );
          }else{ /* squared exponential or seasonal */
                    if( gp_mod(j) == 2 )
                    {
                         C_i  += (1/theta_i((0+n_parms(j)))) * exp(-(1/theta_i((1+n_parms(j))))
                                                                 *Omega_t);
                    }else{
                         if( gp_mod(j) == 3 ) /*seasonal*/
                         {
                              // p=3 := kappa_s; p=4 := rho_s; p=5 := rho_t
                              // memo: each seasonal term, j, has its own cov mat, */
                              /* Omega_s.slice(j) */
                              C_i  += (1/theta_i(0+n_parms(j))) *
                                        exp(-(1/theta_i(1+n_parms(j)))*Omega_s.slice(pos_s(j)) -
                                        (1/theta_i(2+n_parms(j)))*Omega_t);
                         } /*end conditional statement on whether j is seasonal covariance */
                    } /* end conditional on whether j is an SE covariance */
                    
          } /* end conditional statement on whether covariance trend term j is rational quadratic */            
     } /* end loop j over n_trend additive trend covariance terms */
     
     return C_i;
} /* end function gen_Casym() */

double logy_like(int i, const mat& y, const colvec& ipr, ucolvec& s, 
                 const cube& U_last)
{
        // Compute log(MVN) density under GP prior
        int T                 = y.n_cols;
        colvec weights        = 1/ipr;
        
        /* gp covariance matrix generation */
        /* u = [U']^-1*y -> U'u = y -> u = solve(U',y) */
        /* y'*C^-1*y = u'u */
        mat U_i               = U_last.slice(s(i));
        /* -0.5 * logdet(C_i) */
        double loglike_i      = -(0.5)*weights(i)*sum(2*log(U_i.diag())); /* for an observation, not a cluster */
        rowvec y_i(T); y_i.zeros();
        colvec u_i(T);
        y_i            = y.row(i);
        u_i            = solve(trimatu(U_i).t(),y_i.t());
        /* -0.5*sum(y_i'C_i^-1y_i) */
        loglike_i      -= 0.5*as_scalar( weights(i) * u_i.t()*u_i );
        
        return loglike_i;
}

mat compute_Um(const colvec& thetastar_m, double tau_e, double jitter,  
                    uvec& gp_mod, uvec& n_parms, uvec& pos_s, const mat& Omega_t, const cube& Omega_s,
                    int noise)
{
     mat C_m                  = gen_C(thetastar_m,tau_e,Omega_t,
                                        Omega_s,jitter,gp_mod,
                                        n_parms, pos_s, noise); /* T x T */
     mat U_m                  = trimatu(chol(symmatl(C_m)));

     return U_m;
}

mat compute_Upm(double thetastar_pm, const mat& theta_star, double tau_e, double jitter, int p, 
                    int m, uvec& gp_mod, uvec& n_parms, uvec& pos_s, const mat& Omega_t, 
                    const cube& Omega_s, int noise)
{
     // set up vector of parameters for cluster m, including
     // the parameter to be sampled, thetastar_pm
     colvec thetastar_m       = theta_star.col(m); /* P x 1 */
     /* place parm for which log-density is needed in thetastar_m */
     thetastar_m(p)           = thetastar_pm;  
     mat C_pm                 = gen_C(thetastar_m,tau_e,Omega_t,
                                        Omega_s,jitter,gp_mod,
                                        n_parms, pos_s, noise); /* T x T */
     mat U_pm                 = trimatu(chol(symmatl(C_pm)));
     
     return U_pm;
}

SEXP compute_U(const mat& theta_star, double tau_e, double jitter, uvec& gp_mod,
                    uvec& n_parms, uvec& pos_s, const mat& Omega_t, const cube& Omega_s,
                    int noise, cube& U)
{
     BEGIN_RCPP
     /* set dimensions */
     int T                 = Omega_t.n_cols;
     int M                 = theta_star.n_cols; /* number of clusters */
     int P                 = theta_star.n_rows; /* number of parms-per-cluster */
     
     /* gp covariance matrix generations */
     int m;
     mat C_star_m(T,T); C_star_m.zeros();
     colvec thetastar_m(P);
     for(m = 0; m < M; m++)
     {
          thetastar_m         = theta_star.col(m);
          C_star_m            = gen_C(thetastar_m,tau_e,Omega_t,
                                        Omega_s,jitter,gp_mod,
                                        n_parms, pos_s, noise); /* T x T */ 
          U.slice(m)          = trimatu(chol(symmatl(C_star_m)));
     }
     END_RCPP
}

double logFm_like(const mat& y, int m, const ucolvec& s, 
                    const mat& U_m, const vec& ipr)
{
        // Compute log(MVN) density under GP prior
        // gp gaussian likelihood
        int T                 = y.n_cols;
        uvec pos_m            = find( s == m ); /* s is vector length N */
        int count_m           = pos_m.n_elem;
        double num_m          = sum( 1/ipr.elem(pos_m) );

        /* u = [U']^-1*y -> U'u = y -> u = solve(U',y) */
        /* y'*C^-1*y = u'u */
        double loglike_m      = -(0.5*num_m)*sum(2*log(U_m.diag()));
        rowvec y_mi(T); y_mi.zeros();
        colvec u_mi(T);
        int i;
        /* form exponential term from set of y_i attached to cluster m */
        for(i = 0; i < count_m; i++)
        {
             y_mi            = y.row(pos_m(i));
             u_mi            = solve(trimatu(U_m).t(),y_mi.t());
             //loglike_m      -= 0.5*as_scalar(y_mi * pinv(symmatl(U_m.t() * U_m)) * y_mi.t());
             loglike_m      -= 0.5*as_scalar( u_mi.t()*u_mi ) / ipr(pos_m(i));
        }
        
        return loglike_m;
}

double logFpm_post(double thetastar_pm, const mat& theta_star, double tau_e, double jitter, int p, 
                    int m, uvec& gp_mod, uvec& n_parms, uvec& pos_s, const mat& Omega_t, 
                    const cube& Omega_s, int noise, const mat& y, const ucolvec& s, double a, double b,
                    const vec& ipr)
{
     /* use this only for uni.slice_pm(), not temper() */
     mat U_pm            = compute_Upm(thetastar_pm, theta_star, tau_e, jitter, p, 
                                        m, gp_mod, n_parms, pos_s, Omega_t,  Omega_s,
                                        noise);
     double loglike_pm   = logFm_like(y, m, s, U_pm, ipr);
     double logpost_pm   = loglike_pm + log_prior(thetastar_pm,a,b);
     return logpost_pm;
} /* end function logFpm_post() to compute log-posterior-density for theta_star(p,m) */


double logFtau_post(const mat& theta_star, double tau_e, double jitter, uvec& gp_mod,
                    uvec& n_parms, uvec& pos_s, const mat& Omega_t, const cube& Omega_s,
                    int noise, const mat& y, const ucolvec& s, double a, double b,
                    const vec& ipr)
{
        // Compute log(MVN) density under GP prior
        int T                 = y.n_cols;
        int M                 = theta_star.n_cols; /* number of clusters */
       
        /* gp covariance matrix generations */
        cube U(T,T,M);
        compute_U(theta_star, tau_e, jitter, gp_mod, n_parms, pos_s,
                    Omega_t, Omega_s, noise, U);
        /* compute GP likelihood */
        double logpost_tau    = logFtau_like(y, s, U, ipr);
        // contribution of log-prior to log-posterior
        logpost_tau          += log_prior(tau_e,a,b);
        
        return logpost_tau;
}

double logFtau_like(const mat& y, const ucolvec& s, const cube& U_last, const vec& ipr)
{
        // Compute log(MVN) density under GP prior
        int T                 = y.n_cols;
        int N                 = y.n_rows;
        int M                 = U_last.n_slices; /* number of clusters */
        // gp gaussian likelihood   
        int i, m;
        double num_m;
        uvec pos_m;

        //cube U_last(T,T,M)
        double loglike_tau = 0;
        for(m = 0; m < M; m++)
        {
             pos_m                 = find( s == m ); /* s is vector length N */
             //num_m                 = pos_m.n_elem;
             num_m                 = sum(1/(ipr.elem(pos_m)));
             loglike_tau           += -(0.5*num_m)*sum(2*log(diagvec(U_last.slice(m))));
        }
            
        /* u = [U']^-1*y -> U'u = y -> u = solve(U',y) */
        /* y'*C^-1*y = u'u */
        colvec u_i(T); 
        mat U_i(T,T);
        for(i = 0; i < N; i++)
        {
             U_i              = U_last.slice(s(i));
             u_i              = solve(trimatu(U_i).t(),y.row(i).t());
             loglike_tau       -= 0.5*as_scalar( u_i.t()*u_i ) / ipr(i);
        }
        
        return loglike_tau;
}

double log_prior(double theta_star_pm, double a, double b)
{
     double log_prior_pm = (a-1)*log(theta_star_pm) - b*theta_star_pm;
     return log_prior_pm;  
} /* end function compute gamma(a,b) log-prior */

SEXP uni_slice_pm(colvec& theta_updown, mat& Pi_n, int pos_ud, int dist,
          const mat& theta_star, double tau_e, int p, int m, int& slice_evals_pm,
          double lower, double upper, int n_slice_iter, const mat& w_tot,
          const ucolvec& s, const mat& Omegat_n, const cube& Omegas_n,
          uvec& gp_mod, uvec& n_parms, uvec& pos_s, int noise, double jitter,
          const mat& y_n, double a, double b, int transition, const vec& ipr)
{
        BEGIN_RCPP
        // univariate slice sampler used for tempered sampling, so 
        // theta_updown and Pi_n are vectors of parameter and density values, respectively
        // in order of sampling during the tempering process.
        // e.g. Employ n = 2 layers:
        // Distributions Pi_1, Pi_2 are progressively more coarse
        // dim(Omega_1) = (T_1,T_1); dim(y_1) = (N,T_1)
        // Pin(0,0) = rowvec(5) - x0^,x1^,0,x1v,x0v - the values on either end
        // Pin(1,0) = rowvec(5) - 0,x1^,x2-,x1v,0 - the values in the middle
        // pos_ud \in (0,..,[2+(2n-1)]-1) is global position for c(x0^,x1^,x2-,x1v,x0v) 
        
        // sampler tuning
        slice_evals_pm   = 0; /* counter for num_evals to get new sampled value */
        /* total number of parameter sets to update. the last is tau_e. */
        int P_aug        = theta_star.n_rows + 1; 
        
        // determine slice level
        double theta_old = theta_updown(pos_ud-1); /* start sampling at x1^ (given x0^) */
        double Pin_old;
        // note: T^ is an up-transition (transition == 1) when pos_ud == 1
        // triggers computation of pi_n for theta_old = x_0^
        if( transition == 0 ) /* within a transition distribution */
        {
             // dist = 0 for fine and 1 for coarse
             Pin_old     = Pi_n(dist,(pos_ud-1));
        }else{ /* moving from one transition distribution to another */
             if( p < (P_aug-1) ) /* index p starts at 0 */
             {    
                    Pin_old     = logFpm_post(theta_old, theta_star, tau_e, jitter, p, m, gp_mod, 
                                             n_parms, pos_s, Omegat_n, Omegas_n, noise,
                                             y_n, s, a, b, ipr);
             }else{ /* update tau_e */
                    Pin_old     = logFtau_post(theta_star, theta_old, jitter, gp_mod, n_parms, pos_s, 
                                                  Omegat_n, Omegas_n, noise,
                                                  y_n, s, a, b, ipr);
             }/* end condition on whether updating theta_star(p,m) or tau_e */
             
             /* save log-density for old value computed at -current- sampling density \in (1,..,n) */
             Pi_n(dist,(pos_ud-1))     = Pin_old;
        } /* end leg of conditional statement where transition == 1, meaning must compute Pin_old */     
        double logy      = Pin_old - rexp(1, 1)[0]; 
        
        // determine initial sampling interval for theta_new
        double u         = runif(1, 0, w_tot(p,m))[0];
        double L         = theta_old - u;
        double R         = theta_old + (w_tot(p,m)-u); /* guarantees theta_old \in (L,R) */
        
        // expand interval until ends are outside the slice
        /* set maximum number of moves for number of steps to expand interval */
        int J            = floor(runif(1,0,n_slice_iter)[0]);
        int K            = (n_slice_iter-1) - J;
         
        /* evaluate widening interval in each step of size w_tot(p) */
        /* lower interval */
        double Pin_L     = 0;  
        while(J > 0)
        {
               /* condition on L to ensure L is in the parameter support */
               if( L <= lower ){break;}
               slice_evals_pm += 1;
               if( p < (P_aug-1) )
               {
                    Pin_L   = logFpm_post(L, theta_star, tau_e, jitter, p, m, gp_mod, 
                                             n_parms, pos_s, Omegat_n, Omegas_n, 
                                             noise, y_n, s, a, b, ipr);
               }else{ /* update tau_e */
                    Pin_L   = logFtau_post(theta_star, L, jitter, gp_mod, 
                                                  n_parms, pos_s, Omegat_n, Omegas_n, 
                                                  noise, y_n, s, a, b, ipr); 
               }/* end condition on whether updating theta_star(p,m) or tau_e */
               /* condition on Fn_L to ensure L is within the slice */
               if( Pin_L <= logy ){break;}
               L -= w_tot(p,m); /* p is passed in.  If p set of P_aug, then uses wtau */
               J -= J-1;
        } /* end condition on maximum number of iterations, J, to check lower bound */
        
        /* evaluate widening interval in each step of size w_tot(p) */
        /* upper interval */
        double Pin_R     = 0;
        while(K > 0)
        {
               /* condition on R to ensure R is in the parameter support */
               if( R >= upper ){break;}
               slice_evals_pm += 1;
               if( p < (P_aug-1) )
               {
                    Pin_R   = logFpm_post(R, theta_star, tau_e, jitter, p, m, gp_mod, 
                                   n_parms, pos_s, Omegat_n, Omegas_n, 
                                   noise, y_n, s, a, b, ipr);
               }else{ /* update tau_e */
                    Pin_R   = logFtau_post(theta_star, R, jitter, gp_mod, 
                                   n_parms, pos_s, Omegat_n, Omegas_n, 
                                   noise, y_n, s, a, b, ipr);
               }/* end condition on whether updating theta_star(p,m) or tau_e */
               /* condition on Fn_R to ensure R is within the slice */
               if( Pin_R <= logy ){break;}
               R += w_tot(p,m); /* p is passed in.  If p set of P_aug, then uses wtau */
               K -= K-1;
        } /* end condition on maximum number of iterations, K, to check upper bound */
        
        // shrink interval to lower and upper bounds
        if(L < lower)
        {
             L = lower;
        }
        
        if(R > upper)
        {
             R = upper;
        }
        
        // sample Fn_new from the interval, shrinking the interval
        // on each rejection
        double theta_new         = 0;
        double Pin_new           = 0;
        /* generate theta_new and accept if Pin_new >= logy */ 
        for(;;) /* infinite loop that contains a break statement */
        {
             /* generate theta_new */
             theta_new   = runif(1, L, R)[0];
             /* evaluate likelihood */
             slice_evals_pm    += 1;
             if( p < (P_aug-1) )
             {
                    Pin_new             = logFpm_post(theta_new, theta_star, tau_e, jitter, p, m, gp_mod, 
                                             n_parms, pos_s, Omegat_n, Omegas_n, 
                                             noise, y_n, s, a, b, ipr);
             }else{ /* update tau_e */
                    Pin_new             = logFtau_post(theta_star, theta_new, jitter, gp_mod, 
                                                  n_parms, pos_s, Omegat_n, Omegas_n, 
                                                  noise, y_n, s, a, b, ipr);
             } /* end condition on whether updating theta_star(p,m) or tau_e */
             
             /* condition on whether new sampled value is outside of slice */
             if( Pin_new >= logy ){break;}
             
             /* shrink interval */
             if(theta_new > theta_old)
             {
                  R      = theta_new;
             }else{ /* theta_new < theta_old */
                  L      = theta_new;
             } /* end condition to shrink upper or lower interval */
             
        } /* end infinite loop on shrinking slice with break when Pin_new >= logy */ 
        
        /* save the NEw values of the parameter and associated log-density */
        theta_updown(pos_ud)            = theta_new;
        Pi_n(dist,pos_ud)               = Pin_new;
        
	   END_RCPP
} /* end unislice() function to produce a slice sample under coarser tempered distributions */

int temper_dist_compute(int i, int n)
{
        // return level of distribution from fine-to-coarse to use for slice update 
        // in a tempered implementation.  selects a distribution \in (1,..,n)
        // "i" indexes the value (\in 1,..,(vals-1)) in the tempered chain, theta_updown, being updated
        // we leave out i = 0 (the first value) because we don't sample it, but bring it inup.
        // "n" denotes the total number of layers
        int dist;
        int vals              = 2 + (2*n-1); /* in_out + (up_down * # layers) - summit  */
        /* number of uni.slice() steps is vals - 1, since bring in_up initial value */
        int steps             = vals - 1; /* number of tempered sampling steps */
        if( i <= n ) /* i starts at "1" as x^_0 brought up from previously sampled value */
        {
               dist = (i-1);
        }else{
               if( i == (n+1) ) /* summit */
               {
                    dist = (n-1); /* if n+1 = 3, dist = 1 */
               }else{ /* i > (n+1) */
                    dist = steps -  i; /* e.g. steps = 4, i = 4 -> dist = 0 */
               }
        } /* end conditional statements on setting distribution to sample */
        
        return(dist);
} /* end function temper_dist_compute() to select the dist for each tempered update */


ucolvec temper_dist_alt(int n)
{
        // return level of distribution from fine-to-coarse to use for slice update 
        // in a tempered implementation.  selects a distribution \in (1,..,n)
        // "i" indexes the value (\in 1,..,(vals-1)) in the tempered chain, theta_updown, being updated
        // we leave out i = 0 (the first value) because we don't sample it, but bring it inup.
        // "n" denotes the total number of layers
        int vals                   = 2 + (2*n-1); /* in_out + (up_down * # layers) - summit  */
        /* number of uni.slice() steps is vals - 1, since bring in_up initial value */
        int steps                  = vals - 1; /* number of tempered sampling steps */
        /* build symmetric - harmonic vector of distribution indicators based on n */
        uvec desc                  = ones<uvec>(n);
        desc                       = cumsum(desc) - 1; /* e.g. 0,1,..,(n-1) */
        uvec climb                 = sort( desc, "descend" ); /* (n-1),(n-2),...,0 */
        ucolvec dist(steps); /* 0,1,...,(n-1),(n-1),...,1,0 */
        dist.subvec(0,(n-1))       = desc; 
        if( n > 1 ) /* if n == 1, dist = 0. */
        {
             dist.subvec(n,(steps-1))   = climb;
        }
        
        return(dist);
} /* end function temper_dist_alt() to build vector of coarse distributions to sample */

double temper_logpmove_compute(mat& Pi_n)
{
        // compute portion of tempered probability of move within the transformed space
        // "i" indexes the value in the tempered chain, theta_updown, being updated
        // "n" denotes the total number of layers
        int n                 = Pi_n.n_rows;
        int vals              = 2 + (2*n-1); /* in_out + (up_down * # layers) - summit  */
        /* number of uni.slice() steps is vals - 1, since bring in_up initial value */
        //int steps             = vals - 1; /* number of tempered sampling steps */
        // "Pi_n" is an n x vals matrix of log-posterior densities used for slice updates
        
        // accept reject step for x0v, starting from x0^ 
        double logf_pm        = 0;
        int i;
        for(i = 0; i < n; i++) /* loop to add in transition probabilities */
        {
               logf_pm      += Pi_n(i,i); /* numerator: from left, incl center */
               logf_pm      -= Pi_n(i,(vals-1)-i);  /* denom: from right, incl center */
               if(i < (n-1))
               {
                    logf_pm   += Pi_n(i,((vals-1)-(i+1))); /* numerator: from right */
                    logf_pm   -= Pi_n(i,(i+1)); /*denom: from left */
               }         
        } /* end loop to add in transition probabilities in expanded space */
        
        return logf_pm;
        
} /* end function temper_logpmove_compute() to compute probability of move in the transformed space */

double temper_logpmove_alt(mat& Pi_n, const ucolvec& dist)
{
        // compute portion of tempered probability of move within the temporary space
        // "i" indexes the value in the tempered chain, theta_updown, being updated
        // "n" denotes the total number of layers
        int n                 = Pi_n.n_rows;
        int vals              = 2 + (2*n-1); /* in_out + (up_down * # layers) - summit  */
        /* number of uni.slice() steps is vals - 1, since bring in_up initial value */
        /* "updates" denotes the number of update steps for each distribution */
        int updates           = vals - 2; /* toss center and bring one up or take one down */
        /* start with harmonic dist indicators, e.g. 0,1,2,2,1,0 and delete entry n */
        /* because don't include log density for summit value. e.g. \bar{x}_{2} for n = 2 */
        ucolvec dist_move = dist;
        dist_move.shed_row(n); /* 0,1...,(n-1), (n-2),...,0 */
        // "Pi_n" is an n x vals matrix of log-posterior densities used for slice updates
        
        /* define objects to select row indices (dist) and column indices (desc, climb, numer, denom) */
        uvec update_set(updates), numer(updates), denom(updates);
        
        /* compose column indices of Pi_n to compute logf_pm */
        mat Pi_n_center            = Pi_n;
        Pi_n_center.shed_col(n); /* get rid of center, which isn't included in log_probs */
        update_set                 = ones<uvec>(updates); 
        denom                      = cumsum( update_set ); /* e.g. = c(1,2,3) for n = 2 */
        numer                      = denom - 1; /* e.g. = c(0,1,2) for n = 2 */ 
        
        // accept reject step for x0v, starting from x0^ 
        double logf_pm             = sum(diagvec( Pi_n_center( dist_move, numer ) )) 
                                        - sum(diagvec( Pi_n_center( dist_move, denom ) ));
     
        return logf_pm;
        
} /* end function temper_logpmove_compute() to compute probability of move in the transformed space */

SEXP temper(cube& U_last, mat& theta_star, double& tau_e, uvec& gp_mod, uvec& n_parms, 
          uvec& pos_s, double& slice_levals_theta, double& slice_levals_tau, 
          double& n_slice_theta, double& n_slice_tau, double& accept_temper, 
          double& n_temper_evals, double lower, double upper, 
          int n_slice_iter, const mat& w_tot, const ucolvec& s, const field<mat>& Omegat_ns, 
          const field<cube>& Omegas_ns, const mat& Omega_t, const cube& Omega_s, 
          const field<mat>& y_ns, const mat& y, int noise, double jitter, double a, double b,
          double atau, double btau, const vec& ipr)
{
        BEGIN_RCPP
        // number of layers
        int n                 = Omegat_ns.n_rows;
        // dimension of GP covariance
        int T                 = Omega_t.n_rows;
        // number of parameters to sample
        int M                 = theta_star.n_cols;
        int P                 = theta_star.n_rows;
        /* need P_aug to trigger update of tau_e in uni.slice_pm() */
        int P_aug             = P + 1; /* last parameter updated is tau_e */
        /* intermediate sampled values for each theta_star(p,m) */
        int vals              = 2 + (2*n-1); /* in_out + (up_down * # layers) - summit  */
        /* number of uni.slice() steps is vals - 1, since bring in_up initial value */
        int steps             = vals - 1; /* number of tempered sampling steps */
        /* tempered sample and density values */
        colvec theta_updown(vals); mat Pi_n(n,vals);  /* note: 1 more than "steps" */
        // every step is a transition, except step n+1 where the most coarse distribution
        // is sampled twice.  The ladder is 1,2,...,n,n,...,2,1
        colvec transition(steps); transition.ones(); transition(n) = 0;
        /* distribution indicator to use for / update in sampling steps */
        ucolvec dist             = temper_dist_alt(n); /* for n = 2, dist order is 0 1 1 0 */
        /* loop counters */
        int i = 0, p = 0, m = 0;
        // tempered transition probabilities and M-H acceptance step for theta_star(p,m)
        double logf_pm = 0, thresh_pm = 0;
        double loglike_m_old = 0, logpi_pmold = 0, loglike_pmnew = 0, logpi_pmnew = 0;
        // same for tau_e
        double logpitau_old = 0, logf_tau = 0, thresh_tau = 0, logpitau_new = 0;
        // update U_last when proposing new tau_e
        mat Upm_new(T,T); cube U_new(T,T,M);
        // capture average number of slice sampling evals per tempered step for each theta_star(p,m)
        int slice_evals_pm         = 0; /* counter for each (p,m) + tau_e */
        slice_levals_theta         = 0; /* numberof slice sampling like evals across (p,m) */
        slice_levals_tau           = 0; /* separate counter for tau_e */
        n_slice_theta              = (P*M*steps); /* number of (p,m) slice steps */
        n_slice_tau                = steps; /* average num of slice lik updates per slice sampling */
        accept_temper              = 0; /* global step-up acceptance yes = 1, no = 0 */ 
        n_temper_evals             = 0; /* total number of temper() proposals */
        // n_slice_iter = 10000; /* number of steps to compose (widen) interval / slice */
                                                 
        // sample theta_star[p,m], each in an expanded space
        for(m = 0; m < M; m++) /* sample clusters, m */
        {
               // update log-like for theta_star(,m) since share same observations
               // within cluster can avoid computation of old likelihood
               // log-priors are different for P parameters in cluster m, so have to 
               // compute those.
               /* theta_star[,m] \perp theta_star[,m'] */
               loglike_m_old  = logFm_like(y, m, s, U_last.slice(m), ipr); 
               for(p = 0; p < P; p++) /* sample parameter types, p */
               {
                    /* tick up temper_evals */
                    n_temper_evals   += 1;
                    // perform tempered sampling for theta_star[p,m]
                    theta_updown.zeros(); /* for n = 2, theta_updown is of length 5 */
                    theta_updown(0) = theta_star(p,m);
                    Pi_n.zeros();
                    /* note: updates in theta_updown start at 1 bc bring in previous theta_star(p,m) */
                    for(i = 1; i < vals; i++) /* emphasize that update uni.slice from 1,..,steps */
                    {
                         /* note: theta_updown and Pi_n are of length vals = steps + 1 = 5 */
                         // univariate slice sampler
                         // returns vector, theta_updown and matrix, Pi_n
                         uni_slice_pm(theta_updown, Pi_n, i, dist(i-1),
                              theta_star, tau_e, p, m, slice_evals_pm, lower, upper, n_slice_iter, 
                              w_tot, s, Omegat_ns(dist(i-1),0),  Omegas_ns(dist(i-1),0),
                              gp_mod, n_parms, pos_s, noise, jitter,
                              y_ns(dist(i-1),0), 
                              a, b, transition(i-1), ipr); 
                         /* memo; the "i-1" in transition reflects c++ start count at 0 for vector */
                         slice_levals_theta  += slice_evals_pm;
           
                    } /* end loop over tempered (lower rank) dists to sample theta_star[p,m] */
                  
                    // accept reject step for x0v, starting from x0^ 
                    // compute portion of probability of move only relating to transformed space */
                    logf_pm     = temper_logpmove_alt(Pi_n, dist);
                  
                    // compute probability of move
                    /* add in log-posterior under fine (full) density in original space */
                    /* full data with (Omega,y) of size T */
                    logpi_pmold         = loglike_m_old + log_prior(theta_star(p,m),a,b);/* old value */
                    Upm_new             = compute_Upm(theta_updown((vals-1)), theta_star, tau_e, jitter, 
                                             p, m, gp_mod, n_parms, pos_s,
                                             Omega_t, Omega_s, noise);
                    loglike_pmnew       = logFm_like(y, m, s, Upm_new, ipr);
                    logpi_pmnew         = loglike_pmnew + log_prior(theta_updown((vals-1)),a,b);
                    logf_pm             += logpi_pmnew; /* proposed value */
                    logf_pm             -= logpi_pmold;
                  
                    /* acceptance step */
                    thresh_pm          = runif(1,0,1)[0];
                    if( log(thresh_pm) < logf_pm )
                    {
                         theta_star(p,m)         = theta_updown(vals-1); /* end value of tempered steps */
                         loglike_m_old           = loglike_pmnew;
                         U_last.slice(m)         = Upm_new;
                         accept_temper           += 1;
                    } /* end M-H acceptance step */        
               } /* end sampling loop p over GP parameters within clusters */       
          }/* end sampling loop m over clusters */
          /* end updating of theta_star(p,m) */ 
             
          /* update TAU_E */
          /* tick up temper_evals */
          n_temper_evals   += 1;
          // perform tempered sampling for tau_e
          theta_updown.zeros(); /* for n = 2, theta_updown is of length 5 */
          theta_updown[0] = tau_e;
          Pi_n.zeros();
          for(i = 1; i < vals; i++) /* emphasize that update uni.slice from 1,..,steps */
          {
               /* note: theta_updown and Pi_n are of length vals = steps + 1 = 5 */
               // univariate slice sampler
               // returns vector, theta_updown and matrix, Pi_n
               slice_evals_pm = 0; /* start over for tau_e */
               /* p = P_aug triggers update for tau_e */
               uni_slice_pm(theta_updown, Pi_n, i, dist(i-1),
                    theta_star, tau_e, (P_aug-1), 0, slice_evals_pm, lower, upper, n_slice_iter, 
                    w_tot, s, Omegat_ns(dist(i-1),0),  Omegas_ns(dist(i-1),0),
                    gp_mod, n_parms, pos_s, noise, jitter,
                    y_ns(dist(i-1),0), 
                    atau, btau, transition(i-1), ipr); /* m value not used - set to 0 */
               slice_levals_tau  += slice_evals_pm;
           
          } /* end loop over tempered (lower rank) dists to sample theta_star[p,m] */
          /* compute average number of slice.evals per run */
                  
          // accept reject step for x0v, starting from x0^ 
          // compute portion of probability of move only relating to transformed space */
          logf_tau     = temper_logpmove_alt(Pi_n, dist);
                  
          // compute probability of move
          /* add in log-posterior under fine (full) density in original space */
          /* full data with (Omega,y) of size T */          
          /* old value */
          /* U_last is carried over from by-cluster updates of U_last.slice(m) */
          logpitau_old       = logFtau_like(y, s, U_last, ipr) + log_prior(tau_e,atau,btau);
          /* proposed value */
          compute_U(theta_star, theta_updown((vals-1)), jitter, gp_mod, n_parms, pos_s,
                    Omega_t, Omega_s, noise, U_new);
          logpitau_new       = logFtau_like(y, s, U_new, ipr) 
                                   + log_prior(theta_updown((vals-1)),atau,btau);
          logf_tau           += logpitau_new; /* proposed value */
          /* using last updated log_pi_old from sampling theta_star */
          logf_tau           -= logpitau_old; /* old value */
                  
          /* acceptance step */
          thresh_tau          = runif(1,0,1)[0];
          if( log(thresh_tau) < logf_tau )
          {
               tau_e         = theta_updown(vals-1); /* end value of tempered steps */
               /* save fields C & U */
               /* for likelihood computations on cluster assignments */
               U_last         = U_new;
               /* track global number of tempered acceptances */
               accept_temper   += 1;
          } /* end M-H acceptance step for updating tau_e */ 
        END_RCPP
} /* end function temper() to sample theta_star(p,m) through tempered transitions in expanded space */

SEXP temper_b(cube& U_last, mat& theta_star, double tau_e, uvec& gp_mod, uvec& n_parms, 
          uvec& pos_s, double& slice_levals_theta,
          double& n_slice_theta, double& accept_temper, 
          double& n_temper_evals, double lower, double upper, 
          int n_slice_iter, const mat& w_tot, const ucolvec& s, const field<mat>& Omegat_ns, 
          const field<cube>& Omegas_ns, const mat& Omega_t, const cube& Omega_s, 
          const field<mat>& bb_ns, const mat& bb, double jitter, double a, double b,
          const vec& ipr)
{
        BEGIN_RCPP
        // number of layers
        int n                 = Omegat_ns.n_rows;
        // dimension of GP covariance
        int T                 = Omega_t.n_rows;
        // number of parameters to sample
        int M                 = theta_star.n_cols;
        int P                 = theta_star.n_rows;
        /* intermediate sampled values for each theta_star(p,m) */
        int vals              = 2 + (2*n-1); /* in_out + (up_down * # layers) - summit  */
        /* number of uni.slice() steps is vals - 1, since bring in_up initial value */
        int steps             = vals - 1; /* number of tempered sampling steps */
        /* tempered sample and density values */
        colvec theta_updown(vals); mat Pi_n(n,vals);  /* note: 1 more than "steps" */
        // every step is a transition, except step n+1 where the most coarse distribution
        // is sampled twice.  The ladder is 1,2,...,n,n,...,2,1
        colvec transition(steps); transition.ones(); transition(n) = 0;
        /* distribution indicator to use for / update in sampling steps */
        ucolvec dist             = temper_dist_alt(n); /* for n = 2, dist order is 0 1 1 0 */
        /* loop counters */
        int i = 0, p = 0, m = 0;
        /* since sampling T x 1, b_i, noise is always 0 */
        int noise = -1; /* addition of 1/tau_e to C based on noise > 0 */
        /* reduce jitter to use for slice moves under lower dimensional distributions */
        double jitter_slice = 0.25*jitter;
        // tempered transition probabilities and M-H acceptance step for theta_star(p,m)
        double logf_pm = 0, thresh_pm = 0;
        double loglike_m_old = 0, logpi_pmold = 0, loglike_pmnew = 0, logpi_pmnew = 0;
        mat Upm_new(T,T); 
        // capture average number of slice sampling evals per tempered step for each theta_star(p,m)
        int slice_evals_pm         = 0; /* counter for each (p,m) + tau_e */
        slice_levals_theta         = 0; /* numberof slice sampling like evals across (p,m) */
        n_slice_theta              = (P*M*steps); /* number of (p,m) slice steps */
        accept_temper              = 0; /* global step-up acceptance yes = 1, no = 0 */ 
        n_temper_evals             = 0; /* total number of temper() proposals */
        // n_slice_iter = 10000; /* number of steps to compose (widen) interval / slice */
        
        // sample theta_star[p,m], each in an expanded space
        for(m = 0; m < M; m++) /* sample clusters, m */
        {
               // update log-like for theta_star(,m) since share same observations
               // within cluster can avoid computation of old likelihood
               // log-priors are different for P parameters in cluster m, so have to 
               // compute those.
               /* theta_star[,m] \perp theta_star[,m'] */
               loglike_m_old  = logFm_like(bb, m, s, U_last.slice(m), ipr); 
               for(p = 0; p < P; p++) /* sample parameter types, p */
               {
                    /* tick up temper_evals */
                    n_temper_evals   += 1;
                    // perform tempered sampling for theta_star[p,m]
                    theta_updown.zeros(); /* for n = 2, theta_updown is of length 5 */
                    theta_updown(0) = theta_star(p,m);
                    Pi_n.zeros();
                    /* note: updates in theta_updown start at 1 bc bring in previous theta_star(p,m) */
                    for(i = 1; i < vals; i++) /* emphasize that update uni.slice from 1,..,steps */
                    {
                         /* note: theta_updown and Pi_n are of length vals = steps + 1 = 5 */
                         // univariate slice sampler
                         // returns vector, theta_updown and matrix, Pi_n
                         uni_slice_pm(theta_updown, Pi_n, i, dist(i-1),
                              theta_star, tau_e, p, m, slice_evals_pm, lower, upper, n_slice_iter, 
                              w_tot, s, Omegat_ns(dist(i-1),0),  Omegas_ns(dist(i-1),0),
                              gp_mod, n_parms, pos_s, noise, jitter_slice,
                              bb_ns(dist(i-1),0), 
                              a, b, transition(i-1), ipr); 
                         /* memo; the "i-1" in transition reflects c++ start count at 0 for vector */
                         slice_levals_theta  += slice_evals_pm;
           
                    } /* end loop over tempered (lower rank) dists to sample theta_star[p,m] */
                  
                    // accept reject step for x0v, starting from x0^ 
                    // compute portion of probability of move only relating to temporary space */
                    logf_pm     = temper_logpmove_alt(Pi_n, dist);
                  
                    // compute probability of move
                    /* add in log-posterior under fine (full) density in original space */
                    /* full data with (Omega,y) of size T */
                    logpi_pmold         = loglike_m_old + log_prior(theta_star(p,m),a,b);/* old value */
                    Upm_new             = compute_Upm(theta_updown((vals-1)), theta_star, tau_e, jitter, 
                                             p, m, gp_mod, n_parms, pos_s,
                                             Omega_t, Omega_s, noise);
                    loglike_pmnew       = logFm_like(bb, m, s, Upm_new, ipr);
                    logpi_pmnew         = loglike_pmnew + log_prior(theta_updown((vals-1)),a,b);
                    logf_pm             += logpi_pmnew; /* proposed value */
                    logf_pm             -= logpi_pmold;
                  
                    /* acceptance step */
                    thresh_pm          = runif(1,0,1)[0];
                    if( log(thresh_pm) < logf_pm )
                    {
                         theta_star(p,m)         = theta_updown(vals-1); /* end value of tempered steps */
                         loglike_m_old           = loglike_pmnew;
                         U_last.slice(m)         = Upm_new;
                         accept_temper           += 1;
                    } /* end M-H acceptance step */        
               } /* end sampling loop p over GP parameters within clusters */  
          }/* end sampling loop m over clusters */
          /* end updating of theta_star(p,m) */          
        END_RCPP
} /* end function temper_b() to sample theta_star(p,m) through tempered transitions in expanded space */

// update vector of cluster membership indicators, s(i),....,s(N)
SEXP auxclusterstep(mat& theta_star, mat& wpm, cube& U_last, const mat& Omega_t, const cube& Omega_s,
            const mat& y, double tau_e, int noise, double jitter, uvec& gp_mod, uvec& n_parms, 
            uvec& pos_s, ucolvec& s, ucolvec& num, unsigned int& M, const int& w_star, 
            double& conc, double a, double b, const vec& ipr, colvec& Num)
{
          BEGIN_RCPP
          // sample cluster assignments, s(1), ..., s(N)
          // fixing all other parameters, including B_M
          unsigned int auxSize, startPos, h, ell;
          int N     = s.n_elem;
	        int P     = theta_star.n_rows;
          int T     = Omega_t.n_rows;
          int i;
          mat theta_aux_h; /* added auxSize = h - M auxiliary locations */
          mat wpm_aux_h; /* added slice sampling tuning parameters for auxiliary locations */
          cube U_aux; /* added auxSize auxiliary covariance cholesky locations */
          cube U_keep; /* final set of by-cluster choleskys with observations */
          uvec loc_keep; /* final set of cluster location labels linked to observations */
          unsigned int num_keep = 0; /* number of cluster location labels linked to observations */
          //colvec Num; /* Horvitz-Thompson scaled up num vector with inverse probability weighting */
          for(i = 0; i < N; i++)
          {
               /* set starting index for w_star new cluster locations */
	          /* remember, this is c++, so starting index is 0 for all vectors */
               /*, including s = (0,...,M-1) */
	          if( num(s(i)) > 1 ) /* set start at M */
	          {
		          auxSize 	= w_star;
		          num(s(i))	-= 1; /* decrement cluster count of subject "i" */
              /* scale up num to population totals, Num, based on H-T inverse probability estimator */
              Num(s(i)) -= 1/ipr(i);
		          startPos	= M;
	          }else{ /* num(s(i)) == 1 so draw new clusters staring at M-1 */
		          auxSize 	= w_star - 1;
		          /* re-assign the current location for subject i to the last location, M-1 */
		          theta_star.insert_cols(M,1);
                    wpm.insert_cols(M,1);
                    /* carrying along logfm_old vector to maintain current cluster */
                    /* count and location. Recall, logfm_old is used for Metropolis */
                    /* steps to update rho_star */
		          theta_star.col(M)		= theta_star.col(s(i));
                    wpm.col(M)               = wpm.col(s(i));
                    U_last.insert_slices(M,1);
                    U_last.slice(M)          = U_last.slice(s(i));
                    num.insert_rows(M,1);
		          num(M)			     = 0; /* insert new location ahead of observations */
                    /* scale up num to population totals, Num, based on H-T inverse probability estimator */
                    Num.insert_rows(M,1);
     	          Num(M)			     = 0; /* insert new location ahead of observations */
		          theta_star.shed_col(s(i));
                    wpm.shed_col(s(i));
                    U_last.shed_slice(s(i));
		          num.shed_row(s(i));
                    Num.shed_row(s(i));
		          s( find(s > s(i)) ) 	-= 1;
		          s(i)			          = M-1;
		          startPos		          = M-1;
	          }/* end setting starting index for w_star new cluster locations */
               /* set total number of existing and auxiliary cluster locations, h */
	          h 		= M + auxSize;
	          colvec weights(h); weights.zeros();
               
               /* scale up num to population totals, Num, based on H-T inverse probability estimator */
               //pop_Num(s, ipr, Num); /* computed after set s, num and h */

	          /* sample w_star new location values from prior (F_0) */
	          /* since ahead of observations */
	          /* compute weights used to draw s(i) */
	          if( h > M ) /* num(s(i)) > 1 | (num(s(i)) == 1 & w_star > 1) */
	          {
                    /* create auxSize new locations to insert to theta_star */
                    NumericVector _theta_vec            = rgamma( (P*auxSize), a, (1/b) );
                    vec theta_vec                       = as<vec>(_theta_vec);
                    theta_aux_h                         = reshape( theta_vec, P, auxSize );
		                theta_star.insert_cols(M,theta_aux_h); /* theta_star now has h columns */
                    /* set w_aux_h(p,m) to [0.1,0.9] quantiles of gamma distribution with */
                    /* shape = theta_star(p,m), rate = 1 */
                    wpm_aux(wpm_aux_h, theta_aux_h);
                    wpm.insert_cols(M, wpm_aux_h);
                    /* compute T x T x auxSize, U_aux */
                    U_aux.set_size(T,T,auxSize);
                    /* memo: use theta_aux_h to compute U_aux; doesn't appear directly in like */
                    compute_U(theta_aux_h, tau_e, jitter, gp_mod, n_parms, pos_s,
                                   Omega_t, Omega_s, noise, U_aux);       
                    U_last.insert_slices(M,U_aux); /* U_last now has h slices */
                    /* memo: no observations, yet */
                    /* adding entries for aux vars and set to 0 */
                    /* by default, new rows/cols/slices set to 0 */
		                num.insert_rows(M,auxSize,true);
                    Num.insert_rows(M,auxSize,true); /* population uplifted */
                    /* inserts 0 for M -> h-l positions */
	     	      weights			     = Num / (double(N)- (1/ipr(i)) +conc); 
		          weights.subvec(M,(h-1))	= ( (conc/w_star) / (double(N) - (1/ipr(i)) + conc) )*ones(auxSize);
	          }else{ /* num(s(i)) == 1 & w_star == 1, so new location already sampled at M-1 */
		          weights			     = Num / (double(N)- (1/ipr(i)) +conc);
		          weights(M-1)		     = (conc/w_star) / (double(N) - (1/ipr(i)) + conc);
               } /* end sampling of new locations and computing weights */

	          /* draw value for cluster index, s(i) */
	          /* note that rdrawone() function has minimum value of 0 */
	          for( ell = 0; ell < h; ell++)
	          {
		          s(i)           = ell; /* temporary assignment */
		          weights(ell)	= weights(ell) * 
                                        exp( logy_like(i, y, ipr, s, U_last) );
	          }
	          weights.elem( find_nonfinite(weights) ).zeros();
	          weights		     = weights / sum(weights);
	          // conduct discrete posterior draw for s(i), where min(s(i)) = 0
               s(i) 		           = rdrawone(weights, h);
	             num(s(i))		       += 1;
               Num(s(i))           += 1/ipr(i);
	
	          /* remove clusters with no observations  - must be among the w_star new clusters */
               loc_keep            = find(num != 0);
               num_keep            = loc_keep.n_elem;
               num                 = num( loc_keep );
               Num                 = Num( loc_keep );
	          theta_star		= theta_star.cols( loc_keep );
               wpm                 = wpm.cols( loc_keep );
               /* Armadillo doesn't have a non-contiguous operation on slices */
               /* must update slice-by-slice */
               U_keep.set_size(T,T,num_keep);
               for( ell = 0; ell < num_keep; ell++ )
               {
                    U_keep.slice(ell)   = U_last.slice( loc_keep(ell) );
               }
               
               /* replace U_last with U_keep, which holds choleskys linked to observations */
               U_last.set_size(T,T,num_keep);
               U_last              = U_keep;
	          /* if s(i) > startPos-1, then the drawn cluster is from the w_star new */
	          /* so we need to make it contiguous through shifting it to the last cluster */
	          if( s(i) > (startPos-1) ) /* we're replacing the original startPos position */
	          {
		          s(i)	          = startPos;
                    // since sampling one observation, i = 1,...N, at a time,
                    // if one of the w_star new clusters is picked, it can only have
                    // a single observation (for a single new cluster at location, startPos)
                    // so we compute a double value for the log-likelihood of 
                    // rho_star(startPos) to later use for sampling.
	          }
	   	
	          /* reset M */
	          M = theta_star.n_cols;
         
        } /* end loop i for cluster assignment */
        END_RCPP
} /* end function auxclusterstep for cluster assignments, s  */
	          
SEXP move_s(ucolvec& s, ucolvec& num, const mat& y, const cube& U_last,
            const colvec& p, const double& conc, const colvec& ipr)
{
  BEGIN_RCPP
  int M     = num.n_elem; /* also length of p */
  int n     = s.n_elem;
  int i     = 0, ell = 0;
  
  colvec weights = 1/ipr;
  colvec p_w(M); p_w.zeros();
  
  for( i = 0; i < n; i++ )
  {
    for( ell = 0; ell < M; ell++)
    {
      s(i)      = ell; /* temporary assignment */
      p_w(ell)	= pow(p(ell),weights(i)) * 
      exp( logy_like(i, y, ipr, s, U_last) );
    }
    p_w		          = p_w / sum(p_w);
    // conduct discrete posterior draw for s(i), where min(s(i)) = 0
    s(i) 		          = rdrawone(p_w, M);
    num(s(i))		      += 1;
  } /* end loop i over n sample observations */
  
  END_RCPP
} /* end function move_s() to sample cluster indicators, s(1), .... , s(n) */

// // update vector of cluster membership indicators, s(i),....,s(N)
// SEXP auxclusterstep2(mat& theta_star, mat& wpm, cube& U_last, const mat& Omega_t, const cube& Omega_s,
//                    const mat& y, double tau_e, int noise, double jitter, uvec& gp_mod, uvec& n_parms, 
//                    uvec& pos_s, ucolvec& s, ucolvec& num, unsigned int& M, const int& w_star, 
//                    double& conc, double a, double b, const vec& ipr, colvec& Num)
// {
//   BEGIN_RCPP
//   // sample cluster assignments, s(1), ..., s(N)
//   // fixing all other parameters, including B_M
//   unsigned int auxSize, startPos, h, ell;
//   int N     = s.n_elem;
//   int P     = theta_star.n_rows;
//   int T     = Omega_t.n_rows;
//   int i;
//   mat theta_aux_h; /* added auxSize = h - M auxiliary locations */
//   mat wpm_aux_h; /* added slice sampling tuning parameters for auxiliary locations */
//   cube U_aux; /* added auxSize auxiliary covariance cholesky locations */
//   cube U_keep; /* final set of by-cluster choleskys with observations */
//   uvec loc_keep; /* final set of cluster location labels linked to observations */
//   unsigned int num_keep = 0; /* number of cluster location labels linked to observations */
//   //colvec Num; /* Horvitz-Thompson scaled up num vector with inverse probability weighting */
//   for(i = 0; i < N; i++)
//   {
//     /* set starting index for w_star new cluster locations */
//     /* remember, this is c++, so starting index is 0 for all vectors */
//     /*, including s = (0,...,M-1) */
//     if( num(s(i)) > 1 ) /* set start at M */
//     {
//       auxSize 	= w_star;
//       num(s(i))	-= 1; /* decrement cluster count of subject "i" */
//       /* scale up num to population totals, Num, based on H-T inverse probability estimator */
//       Num(s(i)) -= 1/ipr(i);
//       startPos	= M;
//       }else{ /* num(s(i)) == 1 so draw new clusters staring at M-1 */
//       auxSize 	= w_star - 1;
//       /* re-assign the current location for subject i to the last location, M-1 */
//       theta_star.insert_cols(M,1);
//       wpm.insert_cols(M,1);
//       /* carrying along logfm_old vector to maintain current cluster */
//       /* count and location. Recall, logfm_old is used for Metropolis */
//       /* steps to update rho_star */
//       theta_star.col(M)		= theta_star.col(s(i));
//       wpm.col(M)               = wpm.col(s(i));
//       U_last.insert_slices(M,1);
//       U_last.slice(M)          = U_last.slice(s(i));
//       num.insert_rows(M,1);
//       num(M)			     = 0; /* insert new location ahead of observations */
//       /* scale up num to population totals, Num, based on H-T inverse probability estimator */
//       Num.insert_rows(M,1);
//       Num(M)			     = 0; /* insert new location ahead of observations */
//       theta_star.shed_col(s(i));
//       wpm.shed_col(s(i));
//       U_last.shed_slice(s(i));
//       num.shed_row(s(i));
//       Num.shed_row(s(i));
//       s( find(s > s(i)) ) 	-= 1;
//       s(i)			            = M-1;
//       startPos		          = M-1;
//     }/* end setting starting index for w_star new cluster locations */
//     /* set total number of existing and auxiliary cluster locations, h */
//     h 		= M + auxSize;
//     colvec weights(h); weights.zeros();
//     
//     /* scale up num to population totals, Num, based on H-T inverse probability estimator */
//     //pop_Num(s, ipr, Num); /* computed after set s, num and h */
//     
//     /* sample w_star new location values from prior (F_0) */
//     /* since ahead of observations */
//     /* compute weights used to draw s(i) */
//     if( h > M ) /* num(s(i)) > 1 | (num(s(i)) == 1 & w_star > 1) */
//     {
//       /* create auxSize new locations to insert to theta_star */
//       NumericVector _theta_vec            = rgamma( (P*auxSize), a, (1/b) );
//       vec theta_vec                       = as<vec>(_theta_vec);
//       theta_aux_h                         = reshape( theta_vec, P, auxSize );
//       theta_star.insert_cols(M,theta_aux_h); /* theta_star now has h columns */
//       /* set w_aux_h(p,m) to [0.1,0.9] quantiles of gamma distribution with */
//       /* shape = theta_star(p,m), rate = 1 */
//       wpm_aux(wpm_aux_h, theta_aux_h);
//       wpm.insert_cols(M, wpm_aux_h);
//       /* compute T x T x auxSize, U_aux */
//       U_aux.set_size(T,T,auxSize);
//       /* memo: use theta_aux_h to compute U_aux; doesn't appear directly in like */
//       compute_U(theta_aux_h, tau_e, jitter, gp_mod, n_parms, pos_s,
//                Omega_t, Omega_s, noise, U_aux);       
//       U_last.insert_slices(M,U_aux); /* U_last now has h slices */
//       /* memo: no observations, yet */
//       /* adding entries for aux vars and set to 0 */
//       /* by default, new rows/cols/slices set to 0 */
//       num.insert_rows(M,auxSize,true);
//       Num.insert_rows(M,auxSize,true); /* population uplifted */
//       /* inserts 0 for M -> h-l positions */
//       weights			     = num / (double(N)- 1 +conc); 
//       weights.subvec(M,(h-1))	= ( (conc/w_star) / (double(N) - 1 + conc) )*ones(auxSize);
//     }else{ /* num(s(i)) == 1 & w_star == 1, so new location already sampled at M-1 */
//       weights			     = num / (double(N)- 1 +conc);
//       weights(M-1)		     = (conc/w_star) / (double(N) - 1 + conc);
//     } /* end sampling of new locations and computing weights */
//       
//     /* draw value for cluster index, s(i) */
//     /* note that rdrawone() function has minimum value of 0 */
//     for( ell = 0; ell < h; ell++)
//     {
//       s(i)           = ell; /* temporary assignment */
//       weights(ell)	= weights(ell) * 
//                             exp( logy_like(i, y, ipr, s, U_last) );
//     }
//     weights		     = weights / sum(weights);
//     // conduct discrete posterior draw for s(i), where min(s(i)) = 0
//     s(i) 		           = rdrawone(weights, h);
//     num(s(i))		       += 1;
//     Num(s(i))          += 1/ipr(i);
//     
//     /* remove clusters with no observations  - must be among the w_star new clusters */
//     loc_keep            = find(num != 0);
//     num_keep            = loc_keep.n_elem;
//     num                 = num( loc_keep );
//     Num                 = Num( loc_keep );
//     theta_star		      = theta_star.cols( loc_keep );
//     wpm                 = wpm.cols( loc_keep );
//     /* Armadillo doesn't have a non-contiguous operation on slices */
//     /* must update slice-by-slice */
//     U_keep.set_size(T,T,num_keep);
//     for( ell = 0; ell < num_keep; ell++ )
//     {
//       U_keep.slice(ell)   = U_last.slice( loc_keep(ell) );
//     }
//     
//     /* replace U_last with U_keep, which holds choleskys linked to observations */
//     U_last.set_size(T,T,num_keep);
//     U_last              = U_keep;
//     /* if s(i) > startPos-1, then the drawn cluster is from the w_star new */
//     /* so we need to make it contiguous through shifting it to the last cluster */
//     if( s(i) > (startPos-1) ) /* we're replacing the original startPos position */
//     {
//       s(i)	          = startPos;
//       // since sampling one observation, i = 1,...N, at a time,
//       // if one of the w_star new clusters is picked, it can only have
//       // a single observation (for a single new cluster at location, startPos)
//       // so we compute a double value for the log-likelihood of 
//       // rho_star(startPos) to later use for sampling.
//     }
//     
//     /* reset M */
//     M = theta_star.n_cols;
//     
//   } /* end loop i for cluster assignment */
//   END_RCPP
// } /* end function auxclusterstep for cluster assignments, s  */	          
	          
SEXP pop_Num(const ucolvec& s, const vec& ipr, colvec& Num)
{
     BEGIN_RCPP 
     ucolvec clusters    = unique( s );
     int M               = clusters.n_elem;
     int m;
     uvec pos_m;
     Num.set_size(M);
     
     for( m = 0; m < M; m++ )
     {
         pos_m      = find(s == m);
         Num(m)     = sum( 1 / ipr.elem(pos_m) ); /* element-wise division */
     }
     END_RCPP
} /* end function to create Horvitz-Thompson summary for cluster membership counts */

SEXP concstep(double& conc, int M, int N,  double a6, double b6)
{
        BEGIN_RCPP
        // sample DP concentration parameter, c
        // c|s independent of sample, y
        // see page 585 of Escobar and West (1995)
        double eta = rbeta( 1,conc+1,double(N) )[0];
        double drawP;
        drawP = (a6+double(M)-1)/( a6+double(M)-1 + double(N)*(b6-log(eta)) );
        int mixcomp = rbinom(1,1,drawP)[0];
        if( mixcomp == 1 )
        {
            conc = rgamma( 1, a6 + M, 1/(b6-log(eta)) )[0];
        }
        else
        {
            conc = rgamma( 1, a6 + M - 1, 1/(b6-log(eta)) )[0];
        }
        END_RCPP
} /* end function to sample concentration parameter, conc */

SEXP wp_tune(cube& Theta_tune, colvec& wtune)
{
        BEGIN_RCPP
        // Theta_tune is a K x (N*P) matrix of theta returned from 
        // posterior sampling during a tuning run (of K iterations).
        int M  = Theta_tune.n_cols;
        int K  = Theta_tune.n_slices;
        int P  = Theta_tune.n_rows;
        /* P x M, w, holds tuned step width for each (p,m) */
        mat w(P,M);
        /* K x 1, thetastar_pm, holds K samples used to tune for each (p,m) */
        colvec thetastar_pm(K);
        /* compute the (0.05,0.95) sampled quantiles to reset wtune */
        int plow_ind, phigh_ind;
        int m, p;
        for( p = 0; p < P; p++ )
        {    
             for( m = 0; m < M; m++ )
             {
                    thetastar_pm     = Theta_tune.tube(p,m); /* K x 1 */
                    thetastar_pm     = sort(thetastar_pm);
                    plow_ind         = floor( 0.05*K ) - 1; /* - 1 since index starts at 0 */
                    phigh_ind        = floor( 0.95*K ) - 1;
                    /* update slice sampling step width for parameters (p,1:M) */
                    w(p,m)           = thetastar_pm(phigh_ind) - thetastar_pm(plow_ind);
             } /* end loop m over (fixed number of) clusters */
             /* select each element, p, of wtune by choosing the minimum over the M clusters for each p */
             wtune(p)             = max(w.row(p));
        } /* end loop p over GP paramters */
        
        END_RCPP
} /* end function wp_tune() to tune the slice sampling step width for theta_star */

SEXP wpm_tune(cube& Theta_tune, mat& wtune)
{
        BEGIN_RCPP
        // Theta_tune is a K x (N*P) matrix of theta returned from 
        // posterior sampling during a tuning run (of K iterations).
        int M  = Theta_tune.n_cols;
        int K  = Theta_tune.n_slices;
        int P  = Theta_tune.n_rows;
        /* P x M, w, holds tuned step width for each (p,m) */
        mat w(P,M);
        /* K x 1, thetastar_pm, holds K samples used to tune for each (p,m) */
        colvec thetastar_pm(K);
        /* compute the (0.05,0.95) sampled quantiles to reset wtune */
        int plow_ind, phigh_ind;
        int m, p;
        double wpm_sample;
        for( p = 0; p < P; p++ )
        {    
             for( m = 0; m < M; m++ )
             {
                    thetastar_pm     = Theta_tune.tube(p,m); /* K x 1 */
                    thetastar_pm     = sort(thetastar_pm);
                    plow_ind         = floor( 0.05*K ) - 1; /* - 1 since index starts at 0 */
                    phigh_ind        = floor( 0.95*K ) - 1;
                    /* update slice sampling step width for parameters (p,1:M) */
                    wpm_sample           = thetastar_pm(phigh_ind) - thetastar_pm(plow_ind);
                    w(p,m)               = max(wpm_sample,0.05);
             } /* end loop m over (fixed number of) clusters */
        } /* end loop p over GP paramters */
        wtune            = w; /* now dimensioned as a (p,m) matrix */
        
        END_RCPP
} /* end function wp_tune() to tune the slice sampling step width for theta_star */


SEXP wtau_tune(colvec& Taue_tune, double& wtune)
{
        BEGIN_RCPP
        int K                 = Taue_tune.n_elem; /* number of tuning iterations */
        /* compute sample (0.05,0.95) quantiles to set wtune */
        colvec Taue_ts     = sort( Taue_tune );
        int low_ind        = floor( 0.05*K ) - 1; /* - 1 since index starts at 0 */
        int high_ind       = floor( 0.95*K ) - 1;
        double wsample     = Taue_ts(high_ind) - Taue_ts(low_ind);
        wtune              = max(wsample,1.0);
        END_RCPP
} /* end function wtau_tune() to tune the slice sampling step width for tau_e */

SEXP wpm_aux(mat& wpm_aux_h, const mat& theta_aux_h)
{
        BEGIN_RCPP
        // Theta_tune is a K x (N*P) matrix of theta returned from 
        // posterior sampling during a tuning run (of K iterations).
        int auxSize        = theta_aux_h.n_cols;
        int P              = theta_aux_h.n_rows;
        wpm_aux_h.set_size(P,auxSize);
        int m, p;
        double shape_mean_pm;
        vec range_pm(2);
        /* P x M, w, holds tuned step width for each (p,m) */
       
        /* fill values for new columns (clusters) of theta_star */
        for( p = 0; p < P; p++ )
        {
             for( m = 0; m < auxSize; m++ )
             {
                  /* score_pm is M_del x 1 */
                  shape_mean_pm    =  double(theta_aux_h(p,m));
                  range_pm(0)      = R::qgamma(0.9,shape_mean_pm, 1.0, 1, 0) - 
                                        R::qgamma(0.1,shape_mean_pm, 1.0, 1, 0);
                  range_pm(1)      = 0.05; /* this is minimum range to avoid degeneracy */
                  wpm_aux_h(p,m)   = max( range_pm );
             } /* end loop m over M_del new clusters */
        } /* end loop p over P GP covariance parameters. */
        END_RCPP
} /* end function wp_tune() to tune the slice sampling step width for theta_star */

mat update_w(const mat& wpm, double wtau)
{
     int M                                        = wpm.n_cols;
     int P                                        = wpm.n_rows;
     mat w_tot((P+1),M);
     w_tot.submat(span(0,(P-1)),span::all)        = wpm; 
     w_tot.row(P)                                 = wtau*ones<rowvec>(M);
     return w_tot;
} /* end function update_w to update slice sampling (p,m) step size hyperparameters */

SEXP move_b(mat& bb, cube& invG_star, const cube& U_last, const ucolvec& s, 
               const mat& y, double tau_e)
{
     BEGIN_RCPP
     // N x T matrix, bb, includes de-noised functions
     // N x T matrix, y, holds the data for observation i at time t
     // return T x T x M cube, invG_star, of (no noise) GP covariance matrices 
     int N          = y.n_rows;
     int M          = U_last.n_slices;
     int T          = y.n_cols;
     int i = 0, m = 0;
     invG_star.set_size(T,T,M); /* M varies iteration-by-iteration for returned invG_star */
     mat eye_T = eye(T,T);
     mat inv_Um(T,T);
     colvec e_i(T), b_i(T); /* posterior sample of bb_i */
     mat phi_i(T,T); /* posterior precision for bb_i */
     /* compute GP precision matrix for each cluster location */
     for( m = 0; m < M; m++ )
     {
         inv_Um                 = inv(trimatu(U_last.slice(m)));
         invG_star.slice(m)     = inv_Um * inv_Um.t();
     } /* end loop m over cluster locations */
     
     /* sample T x 1, bb_i */
     for( i = 0; i < N; i++ )
     {
          e_i            = tau_e*y.row(i).t(); /* T x 1 */
          phi_i          = tau_e*eye_T + invG_star.slice(s(i)); /* T x T precision matrix */
          rmvnbasic(phi_i, e_i, b_i);
          bb.row(i)      = b_i.t();
     }
     END_RCPP
} /* end function miss_ystep() to sample missing data values */

SEXP move_bslice(mat& bb, const cube& U_last, const ucolvec& s,  
               const mat& y, double tau_e, int R)
{
     BEGIN_RCPP
     // N x T matrix, bb, includes de-noised functions
     // N x T matrix, y, holds the data for observation i at time t
     // return T x T x M cube, invG_star, of (no noise) GP covariance matrices 
     int N          = y.n_rows;
     int T          = y.n_cols;
     int i = 0, r = 0;
     colvec h(T), vi(T); /* posterior sample of bb_i */
     h.zeros();
     rowvec yi_vec(T), bbi_vec(T), bbi_prime(T);
     double logLike_starti, logLike_candi, u, logy;
     double theta, theta_min, theta_max;
     
     for( r = 0; r < R; r++ ) /* R re-samples takes advantage of easy generation from prior */
     {
          /* sample T x 1, bb_i */
          for( i = 0; i < N; i++ )
          {
               /* draw ellipse, v_i */
               rmvnchol(U_last.slice(s(i)), h, vi); /* v_i ~ N_T(0,C*(s(i))) */
               /* initial evaluation of  L(bb_i) */
               bbi_vec        = bb.row(i); /* T x 1 row i of bb_mat */
               yi_vec         = y.row(i);
               /* need to employ Rcpp containers to get vectorized dnorm */
               logLike_starti = log_dnorm_vec(yi_vec, bbi_vec, tau_e); /* T independent sums */
               /* compute likelihood threshold for slice */
               u              = runif(1, 0, 1)[0];
               logy           = logLike_starti + log(u); /* slice threshold */
               /* draw an initial proposal for mixture element, theta */
               theta          = runif(1, 0, 2*M_PI)[0]; /* in radians */
               theta_min      = theta - 2*M_PI;
               theta_max      = theta;
               /* generate a proposal for bb_i */
               bbi_prime      = bbi_vec * cos(theta) + vi.t() * sin(theta);
               /* evaluate L(bbi_prime) and accept-break if L(bbi_prime) >= logy */
               /* else, SHRINK interval for drawing theta and re-compute bbi_prime */
               for(;;) /* infinite loop that contains a break statement */
               {                 
                    logLike_candi    = log_dnorm_vec(yi_vec, bbi_prime, tau_e); 
                         /* T independent sums */
             
                    /* condition on whether new sampled value is outside of slice */
                    if( logLike_candi > logy ){break;}
             
                    /* shrink interval */
                    if(theta < 0)
                    {
                         theta_min      = theta;
                    }else{ /* theta > 0  */
                         theta_max      = theta;
                    } /* end condition to shrink upper or lower interval */
             
                    theta            = runif(1,theta_min,theta_max)[0];
                    bbi_prime        = bbi_vec * cos(theta) + vi.t() * sin(theta);
               } /* end infinite loop on shrinking slice with break when Pin_new >= logy */
               /* set row i of bb to the accepted candidate */
               bb.row(i)        = bbi_prime;   
          } /* end loop over rows i of N x T, bb */  
     } /* end duplicate repeated sampling loops for bb within an MCMC sampling iteration */
     END_RCPP
} /* end function move_bslice to sample the N x T matrix of GP functions, bb */

SEXP gen_bb_ns(const mat& bb, const field<uvec>& index, 
               field<mat>& bb_ns)
{
     BEGIN_RCPP
     int N     = bb.n_rows;
     int n     = index.n_rows; /* number of progressively coarser distributions */
     int T_i; /* number of sub-sampled time points per distribution */
     int i = 0;
     for(i = 0; i < n; i++)
     {
         T_i                  = index(i,0).n_elem;
         /* set size for element i in field<mat> Omegat_ns, field<cube> Omegas_ns, field<colvec> y_ns */
         bb_ns(i,0).set_size(N, T_i);
         bb_ns(i,0)            = bb.cols( index(i,0) ); 
     } /* end loop i over number of tempered distributions - steps */

     END_RCPP
} /* end function gen_bb_ns */
