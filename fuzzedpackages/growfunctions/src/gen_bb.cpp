#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

SEXP gen_bb(const mat& y, const mat& Tau_e, umat& S, const field<mat>& Theta_star, const mat& Omega_t, 
               const cube& Omega_s, uvec& gp_mod, uvec& n_parms, uvec& pos_s, uvec& P_vec,
               mat& bb, field<mat>& f, field<cube>& invG_star, double jitter)
{
        BEGIN_RCPP
        // rhotrend_star and ksi_star is of length, L, equal to total number of terms
        // only contains trend components
        // rhoseas_star, Omega_seas is of length equal to number of seasonal components
        // "seasonal" is of length L
        // seasonal = [0 0 1 0 2]' for L = 5 where non-zero entries index 
        // which seasonal entry in rhoseas_star, Omega_seas is associated to that trend term
        // int M denotes the total number of clusters
        /* indices and loop counters */
        int n_iter  = S.n_rows;
        int N       = S.n_cols;
        int L       = gp_mod.n_elem;
        int T       = Omega_t.n_rows;
        int P       = Theta_star(0,0).n_rows; /* all theta_star matrices have P rows */
        int i, j, l, k, m, M_k;
        
        //create objects
        /* objects need to compose total covariances, G_star */
        cube G_star_lk; /* l indexes term and k, mcmc iteration */
        cube G_star_k; /* sum G_star_l over l = 1,..,L terms */
        cube invG_star_k; /* returning this object b/c use for prediction */
        field<cube> G_star_term(n_iter,L); /* by term T x T x M_k covariances */
        for( k = 0; k < n_iter; k++ )
        {
             M_k       = Theta_star(k,0).n_cols;
             for( l = 0; l < L; l++ )
             {
                  G_star_term(k,l).set_size(T,T,M_k);
             } /* end loop L over terms */
        } /* end loop k over sampling iterations */
        /* parameters used to compute G_star_lkm */
        colvec theta_star_kml; /* contains P(l) parameters */
        colvec theta_star_km(P); /* across all terms */
        /* mean functions to draw GP posterior for bb */
        colvec mu_fkli(T); mu_fkli.zeros(); 
        cube musamples_fl(n_iter,T,N); musamples_fl.zeros();
        mat mu_fl(N,T); mu_fl.zeros();
        /* variance functions to draw GP posterior for f(l,0) that sum to bb */
        cube var_lk(T,T,N), var_l(T,T,N);
        var_lk.zeros(); var_l.zeros(); 
        /* object to capture posterior draws */
        colvec fli_draw(T);
        /* by term indexing for covariance computations*/
        int start_l, end_l;
        /* cube to capture sampled values of f(l,0) */
        cube fl(N,T,n_iter);  /* used to build K x (N*T) matrices, f(l,0) */
        //field<mat> f(L,1);
        /* sum of all functions */
        //mat bb(n_iter,(N*T)); 
        bb.zeros(); /* N is fast moving */
        int noise        = 1; /* noise term for covariance matrix of the observations, y */
        int noise_term   = 0; /* noise term for covariance matrix of a functional term, f */
        
        /* FIRST, create total G_star of T x T cov mats for all (k,m) */
        for( k = 0; k < n_iter; k++ ) /* k denotes MCMC sampling iteration */
        {
               M_k       = Theta_star(k,0).n_cols;
               G_star_lk.set_size(T,T,M_k);
               G_star_k.set_size(T,T,M_k);
               invG_star_k.set_size(T,T,M_k);
               for( m = 0; m < M_k; m++ ) /* m denotes cluster */
               {
                    G_star_k.slice(m)        = gen_C(Theta_star(k,0).col(m), Tau_e(k), Omega_t, Omega_s,
                                                  jitter, gp_mod, n_parms, pos_s, noise);
                    invG_star_k.slice(m)     = inv_sympd(symmatl(G_star_k.slice(m)));
                    for( l = 0; l < L; l++ )  /* l denotes functional term */
                    {
                         start_l             = n_parms(l); /* indexing starts at 0 */
                         end_l               = start_l + P_vec(l) - 1;
                         /* number of covariance parameters for term l */
                         theta_star_kml.set_size( P_vec(l) );
                         theta_star_km                 = Theta_star(k,0).col(m); /* P x 1 */
                         theta_star_kml                = theta_star_km.subvec( start_l, end_l ); /*P(l)*/
                         G_star_lk.slice(m)            = gen_Cterm(theta_star_kml, Tau_e(k), Omega_t, 
                                                       Omega_s.slice(pos_s(l)), jitter, 
                                                       gp_mod(l), noise_term); /* no noise */
                         G_star_term(k,l).slice(m)     =  G_star_lk.slice(m);
                    } /* end loop l over function terms */                 
               } /* end loop m over clusters of G_star */
               invG_star(k,0)                = invG_star_k; /* each k cube is defined with M_k slices */     
          } /* end loop k over MCMC samples */
               
        /* SECOND, compose n_iter x (N*T) functions, f_l */ 
        for( l = 0; l < L; l++ )
        {
               /* set T x T x N cube that sums variances over K MCMC iterations */
               /* to zeros before filling */
               var_l.zeros();
               f(l,0).set_size(n_iter,(N*T));
               
               for( k = 0; k < n_iter; k++ )
               {
                    /* compute E{f(l,0)} */
                    for( i = 0; i < N; i++ )
                    {
                         /* conditional expectation for T x 1, f_{lki} */
                         mu_fkli                       = G_star_term(k,l).slice(S(k,i)) *
                                                            invG_star(k,0).slice(S(k,i)) *
                                                            y.row(i).t();
                         musamples_fl.slice(i).row(k)  = mu_fkli.t(); /* K x T x N */
                         /* T x T predictive variance for the test function of unit i */
                         var_lk.slice(i)      = G_star_term(k,l).slice(S(k,i))  - 
                                                    G_star_term(k,l).slice(S(k,i))*
                                                       invG_star(k,0).slice(S(k,i)) *
                                                       G_star_term(k,l).slice(S(k,i)).t(); /* T x t */
                    } /* end loop i over sample observations */
                    var_l     += var_lk;
                    var_lk.zeros();
               } /* end loop k over sampling iterations */
               var_l          /= n_iter; /* T x T x N */
 
               for( i = 0; i < N; i++ )
               {
                    /* compose N x T matrix of GP means, mu_fl */
                    mu_fl.row(i)    = mean( musamples_fl.slice(i), 0 ); /* by column */
                    for( j = 0; j < n_iter; j++ )
                    {
                         rmvncov(var_l.slice(i), mu_fl.row(i).t(), fli_draw);
                         fl.slice(j).row(i) = fli_draw.t(); /* N x T x K, fl */
                    } /* end loop over J samples generated from posterior distribution over f(l,0) */
               } /* end loop over N units for sampling K x (N*T),  f(l,0) */
               
               /* build functions, f(l,0) as K x (N*T) matrices, with N fast-moving */
               for( k = 0; k < n_iter; k++ )
               {
                    /* each fl.slice(k) is an N X T matrix */
                    f(l,0).row(k)  = vectorise( fl.slice(k) ).t(); /* by column makes N fast-moving */
               }
               bb                  += f(l,0); /* K x (N*T) */
        }/* end loop l over functions or terms */
        
        END_RCPP
} /* end function gen_bb */


SEXP gen_f(const mat& BB, const mat& Tau_e, umat& S, const field<mat>& Theta_star, const mat& Omega_t, 
               const cube& Omega_s, uvec& gp_mod, uvec& n_parms, uvec& pos_s, uvec& P_vec,
               field<mat>& f, const field<cube>& invG_star, double jitter)
{
        BEGIN_RCPP
        // rhotrend_star and ksi_star is of length, L, equal to total number of terms
        // only contains trend components
        // rhoseas_star, Omega_seas is of length equal to number of seasonal components
        // "seasonal" is of length L
        // seasonal = [0 0 1 0 2]' for L = 5 where non-zero entries index 
        // which seasonal entry in rhoseas_star, Omega_seas is associated to that trend term
        // int M denotes the total number of clusters
        /* indices and loop counters */
        int n_iter  = S.n_rows;
        int N       = S.n_cols;
        int L       = gp_mod.n_elem;
        int T       = Omega_t.n_rows;
        int P       = Theta_star(0,0).n_rows; /* all theta_star matrices have P rows */
        int i, l, k, m, M_k;
        
        //create objects
        field<cube> G_star_term(n_iter,L); /* by term T x T x M_k covariances */
        for( k = 0; k < n_iter; k++ )
        {
             M_k       = Theta_star(k,0).n_cols;
             for( l = 0; l < L; l++ )
             {
                  G_star_term(k,l).set_size(T,T,M_k);
             } /* end loop L over terms */
        } /* end loop k over sampling iterations */
        /* parameters used to compute G_star_lkm */
        colvec theta_star_kml; /* contains P(l) parameters */
        colvec theta_star_km(P); /* across all terms */
        /* functions drawn from GP */
        colvec h(T); h.zeros();
        colvec fkli(T); fkli.zeros(); 
        mat f_kl(N,T); f_kl.zeros(); 
        mat fl(n_iter,(N*T)); fl.zeros();
        /* by term indexing */
        int start_l, end_l;
        //field<mat> f(L,1);
        int noise_term   = 0; /* noise term for covariance matrix of a functional term, f */
        
        /* FIRST, create total G_star of T x T cov mats for all (k,m) */
        for( k = 0; k < n_iter; k++ ) /* k denotes MCMC sampling iteration */
        {
               M_k       = Theta_star(k,0).n_cols;
               for( m = 0; m < M_k; m++ ) /* m denotes cluster */
               {
                    for( l = 0; l < L; l++ )  /* l denotes functional term */
                    {
                         start_l             = n_parms(l); /* indexing starts at 0 */
                         end_l               = start_l + P_vec(l) - 1;
                         /* number of covariance parameters for term l */
                         theta_star_kml.set_size( P_vec(l) );
                         theta_star_km                 = Theta_star(k,0).col(m); /* P x 1 */
                         theta_star_kml                = theta_star_km.subvec( start_l, end_l ); /*P(l)*/
                         G_star_term(k,l).slice(m)     = gen_Cterm(theta_star_kml, Tau_e(k), Omega_t, 
                                                       Omega_s.slice(pos_s(l)), jitter, 
                                                       gp_mod(l), noise_term); /* no noise */
                    } /* end loop l over function terms */                 
               } /* end loop m over clusters of G_star */     
          } /* end loop k over MCMC samples */
               
        /* SECOND, compose n_iter x (N*T) functions, f_l */ 
        for( l = 0; l < L; l++ )
        {
               for( k = 0; k < n_iter; k++ )
               {
                    /* compute E{f(l,0)} */
                    for( i = 0; i < N; i++ )
                    {
                         /* conditional expectation for T x 1, f_{lki} */
                         fkli           = G_star_term(k,l).slice(S(k,i)) *
                                                  invG_star(k,0).slice(S(k,i)) *
                                                  BB.row(i).t();
                         f_kl.row(i)    = fkli.t(); /* N x T matrix */
                    } /* end loop i over sample observations */
                    
                         /* stack columns - fast-moving index is N, slow T */
                    fl.row(k)           = vectorise(f_kl,0).t(); /* stack columns */
               } /* end loop k over sampling iterations */
                    
               f(l,0)              = fl; /* fl is n_iter x (N*T) */
        }/* end loop l over functions or terms */
        
        END_RCPP
} /* end function gen_f */

SEXP predict_bb(SEXP res, SEXP o_Omegas_tetr, SEXP o_Omegas_tete, SEXP o_Omegat_tetr,
                SEXP o_Omegat_tete, SEXP o_J)
{
     BEGIN_RCPP
     // "res" is the List object returned by gpdpmix.cpp
     // we next unwrap the saved parameter objects needed for prediction
     List res_r(res);
     // loop variables
     int i, j, k;
     // read in invGstar
     List invGstar_r = res_r["invG_star"];
     int K          = invGstar_r.size(); /* Number of retained sampling iterations */
     field<cube> invG_star(K,1);
     for( k = 0; k < K; k++ )
     {
          NumericVector tmp            = invGstar_r[k];
          IntegerVector array_dims_k   = tmp.attr("dim");
          cube cubeArray_k(tmp.begin(), array_dims_k[0], array_dims_k[1], array_dims_k[2], false);
          invG_star(k,0) = cubeArray_k;
     }
     
     // read in cube of seasonal Euclidean distance cube, Omegas_tetr - of L_s slices of T_pred x T_train 
     NumericVector tmp_tetr_r(o_Omegas_tetr);
     IntegerVector array_dims_tetr = tmp_tetr_r.attr("dim");
     cube Omegas_tetr(tmp_tetr_r.begin(), array_dims_tetr[0], array_dims_tetr[1], 
                                   array_dims_tetr[2], false);
     
     // read in cube of seasonal Euclidean distance cube, Omegas_tete - of L_s slices of T_pred x T_pred 
     NumericVector tmp_tete_r(o_Omegas_tete);
     IntegerVector array_dims_tete = tmp_tete_r.attr("dim");
     cube Omegas_tete(tmp_tete_r.begin(), array_dims_tete[0], array_dims_tete[1], 
                                   array_dims_tete[2], false);
     
     // read in T_pred x T_train trend Euclidean distance matrix
     NumericMatrix Otr(o_Omegat_tetr);
     int T_pred          = Otr.nrow(); /* number of predicted time points */
     int T               = Otr.ncol(); /* n_cols = T, number of training cases */
     mat Omegat_tetr(Otr.begin(),T_pred,T,false);  
     
     // read in T_test x T_train trend Euclidean distance matrix
     NumericMatrix Ote(o_Omegat_tete);
     mat Omegat_tete(Ote.begin(),T_pred,T_pred,false); 

     // read in Theta(K,(P*N)) - P is fast-moving
     mat Theta           = as<mat>(res_r["Theta"]);
     
     // read covariance function uvecs
     List gp_r           = res_r["gp_indices"];
     uvec gp_mod         = as<uvec>(gp_r["gp_mod"]);
     /* vector of length L where positive values index count for seasonal covariance matrices */
     uvec pos_s          = as<uvec>(gp_r["pos_s"]);
     /* cumulative sum of covariance parameters counts, offset by 1 */
     uvec n_parms        = as<uvec>(gp_r["n_parms"]);
     
     // read in K x N cluster assignment matrix, S
     List op_r           = res_r["optpartition"];
     umat S              = as<umat>(op_r["S"]);
     mat y               = as<mat>(op_r["y"]); /* read in the N x T training data, y */
     int N               = S.n_cols; /* Number of experimental units */
     int P               = Theta.n_cols / N; /* number of covariance parameters */
     int J               = as<int>(o_J); /* number of draws from posterior predictive distributions */
     
     // extract K x (T*N) bb to kit together with prediction
     // later used for plotting predicted output
     mat bb              = as<mat>(res_r["bb"]);
     
     /* Code to fill in predicted function values, bb_pred(K,T,N) */
     /* where K denotes number of sampling iterations, T the number of time points */
     /* and N, the number of units */
     
     /* temporary objects for sampling draw, k, and unit i */
     cube var_tete_k(T_pred,T_pred,N), var_tete_sumk(T_pred,T_pred,N);
     mat G_tetr_ki(T_pred,T), G_tete_ki(T_pred,T_pred); 
     colvec f_ki(T_pred), bi_draw(T_pred);
     mat bb_k(N,T); /* estimated functions from training data for MCMC iteration k */
     /* return objects */
     cube bb_pred(K,T_pred,N); /* each slice has the set of K sampled curves (of length T_pred) for unit i */
     cube bbpred_draws(J,T_pred,N); /* each slice has the set of J sampled curves (of length T_pred) */
     /* from the posterior predictive distribution for unit i */
     cube var_tete(T_pred,T_pred,N); /* compute predictive distribution variance for each MCMC draw, k */
     mat E_bb_pred(N,T_pred); /* average over the K samples at each time point, t, for the N units */
     mat theta_k(P,N); /* matrix of covariance parameters for sampling draw k */
     
     /* set T_pred x T_pred x N cube that sums variances over K MCMC iterations to zeros before filling */
     var_tete_sumk.zeros(); var_tete_k.zeros();
     for( k = 0; k < K; k++ )
     {
          /* snap 1 x P*N vector into a P x N matrix - P is fast-moving */
          theta_k       = reshape( Theta.row(k), P, N ); /* by column */
          bb_k          = reshape( bb.row(k), N, T ); /* N is fast-moving */
          /* compute f(K,T).slice(i) */
          for( i = 0; i < N; i++ )
          {
               /* memo: Omegat_tetr is T_pred x T */
               G_tetr_ki           = gen_Casym(theta_k.col(i), Omegat_tetr, Omegas_tetr,
                                                  gp_mod, n_parms, pos_s); 
               /* tete with no noise or jitter */
               G_tete_ki           = gen_C(theta_k.col(i), 1e300, Omegat_tete, Omegas_tete,
                                                  0.0, gp_mod, n_parms, pos_s, 0);
               /* conditional expectation for T_pred x 1, f_{lki} */
               if( any(vectorise(y) == -9) ) /* any(is.na(y)) in R */
               { /* sample from de-noised estimates (obtained from a Gibbs sampler) */
                 /* estimates f_test | y by marginalizing over posterior distributions, given y */
                    f_ki                          = G_tetr_ki * invG_star(k,0).slice(S(k,i)) *
                                                  bb_k.row(i).t(); /* T x 1 */
               }else{ /* sample from noisy observations, so invG_star has noise, tau_e */
                    f_ki                          = G_tetr_ki * invG_star(k,0).slice(S(k,i)) *
                                                  y.row(i).t(); /* T x 1 */
               } /* end conditional statement on whether bb sampled or marginalized out of sampler */
               
               /* T_pred x T_pred predictive variance for the test function of unit i */
               var_tete_k.slice(i)           = symmatl(G_tete_ki  - G_tetr_ki*invG_star(k,0).slice(S(k,i)) *
                                                             G_tetr_ki.t()); /* T_pred x T_pred */
               bb_pred.slice(i).row(k)        = f_ki.t(); /* K x T_pred matrix for each i */
          } /* end loop i over sample observations */
          var_tete_sumk       += var_tete_k;
          /* re-set var_tete_k to zeros */
          var_tete_k.zeros();
     }/* end loop K over MCMC draws */
     /* divide sumk by K to get T_pred x T_pred x N rao-blackwellized posterior predictive var matrices */
     var_tete       = var_tete_sumk / K; /* T_pred x T_pred x N */
     /* compute N x T_pred, E_f */
     for( i = 0; i < N; i++ )
     {
          /* bb_pred.slice(i) is a K x T_pred matrix */
          E_bb_pred.row(i)     = mean( bb_pred.slice(i), 0 ); /* 1 x T_pred column means */
          for( j = 0; j < J; j++ )
          {
               rmvncov(var_tete.slice(i), E_bb_pred.row(i).t(), bi_draw);
               bbpred_draws.slice(i).row(j) = bi_draw.t();
          }
     }
     
     return Rcpp::List::create(Rcpp::Named("bbpred_draws")       = bbpred_draws,
                               Rcpp::Named("var_tete")           = var_tete,
                               Rcpp::Named("E_bb_pred")          = E_bb_pred,
                               Rcpp::Named("bb_pred")            = bb_pred,
                               Rcpp::Named("bb")                 = bb
                               );
     END_RCPP
}


SEXP predict_gmrf_bb(SEXP res, SEXP o_R, SEXP o_J)
{
     BEGIN_RCPP
     // loop variables
     int i, k, m, j, p;
     
     // "res" is the List object returned by gpdpmix.cpp
     // we next unwrap the saved parameter objects needed for prediction
     List res_r(res);
     
     // read in P x N cluster assignment matrix, S
     List op_r           = res_r["optpartition"];
     umat S              = as<umat>(op_r["S"]);
     mat y               = as<mat>(op_r["y"]); /* read in the N x T training data, y */
     int N               = S.n_cols; /* Number of experimental units */
     int T_train         = y.n_cols; /* Number of training time points */
     int J               = as<int>(o_J);
       
     // extract P x (T_train*N) estimated functions, bb, to kit together with prediction
     // later used for plotting predicted output. N is fast-moving
     mat bb               = as<mat>(res_r["bb"]); /* read in the P x N*T_train estimated functions */
     
     // read in P length list, Kappa_star, where each entry is K x M_p
     List Kappastar_r    = res_r["Kappa_star"];
     int P               = Kappastar_r.size(); /* Number of retained sampling iterations */
     field<mat> Kappa_star(P,1);
     for( p = 0; p < P; p++ )
     {
          Kappa_star(p,0)     = as<mat>(Kappastar_r[p]); /* K x M_p */
     }
     
     /* GMRF precision matrices */
     List R_r(o_R); /* K, T_tot x T_tot Precision matrices.  Order aligned with rows of Kappa_star(p) */
     int K                    = R_r.size(); /* Number of iGMRF terms */
     int T_tot                = as<mat>(R_r(0)).n_rows; /* T_tot = T_train + T_test */
     cube R(T_tot,T_tot,K); /* each entry is a T_tot x T_tot precision matrix for term k \in (1,..,K) */
     /* estimated functions, f[[k]] */
     List f_r    = res_r["f"];
     field<mat> f(K,1), f_p(K,1); /* f(K,1) is filled with P x (N*T_train) matrices */
                                /* f_p(K,1) is filled with N x T_train matrices for each p \in (1,..,P) */
     for( k = 0; k < K; k++ )
     {
          R.slice(k)          = as<mat>(R_r[k]);
          f(k,0)              = as<mat>(f_r[k]); /* each f(k,0) is P x (N*T_train), N is fast-moving */
     }
     int T_pred               = T_tot - T_train;
     
     colvec Tau_e             = as<colvec>(res_r["Tau_e"]);
     

     /* Code to fill in predicted function values, bb_pred(P,T_pred,N) */
     /* where P denotes number of sampling iterations, T_test the number of prediction time points */
     /* and N, the number of units */
     
     /* temporary objects for sampling draw, p, and unit i */
     field<mat> Qstar_p, Ustar_tete_p, Qstar_tete_p, Qstar_tetr_p; /* for each (p,m) */
     cube Qtete_p(T_pred,T_pred,(K*N)); /* T_pred x T_pred terms for each (p,k*i) */
     cube Q_tete(T_pred,T_pred,(K*N)); /* average over P, T_pred x T_pred terms for each (k*i) */
     mat Qtetr_pki(T_pred,T_train);
     mat kappastar_p; /*kappastar_p is  K x M_p */
     /* T_pred x 1 vectors where v_pki is an intermediate product and mu_pki is the predicted mean */
     colvec mu_pki(T_pred), v_pki(T_pred), bki_draw(T_pred);
     /* return objects */
     /* Predictive mean, E_bb_pred */
     cube f_pred(P,T_pred,(K*N)); /* each slice has the set of P sampled curves (of length T_pred) */
                                 /* for unit i and GMRF term, k */
     cube bb_pred(P,T_pred,N); /* each slice has the set of P sampled curves (of length T_pred) for unit i */
     bb_pred.zeros(); /* will sum over K GMRF terms to get total prediction */
     /* from the posterior predictive distribution for unit i */
     mat E_f_pred((K*N),T_pred); /* average over the P samples at each time point, t, for the N units */
                                 /* for each of K GMRF terms */
     mat E_bb_pred(N,T_pred); /* average over the P samples at each time point, t, for the N units */
     /* credible intervals - draws */
     cube fpred_draws(J,T_pred,(K*N)); /* each (k,i) slice has the set of J sampled curves */
                                        /* (of length T_pred) */
     fpred_draws.zeros();
     cube bbpred_draws(J,T_pred,N); /* each slice has the set of J sampled curves (of length T_pred) */
     bbpred_draws.zeros(); /* since filling by summation over K terms, not setting exactly equal to */
     
     int M_p; /* number of clusters on MCMC iteration p */
     /* set T_pred x T_pred x (K*N) cube that sums precisions over P MCMC iterations */
     /* to zeros before filling */
     Q_tete.zeros(); Qtete_p.zeros(); /* T_pred x T_pred x (K*N) */
     field<mat> f_bar(K,1);
     for( k = 0; k < K; k++ ) /* loop over K GMRF terms */
     {
          f_bar(k,0)                      = reshape( mean(f(k,0),0), N, T_train ); /* by column */
     }
     
     for( p = 0; p < P; p++ ) /* loop over P MCMC sampler iterations */
     {
          kappastar_p       = Kappa_star(p,0); /* K x M_p matrix */
          M_p               = kappastar_p.n_cols;
          Qstar_p.set_size(M_p,K); Ustar_tete_p.set_size(M_p,K); Qstar_tete_p.set_size(M_p,K);
          Qstar_tetr_p.set_size(M_p,K);
          for( k = 0; k < K; k++ ) /* loop over iGMRF terms */
          {
               /* snap 1 x T_train*N vector into a N x T_train matrix - N is fast-moving */
               f_p(k,0)                      = reshape( f(k,0).row(p), N, T_train ); /* by column */
               for( m = 0; m < M_p; m++ ) /* iterate over K precision terms */
               {
                    Qstar_p(m,k).set_size(T_tot,T_tot);
                    Qstar_p(m,k)             = R.slice(k) * kappastar_p(k,m);
                    Qstar_tete_p(m,k)        = Qstar_p(m,k).submat(span(T_train,(T_tot-1)),
                                                       span(T_train,(T_tot-1)));
                    Qstar_tetr_p(m,k)        = Qstar_p(m,k).submat(span(T_train,(T_tot-1)),
                                                       span(0,(T_train-1)));
                    Ustar_tete_p(m,k)        = trimatu(chol(symmatl(Qstar_tete_p(m,k))));
               }  /* end loop over M_p clusters */      
          } /* end loop over K iGMRF terms */ 
          
          /* compute mu(P,T_pred).slice(i) */
          for( i = 0; i < N; i++ )
          {
               for( k = 0; k < K; k++ ) /* iterate over K precision terms */
               {
                    /* memo: Qtetr_pi is T_pred x T_train */
                    Qtetr_pki                     = Qstar_tetr_p(S(p,i),k);
                    /* T_pred x T_pred matrix for each (p,k,i) */
                    /* fill the T_pred x T_pred x (N*K) cube for each p */
                    Qtete_p.slice((K*i+k))       = Qstar_tete_p(S(p,i),k); /* k is fast-moving */
                    /* conditional expectation for T_pred x 1, mu_pi */
                    //v_pki           = solve( trimatu(Ustar_tete_p(S(p,i),k)).t() , 
                    //                         (-Qtetr_pki*f_bar(k,0).row(i).t()) );
                    v_pki           = solve( trimatu(Ustar_tete_p(S(p,i),k)).t() , 
                                             (-Qtetr_pki*f_p(k,0).row(i).t()) );
                    //v_pki           = solve( trimatu(Ustar_tete_p(S(p,i),k)).t() , 
                    //                         (-Qtetr_pki*y.row(i).t()) );
                    mu_pki          = solve( trimatu(Ustar_tete_p(S(p,i),k)), v_pki ); /* T_pred x 1 */
                    f_pred.slice((K*i+k)).row(p) = mu_pki.t(); /* row p of P x T_pred matrix for each k,i */
                    /* summed over K iGMRF terms for each (p,i) */
                    /* row p of P x T_pred matrix for each i, summing over K iGMRF terms */
                    /* not summing over P MCMC iterations */
                    bb_pred.slice(i).row(p)  += f_pred.slice((K*i+k)).row(p); /* (i,p) fixed */ 
               } /* end loop over K iGMRF terms */
          } /* end loop i over sample observations */
          Q_tete     += Qtete_p; /* sum over P, T_pred x T_pred x (K*N) cubes */
          Qtete_p.zeros(); /* T_test x T_test x (K*N) cube for each MCMC iteration, p. Refill it. */
     }/* end loop p over MCMC draws */
     /* compute ergodic average for T_pred x T_pred x (K*N) cube, Q_tete */
     /* for the summation over the P draws */
     Q_tete /= P;
     
     /* compute N x T, E_f */
     for( i = 0; i < N; i++ )
     {
          /* bb_pred.slice(i) is a P x T_pred matrix */
          E_bb_pred.row(i)     = mean( bb_pred.slice(i), 0 ); /* 1 x T_pred column means */
          
          for( k = 0; k < K; k++ ) /* iterate over K iGMRF terms */
          {    /* compose rao-blackwellized expectation over the P MCMC draws */
               E_f_pred.row((K*i+k)) = mean( f_pred.slice((K*i+k)), 0 ); /* 1 x T_pred */
               for( j = 0; j < J; j++ )
               {
                    rmvnsample(Q_tete.slice((K*i+k)), E_f_pred.row((K*i+k)).t(), bki_draw);
                    //bki_draw += randn<colvec>(T_pred)*sqrt(1/mean(Tau_e));
                    fpred_draws.slice((K*i+k)).row(j) = bki_draw.t();
               } /* end loop j over J sample draws from posterior predictive */
               /* sum over K individual J x T_pred draws for unit i */
               bbpred_draws.slice(i)    += fpred_draws.slice((K*i+k)); /* J x T_pred for eack (k,i) */
          } /* end loop over K iGMRF terms */       
     } /* end loop i over N units */
     
     return Rcpp::List::create(Rcpp::Named("bbpred_draws")       = bbpred_draws,
                               Rcpp::Named("Q_tete")             = Q_tete,
                               Rcpp::Named("E_bb_pred")          = E_bb_pred,
                               Rcpp::Named("E_f_pred")           = E_f_pred,
                               Rcpp::Named("bb")                 = bb
                               );
     END_RCPP
}

