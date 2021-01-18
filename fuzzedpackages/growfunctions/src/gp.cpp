#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

SEXP GP(SEXP Ymat, SEXP Otrend, SEXP Oseas, SEXP o_gp_mod, SEXP o_jitter,
        SEXP o_a, SEXP o_b, SEXP o_atau, SEXP o_btau, SEXP o_lower, SEXP o_upper,
        SEXP o_w, SEXP o_n_slice_iter, 
        SEXP o_y_index, SEXP niterInt, SEXP nburnInt, SEXP nthinInt, 
        SEXP ntuneInt, SEXP o_progress, SEXP o_ipr)
{
    BEGIN_RCPP
    // Time run
    // wall_clock timer;
    // timer.tic();
    // Wrap SEXP objects in Rcpp Vector and Matrix objects.  Not copies,
    // but wrappers.
    NumericMatrix Yr(Ymat); /* N x T data matrix */
    NumericMatrix Otr(Otrend); /* T x T trend covariance matrix */
    List Osr(Oseas); /* L_s T x T seasonal covariance matrices */
    List yr_index(o_y_index); /* of length n, number of tempered steps */
    NumericVector gpr_mod(o_gp_mod); /* covariance kernel indicators for L covariance terms */
    NumericVector ipr_r(o_ipr);
    vec ipr   = as<vec>(ipr_r); /* Inclusion probabilities to produce Horvitz-Thompson plug-in */
    int ntune = as<int>(ntuneInt); /* tuning iterations for slice widths, (wp,wtau) */
    int niter = as<int>(niterInt); /* sampler iterations */
    int nburn = as<int>(nburnInt); /* sampler burnin */
    int nthin = as<int>(nthinInt); /* sampler thinning */
    int nkeep = (niter -  nburn)/nthin; 
    unsigned int M     = 1; /* fix number of clusters to M = 1 */
    double a = as<double>(o_a); /* GP parameters (theta_star, tau) hyperparameters */
    double b = as<double>(o_b); /* GP parameters (theta_star, tau) hyperparameters */
    double atau = as<double>(o_atau); /* Noise precision, tau_e, hyperparameters */
    double btau = as<double>(o_btau); /* Noise precision, tau_e, hyperparameters */
    double lower = as<double>(o_lower); /* GP parameters (theta_star, tau) range */
    double upper = as<double>(o_upper); /* GP parameters (theta_star, tau) range */
    double jitter = as<double>(o_jitter); /* GP cov matrix jitter */
    int n_slice_iter = as<int>(o_n_slice_iter); /* limit on num steps to increase slice width */
    double w = as<double>(o_w); /* step size for increasing slice sampling width */
    int progress = as<int>(o_progress); /* indicator for whether to display a progress bar */
    
    // Extract key row and column dimensions we will need to sample parameters
    int N      = Yr.nrow(), T = Yr.ncol(); /* data dimensions */
    int n      = yr_index.size(); /* number of progressively coarser tempered distributions */
    int L_s    = Osr.size(); /* number of seasonal covariance terms */
    int L      = gpr_mod.length(); /* number of covariance terms */
    
    // define (y,Omega_t,Omega_s) for full data 
    mat y(Yr.begin(), N, T, false);
    mat Omega_t(Otr.begin(), T, T, false);
    cube Omega_s(T,T,L_s);
    int i, j, k; /* initialize loop variables */
    for(j = 0; j < L_s; j++)
    {
          Omega_s.slice(j) = as<mat>(Osr[j]);
    }
    
    // define field objects, (y_ns, Omegas_ns, Omegat_ns), for n coarse tempered steps 
    field<uvec> index(n,1); /* indexes columns (1:T) in y for sub-sampling time points */
    field<mat> Omegat_ns(n,1); /* each entry is a progressively coarser, T_i x T_i trend covariance */
    field<cube> Omegas_ns(n,1); /* each entry is a progressively coarser, T_i x T_i x L_s seas cube */
    field<mat> y_ns(n,1); /* each entry is a progressively coarse, N x T_i data matrix */
    int T_i;
    
    for(i = 0; i < n; i++)
    {
         index(i,0)           = as<uvec>(yr_index[i]);
         T_i                  = index(i,0).n_elem;
         /* set size for element i in field<mat> Omegat_ns, field<cube> Omegas_ns, field<colvec> y_ns */
         Omegat_ns(i,0).set_size( T_i, T_i );
         Omegas_ns(i,0).set_size( T_i, T_i, L_s );
         y_ns(i,0).set_size(N, T_i);
         Omegat_ns(i,0)       = Omega_t.submat( index(i,0), index(i,0) );
         y_ns(i,0)            = y.cols( index(i,0) );
         for(j = 0; j < L_s; j++)
         {
              Omegas_ns(i,0).slice(j)   = Omega_s.slice(j).submat( index(i,0), index(i,0) );
         } /* end loop j over seasonal covariance terms */  
    } /* end loop i over number of tempered distributions - steps */
    
    // generate total number of parameters, P, across all covariance terms 
    /* P_vec of length L contains the number of parameters for each term */
    /* based on the term identified for each entry in gpr_mod */
    /* memo: since clustering units, i = 1,...,N, will have M clusters, each of P parameters */
    //uvec gp_mod(gpr_mod.begin(),L,false); /* term of length L identifies cov formulations */
    uvec gp_mod	= as<uvec>(gpr_mod);
    uvec P_vec(L); P_vec.zeros();
    gen_P(gp_mod, P_vec); /* returns length L, P_vec, of parameter counts per covariance term */
    int P = sum( P_vec );
    
    //compute starting position in 1:P for each GP covariance term in C(theta_star_m)
    /* shift cumsum length up 1 so that use last value in current loop */
    /* e.g. P_vec = c(2,2,3,2) -> n_parms = c(0,2,4,7,9) */
    uvec n_parms(L+1); n_parms.zeros();
    n_parms.subvec(1,L)    = cumsum( P_vec ); /* 1st element is 0, last element unused */
    
    // index seasonal covariance matrices, Omega_s.slices(1,L_s)
    /* index of seasonal covariances, pos_s, of length L */
    /* e.g. j = 0,1,2,3,4.  j = 1,3,4 are seasonal */
    /* then pos_s = 0,0,0,1,2 and when gp_mod(j) == 3 -> pos_s | gp_mod == 3  = (0,1,2) */
    uvec ind_seas        = find( gp_mod == 3 ); /* e.g. (1,3,4) */
    int n_seas           = ind_seas.n_elem;
    int count_s          = 0; /* index increment of seasonal terms in L x 1, pos_s */
    uvec pos_s(L); pos_s.zeros();  /* non-seasonal entries hold 0 */
    if( n_seas > 0 ) /* there are seasonal terms */
    {
          for( i = 0; i < n_seas; i++ )
          {
               pos_s( ind_seas(i) ) = count_s; /* count starts at 0 */
               count_s              += 1;
          }
    } 
    // okay that pos_s has duplicate 0's b/c only Omega_s.slice(pos_s(j)) 
    // only used if gp_mod == 3; else, Omega_s.slice(0) is input to functions (such as gen_C())
    // but it is never actually used to compute C.
    
 
    // Set random number generator state
    RNGScope scope; /* Rcpp */
    //arma_rng::set_seed_random(); /* arma */

    // Initialize SAMPLED parameter values   
    /* cluster capture variables */
    int reps                  = N/M; /* truncato resulto */
    IntegerVector clust       = seq_len(M); /* 1:M */
    IntegerVector s_          = rep_each(clust,reps);  /* 111222...MMM */
    ucolvec s                 = as<ucolvec>(s_);
    /* ensure there are N objects slotted into M clusters */
    /* if N % M = R, then the R units are allocated cluster M */
    ucolvec srest((N-M*reps)); srest.ones(); srest *= M;
    s	                    = join_cols(s,srest); s -= 1; /* use s as position selector */
    s 	                    = shuffle(s,0); 
  
    /* cluster locations, theta_star */
    NumericVector _theta_vec            = rgamma( (P*M), a, (1/b) );
    vec theta_vec                       = as<vec>(_theta_vec);
    mat theta_star                      = reshape( theta_vec, P, M );
    /* since M changes iteration-to-iteration, use theta to capture parms */
    mat theta(P,N); theta.zeros();
    /* associated cholesky decompositions of C_star using theta_star */
    /* U_last is repeatedly updated in temper() and then used to */
    /* compute likelihoods for s(i) under auxclusterstep() */
    double tau_e    = rgamma(1, a, (1/b))[0];
    cube U_last(T,T,M);
    int noise       = 1; /* add noise since integrated over {bb_i}, GP functional draws */
    compute_U(theta_star, tau_e, jitter, gp_mod, n_parms, pos_s, Omega_t, Omega_s,
                    noise, U_last);
                    
    // initialize sampler tuning parameters for (theta_star, tau_e)
    /* set (wp,wtau) = w */
    mat wpm(P,M); wpm.ones(); wpm *= w;
    double wtau = w;
    mat w_tot((P+1),M); w_tot = update_w(wpm, wtau);
    
    // initialize counter variables for performance of uni.slice() and temper() 
    double slice_levals_theta = 0; /* number of theta_star(p,m) slice sample like computes */
    double slice_levals_tau   = 0; /* number of tau slice sample like computes */
    double n_slice_theta      = 0; /* number of slice sampling steps for theta_star(p,m) */
    double n_slice_tau        = 0; /* number of slice sampling steps for tau */
    double accept_temper      = 0; /* number of temper() acceptances */
    double n_temper_evals     = 0; /* number of temper() steps */
    
    /* initialize vector of residuals */
    colvec resid(N*T); resid.zeros();
    /* fit assessment - related measures */
    double deviance = 0;  ucolvec ordscore(nkeep); ordscore.zeros();
    mat phat(N,N); phat.zeros(); /* empirical co-clustering probability matrix */
    rowvec devmarg(N); colvec devres(4); devres.zeros();
    rowvec logcpo(N); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0,  kk = 0;
    /* samples to update tuning parameters */
    int nn = 0;
   
    // Armadillo structures to capture results
    // Will be wrapped into a list element on RETURN
    // DP return objects
    ucolvec numM(nkeep); /* M */
    umat S(nkeep,N); /* cluster indicator  vectors of length N */
    field<ucolvec> bigS(1,1); /* create as dummy to support plotting */
    bigS(0,0) = cumsum(ones<ucolvec>(N)); /* all the same cluster */
    /* locations */
    mat Theta(nkeep,P); /* P parameters global for all of N units */
    field<mat> Theta_star(nkeep,1); /* each will hold a P x M_k matrix used to compute invG_star */
    field<cube> invG_star(nkeep,1); /* a T x T x M_k cube of sampled inverse C_star locations */
    /* use invG_star for prediction */
    /* return indices describing the size characteristics of the P covariance parameters (per cluster) */
    /* non-cluster parameters */
    colvec Tau_e(nkeep);
    /* draws for GP function using theta */
    mat bb(nkeep,(N*T)); /* filled by gen_bb() */
    /* draws for each additive GP function, (f_1,...,f_L) */
    field<mat> f(L,1); /* filled by gen_bb() */
    mat Resid(nkeep,(N*T));
    /* fit indicators */
    colvec Deviance(nkeep); Deviance.zeros();  
    mat Devmarg(nkeep,(N*T));     
    
    // TUNING sample run
    int n_blocks    = 1;
    int block       = floor( ntune / n_blocks);
    cube Theta_tune(P,M,ntune); /* set of sampled parameters used to tune slice step width */
    cube Theta_tune_block(P,M,block); /* sub-block of Theta_tune for each of n_block tuning updates */
    colvec Taue_tune(ntune); /* set of parameters used to tune slice sampling step width */
    colvec Taue_tune_block(block); /* sub-block of Taue_tune for each of n_block tuning updates */
    for( k = 0; k < ntune; k++ )
    {
         if(progress == 1)
         {
              if( (k % 90) == 0 ) Rcout << "Tuning Interation: " << k << endl;
         }
        /* M changes iteration-to-iteration, so dimension of theta_star will update */
          temper(U_last, theta_star, tau_e, gp_mod, n_parms, pos_s, slice_levals_theta,
                    slice_levals_tau, n_slice_theta, n_slice_tau, accept_temper, n_temper_evals, 
                    lower, upper, n_slice_iter, w_tot, s, Omegat_ns, Omegas_ns, 
                    Omega_t, Omega_s, y_ns, y, noise, jitter, a, b, atau, btau, ipr); /* noise set to 1 */

          /* capture parameter updates for tuning (wp,w) */
          Theta_tune.slice(k)          = theta_star;
          Taue_tune(k)                 = tau_e;    
          if( k ==  (nn+1)*block-1 ) /* update (wp,w) n_blocks times during tuning */
          {
               // compute p x 1 tuned slice sampling step width parameters, (wp,wtau)
               Theta_tune_block       = Theta_tune.slices( block*nn, (block*(nn+1)-1) );
               Taue_tune_block        = Taue_tune.subvec( block*nn, (block*(nn+1)-1) );
               wpm_tune(Theta_tune_block, wpm); /* min for p of avg for K samples for each (p,m) */
               wtau_tune(Taue_tune_block, wtau);
               /* reset P x 1, w_tot */
               w_tot = update_w(wpm, wtau);
               nn      += 1;     
          }  /* end a block update of (wp,wtau) */     
        
    } /* end loop k on tuning sampler runs to update (wp,wtau) */
    
    
    // POSTERIOR samples run with tuned and fixed (wp, wtau)
    /*  initialize counter variables for performance of uni.slice() and temper() */
    slice_levals_theta = 0; /* number of theta_star(p,m) slice sample like computes */
    slice_levals_tau   = 0; /* number of tau slice sample like computes */
    n_slice_theta      = 0; /* number of slice sampling steps for theta_star(p,m) */
    n_slice_tau        = 0; /* number of slice sampling steps for tau */
    accept_temper      = 0; /* number of temper() acceptances */
    n_temper_evals     = 0; /* number of temper() steps */
    /* capture results over all saved iterations */
    double accept_temper_global = 0;
    double n_temper_evals_global = 0;
    double avg_slice_levals_theta = 0;
    double avg_slice_levals_tau = 0;
    double n_slice_theta_global = 0;
    double n_slice_tau_global = 0;
    
    kk = 0;
    /* sampling iterations */
    for( k = 0; k < niter; k++ )
    {
        if(progress == 1)
        {
             if( (k % 450) == 0 ) Rcout << "Production Interation: " << k << endl;
        }
        /* M changes iteration-to-iteration, so dimension of theta_star will update */
        temper(U_last, theta_star, tau_e, gp_mod, n_parms, pos_s, slice_levals_theta,
                    slice_levals_tau, n_slice_theta, n_slice_tau, accept_temper, n_temper_evals, 
                    lower, upper, n_slice_iter, w_tot, s, Omegat_ns, Omegas_ns, 
                    Omega_t, Omega_s, y_ns, y, noise, jitter, a, b, atau, btau, ipr);
        
        if(k >= nburn)
        {
            kk = k - nburn;
            if( kk == ((oo+1)*nthin - 1) )
            {
               /* monitor average chain acceptance rate */
     
               /* capture P x N parameter matrix for sampling iteration k */
               /* p is the fast-moving index, N is the slow-moving index */
               Theta.row(oo)                 = vectorise( theta_star, 0 ).t(); /* global P x 1 GP parms */
               Theta_star(oo,0)              = theta_star;
               Tau_e(oo)                     = tau_e;
               numM(oo)                      = M;
               S.row(oo)                     = s.t();
               /* compute chain acceptance statistics */ 
               accept_temper_global          += accept_temper;
               n_temper_evals_global         += n_temper_evals;
               avg_slice_levals_theta        += slice_levals_theta;
               avg_slice_levals_tau          += slice_levals_tau;
               n_slice_theta_global          += n_slice_theta;
               n_slice_tau_global            += n_slice_tau;
               oo += 1; /* increment sample return counter */
            } /* end conditional statement on returning sample */
        } /* end if k > burnin, record results for return */
    } /* end MCMC loop over k */
    
    // compose chain diagnostic return objects
    accept_temper_global      /= (n_temper_evals_global);
    avg_slice_levals_theta    /= (n_slice_theta_global);
    avg_slice_levals_tau      /= (n_slice_tau_global);
    List chaintune = Rcpp::List::create(Rcpp::Named("accept_temper_global")    
                                                            = accept_temper_global,
                                        Rcpp::Named("avg_slice_levals_theta")    
                                                            = avg_slice_levals_theta,
                                        Rcpp::Named("avg_slice_levals_tau")    
                                                            = avg_slice_levals_tau
                                           );
                         /* L x 1 vector that counts number of cov parms in each cov term */
    List gp_indices = Rcpp::List::create(Rcpp::Named("P_vec")         = P_vec,
                         /* (L+1) x 1 vector that composes cumsum of P_vec, shifted right by one */
                                           Rcpp::Named("n_parms")     = n_parms,
                         /* L x 1 vector where positive values denote index of seasonal terms */
                                           Rcpp::Named("pos_s")       = pos_s,
                         /* L x 1 vector \in {1,2,3} indicating form for covariance function */
                                           Rcpp::Named("gp_mod")      = gp_mod
                                           );
    
    // compute nkeep x (N*T) fitted functions, bb, where N is the fast moving indx
    gen_bb(y, Tau_e, S, Theta_star, Omega_t, Omega_s,  
               gp_mod, n_parms, pos_s, P_vec, bb, f, invG_star, jitter);
    if(progress == 1)
    {
         Rcout << "Finished generating estimated GP functions, bb "  << endl;
    }
    
    // compute FIT statistics
    //sample cluster functions, bb[i,1:T], and use to compute LPML
    for( k = 0; k < nkeep; k++ )
    {
         /* N is the fast-moving index, T is the slow-moving index */
         resid           = vectorise( y ); /* y is an N x T mat */
         resid           -= bb.row(k).t(); /* resid is a colvec */
         deviance        =  dev(resid, tau_e); /* scalar double deviance */
         dmarg(resid, tau_e, devmarg); /* 1 x N vector of densities*/
         Deviance(k)     = deviance;
         Devmarg.row(k)  = devmarg;
         Resid.row(k)    = resid.t();
    }
    
    // DIC
    dic3comp(Deviance, Devmarg, devres); /* devres = c(dic,dbar,dhat,pd) */
    cpo(Devmarg, logcpo, lpml);

    List optpartition = Rcpp::List::create(Rcpp::Named("S")           = S,
                                           /* N x T data matrix need for plotting fit */
                                           Rcpp::Named("y")           = y
                                           );

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance")            = Deviance,
                                  Rcpp::Named("Devmarg")         = Devmarg,
                                  Rcpp::Named("devres")          = devres,
                                  Rcpp::Named("logcpo")          = logcpo,
                                  Rcpp::Named("lpml")            = lpml,
               	              Rcpp::Named("Theta")           = Theta,
                                  Rcpp::Named("invG_star")       = invG_star,
				              Rcpp::Named("Chain_tune")      = chaintune,
                                  Rcpp::Named("gp_indices")      = gp_indices,
                                  Rcpp::Named("w_tot")           = w_tot,
                                  Rcpp::Named("bb")              = bb,
                                  Rcpp::Named("f")               = f,
                                  Rcpp::Named("Tau_e")           = Tau_e,
                                  Rcpp::Named("Residuals")       = Resid,
                                  Rcpp::Named("M")               = numM,
                                  Rcpp::Named("optpartition")    = optpartition,
                                  Rcpp::Named("bigSmin")         = bigS
				  );

    // Print CPU runtime
     //double n_secs = timer.toc();
     //cout << "took " << n_secs << " seconds" << endl;
END_RCPP
} /* end MCMC function returning SEXP */
