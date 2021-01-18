#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

SEXP IGMRFDPMIXCOUNT(SEXP Ymat, SEXP o_E, SEXP o_ksi, SEXP o_C, SEXP o_D, SEXP o_order,  
        SEXP niterInt, SEXP nburnInt, SEXP nthinInt,
        SEXP Minit, SEXP o_w_star, SEXP o_a, SEXP o_b, SEXP o_a_tau, SEXP o_b_tau,
        SEXP shapealph, SEXP ratebeta, SEXP o_nu, SEXP o_Rep, 
        SEXP o_progress, SEXP o_jitter, SEXP o_kappa_fast, SEXP o_stable_launch,
        SEXP o_ipr)
{
    BEGIN_RCPP
    // Time run
    // wall_clock timer;
    // timer.tic();
    // Wrap SEXP objects in Rcpp Vector and Matrix objects.  Not copies,
    // but wrappers.
    NumericMatrix Yr(Ymat); /* N x T data matrix */
    /* N x R predictor matrix to use for co-clustering */
    /* will be input in R as NULL object if user does not invoke
     * predictor-depedent cluster assignments.
     */
    mat ksi;
    if( o_ksi == R_NilValue )
    {
      ksi = datum::nan;
    }else{
      ksi = as<mat>(o_ksi);
    }
    /* Normalizing offset, N x T, E.
     * 
     */
    mat E   = as<mat>(o_E);
    E.elem(find(E == 0)).ones(); /* set 0 values in county-month population equal to 1 */
    /* memo: on R side: Q_s12  = as.matrix(precmat.season(T, season=12)) */
    /*                  Q_s12  <- as(Q_s12, "dgCMatrix") */
    List Cr(o_C); /* K, T x T iGMRF normalized adjacency matrices - defined from dgCMatrix in R */
    NumericMatrix Dr(o_D); /* K x T matrix, where each row holds the diagonal elements of Q_k */
    IntegerVector o_r(o_order); /* covariance kernel indicators for L covariance terms */
    NumericVector ipr_r(o_ipr);
    vec ipr   = as<vec>(ipr_r); /* Inclusion probabilities to produce Horvitz-Thompson plug-in */
    vec ipr_dummy(ipr.n_elem); ipr_dummy.ones(); /* Use to make no direct adjustments to weights in clustering algorithm */
    int nu    = as<int>(o_nu); /* df parameter for H&W prior on precision locs, Lambda_star */
    int w_star= as<int>(o_w_star); /* tuning parameter for Algorithm 8 cluster assignment */
    int Rep   = as<int>(o_Rep); /* Number of draws per sampling iteration for Psi */
    int niter = as<int>(niterInt); /* sampler iterations */
    int nburn = as<int>(nburnInt); /* sampler burnin */
    int nthin = as<int>(nthinInt); /* sampler thinning */
    int kappa_fast = as<int>(o_kappa_fast); /* use kappa generated from full conditionals vs joint */
    int nkeep = (niter -  nburn)/nthin; 
    unsigned int M     = as<int>(Minit); /* initial number of clusters for sampling */
    double ac = as<double>(shapealph); /* DP conc hyperparameter */
    double bc = as<double>(ratebeta);  /* DP conc hyperparameter */
    double a = as<double>(o_a); /* GP parameters (theta_star) hyperparameters */
    double b = as<double>(o_b); /* GP parameters (theta_star) hyperparameters */
    double a_tau = as<double>(o_a_tau); /* GP parameters (tau) hyperparameters */
    double b_tau = as<double>(o_b_tau); /* GP parameters (tau) hyperparameters */
    int progress = as<int>(o_progress); /* indicator for whether to display a progress bar */
    double jitter = as<double>(o_jitter); /* jitter on posterior rate of kappa_star to stabilize */
    bool stable_launch = as<bool>(o_stable_launch); /* how to initialize log-mean, Psi */
    
    // Extract key row and column dimensions we will need to sample parameters
    int N      = Yr.nrow(), T = Yr.ncol(); /* data dimensions */
    int K      = Cr.size(); /* number of iGMRF terms */
    int i, j, k, p; /* initialize loop variables */
    unsigned int m;
    
    /* data */
    mat Y(Yr.begin(), N, T, false);
    // arma objects
    mat D(Dr.begin(),K,T,false); /* dense matrix where each row holds the T diagonal for Q_k */
    // field<sp_mat> C(K,1); /* list of normalized adjacency matrices */
    field<mat> C(K,1); /* list of normalized adjacency matrices */
    cube Q(T,T,K); /* iGMRF precision matrices for K terms, computed from C */
    /* uvec, o, to capture order for each of the K T x T precision matrices, Q_k */
    // uvec o(o_r.begin(),K,false); /* term of length K identifies order of cov formulations */
    uvec o = as<uvec>(o_r);
    //field<vec> eigraw(K,1);
    //vec eigval_k;
    mat D_k(T,T); /* diagonal matrix with values the diagonals of Q.slice(k) */
    for(k = 0; k < K; k++)
    {
          //C(k,0)          = as<sp_mat>(Cr[k]); 
          C(k,0)          = as<mat>(Cr[k]);
          D_k             = eye(T,T);
          D_k.diag()      = D.row(k);
          Q.slice(k)      = D_k * (eye(T,T) - C(k,0));
          //Q.slice(k)      = D_k - (D_k*C(k,0));
          /* compute raw eigenvalues of Q (not uplifed by kappa_star) */
          //eig_sym(eigval_k, symmatl(Q.slice(k)));
          /* use non-zero eigenvalues to compute determinant */
          //eigraw(k,0)    = sort( eigval_k, "descend" );
          //eigraw(k,0)    = eigraw(k,0).subvec(0,(T-o(k)-1));/* Q is rank-deficient */
    } /* end loop k over iGMRF terms */
 
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
    s                         = join_cols(s,srest); s -= 1; /* use s as position selector */
    s 	                      = shuffle(s,0); 
    ucolvec num(M); uvec pos;
    /* Horvitz-Thompson scaled up num vector with inverse probability weighting */
    colvec num_ht(M); num_ht.zeros();
    for(m = 0; m < M; m++)
    {
          pos	                = find(s == m);
	        num(m)              = pos.n_elem;
          //num_ht(m)           = sum(1/ipr(pos));
          num_ht(m)           = sum(1/ipr_dummy(pos));
    }   
    /* cluster locations, kappa_star */
    NumericVector _kappa_vec            = rgamma( (K*M), a, double(1/b) );
    //NumericVector _kappa_vec            = rgamma( (K*M), 1, (1/1) );
    vec kappa_vec                       = as<vec>(_kappa_vec);
    mat kappa_star                      = reshape( kappa_vec, K, M );
    /* since M changes iteration-to-iteration, use kappa to capture parms */
    mat kappa(K,N); kappa.zeros();
    /* estimated functions, B(N,T,K) */
    mat gamma(N,T); gamma.zeros(); /* gamma = B_1 + B_2 + ... + B_K */
    colvec bki(T); bki.zeros();/* T draws for term k, row i for function B_k */
    cube B(N,T,K); /* each term k holds N functions, each (bki) of length T */
    for( k = 0; k < K; k++ )
    {
        for( i = 0; i < N; i++ )
        {
            for( j = 0; j < T; j++ )
            {
              bki(j)   = R::rnorm(0.0,(1/sqrt(kappa_star(k,s(i))*D(k,j))));
            }
            (B.slice(k)).row(i)  =  bki.t();
        } /* end loop i over rows of B_k */
          /* generate N x T, gamma = sum_k B_k */
          gamma += B.slice(k);
    } /* end loop k over slices of N x T x K, B */
    
    /* create K x N matrix of quadratic products, */
    /* B1(k,i)       = 0.5*dot( D.row(k), pow((bki-bbar_ki),2) ), where bki is T x 1 */
    /* used both compose q0,i for DP and to sample kappa_star(k,m) */
    mat B1(K,N); B1.zeros();
    
    /* if !is.null(ksi), include include predictors, ksi, 
     * in prior for global (s[i])
     * cluster assignments.
     */
    int R = 2;
    if( is_finite(ksi) ) /* grab number of predictors, R. */
    {
      R      = ksi.n_cols; /* Number of spatial locations for each area and harmonic */
    }
    mat u_star(R,M), P_bar(R,R), as_star(R,M); 
    colvec u_bar(R), as_star_m(R);
    cube Lambda_star(R,R,M); 
    
    
    if( is_finite(ksi) )
    {
      mat inv_Vm;
      u_bar.zeros();
      P_bar = eye<mat>(R,R);
      colvec u_star_m;
      for( m = 0; m < M; m++ )
      {
        /* construct precision cluster locations for mixture prior on Ksi */
        as_star.col(m).ones();
        inv_Vm  = 2*nu*diagmat(as_star.col(m));
        wishrnd( Lambda_star.slice(m), inv_sympd(symmatl(inv_Vm)), (R+nu-1) );
        /* construct mean cluster locations for mixture prior on Ksi */
        rmvnsample( P_bar, u_bar, u_star_m );
        u_star.col(m) = u_star_m;
        
      } /* end loop k over K truncated number of global clusters */
        
    } /* end condition on whether employ spatial weights for local cluster assignment */
        
    double conc     = 1; /* DP concentration parameter */
    /* global noise precision parameter */
    double tau_e    = rgamma(1, a_tau, (1/b_tau))[0];
    
    /* define a replicated data matrix to capture sampled missing values */
    mat Y_rep  = Y; /* if no missing values, this will never be updated, which is fine */
    
    /* draw initial values for n x T latent log-mean, Psi */
    mat Psi(N,T); Psi.zeros();
    mat psi_bar     = gamma;
    if( !stable_launch )  
    {
      for( i = 0; i < N; i++ )
      {
        for( j = 0; j < T; j++ )
        {
          Psi(i,j) = R::rnorm(psi_bar(i,j),1/sqrt(tau_e));
        } /* end loop j over time points */
      } /* end loop i over areas */
    }else{ /* stable launch */
      /* initial value for numerical stability */
      /* fill in missing values in Y_rep */
      uvec na_vals                = find(Y == -9);
      int n_na                    = na_vals.n_elem;
      Y_rep.elem( na_vals )       = randu<vec>(n_na) + E.elem( na_vals ); /* always > E */
      Psi                         = log( Y_rep + 1 ) - log( E + 1 );
    } /* end condition on whether to launch Psi in a more stable manner */
    /* rather than by drawing it from the model under initial hyperparameter values */
    
    /* mean, Y_bar, of Y */
    mat Y_mean    = E % exp( Psi );
    mat Y_rate    = Y_rep / E; /* there are no zeros in E */
    
    
    /* initialize vector of residuals */
    colvec resid(N*T), y_vec(N*T), mu_vec(N*T); 
    resid.zeros(); y_vec.zeros(); mu_vec.zeros();
    /* fit assessment - related measures */
    double deviance = 0;  ucolvec ordscore(nkeep); ordscore.zeros();
    mat phat(N,N); phat.zeros(); /* empirical co-clustering probability matrix */
    rowvec devmarg(N); colvec devres(4); devres.zeros();
    rowvec logcpo(N); logcpo.zeros(); double lpml;
    /* capture samples to return */
    int oo = 0, kk = 0;
   
    // Armadillo structures to capture results
    // Will be wrapped into a list element on RETURN
    // DP return objects
    ucolvec numM(nkeep); /* M */
    umat S(nkeep,N); /* cluster indicator  vectors of length N */
    field<ucolvec> Num(nkeep,1); /* number of observations per cluster */
    field<colvec> Num_ht(nkeep,1); /* number of observations per cluster scaled up to finite pop */
    field<ucolvec> bigS; /* best fit clustering - units bucketed by cluster */
    colvec Conc(nkeep); /* DP concentration parameter */
    /* locations */
    mat Kappa(nkeep,(K*N)); /* K parameters for each of N units */
    field<mat> Kappa_star(nkeep,1); /* each will hold a K x M_p matrix used to compute invG_star */
    /* predictor indexed locations for a priori probability of co-clustering */
    field<mat> u_stars(nkeep,1); /* holds R x M matrix of mean locations */
    field<cube> Lambda_stars(nkeep,1); /* holds R x R x M cube of covariance locations */
    field<mat> as_stars(nkeep,1); /* R x M matrix of H&W hyperparms for Lambda_star */
    mat u_bars(nkeep,R); /* R x 1 mean of u_star */
    cube P_bars(R,R,nkeep); /* R x R prior precision of u_star */
    /* use kappa_star for prediction */
    /* sum over K functions, gamma */
    mat bb(nkeep,(N*T)); /* N is fast-moving */
    /* functions, B(N,T,K) */
    field<mat> f(K,1); /* each element holds an nkeep x N*T matrix, with N fast-moving */
    for( k = 0; k < K; k++ )
    {
          f(k,0).set_size(nkeep,(N*T));
    }
    /* capture replicated y values */
    mat Psis(nkeep,(N*T)); /* n is fast, T is slow */
    mat Y_reps(nkeep,(N*T));
    mat Y_means(nkeep,(N*T)); /* (E[i,j]*exp(Psi[i,j])), n is fast, T is slow */
    mat Y_rates(nkeep,(N*T)); /* (Y[i,j]/E[i,j]), n is fast, T is slow */
    /* non-cluster parameters */
    colvec Tau_e(nkeep);
 
    /* fit indicators */
    mat Resid(nkeep,(N*T));
    colvec Deviance(nkeep); Deviance.zeros();  
    mat Devmarg(nkeep,(N*T));     
    
    // POSTERIOR samples 
    for( p = 0; p < niter; p++ )
    {
        if(progress == 1)
        {
             if( (p % 450) == 0 ) Rcout << "Production Interation: " << p << endl;
        }
      
        if( any(vectorise(Y) == -9) )
        {
          miss_ycount(Y_rep, Y, E, Psi);
        }
        /* sample N x T log-mean, Psi, of Poisson likelihood on Y */
        move_Psi_i(Psi, Y_rep, E, gamma, tau_e, Rep);
        /* after filling in missing values in y_rep, update estimated functions in N x T x K, B */
        //move_B(y_rep, B, kappa_star, C, gamma, D, s, tau_e);
        move_B_alt(Psi, B, kappa_star, C, gamma, D, s, tau_e);
        /* M changes iteration-to-iteration, so dimension of theta_star will update */
        /* will adjust size of kappa_star(k,m) as M changes, as well as updating s */
        /* Have to move clusters before kappa because use B1 (quadratic product) in move_kappastar */
        /* use improper joint iGMRF dens */
        //clusterstep_alt(B, kappa_star, B1, o, Q, s, num, M, conc, a, b, ipr_dummy, num_ht); 
        if( !is_finite(ksi) ) /* no dependent PPM */
        {
          clusterstep(B, kappa_star, B1, o, C, D, s, 
                      num, M, conc, a, b, ipr_dummy, num_ht); /* proper iGMRF full conds */
          if(kappa_fast == 1) /* because have already computed B1 in clusterstep */
          {
            move_kappastar(kappa_star, B1, s, o, T, a, b, ipr);
          }else{
            move_kappastar_alt(kappa_star, B, Q, s, o, T, a, b, ipr);
          }
        }else{ /* predictor dependent PPM */
          auxclusterstep_gmrf(B, ksi, kappa_star, o, C, D, 
                         u_star, Lambda_star, as_star, nu,
                         u_bar, P_bar, 
                         s, num, M, w_star, conc, a, b,
                         ipr_dummy, num_ht);
          move_kappastar_alt(kappa_star, B, Q, s, o, T, a, b, ipr);
          move_ustar(u_star, ksi, Lambda_star, u_bar, P_bar, s);
          move_Lambdastar(Lambda_star, as_star, ksi, u_star, s, nu);
          for( m = 0; m < M; m++ )
          {
            //as_star_m = as_star.col(m); /* R x 1 */
            move_as(Lambda_star.slice(m), as_star_m, nu, b);
            as_star.col(m) = as_star_m;
          }
          move_ubar(u_bar, u_star, P_bar);
          move_Pbar(P_bar, u_bar, u_star);
        } /* 
           * end sampling of cluster assignments and locations conditional on whether
           * a non-null input of predictor matrix, ksi, indicates a predictor-dependent
           * product partition model.
           */
            
        concstep(conc, M, N, ac, bc);
        move_taue_jitter(Psi, gamma, tau_e, a_tau, b_tau, jitter, ipr);
        
        Y_mean      = E % exp( Psi ); /* n x T mean of Y, which is our modeled estimate */
        Y_rate      = Y_rep / E ; /* n x T employment rate, Y / E to compare to exp(Psi)*/
        
        if(p >= nburn)
        {
            kk = p - nburn;
            if( kk == ((oo+1)*nthin - 1) )
            {
               /* monitor average chain acceptance rate */
     
               /* capture K x N parameter matrix for sampling iteration k */
               /* M changes on every iteration so need to store by N units clustered */
               kappa                         = kappa_star.cols( s );
               /* K is the fast-moving index, N is the slow-moving index */
               Kappa.row(oo)                 = vectorise( kappa ).t();
               Kappa_star(oo,0)              = kappa_star;
               Lambda_stars(oo,0)            = Lambda_star; /* R x R x M */
               u_stars(oo,0)                 = u_star; /* R x M */
               as_stars(oo,0)                = as_star; /* R x M */
               u_bars.row(oo)                = u_bar.t(); /* 1 x R */
               P_bars.slice(oo)              = P_bar; /* R x R */
               Conc(oo)                      = conc;
               /* vectorize N x T gamma, by columns, so N is fast-moving */
               bb.row(oo)                    = vectorise( gamma ).t();
               Y_reps.row(oo)                = vectorise( Y_rep ).t();
               Y_means.row(oo)               = vectorise( Y_mean ).t();
               Y_rates.row(oo)               = vectorise( Y_rate ).t();
               /* vectorise n x T matrix, Psi, by column - so n is faster than T */
               Psis.row(oo)                  = vectorise( Psi ).t();
               for( k = 0; k < K; k++ ) /* loop over K iGMRF terms */
               {                               /* 1 x N*T */
                  f(k,0).row(oo)             = vectorise( B.slice(k) ).t(); 
               } 
               Tau_e(oo)                     = tau_e;
               numM(oo)                      = M;
               S.row(oo)                     = s.t();
               Num(oo,0)                     = num;
               Num_ht(oo,0)                  = num_ht;
               /* compute chain acceptance statistics */ 
               /* N is the fast-moving index, T is the slow-moving index */
               y_vec                         = vectorise( Y_rep ); /* y is an N x T mat */
               mu_vec                        = exp( vectorise(Psi) );
               /* residual for a Poisson likelihood */
               resid                         = y_vec * log(y_vec / mu_vec)
                                                    - (y_vec - mu_vec); 
               deviance                      =  2 * sum(resid); /* scalar double deviance */
               dmarg_count(y_vec, mu_vec, devmarg); /* 1 x N*T vector of densities*/;
               Deviance(oo)                  = deviance;
               Devmarg.row(oo)               = devmarg;
               Resid.row(oo)                 = resid.t();
               oo += 1; /* increment sample return counter */
            } /* end conditional statement on returning sample */
        } /* end if p > burnin, record results for return */
    } /* end MCMC loop over p */
    
    // compute y_bar based on averaging over MCMC samples that include missing data.
    rowvec y_bar_vec     = mean( Y_rep ); /* by columns */
    mat Y_bar            = reshape( y_bar_vec, N, T ); /* by column since N is fast-moving (for an N x T matrix) */
    
    // compute FIT statistics
    // compute least squares cluster
    
    lsqcluster(S, ordscore, phat, bigS);
    List optpartition = Rcpp::List::create(Rcpp::Named("ordscore")    = ordscore,
                                           Rcpp::Named("S")           = S,
                                           Rcpp::Named("Num")         = Num,
                                           Rcpp::Named("Num_ht")      = Num_ht,
                                           Rcpp::Named("ipr")         = ipr,
                                           Rcpp::Named("phat")        = phat,
                                           //Rcpp::Named("eigraw")      = eigraw,
                                           /* N x T data matrix need for plotting fit */
                                           Rcpp::Named("y")           = Y,
                                           Rcpp::Named("y_bar")       = Y_bar,  
                                           Rcpp::Named("Y_reps")      = Y_reps, 
                                           Rcpp::Named("Y_rates")     = Y_rates,
                                           Rcpp::Named("Y_means")     = Y_means,
                                           Rcpp::Named("Psis")        = Psis,
                                           Rcpp::Named("Q")           = Q,
                                           Rcpp::Named("q_order")     = o,
                                           Rcpp::Named("u_bars")      = u_bars,
                                           Rcpp::Named("P_bars")      = P_bars
                                           );
    // DIC
    dic3comp(Deviance, Devmarg, devres); /* devres = c(dic,dbar,dhat,pd) */
    cpo(Devmarg, logcpo, lpml);

    // Return results
    return Rcpp::List::create(Rcpp::Named("Deviance")            = Deviance,
                                  Rcpp::Named("Devmarg")         = Devmarg,
                                  Rcpp::Named("devres")          = devres,
                                  Rcpp::Named("logcpo")          = logcpo,
                                  Rcpp::Named("lpml")            = lpml,
               	                  Rcpp::Named("Kappa")           = Kappa,
                                  Rcpp::Named("Kappa_star")      = Kappa_star,
                                  Rcpp::Named("Lambda_stars")    = Lambda_stars,
                                  Rcpp::Named("u_stars")         = u_stars,
                                  Rcpp::Named("as_stars")        = as_stars,
                                  Rcpp::Named("Conc")            = Conc,
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
