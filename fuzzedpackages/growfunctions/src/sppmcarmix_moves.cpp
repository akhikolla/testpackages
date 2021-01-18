#include "growfunctions.h"
using namespace Rcpp;
using namespace arma;
using namespace std;

// update vector of cluster membership indicators, s(i),....,s(N)
SEXP auxclusterstep_gmrf(const cube& B, const mat& ksi, mat& kappa_star, const uvec& o,
                 const field<mat>& C, const mat& D, 
                 mat& u_star, cube& Lambda_star,
                 mat& as_star, double nu,
                 const colvec& u_bar, const mat& P_bar, 
                 ucolvec& s, 
                 ucolvec& num, unsigned int& M, const int& w_star,
                 double& conc, int a, int b,
                 const vec& ipr, colvec& Num)
{
  BEGIN_RCPP
  // sample cluster assignments, s(1), ..., s(N)
  // B = (B_1,...,B_K), where B_k is N x T matrix for iGMRF term k
  // Q = (Q_1,...,Q_K), where Q_k is a T x T de-scaled iGMRF precision matrix
  // C = (C_1,...,C_K), where T x T, C_k = D_k^-1 * Omega_k, 
  // where Omega_k is the T x T adjacency matrix for iGMRF term, k
  // D is a K x T matrix where row k contains T diagonal elements of Q_k
  // n x R, ksi is matrix of predictors; e.g., spatial locations (lat, lon)
  //        to use to set a predictor-dependent prior on cluster assignments, s
  // K x M matrix, kappa_star records locations for each iGMRF term
  // R x M matrix, u_star are mean locations for DP mix prior on R x 1, ksi_i
  // R x R x M cube, Lambda_star are R x R precision matrix locations for mix prior on ksi
  // R x M, as_star is matrix of by-cluster precision parameters of Huang & Wand (2013)
  // o = (o_1,...,o_k) is a vector where each entry denotes the order of term K.
  // e.g. RW(1) -> o = 2, RW(2) -> o = 3, seas(3) -> o = 3
  unsigned int auxSize, startPos, h, ell, m;
  int N = B.n_rows;  /* number of observations */
  int T = B.n_cols; /* number of months */
  int K = C.n_rows; /* number of random effect terms */
  int R = Lambda_star.n_rows; /* number of predictors for prior on cluster assignment */
  int i, j, k;
  mat kappa_aux_h; /* added auxSize = h - M auxiliary locations */
  mat u_aux_h; /* added R x 1 auxSize mean locations for DP mix prior on spatial locations */
  colvec e_m    = P_bar * u_bar; /* R x 1 prior mean of u_star.col(m) */ 
  mat phi_m     = P_bar; /* R x R prior precision of u_star.col(m) */
  cube Lambda_aux_h; /*added auxSize R x R precision locations for DP mix on spatial locs */
  /* where r+nu is the prior df for the Wishart on Lambda_star.slice(m) */
  mat as_aux_h; /* cluster-indexed hyperparameters for H&W prior for Lambda_star */
  mat inv_Vm; /* prior precision matrix mean of the Wishart */
  uvec loc_keep; /* final set of cluster location labels linked to observations */
  unsigned int num_keep = 0; /* number of cluster location labels linked to observations */
  //colvec Num; /* Horvitz-Thompson scaled up num vector with inverse probability weighting */
  cube Lambda_keep; /* final set of by-cluster precision matrices linked to observations */
  double llike_iell = 0; /* log_like for (B,ksi).row(i) under assign to cluster ell */
  colvec bki(T), bbar_ki(T); /* T x 1, bbar_ki = D_k^-1*Omega_k*b_ki = C(k,0)*b_ki */
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
      kappa_star.insert_cols(M,1);
      u_star.insert_cols(M,1);
      as_star.insert_cols(M,1);
      Lambda_star.insert_slices(M,1);
      /* carrying along logfm_old vector to maintain current cluster */
      /* count and location. Recall, logfm_old is used for Metropolis */
      /* steps to update rho_star */
      kappa_star.col(M)		      = kappa_star.col(s(i));
      u_star.col(M)             = u_star.col(s(i));
      as_star.col(M)            = as_star.col(s(i));
      Lambda_star.slice(M)      = Lambda_star.slice(s(i));
      num.insert_rows(M,1);
      num(M)			     = 0; /* insert new location ahead of observations */
      /* scale up num to population totals, Num, based on H-T inverse probability estimator */
      Num.insert_rows(M,1);
      Num(M)			     = 0; /* insert new location ahead of observations */
      kappa_star.shed_col(s(i));
      u_star.shed_col(s(i));
      as_star.shed_col(s(i));
      Lambda_star.shed_slice(s(i));
      num.shed_row(s(i));
      Num.shed_row(s(i));
      s( find(s > s(i)) ) 	-= 1;
      s(i)			            = M-1;
      startPos		          = M-1;
    }/* end setting starting index for w_star new cluster locations */
    /* set total number of existing and auxiliary cluster locations, h */
    h 		= M + auxSize;
    colvec weights_i(h); weights_i.zeros();
    
    /* scale up num to population totals, Num, based on H-T inverse probability estimator */
    //pop_Num(s, ipr, Num); /* computed after set s, num and h */
    
    /* sample w_star new location values from prior (F_0) */
    /* since ahead of observations */
    /* compute weights used to draw s(i) */
    if( h > M ) /* num(s(i)) > 1 | (num(s(i)) == 1 & w_star > 1) */
    {
      /* create auxSize new locations to insert to kappa_star */
      NumericVector _kappa_vec            = rgamma( (K*auxSize), a, (1/b) );
      vec kappa_vec                       = as<vec>(_kappa_vec);
      kappa_aux_h                         = reshape( kappa_vec, K, auxSize );
      kappa_star.insert_cols(M,kappa_aux_h); /* kappa_star now has h columns */
      /* create auxSize new locations to insert to Lambda_star */
      /* first augment the cluster-indexed hyperparameters */
      NumericVector _as_vec               = rgamma( (R*auxSize), 0.5, (1/b) );
      vec as_vec                          = as<vec>(_as_vec);
      as_aux_h                            = reshape( as_vec, R, auxSize ); /* R x auxSize */
      as_star.insert_cols(M,as_aux_h); /* as_star now has h columns */
      /* next augment the cube of cluster-indexed precision matrices, Lambda_star */
      Lambda_aux_h.set_size(R,R,auxSize);
      u_aux_h.set_size(R,auxSize);
      colvec u_aux_m;
      for( m = 0; m < auxSize; m++ )
      {
        inv_Vm            = 2*nu*diagmat( as_star.col(m) );
        wishrnd(Lambda_aux_h.slice(m), inv_sympd(symmatl(inv_Vm)) , (R+nu-1));
        rmvnbasic(phi_m,e_m,u_aux_m);
        u_aux_h.col(m) = u_aux_m;
      } /* end loop m over auxiliary clusters under algorithm 8 of Neal */
      Lambda_star.insert_slices(M,Lambda_aux_h);
      u_star.insert_cols(M,u_aux_h);
      /* create auxSize new locations to insert to u_star */
      
      /* memo: no observations, yet */
      /* adding entries for aux vars and set to 0 */
      /* by default, new rows/cols/slices set to 0 */
      num.insert_rows(M,auxSize,true);
      Num.insert_rows(M,auxSize,true); /* population uplifted */
      /* inserts 0 for M -> h-l positions */
      weights_i			     = Num / (double(N)- (1/ipr(i)) +conc); 
      weights_i.subvec(M,(h-1))	= 
            ( (conc/w_star) / (double(N) - (1/ipr(i)) + conc) )*ones(auxSize);
    }else{ /* num(s(i)) == 1 & w_star == 1, so new location already sampled at M-1 */
      weights_i			     = Num / (double(N)- (1/ipr(i)) +conc);
      weights_i(M-1)		 = (conc/w_star) / (double(N) - (1/ipr(i)) + conc);
    } /* end sampling of new locations and computing weights */
      
    /* draw value for cluster index, s(i) */
    /* note that rdrawone() function has minimum value of 0 */
    for( ell = 0; ell < h; ell++)
    {
        s(i)            = ell; /* temporary assignment */
        llike_iell      = 0; /* hold log densities for K computations */
        for(k = 0; k < K; k++)
        {
          bki             = B.slice(k).row(i).t();
          bbar_ki         = C(k,0) * bki; /* T x 1 */
          for( j = 0; j < T; j++ )
          {
            /* effectively making assignment, s(i) = l */
            llike_iell   += (
                              trunc_log(R::dnorm(bki(j),bbar_ki(j),
                                 double(1/sqrt(kappa_star(k,s(i))*D(k,j))),false))
                            );
          } /* end loop j over time index */
        } /* end loop k over iGMRF terms */
        /* Add likelihood for predictors used to index cluster assignment */
        llike_iell      += logmvndens(ksi.row(i).t(), 
                                    u_star.col(s(i)), Lambda_star.slice(s(i)));
        weights_i(ell) *= trunc_exp(llike_iell);
    } /* end loop ell over h cluster assignment evaluations */
      weights_i.elem( find_nonfinite(weights_i) ).zeros();
      weights_i		     /= sum(weights_i);
    // conduct discrete posterior draw for s(i), where min(s(i)) = 0
    s(i) 		           = rdrawone(weights_i, h);
    num(s(i))		       += 1;
    Num(s(i))          += 1/ipr(i);
    
    /* remove clusters with no observations  - must be among the w_star new clusters */
    loc_keep            = find(num != 0);
    num_keep            = loc_keep.n_elem;
    num                 = num( loc_keep );
    Num                 = Num( loc_keep );
    kappa_star		      = kappa_star.cols( loc_keep );
    u_star		          = u_star.cols( loc_keep );
    as_star		          = as_star.cols( loc_keep );
    /* Armadillo doesn't have a non-contiguous operation on slices */
    /* must update slice-by-slice */
    Lambda_keep.set_size(T,T,num_keep);
    for( ell = 0; ell < num_keep; ell++ )
    {
      Lambda_keep.slice(ell)   = Lambda_star.slice( loc_keep(ell) );
    } /* end loop ell over clusters linked to >= 1 observations */
    
    /* replace U_last with U_keep, which holds choleskys linked to observations */
    Lambda_star.set_size(R,R,num_keep);
    Lambda_star         = Lambda_keep;

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
    M = kappa_star.n_cols;
    
  } /* end loop i for cluster assignment */
  END_RCPP
} /* end function auxclusterstep for cluster assignments, s  */

         
SEXP move_ustar(mat& u_star, const mat& ksi,
                const cube& Lambda_star, const colvec& u_bar,
                const mat& P_bar, const ucolvec& s)
{ 
  BEGIN_RCPP
  // sample set of K, r x r precision matrices, Lambda_star, under Huang & Wand (2013)
  // Lambda_star is r x r x K collection of cluster-indexed precision matrices
  // r x K, as_star is matrix of by-cluster precision parameters of Huang & Wand (2013)
  // r x n_ell x n, Ksi, is cube of (lat,long,resolution) spatial locations.
  // r x K, u_star is matrix of cluster-indexed means for each Ksi(i,ell) linked to k 
  // n x ell, s is matrix of global cluster assignments for each (area,harmonic) = (i,ell)
  int K     = Lambda_star.n_slices; /* number of predictors for each time point */
  int r     = Lambda_star.n_rows;

  /* where r+nu is the prior df for the Wishart on Lambda_star.slice(k) */
  colvec e_k(r);
  mat phi_k(r,r);
  uvec pos_k; /* find areas i linked to global cluster k */
  int n_k = 0;
  colvec u_star_k;
  
  int k = 0, i = 0;
  for( k = 0; k < K; k++ )
  {
    /* initialize objects formed from each (i,ell) linked to k */
    e_k   = P_bar * u_bar; /* r x r, prior mean of u_star, P_bar, times r x 1, u_bar */ 
    phi_k = P_bar;
    /* areas i linked to global cluster k */
    pos_k       = find( s == k );
    n_k         = pos_k.n_elem;
    // if n_ik = 0, meaning there are no harmonics in area i linked to global cluster k
    // then quad_k is not incremented.
    for( i = 0; i < n_k; i++ )
    {
      e_k     += Lambda_star.slice(k) * ksi.row(pos_k(i)).t();
    } /* end loop ell over n_ik harmonics in area i linked to glob cluster k */
    /* posterior precision for u_star.col(k) */
    phi_k += n_k * Lambda_star.slice(k);  
    /* note V is a scale matrix */
    rmvnbasic(phi_k,e_k,u_star_k);
    u_star.col(k) = u_star_k;
  } /* end loop k over K clusters */
  END_RCPP
} // end function move_ustar() to sample cluster-indexed means of spatial locations,
// Ksi(i,ell)

          
SEXP move_Lambdastar(cube& Lambda_star, const mat& as_star, const mat& ksi,
                     const mat& u_star, const ucolvec& s, double nu)
{
  BEGIN_RCPP
  // sample set of K, r x r precision matrices, Lambda_star, under Huang & Wand (2013)
  // Lambda_star is r x r x K collection of cluster-indexed precision matrices
  // r x K, as_star is matrix of by-cluster precision parameters of Huang & Wand (2013)
  // r x n_ell x n, Ksi, is cube of (lat,long,resolution) spatial locations.
  // r x K, u_star is matrix of cluster-indexed means for each Ksi(i,ell) linked to k 
  // n x ell, s is matrix of global cluster assignments for each (area,harmonic) = (i,ell)
  int K     = Lambda_star.n_slices; /* number of predictors for each time point */
  int r     = Lambda_star.n_rows;
  
  /* where r+nu is the prior df for the Wishart on Lambda_star.slice(k) */
  mat inv_Vk; /* prior precision matrix mean of the Wishart */
  uvec pos_k; /* find observations i linked to global cluster k */
  int n_k = 0;
  colvec ksitilde_ik(r); /* ksi.row(s(i)).t() - u_star(k) */
  /* form rank 1 matrices for each i */
  mat quad_ik(r,r), quad_k(r,r);
  quad_k.zeros();
  
  int k = 0, i = 0;
  for( k = 0; k < K; k++ )
  {
    /* initialize objects formed from each (i,ell) linked to k */
    quad_k.zeros();
    inv_Vk          = 2*nu*diagmat( as_star.col(k) );
    /* harmonics, ell, in area i linked to global cluster k */
    pos_k       = find( s == k );
    n_k         = pos_k.n_elem;
    for( i = 0; i < n_k; i++ )
    {
      ksitilde_ik    = ksi.row(pos_k(i)).t() - u_star.col(k);
      quad_ik        = ksitilde_ik * ksitilde_ik.t(); /* r x r */
      quad_k         += quad_ik;
    } /* end loop ell over n_ik harmonics in area i linked to glob cluster k */
      
    inv_Vk  += quad_k;
    /* note V is a scale matrix */
    wishrnd(Lambda_star.slice(k), inv_sympd(symmatl(inv_Vk)) , (n_k + (r+nu-1))); 
  } /* end loop k over K clusters */
  END_RCPP
} /* end function move_Lambdastar() to sample rxr cluster-indexed precisions for Ksi */
      


SEXP move_ubar(colvec& u_bar, const mat& u_star,
               const mat& P_bar)
{ 
  BEGIN_RCPP
  int K         = u_star.n_cols;
  int r         = u_star.n_rows;
  colvec e_bar(r); e_bar.zeros();
  mat phi_bar   = K * P_bar + eye<mat>(r,r);
  int k         = 0;
  for( k = 0; k < K; k++ )
  {
    e_bar += P_bar*u_star.col(k);
  } /* end loop k over populated global clusters */
  rmvnbasic(phi_bar,e_bar,u_bar);
  
  END_RCPP
} /* end function move_ubar() to sample R x 1 mean of u_stars */
      
      
SEXP move_Pbar(mat& P_bar, const colvec& u_bar, const mat& u_star)
{ 
  BEGIN_RCPP
  int K         = u_star.n_cols;
  int r         = u_star.n_rows;
  mat inv_Vbar  = eye<mat>(r,r);
  int k         = 0;
  for( k = 0; k < K; k++ )
  {
    inv_Vbar +=  u_star.col(k) * u_star.col(k).t();
  } /* end loop k over populated global clusters */
  /* note V is a scale matrix */
  wishrnd(P_bar, inv_sympd(symmatl(inv_Vbar)) , (K + (r+1)));
  
  END_RCPP
} /* end function move_Pbar() to sample R x R prior precision of u_stars */
  
       
SEXP move_as(const mat& P_8, vec& as, double nu, double b)
{
  BEGIN_RCPP
  // Sample P_8, Q x Q precision matrix drawn under Huang & Wand (2013) prior
  int Q           = P_8.n_cols;
  double a_q      = 0.5*(nu + Q);
  
  /* memo: weights of length n, weights(est) of length, n_c */
  /* sum of weights over all n_c est-time cases */
  double b_1q = 0;
  int q = 0;
  for( q = 0; q < Q; q++ )
  {
    b_1q      = nu * P_8(q,q) + b;
    as(q)     = R::rgamma(a_q, (1/b_1q));
  } /* end loop q over elements in taus_q */
  
  END_RCPP
} /* end function move_as() to sample by-variable over-dispersion parameter */
  
/* function to compute likelihood for T x 1 latent, log-mean of Poisson likelihood */
double loglike_psi_i(const colvec& psi_i, const mat& Y, const mat& E,
                     int i)
{
  double y_term     = dot( Y.row(i), (log(E.row(i)) + psi_i.t()) ) - 
    dot( E.row(i), exp(psi_i) );
  //double sig_term = 0.5 * tau_y * sum( pow( (logSig2.row(i)- (log(E.row(i)) + psi_i.t())), 2 ) );
  double logl_psi_i = y_term;// - sig_term;
  return logl_psi_i;
} /* end function loglike_psi_i() to compute like of T x 1, Psi.row(i) */

/* function to move T x 1, Psi.row(i), under ESS */
SEXP move_Psi_i(mat& Psi, const mat& Y, const mat& E,  
                const mat& gamma, double tau_e, int R)
{
  BEGIN_RCPP
  // n x T matrix of log-mean, Psi
  // n x T matrix E is offset (e.g., population for each area i in month j)
  // n x P x T, X, contains n x P matrix of predictors over T time points
  // H^* = (eta^*_1,...,eta^*_K), where eta^*_k is T x 1 vector of global cluster locations
  // n x n_ell spatial basis matrix, M
  // n x n_ell, s, holds global cluster assignments, k in 1,...K for each area i and harmonic ell
  // B is P x T matrix of time-indexed regression coefficients
  int n           = Y.n_rows;
  int T           = Y.n_cols;
  int i = 0, j = 0, r = 0;
  colvec v_i(T);
  colvec psi_i(T), psibar_i(T);
  colvec psi_iprime(T); /* new ESS proposal for T x 1, psi_i */
  double loglike_start_i, loglike_cand_i, u, logy;
  double theta, theta_min, theta_max;
  
  
  for( r = 0; r < R; r++ ) /* R re-samples takes advantage of easy generation from prior */
  {
    // int iteration = 0;
    for( i = 0; i < n; i++ )
    {
      /* initialize variable to be sampled */
      psi_i             = (Psi.row(i)).t();
      psibar_i          = gamma.row(i).t();
      
      /* sample T x 1, v_i, for each j in 1,...,T from independent normal prior */
      for( j = 0; j < T; j++ )
      {
        v_i(j) = R::rnorm(gamma(i,j),1/sqrt(tau_e));
      } /* end loop j over time points */
      
      /* compute likelihood threshold for slice */
      loglike_start_i   = loglike_psi_i(psi_i, Y, E, i);//, tau_y, logSig2);
      // Rcout << "loglike_start_c = " << loglike_start_c << endl;
      u                 = runif(1, 0, 1)[0];
      logy              = loglike_start_i + log(u); /* slice threshold */
      /* draw an initial proposal for mixture element, theta */
      theta             = runif(1, 0, 2*M_PI)[0]; /* in radians */
      theta_min         = theta - 2*M_PI;
      theta_max         = theta;
      /* generate a proposal for psi_i */ 
      psi_iprime         = (psi_i - psibar_i) * cos(theta) 
        + (v_i - psibar_i) * sin(theta) 
        + psibar_i;
        /* evaluate L(Psi_prime) and accept-break if L(Psi_prime) >= logy */
        /* else, SHRINK interval for drawing theta and re-compute Psi_prime */
        for(;;) /* infinite loop that contains a break statement */
        {  
          // iteration += 1;
          // Rcout << "iteration: " << iteration << endl;
          loglike_cand_i     = loglike_psi_i(psi_iprime, Y, E, i);//, tau_y, logSig2); 
          // Rcout << "loglike_cand_c = " << loglike_cand_c << endl;
          /* condition on whether new sampled value is outside of slice */
          if( loglike_cand_i > logy ){break;}
          
          /* shrink interval */
          if(theta < 0)
          {
            theta_min      = theta;
          }else{ /* theta > 0  */
          theta_max      = theta;
          } /* end condition to shrink upper or lower interval */
          
          theta           = runif(1,theta_min,theta_max)[0];
          // Rcout << "theta_min = " << theta_min << endl;
          // Rcout << "theta_max = " << theta_max << endl;
          psi_iprime      = (psi_i - psibar_i) * cos(theta) 
            + (v_i - psibar_i) * sin(theta) 
            + psibar_i;
            //             if( i == 6 ){
            // Rcout << "psi_c = " << psi_c << endl;
            // Rcout << "v_c = " << v_c << endl;
            // Rcout << "psibar_c = " << psibar_c << endl;
            // Rcout << "theta = " << theta << endl;
            //             }
        } /* end infinite loop on shrinking slice with break when Pin_new >= logy */
          Psi.row(i)      = psi_iprime.t();
    } /* end loop over T x 1 observations for sampling from prior for Psi.row(i) */
          //Rcout << "Psi.row(6) = " << Psi.row(6) << endl;     
  } /* end duplicate repeated sampling loops for Pis within an MCMC sampling iteration */ 
          
  END_RCPP
} /* end function move_Psi_i() to sample n x T, latent log-mean, Psi, by row */
          
SEXP miss_ycount(mat& Y_rep, const mat& Y, const mat& E, const mat& Psi)
{
  BEGIN_RCPP
  // sample missing elements of n x Q, T, from its posterior
  // predictive distribution, given Psi.  
  // if( !is_finite(Psi) )
  // {
  //   uvec pos  = find_nonfinite( Psi );
  //   int num   = pos.n_elem;
  //   Rcout << "num_non_finite =  " << num << endl;
  // }
  int n                = Y_rep.n_rows;
  
  int i = 0, j = 0, n_posi = 0;
  uvec pos_i; // find non-finite (NA) values in row i of Y
  for( i = 0; i < n; i++ )
  {
    /* find missing elements in row(i) of y */
    /* Y_rep are the data, incl missing vals */
    pos_i               = find( Y.row(i) == -9 ); 
    n_posi              = pos_i.n_elem; /* number of variables with missing vals in case i */
    for( j = 0; j < n_posi; j++ )
    {
      /* y_rep is the completed data matrix, with missing vals of y filled in */
      /* note that, conditioned on model parameters, each observation */
      /* in row c is independent */
      Y_rep(i,pos_i(j))  = R::rpois( E(i,pos_i(j)) * exp(Psi(i,pos_i(j))) );
    }
  }
  
  END_RCPP
} /* end function miss_ystep() to sample missing data values */  
      
SEXP dmarg_count(const colvec& y_vec, const colvec& mu_vec, rowvec& devmarg)
{
  BEGIN_RCPP
  // devmarg is an nc x 1 vector
  int N = y_vec.n_elem;
  double y_i, mean_i;
  int i; /* loop variable */
  for( i = 0; i < N; i++ )
  {
    y_i           = y_vec(i);
    mean_i        = mu_vec(i);
    devmarg(i)    = R::dpois( y_i, mean_i, 0 ); /* nc x 1, f(y|mu)*/
  } /* end loop i over N elements in y_vec */
  END_RCPP
}
