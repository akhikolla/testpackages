/**
 * GUTS: Fast Calculation of the Likelihood of a Stochastic Survival Model.
 * Function guts_engine(GUTS-Object, Parameters).
 * soeren.vogel@posteo.ch, carlo.albert@eawag.ch, alexander singer@rifcon.de, oliver.jakoby@rifcon.de, dirk.nickisch@rifcon.de
 * License GPL-2
 * 2017-10-09 (last changes 2019-01-29)
 */

#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
void guts_engine( Rcpp::List gobj, Rcpp::NumericVector par, Rcpp::Nullable<Rcpp::NumericVector > z_dist = R_NilValue) { 
  
  /*
   * Check if this is an object of class GUTS.
   */
  if ( !gobj.inherits("GUTS") ) {
    Rcpp::stop( "No GUTS object. Use `guts_setup()` to create or modify objects." );
  }
  
  /*
   * Get object elements that are necessary until here.
   */
  NumericVector       par_pos  =  gobj.attr("par_pos");
  NumericVector       wpar     =  gobj.attr("wpar");
  NumericVector       opar     =  gobj["par"];
  NumericVector       oS       =  gobj["S"];
  NumericVector       oD       =  gobj["D"];
  std::vector<double> par_NA(opar.size(), NA_REAL);
  std::vector<double> S_NA(oS.size(), NA_REAL);
  std::vector<double> D_NA(oD.size(), NA_REAL);
  
  /*
   * Check parameters.
   * Error: Stop and do nothing. Warning: Assign NAs.
   * Assign working parameters.
   * Assign parameters to object.
   */
  if ( par.size() != par_pos.size() ) {
    // Error.
    Rcpp::stop( "Vector of parameters has wrong length." );
  } else if ( *std::min_element( par.begin(), par.end() ) < 0.0 ) {
    // Warning.
    Rcpp::warning( "Parameters must be non-negative values -- ignored." );
    gobj["par"] = par_NA;
    gobj["S"]   = S_NA;
    gobj["D"]   = D_NA;
    gobj["LL"]  =  -std::numeric_limits<double>::infinity();
    return;
  } else {
    // assign
    for ( int i = 0; i < par_pos.size(); ++i ) {
      wpar[(par_pos[i] - 1)] = par[i]; // R positions to C++ position.
    }
    // adjust wpar at 3.
    if ( std::isinf(wpar[2]) ) {
      wpar[2] = std::numeric_limits<double>::max();
    }
    // assign if not returned earlier
    gobj["par"] = par;
    gobj.attr("wpar") = wpar;
  }
  
  /*
   * Experiment type.
   * Vectors.
   * S.
   * LL.
   */
  std::vector<double> C  = gobj["C"];
  std::vector<double> Ct = gobj["Ct"];
  std::vector<int>    y  = gobj["y"];
  std::vector<double> yt = gobj["yt"];
  unsigned N = gobj["N"];        
  int M = gobj["M"];         
  int experiment = gobj.attr("experiment");
  std::vector<double> S(yt.size(), 0.0);
  std::vector<double> D( M, 0.0 );
  double LL = 0.0;
  
  std::vector<double> z(N);
  std::vector<double> zw(N);
  if ( experiment < 20 && wpar[4] != 0 ) {
    // lognormal model
    if ( wpar[3] == 0.0) {
      Rcpp::warning( "mn = 0 and sd != 0 -- incomplete lognormal model ignored." );
      gobj["D"]   = D_NA;
      gobj["S"]   = S_NA;
      gobj["LL"]  =  -std::numeric_limits<double>::infinity();
      return;
    }
    const double R = 4.0;
    double sigma2   =  log(   1.0  +  pow( (wpar[4] / wpar[3]), 2.0 )   );
    double mu       =  log(wpar[3])  -  (0.5 * sigma2);
    double sigmaD   =  sqrt(sigma2) * R;
    double ztmp;
    
    if (sigmaD + mu > 700) {
      Rcpp::warning( "Approximating lognormal distribution: infinite variates. Please check parameter values." );
    }
    
    for ( unsigned i = 0; i < N; ++i ) {  
      ztmp = (2.0 * i - N + 1.0) / (N-1.0);
      z.at(i)  =  exp( ztmp * sigmaD + mu );   
      zw.at(i) =  -0.5 * ztmp * ztmp * R * R;
    }
  } else if ( experiment >20 && experiment < 30 ) {
    // Delta model (mode Delta IT is possible but undocumented)
    z.assign(N, wpar[3]);
    zw.assign(N, 0.0);
  } else if ( experiment > 30 && experiment < 40 && wpar[4] != 0 ) {
    const double R      =  50.0; // sampling range importance sampling
    
    // parameters are given as alpha = scale and beta = shape
    // transform parameters to mu and s
    double mu       =  log(wpar[3]);
    double s        =  1 / wpar[4];
    double ztmp;
    
    
    /* if scale (wpar3]) <= 0 or shape (wpar[4]) <= 0: 
     *  the loglogistic distribution is undefined.
     *  These cases are excluded.
     */
    if (wpar[3] <= 0) {
      Rcpp::warning( "Loglogistic distribution undefined for scale parameter <= 0. \nPlease check parameter values." );
      gobj["D"]   = D_NA;
      gobj["S"]   = S_NA;
      //gobj["LL"]  =  -std::numeric_limits<double>::infinity();
      gobj["SPPE"]  =  std::numeric_limits<double>::infinity();
      gobj["squares"]  =  std::numeric_limits<double>::infinity();
      return;
    }
    if (wpar[4] <= 0) {
      Rcpp::warning( "Loglogistic distribution undefined for shape parameter <= 0. \nPlease check parameter values." );
      gobj["D"]   = D_NA;
      gobj["S"]   = S_NA;
      gobj["LL"]  =  -std::numeric_limits<double>::infinity();
      gobj["squares"]  =  std::numeric_limits<double>::infinity();
      return;
    } else {
      /* if shape (wpar4]) <=1: the loglogistic mode = 0 and mean undefined
       * To avoid loglogistic distributions that peak at 0, wpar[4] <= 1 throws a warning
       * Excluding this distribution shape still allows approximation of a concentration threshold of 0,
       * by setting scale \approx 0
       */
      if (wpar[4] <= 1) {
        Rcpp::warning( "Approximating loglogistic distribution: \nShape parameter should be above 1 to avoid an unrealistic concentration threshold distribution that peaks at 0. A concentration threshold close to 0 is better described by a scale parameter that approximates 0. \nNummeric approximation might be wrong. Please check parameter values." );
      }
    }
    // if s * R + mu is above 700, z(N-1) -> infty; returning nan for S and LL   
    if (s * R + mu > 700) {
      Rcpp::warning( "Approximating loglogistic distribution: infinite variates. \nPlease check parameter values." );
    }
    
    /* loglogistic weights formula
     * zw.at(i) =  log(R / 2.0)  - 2.0 *  log( cosh( (log(z.at(i)) - mu) / 2.0 / s ) );
     */
    for ( unsigned i = 0; i < N; ++i ) {
      ztmp = (2.0 * i - N + 1.0) / (N-1.0);
      z.at(i)  =  exp( ztmp * s * R + mu );   
      zw.at(i) =  - 2.0 *  log( cosh( ztmp * R / 2.0 ) );
    }
  } else if ( experiment > 40 && experiment < 50 && z_dist.isNotNull() ) {
    // generic external distribution
    z = as<std::vector<double > >(NumericVector(z_dist.get()));
    std::sort(z.begin(), z.end());
    if (N != z.size()) {
      Rcpp::warning("Sample length is reset to the length of the external distribution.");
    }
    N = z.size();
    zw.assign(N, 0.0);
    gobj["N"] = N;
  } else {
    Rcpp::warning("Error in GUTS model specification. \n -- model ignored");
    gobj["D"]   = D_NA;
    gobj["S"]   = S_NA;
    gobj["LL"]  =  -std::numeric_limits<double>::infinity();
    return;
  }
  
  
  // Model evaluation
  /*
   * Temporary variables.
   */
  std::vector<double> ee( N, 0.0 );
  std::vector<int>    ff( N, 0 );
  double Scale = 1.0;          
  double tmp;
  double summand3;
  double E;
  int    F;
  
  /*
  * tau, dtau, tauit, wpar[2] * dtau.
  */
  double tau     = 0.0;
  double dtau    = (yt.back() - yt.front()) / M;
  int    tauit   = 0;
  double kkXdtau = wpar[2] * dtau;
  
  /* 
   * Iterators, positions.
   */
  int dpos    = 0;
  int ii      = 0;
  unsigned zpos    = 0; 
  int k       = 0;
  
  /*
   * Diffs.
   */
  std::vector<double>   diffS(yt.size());
  std::vector<int>      diffy(yt.size());
  std::vector<double>   diffC(C.size());
  std::vector<double>  diffCt(Ct.size());
  std::vector<double> diffCCt(Ct.size());
  
  for ( unsigned int i = 1; i < C.size(); ++i ) {
    diffC.at(i-1) =  C.at(i) -  C.at(i-1);
    diffCt.at(i-1) = Ct.at(i) - Ct.at(i-1);
    diffCCt.at(i-1) = diffC.at(i-1) / diffCt.at(i-1);
  }
  
  /*
   * Loop over yt.
   */
  for ( unsigned int ytpos = 0; ytpos < yt.size(); ++ytpos ) {
    
    while ( tau < yt.at(ytpos) && dpos < M ) {
      
      tmp = exp( -wpar[1] * (tau - Ct.at(k)) );
      if ( wpar[1] > 0.0 ) {
        summand3 = (tau - Ct.at(k) - (1.0-tmp)/wpar[1])  *  diffCCt.at(k);
      } else {
        summand3 = 0.0;
      }
      D.at(dpos) =  tmp * (D.at(ii) - C.at(k))  +  C.at(k)  +  summand3;
      
      /*
       * Find zpos from D.at(i).
       * Update ee and ff.
       */
      while ( D.at(dpos) < z.at(zpos) && zpos > 0 ) {
        --zpos;
      }
      while ( D.at(dpos) > z.at(zpos) && zpos < (N-1) ) {
        ++zpos;
      }
      if ( D.at(dpos) > z.at(N-1) ) {
        ee.at(N-1) += D.at(dpos);
        ff.at(N-1)++;
      } else if ( D.at(dpos) > z.at(0) ) {
        ee.at(zpos-1) += D.at(dpos);
        ff.at(zpos-1)++;
      }
      /*
       * Increment or decrement dpos, tau, k.
       */
      ++dpos;
      tau = dtau * (++tauit);
      if ( tau > Ct[k+1] ) {
        ++k;
        ii = dpos-1;
      }
      
    } // End of while ( tau < yt[ytpos] && dpos < M )
    
    /*
     * Write S at ytpos.
     */
    E = 0.0;
    F = 0;
    
    for ( int u=N-1; u >= 0; --u ) {
      E += ee.at(u);
      F += ff.at(u);
      S.at(ytpos) += exp(   (kkXdtau * (z.at(u) * F - E))  +  zw.at(u)   );
    }
    S.at(ytpos) *= exp( -wpar[0] * yt[ytpos] ) / N;
    /*
     * Scale S.
     * Calculate diffy, diffS.
     */
    if ( ytpos == 0 ) {
      if ( S.at(0) > 0.0 ) {
        Scale = 1.0 / S.at(0);
      } else {
        Rcpp::warning("Numeric overflow: Survival cannot be calculated for given parameter values. Returning NaN.");
      }
    } else {
      diffS.at(ytpos-1) = S.at(ytpos-1) - Scale * S.at(ytpos);
      diffy.at(ytpos-1) = y.at(ytpos-1) - y.at(ytpos);
    }
    S.at(ytpos) *= Scale;
    
  } // End of for ( unsigned int ytpos = 0; ytpos < yt.size(); ytpos++ )
  
  /*
   * Last diffy, diffS.
   */
  diffS.back() = S.back();
  diffy.back() = y.back();
  
  /*
   * Update object fields.
   */
  gobj["S"]   = S;
  gobj["D"]   = D;
  
  /*
   * Calculate Loglikelihood.
   */
  for ( unsigned int i=0; i < diffy.size(); ++i ) {
    if ( diffy.at(i) > 0 ) {
      if ( diffS.at(i) == 0.0 ) {
        gobj["LL"] = -std::numeric_limits<double>::infinity();
        return;
      }
      LL += diffy.at(i) * log(diffS.at(i));
    }
  } // End of for ( unsigned int i=0; i < diffy.size(); ++i ).
  

  //Calculate SPPE
  double SPPE = 0;
  int tend = 0;
  tend = diffy.size()-1;
  SPPE =(y.at(tend)-y.at(0)*S.at(tend))*100 / y.at(0);

  // End of SPPE

  /*
  * Calculate Sum of squares
  */
  double squares = 0.00;
  for ( unsigned int i=0; i < diffy.size(); ++i ) {
    squares += (y.at(i)-y.at(0)*S.at(i))*(y.at(i)-y.at(0)*S.at(i));
  }

  /*
   * Update object fields.
   */
  gobj["squares"] = squares;
  gobj["SPPE"] = SPPE;
  gobj["LL"]  = LL;
  gobj["zt"]   = z;
  
  
} // End of void guts_engine( Rcpp::List gobj, Rcpp::NumericVector par, Rcpp::Nullable<Rcpp::NumericVector > z_dist = R_NilValue).
