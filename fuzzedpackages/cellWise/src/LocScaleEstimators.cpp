
#include "LocScaleEstimators.h"
#include "RcppArmadilloExtensions/sample.h"


arma::uvec LocScaleEstimators::sample(const arma::uvec & sampleFrom,
                       const unsigned int nbSamples,
                       const bool replace) {
  //acess function for RcppArmadillo::sample
  arma::uvec output = Rcpp::RcppArmadillo::sample(sampleFrom, nbSamples, replace);
  return(output);
}

//////////////
// Location //
//////////////

void LocScaleEstimators::locWeightBiweight(arma::vec &x) {
  // Biweight weight function to be used in location M estimators
  // args:
  //   x:
  //   b:
  // Returns: 
  //   w: weights
  //
  x = x *  1.482602218505602 / 3;
  x.transform( [](double val) { return (1 - std::pow(val,2)); } );
  x.transform( [](double val) { return (std::pow((val + std::abs(val)) / 2, 2)); } );
}


void LocScaleEstimators::locWeightHuber15(arma::vec &x) {
  // Huber weight function to be used in location M estimators
  // args:
  // x:
  // b:
  // Returns:
  // w: weights
  // 
  x.transform( [](double val) { 
    if(std::abs(val) < 1.5 ) {
      return(double(1));
    } else{
      return(1.5 / std::abs(val));
    }
  });
}

void LocScaleEstimators::locWeightTanh154(arma::vec &x) {
  // Hyperbolic Tangent weight function to be used in location M estimators
  // args:
  // x:
  // b:
  // c:
  // k:
  // A:
  // B:
  // Returns:
  // w: weights
  // 
  
  double b = 1.5;
  double c = 4;
  double A = 0.7532528;
  double B = 0.8430849;
  double k = 4.1517212;
  
  x.transform( [b, c , k, A, B](double val) { 
    if(std::abs(val) < b ) {
      return(double(1));
    }else if (std::abs(val) > c){
      return(double(0));
    }else{
      return(sqrt(A * (k - 1)) * std::tanh(0.5 * sqrt((k - 1) *
             std::pow(B, 2) / A) * (c - std::abs(val)))/(std::abs(val)));
    }
  });
}

double LocScaleEstimators::loc1StepM(const arma::vec &x, std::function<void (arma::vec&)> weightFunction, 
                 double initLoc, double initScale, double precScale) {
  // Computes the first step of an algorithm for
  //   a location M-estimator using a weight function.
  // The weights are computed on (med,mad)-standardized data.
  // args:
  //   b:
  //   c:
  //   maxit:
  //   precScale:
  //   Returns:
  //   
  
  if (x.is_empty()) {
    return 0.0;
  }
  arma::uvec idx = arma::find_finite(x); //indices of finite elements 
  double m0,s0,mu;
  if(!std::isfinite(initLoc)) { 
    m0 = arma::median(x.elem(idx));
  } else {
    m0 = initLoc;
  }
  if(!std::isfinite(initScale)) { 
    s0 =  1.482602218505602 * arma::median(arma::abs(x.elem(idx) - m0));
  } else {
    s0 = initScale;
  }
  if(s0 > precScale) {
    arma::vec w = (x.elem(idx) - m0) / s0;
    weightFunction(w);
    mu =   arma::sum(x.elem(idx) % w) / arma::sum(w);
  } else {
    mu = m0;
  }
  return(mu);
}


///////////
// SCALE //
///////////



void LocScaleEstimators::psiTanh(arma::vec &x, double b, double c, double k, double A,
                   double B) {
  // Psi function of the hyperbolic tangent estimator. The input values
  // have 2 degrees of freedom, be carefull in changing them!
  //   
  //   args:
  // x: univariate input data
  // b: parameter of the hyperbolic tangent psi function
  // c: parameter of the hyperbolic tangent psi function
  // k: parameter of the hyperbolic tangent psi function
  // A: parameter of the hyperbolic tangent psi function
  // B: parameter of the hyperbolic tangent psi function
  // Returns:
  // u: a vector with location, scale and the transformed input data
  // 
  
  x.transform( [c](double val) { return (std::abs(val) > c ? double(0) : val); } );
  
  x.transform( [b, c , k, A, B](double val) { 
    return (std::abs(val) > b ? sqrt(A * (k - 1)) *
            std::tanh(0.5 * sqrt((k - 1) * std::pow(B,2) / A) * (c - std::abs(val)))
              * ((val > 0) - (val < 0)) : val);
  });
}


void LocScaleEstimators::rhoTanh154(arma::vec &x) {
  // Computes the hyperbolic tangent rho function for b = 1.5 and c = 4
  // args:
  //   x:
  //   Returns:
  //   
  //   integrate(function(x) psiTanh(x)^2*dnorm(x),-Inf,Inf)$value *2
  double b = 1.5;
  double c = 4;
  double A = 0.7532528;
  double B = 0.8430849;
  double k = 4.1517212;
  psiTanh(x, b, c, k, A, B);
  x = arma::pow(x, 2) / 1.506506;
}



void LocScaleEstimators::rhoHuber25(arma::vec &x) {
  // Computes the Huber rho function for b = 2.5*qnorm(0.75) = 1.6862...
  //   args:
  //   x:
  //   Returns:
  //   
  //   integrate(pmin(abs(x^2), c^2)*dnorm(x), -Inf, Inf)$value *2 = 1.688936
  // 2.8434 = (2.5*qnorm(0.75))^2
  x = arma::pow(x, 2);
  x.transform( [](double val) { return (val >   2.8433526444973292 ?   2.8433526444973292 : val);});
  x = x /  1.688942410165249;
}


void LocScaleEstimators::rhoHuber15(arma::vec &x) {
  // Computes the Huber rho function for b 1.5
  //   args:
  //   x:
  //   Returns:
  //   
  //   integrate(pmin(abs(x^2), c^2)*dnorm(x), -Inf, Inf)$value *2
  // 2.25 = 1.5^2
  x = arma::pow(x, 2);
  x.transform( [](double val) { return (val >  2.25 ?  2.25 : val);});
  x = x / 1.556931;
}


double LocScaleEstimators::scale1StepM(const arma::vec &x,
                                       std::function<void (arma::vec&)> rhoFunction,
                   double initScale, double precScale) {
  // Computes the first step of an algorithm for
  //   a scale M-estimator using the given rho function.
  // The scatter is computed relative to zero.
  // args:
  //   x:
  //   rhoFunction:
  //   maxit:
  //   precScale:
  //   Returns:
  //   
  if (x.is_empty()) {
    return 0.0;
  }
  arma::uvec idx = arma::find_finite(x); 
  
  double s0;
  
  if(!std::isfinite(initScale)) { 
    s0 =  1.482602218505602 * arma::median(arma::abs(x.elem(idx)));
  } else {
    s0 = initScale;
  }
  if(s0 < precScale) {
    s0 = 0;
  } else {
    arma::vec rho = x.elem(idx) / s0;
    rhoFunction(rho);
    s0 = s0 * sqrt(sum(rho) / (0.5 * idx.size()));
  }
  return(s0); 
}



//////////////
// LOCSCALE //
//////////////

LocScaleEstimators::locscale LocScaleEstimators::uniMcd(arma::vec y, double alpha) {
  // Computes reweighted univariate MCD estimator
  // Still needs to be finetuned: check for finite values, return 0 instead of NA
  // assumes alpha is based on y.size(), i.e. before checking for non finite values
  
  
  locscale out ={0, 1, 0, 1, 1, 1, arma::zeros<arma::uvec>(y.size())}; //output vector
  
  arma::uvec indfin = arma::find_finite(y);
  
  
  
  //determine quan
  int quan = std::ceil(indfin.size() * alpha);
  y = y.elem(indfin);
  
  int len = y.size() - quan + 1;
  
  if (len == 1) {
    out.loc = arma::mean(y);
    out.scale = sqrt(arma::var(y));
    out.rawloc = out.loc;
    out.rawscale = out.scale;
  } else {
    arma::vec sh(len, arma::fill::zeros);
    arma::vec sh2(len, arma::fill::zeros);
    arma::vec sq(len, arma::fill::zeros);
    
    
    arma::uvec I = arma::sort_index(y);
    y = y(I);
    sh(0) = arma::sum(y.head(quan));
    sh2(0) = std::pow(sh(0), 2) / quan;
    sq(0) = arma::conv_to<double>::from(arma::sum(arma::pow(y.head(quan), 2)) - sh2(0));
    
    for (int i = 1; i < len ; i++) {
      sh(i) = sh(i - 1) - y(i - 1) + y(i + quan - 1);
      sh2(i) = std::pow(sh(i), 2) / quan;
      sq(i) = sq(i - 1) - std::pow(y(i - 1), 2) + std::pow(y(i + quan - 1), 2) - sh2(i) + sh2(i - 1);
    }
    
    double sqmin = sq.min();
    arma::uvec Isq = arma::sort_index(sq);
    unsigned int ntied =  std::count(sq.begin(), sq.end(), sqmin);
    // number of duplicates
    
    arma::uvec idtied = Isq.head(ntied);
    double initmean = arma::mean(sh(idtied))/quan;
    arma::vec sqres = arma::pow((y - initmean), 2);
    arma::vec sortsqres = arma::sort(sqres); // could be replaced by std::partial_sort
    
    if(y.size() < 10){
      arma::vec allc = {4.252,1.970,2.217,1.648,1.767,1.482,1.555};
      out.cfac1 = allc(y.size()-3);
    } else {
      if(y.size() % 2 == 0) {
        out.cfac1 = (double) y.size() / ( (double) y.size() - 3.0);
      } else {
        out.cfac1 = (double) y.size() /( (double) y.size() - 3.4);
      }    
    }
    
    double rawsqscale = std::pow(out.cfac1, 2) * sortsqres(quan - 1) / R::qchisq(((double) quan) / ((double) y.size()), 1, true, false);
    double cutoff = rawsqscale *  R::qchisq(0.975, 1, true, false) ;
    arma::uvec weights = arma::find(sqres <= cutoff);

    
    out.rawloc = initmean;
    out.rawscale = std::sqrt(rawsqscale);
    out.loc = arma::sum(y(weights)) / weights.size();
    
    arma::uvec weightsOut = arma::zeros<arma::uvec>(y.size());
    weightsOut(weights).ones();
    weightsOut(I) = weightsOut; // I allows us to put the weights in the right order
    out.weights(indfin) = weightsOut; 
    
    
    double tempscale = std::sqrt(std::max(0.0, arma::sum(arma::pow(y(weights) - out.loc, 2))) / (weights.size() - 1));
    if(y.size() < 16){ 
      arma::vec alld = {1.475,1.223,1.253,1.180,1.181,1.140, 1.143,1.105,1.114,1.098,1.103,1.094,1.093};
      out.cfac2 = alld(y.size() - 3) ; 
    } else {
      out.cfac2 = y.size()/(y.size()-1.4);
    }    
    out.scale = out.cfac2 * 1.0835 * tempscale;
  }
  
  return(out);
}



/////////////////////////
// Vectorized LOCSCALE //
/////////////////////////

LocScaleEstimators::Xlocscale LocScaleEstimators::estLocScale(const arma::mat &X,
                                                              unsigned int nLocScale,
                                                              int type,
                                                              double precScale,
                                                              const int center,
                                                              double alpha){
  // estimation of (robust) location and scale of each of the columns in X
  // nLocScale is the number of observations used to calculate the loc/scale
  // If it is smaller than X.n_rows, nLocScale (non-NA) observations are sampled from 
  // every col of X
  // quan is only used for mcd
  // center is boolean determining whether or not to center before calculating scale
  //        it is not used for unimcd options
  // type 0 = biweight location + huber 2.5 scale
  // type 1 = huber location + huber scale
  // type 2 = unimcd + wrapping location/scale
  // type 3 = unimcd
  
  LocScaleEstimators::Xlocscale out;
  out.loc = arma::zeros(X.n_cols);
  out.scale =  arma::zeros(X.n_cols);
  
  if (nLocScale == 0) {nLocScale = X.n_rows;}
  const bool nLocScaleFlag = nLocScale < X.n_rows;
  nLocScale = nLocScale > X.n_rows ? X.n_rows : nLocScale;
  
  switch(type) {
  case 0:
    for (unsigned int i = 0; i < X.n_cols; i++)
    {
      arma::vec tempCol = X.col(i);
      if (nLocScaleFlag) {
        arma::uvec sampleFrom = arma::find_finite(tempCol);
        if (sampleFrom.size() > nLocScale){
        tempCol = tempCol(LocScaleEstimators::sample(sampleFrom, nLocScale, false));
        }
      }
      if (center) {
      out.loc(i) = LocScaleEstimators::loc1StepM(tempCol, 
              LocScaleEstimators::locWeightBiweight,
              arma::datum::nan,
              arma::datum::nan, precScale);
      } else {
        out.loc(i) = 0;
      }
      out.scale(i) = LocScaleEstimators::scale1StepM(tempCol-out.loc(i),
                LocScaleEstimators::rhoHuber25,
                arma::datum::nan, precScale);
    }
    break;
  case 1:
    for (unsigned int i = 0; i < X.n_cols; i++)
    {
      arma::vec tempCol = X.col(i);
      if (nLocScaleFlag) {
        arma::uvec sampleFrom = arma::find_finite(tempCol);
        if (sampleFrom.size() > nLocScale){
          tempCol = tempCol(LocScaleEstimators::sample(sampleFrom, nLocScale, false));
        }
      }
      if (center) {
      out.loc(i) = LocScaleEstimators::loc1StepM(tempCol, 
              LocScaleEstimators::locWeightHuber15,
              arma::datum::nan,
              arma::datum::nan, precScale);
      } else {
        out.loc(i) = 0;
      }
      out.scale(i) = LocScaleEstimators::scale1StepM(tempCol-out.loc(i),
                LocScaleEstimators::rhoHuber15,
                arma::datum::nan, precScale);
    }
    break;
  case 2:
    for (unsigned int i = 0; i < X.n_cols; i++) {
      arma::vec tempCol = X.col(i);
      if (nLocScaleFlag) {
        arma::uvec sampleFrom = arma::find_finite(tempCol);
        if (sampleFrom.size() > nLocScale){
          arma::uvec sampledIndices = LocScaleEstimators::sample(sampleFrom, nLocScale, false);
          tempCol = tempCol(sampledIndices);
        }
      }
      
      // initial estimate
      LocScaleEstimators::locscale uni = LocScaleEstimators::uniMcd(tempCol, alpha);
      double m0 = uni.loc;
      double s0 = uni.scale;
      
      // 1 step M estimate for location
      double m1  =  LocScaleEstimators::loc1StepM(tempCol,
                                                  LocScaleEstimators::locWeightTanh154,
                                                  m0, s0, precScale);
      out.loc(i) = m1;
      out.scale(i) = s0;
    }
    break;
  case 3:
    for (unsigned int i = 0; i < X.n_cols; i++) {
      arma::vec tempCol = X.col(i);
      if (nLocScaleFlag) {
        arma::uvec sampleFrom = arma::find_finite(tempCol);
        if (sampleFrom.size() > nLocScale){
          tempCol = tempCol(LocScaleEstimators::sample(sampleFrom, nLocScale, false));
        }
      }
      LocScaleEstimators::locscale uni = LocScaleEstimators::uniMcd(tempCol, alpha);
      out.loc(i) = uni.loc;
      out.scale(i)  = uni.scale;
    }
    break;
  case 4:
    for (unsigned int i = 0; i < X.n_cols; i++) {
      arma::vec tempCol = X.col(i);
      if (nLocScaleFlag) {
        arma::uvec sampleFrom = arma::find_finite(tempCol);
        if (sampleFrom.size() > nLocScale){
          tempCol = tempCol(LocScaleEstimators::sample(sampleFrom, nLocScale, false));
        }
      }
      LocScaleEstimators::locscale uni = LocScaleEstimators::uniMcd(tempCol, alpha);
      out.loc(i) = uni.rawloc;
      out.scale(i)  = uni.rawscale;
    }
    break;
  default:
    for (unsigned int i = 0; i < X.n_cols; i++) {
      arma::vec tempCol = X.col(i);
      if (nLocScaleFlag) {
        arma::uvec sampleFrom = arma::find_finite(tempCol);
        if (sampleFrom.size() > nLocScale){
          tempCol = tempCol(LocScaleEstimators::sample(sampleFrom, nLocScale, false));
        }
      }
      // 1 step M estimate, starting from median & mad
      double m1;
      if (center) {
       m1  =  LocScaleEstimators::loc1StepM(tempCol, LocScaleEstimators::locWeightTanh154,
                                                  arma::datum::nan,
                                                  arma::datum::nan, precScale);
      } else {
        m1 = 0;
      }
      double s1 = LocScaleEstimators::scale1StepM(tempCol - m1, LocScaleEstimators::rhoTanh154,
                                                  arma::datum::nan, precScale);
      
      // finite sample correction  
      double cn = 1;
      if( tempCol.size() < 16) {
        arma::vec allc = {1.728728, 1.329473, 1.391057, 1.241474, 1.222204, 1.165270,
                          1.168463, 1.130316, 1.129584, 1.107986, 1.107362, 1.094637,
                          1.090304 };
        cn = allc(tempCol.size() - 3);
      } else {
        cn = tempCol.size() / (tempCol.size() - 1.208);
      }
      s1 = s1 * cn;
      out.loc(i) = m1;
      out.scale(i) = s1;
    }
  }
  return(out);
}


///////////
// Ranks //
///////////

arma::vec LocScaleEstimators::rank(arma::vec& v) {
  // rank with averaging ties
  //https://stackoverflow.com/questions/30822729/create-ranking-for-vector-of-double/30827731#30827731
  // does not handle NAs
  arma::vec w(v.size());
  std::iota(w.begin(), w.end(), 0);
  std::sort(w.begin(), w.end(), 
            [&v](std::size_t i, std::size_t j) { return v(i) < v(j); });
  
  arma::vec r(w.size());
  for (std::size_t n, i = 0; i < w.size(); i += n)
  {
    n = 1;
    while (i + n < w.size() && v(w(i)) == v(w(i+n))) ++n;
    for (std::size_t k = 0; k < n; ++k)
    {
      r(w(i+k)) = i + (n + 1) / 2.0; // average rank of n tied values
      // r[w[i+k]] = i + 1;          // min 
      // r[w[i+k]] = i + n;          // max
      // r[w[i+k]] = i + k + 1;      // random order
    }
  }
  return r;
}
