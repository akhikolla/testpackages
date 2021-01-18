#include "synlik.h"

Rcpp::NumericMatrix rickerSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                const int & nBurn, const bool & randInit, const double & initVal);

Rcpp::NumericMatrix genRickerSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                   const int & nBurn, const bool & randInit, const double & initVal);

Rcpp::NumericMatrix pennySimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                               const int & nBurn, const bool & randInit, const double & initVal);

Rcpp::NumericMatrix hassellSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                 const int & nBurn, const bool & randInit, const double & initVal);

Rcpp::NumericMatrix maynardSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                 const int & nBurn, const bool & randInit, const double & initVal);

Rcpp::NumericMatrix varleySimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                const int & nBurn, const bool & randInit, const double & initVal);

SEXP simpleModelsWrap(SEXP model, SEXP days, SEXP nSimul, SEXP param, SEXP nBurn, SEXP randInit, SEXP initVal)
{
  using namespace Rcpp;
  
  try{
    std::string model_ = as<std::string>(model);
    int days_ = as<int>(days);
    int nSimul_ = as<int>(nSimul);
    NumericMatrix param_ = as<NumericMatrix>(param);
    int nBurn_ = as<int>(nBurn);
    bool randInit_ = as<bool>(randInit);
    double initVal_ = as<double>(initVal);
    
    if(model_ == "ricker") return rickerSimul(days_, nSimul_,  param_,  nBurn_, randInit_,  initVal_);
    if(model_ == "genRicker") return genRickerSimul(days_, nSimul_,  param_,  nBurn_, randInit_,  initVal_);
    if(model_ == "penny") return pennySimul(days_, nSimul_,  param_,  nBurn_, randInit_,  initVal_);
    if(model_ == "hassell") return hassellSimul(days_, nSimul_,  param_,  nBurn_, randInit_,  initVal_);
    if(model_ == "maynard") return maynardSimul(days_, nSimul_,  param_,  nBurn_, randInit_,  initVal_);
    if(model_ == "varley") return varleySimul(days_, nSimul_,  param_,  nBurn_, randInit_,  initVal_);
    stop("Model name should be one of \"ricker\", \"genRicker\", \"penny\", \"hassell\", \"maynard\" or \"varley\" ");
    
  } catch( std::exception& __ex__){
    forward_exception_to_r(__ex__);
  } catch(...){
    ::Rf_error( "c++ exception (unknown reason)" );
  }
  return Rcpp::wrap(NA_REAL);
}

/*
  *  Ricker
*/ 
  Rcpp::NumericMatrix rickerSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                  const int & nBurn, const bool & randInit, const double & initVal)
{
    using namespace Rcpp;
    
    RNGScope scope;
    
    int nparam = param.ncol(); 
    int totDays = nBurn + days;
    bool multiparam = false;
    
    if(nparam != 3) stop("Wrong number of parameters");
    if(param.nrow() > 1) { multiparam = true; }
    if(multiparam == true && param.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(param(0, 0));
    double sigma = exp(param(0, 1));
    double phi = exp(param(0, 2));
    
    NumericVector procNoise( rnorm( totDays * nSimul ) );
    NumericVector initState(nSimul);
    if(randInit){ initState = runif(nSimul); } else { initState = initState + initVal; }
    NumericMatrix output( nSimul, days );
    
    NumericVector::iterator noiseIter = procNoise.begin();
    NumericVector::iterator initIter = initState.begin();
    
    double currState;
    
    for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
    {
      
      if( multiparam == true )
      {
        r = exp(param(iRow, 0));
        sigma = exp(param(iRow, 1));
        phi = exp(param(iRow, 2));
      }
      
      currState = *initIter;
      
      for(int iCol = 1; iCol <= nBurn; iCol++, noiseIter++){
        currState = r * currState * exp( - currState + *noiseIter * sigma );
      }
      
      output(iRow, 0) = R::rpois(phi * currState);
      
      for(int iCol = 1; iCol < days; iCol++, noiseIter++){
        currState = r * currState * exp( - currState + *noiseIter * sigma );
        output(iRow, iCol) = R::rpois(phi * currState);
      }
      
    }
    
    return output;
    
    
  }


/*
  * Generalized Ricker 
*/
  Rcpp::NumericMatrix genRickerSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                     const int & nBurn, const bool & randInit, const double & initVal)
{
    using namespace Rcpp;
    
    RNGScope scope;
    
    int nparam = param.ncol(); 
    int totDays = nBurn + days;
    bool multiparam = false;
    
    if(nparam != 4) stop("Wrong number of parameters");
    if(param.nrow() > 1) { multiparam = true; }
    if(multiparam == true && param.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(param(0, 0));
    double theta = exp(param(0, 1));
    double sigma = exp(param(0, 2));
    double phi = exp(param(0, 3));
    
    NumericVector procNoise( rnorm( totDays * nSimul ) );
    
    NumericVector initState(nSimul);
    if(randInit){ initState = runif(nSimul); } else { initState = initState + initVal; }  
    
    NumericMatrix output( nSimul, days );
    
    NumericVector::iterator noiseIter = procNoise.begin();
    NumericVector::iterator initIter = initState.begin();
    
    double currState;
    
    for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
    {
      
      if( multiparam == true )
      {
        r = exp(param(iRow, 0));
        theta = exp(param(iRow, 1));
        sigma = exp(param(iRow, 2));
        phi = exp(param(iRow, 3));
      }
      
      currState = *initIter;
      
      for(int iCol = 1; iCol <= nBurn; iCol++, noiseIter++){
        currState = r * currState * exp( - pow( currState, theta ) + *noiseIter * sigma );
      }
      
      output(iRow, 0) = R::rpois(phi * currState);
      
      for(int iCol = 1; iCol < days; iCol++, noiseIter++){
        currState = r * currState * exp( - pow( currState, theta ) + *noiseIter * sigma );
        output(iRow, iCol) = R::rpois(phi * currState);
      }
      
    }
    
    return output;
    
    
  }




/*
  * Pennycuick 
*/
  
  Rcpp::NumericMatrix pennySimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                 const int & nBurn, const bool & randInit, const double & initVal)
{
    using namespace Rcpp;
    
    
    RNGScope scope;
    
    int nparam = param.ncol(); 
    int totDays = nBurn + days;
    bool multiparam = false;
    
    if(nparam != 4) stop("Wrong number of parameters");
    if(param.nrow() > 1) { multiparam = true; }
    if(multiparam == true && param.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(param(0, 0));
    double a = exp(param(0, 1));
    double sigma = exp(param(0, 2));
    double phi = exp(param(0, 3));
    
    NumericVector procNoise( rnorm( totDays * nSimul ) );
    NumericVector initState(nSimul);
    if(randInit){ initState = runif(nSimul); } else { initState = initState + initVal; } 
    NumericMatrix output( nSimul, days );
    
    NumericVector::iterator noiseIter = procNoise.begin();
    NumericVector::iterator initIter = initState.begin();
    
    double currState;
    
    for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
    {
      
      if( multiparam == true )
      {
        r = exp(param(iRow, 0));
        a = exp(param(iRow, 1));
        sigma = exp(param(iRow, 2));
        phi = exp(param(iRow, 3));
      }
      
      currState = *initIter;
      
      for(int iCol = 1; iCol <= nBurn; iCol++, noiseIter++){
        currState = ( r*currState/( 1 + exp(-a*(1 - currState))) ) * exp( *noiseIter * sigma);
      }
      
      output(iRow, 0) = R::rpois(phi * currState);
      
      for(int iCol = 1; iCol < days; iCol++, noiseIter++){
        currState = ( r*currState/( 1 + exp(-a*(1 - currState))) ) * exp( *noiseIter * sigma);
        output(iRow, iCol) = R::rpois(phi * currState);
      }
      
    }
    
    return output;
    
    
  }


/*
  * Hassell 
*/
  
  Rcpp::NumericMatrix hassellSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                   const int & nBurn, const bool & randInit, const double & initVal)
{
    using namespace Rcpp;
    
    RNGScope scope;
    
    int nparam = param.ncol(); 
    int totDays = nBurn + days;
    bool multiparam = false;
    
    if(nparam != 4) stop("Wrong number of parameters");
    if(param.nrow() > 1) { multiparam = true; }
    if(multiparam == true && param.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(param(0, 0));
    double b = exp(param(0, 1));
    double sigma = exp(param(0, 2));
    double phi = exp(param(0, 3));
    
    NumericVector procNoise( rnorm( totDays * nSimul ) );
    NumericVector initState(nSimul);
    if(randInit){ initState = runif(nSimul); } else { initState = initState + initVal; } 
    NumericMatrix output( nSimul, days );
    
    NumericVector::iterator noiseIter = procNoise.begin();
    NumericVector::iterator initIter = initState.begin();
    
    double currState;
    
    for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
    {
      
      if( multiparam == true )
      {
        r = exp(param(iRow, 0));
        b = exp(param(iRow, 1));
        sigma = exp(param(iRow, 2));
        phi = exp(param(iRow, 3));
      }
      
      currState = *initIter;
      
      for(int iCol = 1; iCol <= nBurn; iCol++, noiseIter++){
        currState = ( r*currState/pow(1.0 + currState, b) ) * exp( *noiseIter * sigma);
      }
      
      output(iRow, 0) = R::rpois(phi * currState);
      
      for(int iCol = 1; iCol < days; iCol++, noiseIter++){
        currState = ( r*currState/pow(1.0 + currState, b) ) * exp( *noiseIter * sigma);
        output(iRow, iCol) = R::rpois(phi * currState);
      }
      
    }
    
    return output;
    
    
  }



/*
  * Maynard Smith 
*/
  Rcpp::NumericMatrix maynardSimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                   const int & nBurn, const bool & randInit, const double & initVal)
{
    using namespace Rcpp;
    
    
    RNGScope scope;
    
    int nparam = param.ncol(); 
    int totDays = nBurn + days;
    bool multiparam = false;
    
    if(nparam != 4) stop("Wrong number of parameters");
    if(param.nrow() > 1) { multiparam = true; }
    if(multiparam == true && param.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(param(0, 0));
    double b = exp(param(0, 1));
    double sigma = exp(param(0, 2));
    double phi = exp(param(0, 3));
    
    NumericVector procNoise( rnorm( totDays * nSimul ) );
    NumericVector initState(nSimul);
    if(randInit){ initState = runif(nSimul); } else { initState = initState + initVal; } 
    NumericMatrix output( nSimul, days );
    
    NumericVector::iterator noiseIter = procNoise.begin();
    NumericVector::iterator initIter = initState.begin();
    
    double currState;
    
    for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
    {
      
      if( multiparam == true )
      {
        r = exp(param(iRow, 0));
        b = exp(param(iRow, 1));
        sigma = exp(param(iRow, 2));
        phi = exp(param(iRow, 3));
      }
      
      currState = *initIter;
      
      for(int iCol = 1; iCol <= nBurn; iCol++, noiseIter++){
        currState = ( r*currState/( 1.0 + pow(currState, b)) ) * exp( *noiseIter * sigma );
      }
      
      output(iRow, 0) = R::rpois(phi * currState);
      
      for(int iCol = 1; iCol < days; iCol++, noiseIter++){
        currState = ( r*currState/( 1.0 + pow(currState, b)) ) * exp( *noiseIter * sigma );
        output(iRow, iCol) = R::rpois(phi * currState);
      }
      
    }
    
    return output;
    
    
  }



/*
  * Varley
*/
  Rcpp::NumericMatrix varleySimul(const int & days, const int & nSimul, const Rcpp::NumericMatrix & param, 
                                  const int & nBurn, const bool & randInit, const double & initVal)
{
    using namespace Rcpp;
    
    
    RNGScope scope;
    
    int nparam = param.ncol(); 
    int totDays = nBurn + days;
    bool multiparam = false;
    
    if(nparam != 5) stop("Wrong number of parameters");
    if(param.nrow() > 1) { multiparam = true; }
    if(multiparam == true && param.nrow() != nSimul) 
      stop("Number of parameters vectors is different from the number of simulations");
    
    double r = exp(param(0, 0));
    double b = exp(param(0, 1));
    double C = exp(param(0, 2));
    double sigma = exp(param(0, 3));
    double phi = exp(param(0, 4));
    
    NumericVector procNoise( rnorm( totDays * nSimul ) );
    NumericVector initState(nSimul);
    if(randInit){ initState = runif(nSimul); } else { initState = initState + initVal; } 
    NumericMatrix output( nSimul, days );
    
    NumericVector::iterator noiseIter = procNoise.begin();
    NumericVector::iterator initIter = initState.begin();
    
    double currState;
    
    for(int iRow = 0; iRow < nSimul; iRow++, initIter++)
    {
      
      if( multiparam == true )
      {
        r = exp(param(iRow, 0));
        b = exp(param(iRow, 1));
        C = exp(param(iRow, 2));
        sigma = exp(param(iRow, 3));
        phi = exp(param(iRow, 4));
      }
      
      currState = *initIter;
      
      for(int iCol = 1; iCol <= nBurn; iCol++, noiseIter++){
        currState = ( currState <= C ? r*currState : r*pow(currState, 1-b) )*exp( *noiseIter * sigma);
      }
      
      output(iRow, 0) = R::rpois(phi * currState);
      
      for(int iCol = 1; iCol < days; iCol++, noiseIter++){
        currState = ( currState <= C ? r*currState : r*pow(currState, 1-b) )*exp( *noiseIter * sigma);
        output(iRow, iCol) = R::rpois(phi * currState);
      }
      
    }
    
    return output;
    
}
