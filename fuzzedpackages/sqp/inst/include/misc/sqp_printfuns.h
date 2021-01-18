#ifndef __sqp_misc_printfuns_included__        // if include guard for 'misc/sqp_printfuns.h' is undefined
#define __sqp_misc_printfuns_included__        // define include guard for 'misc/sqp_printfuns.h'

#include "sqp.h"

namespace sqp {
namespace misc{



inline void printfun(unsigned x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << x << "\n";
}

inline void printfun(int x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << x << "\n";
}

inline void printfun(double x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << x << "\n";
}

inline void printfun(arma::mat x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << arma::size(x);
  
  if(x.n_elem > 0)
  {
    Rcpp::Rcout 
    << ": " << x.min() << "  -  "
    << arma::abs(x).min() << "  -  "
    << x.max();
  }
  
  Rcpp::Rcout << "\n";
}

inline void printfun(arma::cube x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << arma::size(x);
  
  if(x.n_elem > 0)
  {
    Rcpp::Rcout 
    << ": " << x.min() << "  -  "
    << arma::abs(x).min() << "  -  "
    << x.max();
  }
  
  Rcpp::Rcout << "\n";
}

inline void printfun(arma::umat x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << arma::size(x);
  
  if(x.n_elem > 0)
  {
    Rcpp::Rcout 
    << ": " << x.min() << "  -  "
    << arma::abs(x).min() << "  -  "
    << x.max();
  }
  
  Rcpp::Rcout << "\n";
}


inline void printfun(arma::sp_mat x,
                     std::string lab)
{
  Rcpp::Rcout << lab << ": " << arma::size(x);
  
  if(x.n_elem > 0)
  {
    Rcpp::Rcout 
    << ": " << x.min() << "  -  "
    << arma::abs(x).min() << "  -  "
    << x.max();
  }
  
  Rcpp::Rcout << "\n";
} 



}
}

#endif                                        // end of include guard for 'misc/sqp_printfuns.h'
