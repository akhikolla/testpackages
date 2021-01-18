#ifndef __sqp_misc_section_included__        // if include guard for 'misc/sqp_section.h' is undefined
#define __sqp_misc_section_included__        // define include guard for 'misc/sqp_section.h'

#include "sqp.h"

namespace sqp {
namespace misc{



inline void section(std::string x,
                    bool print)
{
  if(print)
  {
    Rcpp::Rcout << "\r                                                                                              ";
    Rcpp::Rcout << "\r" << x << "\n";
  }
} 



}
}

#endif                                      // end of include guard for 'misc/sqp_section.h'
