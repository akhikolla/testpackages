// -*- mode: c++ -*-
///////////////////////////////////////////////////////////////////////////
/*
 This file is part of isqg, a R package for in silico quantitative genetics

              Copyright (C) 2018 Fernando H. Toledo CIMMYT
              
 * Filename: isqg.hpp
 
 * Description: C++ infrastructure for seamless integration of isqg package
 
 * Author: Fernando H. Toledo
 
 * Maintainer: Fernando H. Toledo
 
 * Created: Fr Mar 09 2018
 
  This program is free software; you can redistribute it and/or modify 
  it under the terms of the GNU General Public License as published by 
  the Free Software Foundation; either version 2 of the License, or 
  (at your option) any later version.

  This program is distributed in the hope that it will be useful, but 
  WITHOUT ANY WARRANTY; without even the implied warranty of 
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU 
  General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software Foundation, 
  Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
                                                        
  `` Far better an approximate answer to the right question, which is 
  often vague, than the exact answer to the wrong question, which can
  always be made precise ''
                         --John Tukey, Ann. Math. Stat. 33(1):13 1962
*/
///////////////////////////////////////////////////////////////////////////

# ifndef _HEADERS_H_
# define _HEADERS_H_

// standard libraries and others
# include <string>                               // text strings
# include <vector>                               // standard arrays
# include <tuple>                                // data structure
# include <algorithm>                            // standard algorithms
# include <iterator>                             // standard iterators
# include <RcppCommon.h>                         // allows template specializations
# include <Rcpp/XPtr.h>                          // smart external pointers
# include <boost/dynamic_bitset.hpp>             // runtime allocated bitsets
# include <boost/graph/adjacency_list.hpp>       // direct acyclic graphs
# include <boost/graph/breadth_first_search.hpp> // traverse graphs

// borrowed from glmmsr/RcppR6
// namespace scope for interface declaration
namespace isqg {
  namespace seamless {

  template <typename R6> class Trap ;

  }
}

// internal headers -- review estructures
# include <isqg/FwdDefs.hpp>
# include <isqg/FwdFuncs.hpp>
# include <isqg/Genetics.hpp>

// took from Extending Rcpp 
// non-intrusive extensions via template specialization declarations
namespace Rcpp {

  template <typename R6> SEXP wrap(const isqg::seamless::Trap<R6> & ) ;

  // for populations [WHAT]
  namespace traits {

    template <typename R6> class Exporter<isqg::seamless::Trap<R6> > ;

  }

  template <> Specie as(SEXP) ;
  template <> SEXP wrap(const Specie &) ;

  template <> Specimen as(SEXP) ;
  template <> SEXP wrap(const Specimen &) ;
  
  template <> Trait as(SEXP) ;
  template <> SEXP wrap(const Trait &) ;

}

# include <Rcpp.h> // [CAVEAT] included after declaration extensions

// borrowed from glmmmsr
// class traits specializations declarations
namespace isqg {
  namespace seamless {

    // in house c++ traits
    namespace traits {

      // generics
      template <typename R6> std::string name() { 
        Rcpp::stop("Unknown class") ; 
        return "" ; // do not apply
      }

      template <typename R6> std::string pack() {
        Rcpp::stop("Unknown class") ;
        return "" ; // do not apply
      }

      template <typename R6> std::string make() {
        Rcpp::stop("Unknown class") ;
        return "" ; // do not apply
      }

      // templated specializations
      template <typename R6> std::string name(const R6 &) { return name<R6>() ; }
      template <typename R6> std::string pack(const R6 &) { return pack<R6>() ; }
      template <typename R6> std::string make(const R6 &) { return make<R6>() ; }

    }
    
  }
}

// borrowed from RcppR6
// metaprogramming interface
namespace isqg {
  namespace seamless {

/*  
  -- Rcpp::is<> analogous --
  Throws a _compile time_ error if the type is not a known one.  
  Return false at runtime if _obj_ is not the required type. 
  Doesn't actually check that we inherit from any package... :(
*/
    template <typename R6> bool is(Rcpp::RObject obj) {
      return obj.inherits(traits::name<R6>().c_str()) ; // c string equivalent
    }

    template <typename R6> void check_xptr(Rcpp::XPtr<R6> xptr) {
      R6 * test = xptr ;
      if (test == NULL)
        Rcpp::stop("Pointer is NULL") ;
    }

/*
  -- templated Rcpp::as<> & Rcpp::wrap implementation --  
  Doesn't do anything different than Rcpp::XPtr... differences:
    * check pointer when receives a externalptr;
    * knows the type it holds, checking this at receivement; and
    * allow instantiate R instance on return.
*/
    template <typename R6> class Trap {

    public:
  
      // constructors
      Trap(SEXP       obj) : ptr(R2Cpp(obj))        { }
      Trap(const R6 & obj) : ptr(new R6(obj), true) { }

      // operators
             R6 & operator*   ()  const { return * ptr ;       }    
             R6 * operator->  ()  const { return & (* ptr) ;   }
      inline      operator R6*()        { return (R6 *)(ptr) ; }

      // Convert from an R object from R
      static Rcpp::XPtr<R6> R2Cpp(Rcpp::RObject obj) {
      
        if (is<R6>(obj)) {

          Rcpp::Environment env = Rcpp::as<Rcpp::Environment>(obj) ;
          Rcpp::XPtr<R6>    ptr = Rcpp::as<Rcpp::XPtr<R6>>(env[".ptr"]) ; // [FIXME]

          check_xptr<R6>(ptr) ;
          
          return ptr ;

        } else {

          Rcpp::stop("Expected an object of type " + traits::name<R6>()) ;

          return Rcpp::as<Rcpp::XPtr<R6>>(obj) ; // do not apply

        }
      }

      // Convert from a class instance from C++
      SEXP Cpp2R() const {
       
        Rcpp::Function    catcher = Rcpp::Environment("package:base")["getNamespace"] ;
        Rcpp::Environment package = catcher(traits::pack<R6>()) ;
        Rcpp::Environment factory = package[traits::make<R6>()] ;
        Rcpp::Function    builder = factory["new"] ;

        return builder(ptr) ;

      }

      // solely slot      
      Rcpp::XPtr<R6> ptr ;
      
    } ; // Trap<>
    
  }
}

// borrowed from glmmsr
// class traits inline definitions
namespace isqg {
  namespace seamless {
    namespace traits {

      template <> inline std::string name<Specie>()   { return "Specie" ;           }
      template <> inline std::string pack<Specie>()   { return "isqg" ;             }
      template <> inline std::string make<Specie>()   { return ".R_Specie_ctor" ;   }
      
      template <> inline std::string name<Specimen>() { return "Specimen" ;         }
      template <> inline std::string pack<Specimen>() { return "isqg" ;             }
      template <> inline std::string make<Specimen>() { return ".R_Specimen_ctor" ; }
      
      template <> inline std::string name<Trait>()    { return "Trait" ;            }
      template <> inline std::string pack<Trait>()    { return "isqg" ;             }
      template <> inline std::string make<Trait>()    { return ".R_Trait_ctor" ;    }

    }
  }
}

// took from Extending Rcpp
// non-intrusive extensions via template specialization implementation
namespace Rcpp {

  // --generic--
  template <typename R6> SEXP wrap(const isqg::seamless::Trap<R6> & instance) {
    return instance.Cpp2R() ;
  }

  // --specializations--
  template <> inline SEXP wrap(const Specie & obj) { 
    return wrap(isqg::seamless::Trap<Specie>(obj)) ;
  }
  template <> inline Specie as(SEXP obj) {
    return * (isqg::seamless::Trap<Specie>(obj)) ;
  }

  template <> inline SEXP wrap(const Specimen & obj) {
    return wrap(isqg::seamless::Trap<Specimen>(obj)) ;
  }
  template <> inline Specimen as(SEXP obj) {
    return * (isqg::seamless::Trap<Specimen>(obj)) ;
  }
  
  template <> inline SEXP wrap(const Trait & obj) {
    return wrap(isqg::seamless::Trap<Trait>(obj)) ;
  }
  template <> inline Trait as(SEXP obj) {
    return * (isqg::seamless::Trap<Trait>(obj)) ;
  }

  // --vectorized--
  namespace traits {
 
    template <typename R6> class Exporter<isqg::seamless::Trap<R6> > {

    public:

      Exporter(SEXP instance) : obj(isqg::seamless::Trap<R6>(instance)) {}

      inline isqg::seamless::Trap<R6> get() { return obj ; }

    private:

      isqg::seamless::Trap<R6> obj ;

    } ; // Exporter<>
    
  }

}

# endif // _HEADERS_H_

// \EOF
///////////////////////////////////////////////////////////////////////////
