/*
 * rcppExport.h
 *
 *  Created on: 13.10.2010
 *      Author: daniel
 */

#ifndef RCPPEXPORT_H_
#define RCPPEXPORT_H_

#include <RcppCommon.h>

// load the library defining the templated class
#include <set>

// declaring the partial specialization
namespace Rcpp
{
namespace traits
{
template <typename T> class Exporter<std::set<T> >;
}
}


// this must appear after the specialization,
// otherwise the specialization will not be seen by Rcpp types
#include <Rcpp.h>


// only now the implementation of the new Exporter class:
namespace Rcpp
{
namespace traits
{
template <typename T> class Exporter<std::set<T> >
{
public:
    // we need a constructor taking an SEXP
    Exporter(SEXP s) :
        _x(s)
    {
    }

    // and we need a method called "get" that returns an instance of the wished type:
    inline std::set<T>
    get()
    {
        // allocate return set
        std::set<T> x;

        // insert all elements in _x into the return set x:
        for (R_len_t i = 0; i < _x.size(); ++i)
        {
            x.insert(_x[i]);
        }

        // return the result
        return x;
    }

private:
    Rcpp::IntegerVector _x;
};
}
}


#endif /* RCPPEXPORT_H_ */
