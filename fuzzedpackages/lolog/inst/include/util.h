#ifndef UTILH_
#define UTILH_

#include <Rcpp.h>
#include <RcppCommon.h>
#include <cmath>
#include <memory>
#include<vector>
#include <boost/shared_ptr.hpp>
namespace lolog{

using namespace Rcpp;

const double doubleTolerance = .0000000001;

/*!
 * determine if two numbers are within a tolerance limit of one another
 */
inline bool near(double a, double b){
    return a <= b + doubleTolerance && a >= b - doubleTolerance;
}


/*!
 * Unwraps either an XPtr or Reference Class representation (through an Rcpp Module)
 * of a C++ class
 */
template<class T>
boost::shared_ptr<T> unwrapRobject(const SEXP& s){

    if ( TYPEOF(s) == EXTPTRSXP ) {
        XPtr<T > xp(s);
        return xp->template vShallowCopy<T>();
        //return boost::shared_ptr<T>(new T(*xp));
    }else if ( TYPEOF(s) == S4SXP ) {
        Rcpp::S4 s4obj( s );
        Rcpp::Environment env( s4obj );
        Rcpp::XPtr<T> xp( env.get(".pointer") );
        return xp->template vShallowCopy<T>();
        //return boost::shared_ptr<T>(new T(*xp));
    }
    Rcpp::Rcout << TYPEOF(s);
    ::Rf_error( "unwrapRobject: supplied object is not of correct type." );
}


/*!
 * Creates an R Reference class object from an object.
 * Object must be of a module exported class
 */
template<class T>
Rcpp::RObject wrapInReferenceClass(const T& obj,std::string className){
    XPtr< T > xp = (&obj)->template vShallowCopyXPtr<T>();
    //XPtr< T > xp(new T(obj));
#ifndef INSIDE
    Language call( "new", Symbol( className ),xp);
    return call.eval();
#endif
#ifdef INSIDE
    // when called from RInside, our reference classes are not available.
    return xp;
#endif
}


/*!
 * n choose k optimized for cases where n<k a lot
 */
inline double nchoosek(double n,double k){
    return (n<k) ? 0.0 : Rf_choose(n,k);
}


/*!
 * returns the first index where element occurs in vec
 */
template<class T>
inline int indexOf(const T &element,const std::vector<T> &vec){
    for(int i=0;i<vec.size();i++){
        if(vec[i]==element)
            return i;
    }
    return -1;
}


/*!
 * Coerse a type into a string
 */
template<class T>
std::string asString(const T &item){
    std::ostringstream ss;
    ss << item;
    return ss.str();
}



/*!
 * An enumeriation of the types of edges
 */
enum EdgeDirection {UNDIRECTED, IN, OUT};

}
#endif /* UTILH_ */
