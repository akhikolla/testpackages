// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; tab-width: 4 -*-

# ifndef RCPPTN_H
# define RCPPTN_H

#include <R_ext/Rdynload.h>

namespace RcppTN {
    typedef double(*Fun4)(double, double, double, double) ;
    typedef double(*Fun5)(double, double, double, double, double) ;

    inline double rtn1(double mean, double sd, double low, double high){
        static Fun4 rtn1_ = (Fun4)R_GetCCallable("RcppTN", "RcppTN_rtn1") ;
        return rtn1_(mean, sd, low, high) ;
    }

    inline double etn1(double mean, double sd, double low, double high){
        static Fun4 etn1_ = (Fun4)R_GetCCallable("RcppTN", "RcppTN_etn1") ;
        return etn1_(mean, sd, low, high) ;
    }

    inline double vtn1(double mean, double sd, double low, double high){
        static Fun4 vtn1_ = (Fun4)R_GetCCallable("RcppTN", "RcppTN_vtn1") ;
        return vtn1_(mean, sd, low, high) ;
    }
  
    inline double dtn1(double x, double mean, double sd, double low, double high){
        static Fun5 dtn1_ = (Fun5)R_GetCCallable("RcppTN", "RcppTN_dtn1") ;
        return dtn1_(x, mean, sd, low, high) ;
    }
    
    inline double enttn1(double mean, double sd, double low, double high){
        static Fun4 enttn1_ = (Fun4)R_GetCCallable("RcppTN", "RcppTN_enttn1") ;
        return enttn1_(mean, sd, low, high) ;
    }
}

# endif
