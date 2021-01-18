/* @file GrangerDickeyFullet.h
 *
 * @author Youssef Hmamouche
 *
 * @brief Augmented Dickey-Fullet test od stationary, Granger test of causality
 *
 * @date 01-04-2016
 *
 */
#ifndef TESTS_HEADER
#define TESTS_HEADER

#include <Rcpp.h>

#include "struct.h"

class DickeyFuller
{
    Struct::CVDouble tS;
    double df;
    unsigned lag;
    std::string testResult;
    unsigned nl;
    float Values [3][6] = {{25,50,100,250,500,501}  // Size of the sample, 501 means > 500
        ,{-3.60,-3.50,-3.45,-3.43,-3.42,-3.41} // values for 5% of risk
        ,{-3,-2.93,-2.89,-2.88,-2.87,-2.86 }};  // values for 1% of risk
    
public:
    
    double SBC;
    DickeyFuller () {}
    // First constructor
    DickeyFuller (const Rcpp::NumericVector & tS_, int lag_ = 0); // throw ();
    
    // Second constructor
    DickeyFuller (const Struct::CVDouble & tS_, int lag_ = 0); // throw ();
    
    ~DickeyFuller () {}
    
    unsigned    getLag ()  { return lag; }
    
    double      getDF ()   { return df;  }
    
    double      get5CriticalValue () {
        double res(0);
        
        if (nl > 500)
            return Values [1][5];
        
        for (unsigned i = 0 ; i < 5; i ++)
            if (nl <= Values [0][i]) {
                res = Values [1][i];
                break;
            }
        
        return res;
    }
    
    double      get1CriticalValue () {
        double res(0);
        if (nl > 500)
            return Values [1][5];
        for (unsigned i = 0 ; i < 5; i ++)
            if (nl <= Values [0][i]) {
                res = Values [2][i];
                break;
            }
        return res;
    }
    
    void summary ();
};

int order (const Struct::CVDouble & tS, int p);

#endif
