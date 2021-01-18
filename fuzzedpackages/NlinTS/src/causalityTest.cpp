/**
 *
 * @file    causalityTest.cpp
 *
 * @authors Hmamouche Youssef
 *
 * @date    07/2017
 *
 * @version V1.0
 *
 * @brief  the Granger CausalityT est
 *
 **/

#include <Rcpp.h>
#include "../inst/include/tests.h"
#include "../inst/include/Ftable.h"
#include "../inst/include/fdist.h"
#include "../inst/include/matrixOperations.h"
#include "../inst/include/operateurs.h"
#include "../inst/include/causalityTest.h"


using namespace MatrixOperations;
using namespace Struct;
using namespace std;


causalityTest::causalityTest (Rcpp::NumericVector  ts1_,
                              Rcpp::NumericVector  ts2_, int lag_, bool d /* = false */) //throw (Exception)
{

    // Checking if the lag value is positif
    if (lag_ <= 0)
        throw Exception ("The lag parameter is incorrect!");

    lag  = lag_;

    for (const auto val:ts1_)
        ts1.push_back (val);

    for (const auto val:ts2_)
        ts2.push_back (val);

    if (ts1.size() != ts2.size())
       throw Exception ("The time series have not the same length!");


    // variables
    unsigned nl (ts1.size ());
    double RSS1(0), RSS0(0);

    CVDouble VECT;
    CMatDouble M (1);
    M[0] = ts1;

    /* the VAR model of the system (v1) */
    VECT = VECbivar(M, lag, d);
    RSS0 = VECT[0];

    /* the VAR model of the system (v1,v2) */
    M.push_back(ts2);
    VECT = VECbivar(M, lag, d);
    RSS1 = VECT[0];

    int T = nl - lag;

    // The Granger causality Index
    GCI = log (RSS0 / RSS1);

    // compute the F test
    Ftest = ((RSS0 - RSS1) / lag) / (RSS1 / (T - 2*lag - 1));

    // compute the p-value of the test
    p_value = getPvalue (Ftest , lag , T - 2*lag - 1);

    if (lag <= 20 and T - 2*lag - 1 <= 100)
        criticTest = ftable[T - 2*lag - 2][lag];
    else if (lag > 20 and T - 2*lag - 1 <= 100)
        criticTest = ftable[T - 2*lag - 2][20];
    else if (lag <= 20 and T - 2*lag - 1 > 100)
        criticTest = ftable[100][lag];
    else if (lag > 20 and T - 2*lag - 1 > 100)
        criticTest = ftable[100][20];
}

// The Summary function
void causalityTest::summary ()
{
    Rcpp::Rcout <<  "------------------------------------------------\n";
    Rcpp::Rcout <<  "        the Granger causality test" << "\n";
    Rcpp::Rcout <<  "------------------------------------------------\n";
    Rcpp::Rcout <<  "The lag parameter: p = "<< lag << "\n";
    Rcpp::Rcout <<  "The Granger causality Index: GCI = "<< GCI << "\n";
    Rcpp::Rcout <<  "The value of the F-test: "<< Ftest << "\n";
    Rcpp::Rcout <<  "The p_value of the F-test: "<< p_value << "\n";
    Rcpp::Rcout <<  "The critical value with 5% of risk:: "<< criticTest <<"\n";
}

// Get the p-value of the test
double causalityTest::get_p_value () {
    return p_value;
}

double causalityTest::get_gci () {
    return GCI;
}

// Get the  statistic of the test
double causalityTest::get_F_test () {
    return Ftest;
}
