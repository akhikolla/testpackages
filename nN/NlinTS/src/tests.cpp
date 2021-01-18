/**
 * Class for Augmented Dickey Fuller test
 * @authors Hmamouche Youssef
 **/

#include <cmath>
#include <Rcpp.h>

#include "../inst/include/fdist.h"
#include "../inst/include/operateurs.h"
#include"../inst/include/matrixOperations.h"
#include "../inst/include/tests.h"

using namespace Struct;
using namespace MatrixOperations;
using namespace std;


/* Constructor for the Augmented Dickey Fuller test that takes an R dataFrame as argument */
DickeyFuller::DickeyFuller (const Rcpp::NumericVector & tS_, int lag_)
{
    // Checking if the lag value is positif
    try {
        if (lag_ <= 0)
            throw string("The lag value is incorrect, try strictly positive value.");
    }
    catch(string const& chaine)
    {
        Rcpp::Rcout << chaine << endl;
    }

    lag  = lag_;

    for (const auto val:tS_)
        tS.push_back (val);

    // Initialisation
    Struct::CVDouble one, pBeta, predicted, tStat (tS), errorCorrection, trend;

    Struct::CMatDouble B, target, Var;

    // Differentiation
    for (auto it = tStat.end () -1 ; it != tStat.begin () ; --it)
        *it -= *(it-1) ;
    tStat.erase (tStat.begin ());

    nl = tStat.size ();

    for (unsigned i = lag  ; i < nl ; i++) {
        one.push_back (1);
        errorCorrection.push_back (tS[i]);
        trend.push_back (i);
    }

    B.push_back(one);
    B.push_back(errorCorrection);
    B.push_back(trend);

    P_Part (tStat, target, B, lag);
    regression (B , target[0], pBeta);
    MultCVDouble (Trans(B), B, Var);
    Inverse (Var, Var);
    MultCVDouble (B, pBeta, predicted);

    double Stdv_Value = 0;
    double SSE =0;
    int nl = B[0].size ();

    for (unsigned i = 0  ; i < B[0].size () ; i++)
        Stdv_Value += pow (predicted[i] - target[0][i], 2);

    for (unsigned i = 0  ; i < target[0].size () ; i++)
        SSE += pow (predicted[i] - target[0].Mean(), 2);

    Stdv_Value = sqrt(Stdv_Value * Var[1][1] / (B[0].size () - B.size ()));
    SBC = nl * log(SSE / nl) + (lag + 3) * log (nl);
    df = pBeta[1] / Stdv_Value;
}

/* Constructor for the Augmented Dickey Fuller test that takes an C++ matrix as argument */
DickeyFuller::DickeyFuller (const Struct::CVDouble & tS_, int lag_ ) //throw ()
{
    // Checking if the lag value is positif
    if (lag_ <= 0)
        throw Exception("The lag value is incorrect, try strictly positive value.");

    lag  = lag_;

    // Initialisation
    Struct::CVDouble one, pBeta, predicted, tStat (tS_), errorCorrection, trend;
    Struct::CMatDouble B, target, Var;

    // 0 is the default argument attributed to the lag value length(tS)^1/3
    if (lag == 0)
        lag = pow (tStat.size(), 1/3);

    /* Differentiation */
    for (auto it = tStat.end () -1 ; it != tStat.begin () ; --it)
        *it -= *(it-1) ;
    tStat.erase (tStat.begin ());

    nl = tStat.size ();
    for (unsigned i = lag  ; i < nl ; i++) {
        one.push_back (1);
        errorCorrection.push_back (tS_[i]);
        trend.push_back (i);
    }

    B.push_back(one);
    B.push_back(errorCorrection);
    B.push_back(trend);

    P_Part (tStat, target, B, lag);
    regression (B , target[0], pBeta);
    MultCVDouble (Trans(B), B, Var);
    Inverse (Var, Var);
    MultCVDouble (B, pBeta, predicted);

    double Stdv_Value = 0;
    double SSE =0;
    int nl = B[0].size ();

    for (unsigned i = 0  ; i < B[0].size () ; i++)
        Stdv_Value += pow (predicted[i] - target[0][i], 2);

    for (unsigned i = 0  ; i < target[0].size () ; i++)
        SSE += pow (predicted[i] - target[0].Mean(), 2);

    Stdv_Value = sqrt(Stdv_Value * Var[1][1] / (B[0].size () - B.size ()));

    SBC = nl * log(SSE / nl) + (lag+3) * log(nl);
    //SBC = Stdv_Value;
    df = pBeta[1] / Stdv_Value;
}
/***************************************************/
void DickeyFuller::summary ()
{
    Rcpp::Rcout <<  "------------------------------------------------\n";
    Rcpp::Rcout <<  "         The Augmented Dickey-Fuller test       \n";
    Rcpp::Rcout <<  "------------------------------------------------\n";
    Rcpp::Rcout <<  "The lag parameter: p = "<< lag << "\n";
    Rcpp::Rcout <<  "Critical values: 1% \t 5% \n";
    Rcpp::Rcout <<  "                "<< get1CriticalValue () << "\t" << get5CriticalValue ()<< "\n";
    Rcpp::Rcout <<  "The statistic of the test is: "<< df << "\n";
    Rcpp::Rcout <<  "------------------------------------------------\n";
}

/***************************************************/
int order (const Struct::CVDouble & tS, int p) {
    int res = 0;
    Struct::CVDouble serie (tS);
    DickeyFuller test (serie,p);

    while (test.get5CriticalValue() <= test.getDF())
    {
        res += 1;
        for (auto it = serie.end () -1 ; it != serie.begin () ; --it)
            *it -= *(it-1) ;

        serie.erase (serie.begin ());
        test = DickeyFuller (serie,p);
    }
    return res;
}
