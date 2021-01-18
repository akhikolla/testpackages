/*
 * This file is part of NlinTS package.
 *.
 * See the COPYRIGHT file at the top-level directory of this distribution
 * for details of code ownership.
 *
 */

#include "../inst/include/entropyExport.h"
#include "../inst/include/nsEntropy.h"
#include "../inst/include/exception.h"

using namespace Rcpp;

/************************************************************/
double entropy_disc (IntegerVector  & I, std::string log)
{

  if (I.size() == 0)
            throw string ("Error: the data are empty.");

	std::vector<int>  X;
	for (const auto val : I)
        	X.push_back (val);

    double e = nsEntropy::entropy (X, log);
    return e;
}

/************************************************************/
double joinEntropy_disc (Rcpp::DataFrame & Df, std::string log)
{
    vector<vector<int> > X = Rcpp::as < vector<vector<int> > > (Df);

    if (X.size() == 0)
            throw string ("Error: the data are empty.");

    double je = nsEntropy::joinEntropy (X, log);
    return je;
}

/************************************************************/
double mutualInformation_disc_u (IntegerVector  & I,  IntegerVector  & J, std::string log, bool normalize)
{
    if (I.size() != J.size())
              throw string ("Error: The variables have not the same length.");

    if (I.size() == 0)
              throw string ("Error: the data are empty.");

    std::vector<int>  X, Y;

  	for (const auto val : I)
          	X.push_back (val);

      for (const auto val : J)
          	Y.push_back (val);

      double mi = nsEntropy::mutualInformation (X, Y, log, normalize);
      return mi;
}

/************************************************************/
double mutualInformation_disc (Rcpp::DataFrame & Df, std::string log, bool normalize)
{
    vector<vector<int> > X = Rcpp::as < vector<vector<int> > > (Df);

    if (X.size() == 0)
            throw string ("Error: the data are empty.");

    double mi = nsEntropy::mutualInformation (X, log, normalize);
    return mi;
}

/************************************************************/
// Transfer of information from J to I
double transferEntropy_disc (Rcpp::IntegerVector & I, Rcpp::IntegerVector & J, int p, int q, std::string log, bool normalize)
{

    try {
          if (p <= 0 or q <= 0)
              throw string ("Error: The lag value is incorrect, try strictly positive values.");
      }
      catch(string const& chaine)
      {
          Rcpp::Rcout << chaine << endl;
      }

    if (I.size() != J.size())
            throw string ("Error: The variables have not the same length.");

    if (I.size() == 0)
            throw string ("Error: the data are empty.");

    std::vector<int>  X, Y;

    for (const auto val : I)
            X.push_back (val);

    for (const auto val : J)
            Y.push_back (val);

    double te = nsEntropy::transferEntropy (X, Y, p, q, log, normalize);
    return te;
}


//-------------------------------------------------------------//
double entropy_cont ( Rcpp::NumericVector & I, int k, std::string log)
{
    if (I.size() == 0)
            throw string ("Error: the data are empty.");

    std::vector<double>  X;
    for (const auto val : I)
            X.push_back (val);

    double e = nsEntropy::entropy (X, k, log);
    return e;

}

//-------------------------------------------------------------//
double mutualInformation_cont ( Rcpp::DataFrame & Df, int k, std::string alg, bool normalize)
{
    vector<vector<double> > X = Rcpp::as < vector<vector<double> > > (Df);

    if (X.size() == 0)
            throw string ("Error: the data are empty.");

    double mi = nsEntropy::mutualInformation (X, k, alg, normalize);
    return mi;
}

//-------------------------------------------------------------//
// Transfer of information from J to I
double transferEntropy_cont ( NumericVector & I,  NumericVector & J, int p, int q, int k, bool normalize)
{

    try {
          if (p <= 0 or q <= 0)
              throw string ("Error: The lag value is incorrect, try strictly positive values.");
      }
      catch(string const& chaine)
      {
          Rcpp::Rcout << chaine << endl;
      }

    if (I.size() != J.size())
            throw string ("Error: The variables have not the same length.");

    if (I.size() == 0)
            throw string ("Error: the data are empty.");


    std::vector<double>  X, Y;

    for (const auto val : I)
            X.push_back (val);

    for (const auto val : J)
            Y.push_back (val);

    double te = nsEntropy::transferEntropy (X, Y, p, q, k, normalize);
    return te;
}
