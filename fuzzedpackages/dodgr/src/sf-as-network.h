#include <string>
#include <cmath>

#include <Rcpp.h>
// [[Rcpp::depends(RcppParallel)]]
#include <RcppParallel.h>

const double earth = 6378137.0; // WSG84 definition

constexpr float INFINITE_FLOAT =  std::numeric_limits<float>::max ();
constexpr double INFINITE_DOUBLE =  std::numeric_limits<double>::max ();
constexpr int INFINITE_INT =  std::numeric_limits<int>::max ();

namespace sf{

double haversine (double x1, double y1, double x2, double y2);
void fill_one_row (const unsigned int ngeoms, const Rcpp::NumericMatrix &gi,
        const Rcpp::CharacterVector &rnms,
        const double &d, const double &hw_factor,
        const std::string &hway, const bool &has_names,
        const std::vector <std::string> &way_names,
        const size_t &grownum, const size_t &rownum, const bool &reverse,
        Rcpp::NumericMatrix &nmat, Rcpp::CharacterMatrix &idmat);

}

Rcpp::List rcpp_sf_as_network (const Rcpp::List &sf_lines,
        const Rcpp::DataFrame &pr);
Rcpp::IntegerVector rcpp_points_index (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
Rcpp::IntegerVector rcpp_points_index_par (const Rcpp::DataFrame &xy,
        Rcpp::DataFrame &pts);
