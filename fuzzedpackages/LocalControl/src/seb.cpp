
#include <Rcpp.h>
#include "seb/Seb.h"
using namespace Rcpp;

/*
 Rcpp wrapper function to call the miniball library found at:

    https://github.com/hbf/miniball

 Takes a dataframe containing covariates, returns a double
 which is the maximum distance between any two points.

*/

// [[Rcpp::export]]
double getMaxDist(DataFrame& df) {
  typedef double FT;
  typedef Seb::Point<FT> Point;
  typedef std::vector<Point> PointVector;
  typedef Seb::Smallest_enclosing_ball<FT> Miniball;

  using std::cout;
  using std::endl;
  using std::vector;

  int n = df.nrows();
  int d = df.size();

  //std::vector<NumericVector> clustercols; #bugged, I guess?
  std::vector< std::vector<double> > clustercols;
  clustercols.reserve(d);
  for(int j = 0; j < d; j++){
    clustercols.push_back(df[j]);
  }

  PointVector S;
  S.reserve(n);
  vector<double> coords(d);

  for (int i=0; i<n; ++i) {

    for (int j=0; j<d; ++j) {
      coords[j] = clustercols[j][i];
    }

    S.push_back(Point(d,coords.begin()));
  }

  Miniball mb(d, S);

  // Output
  FT rad = mb.radius();

  return rad*2;
}
