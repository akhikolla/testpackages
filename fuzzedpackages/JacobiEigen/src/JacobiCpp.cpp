#include <Rcpp.h>
using namespace Rcpp;

NumericMatrix Ident(int n) // not exported.
{
    NumericMatrix I(n, n);
    for(int i = 0; i < n; i++) I(i, i) = 1.0;
    return I;
}

// [[Rcpp::export]]
List JacobiCpp(NumericMatrix x, bool only_values = false, double eps = 0.0)
{
    NumericMatrix S(clone(x));
    int nr = S.nrow();
    bool vectors = !only_values;
    NumericMatrix H;

    if(vectors) {
      H = Ident(nr);
    }

    double eps0 =  as<double>((as<List>(Environment::base_env()[".Machine"]))["double.eps"]);
    double tol = eps > eps0 ? eps : eps0;  // i.e. no lower than .Machine$double.eps
    if(only_values & (eps == 0.0)) tol = sqrt(tol); // a lower accuracy is adequate here.

    while(true) {
	    double maxS = 0.0;
	    int i=0, j=0;
	    for(int row = 1; row < nr; row++) {  // find value & position of maximum |off-diagonal|
	        for(int col = 0; col < row; col++) {
	        	double val = fabs(S(row, col));
		        if(maxS < val) {
		          maxS = val;
		          i = row;
		          j = col;
	        	}
	       }
	    }
	    if(maxS <= tol) break;

	    NumericVector Si = S(_, i), Sj = S(_, j);

	    double theta = 0.5*atan2(2.0*Si(j), Sj(j) - Si(i));
	    double s = sin(theta), c = cos(theta);

	    S(i, _) = S(_, i) = c*Si - s*Sj;
	    S(j, _) = S(_, j) = s*Si + c*Sj;
	    S(i, j) = S(j, i) = 0.0;
	    S(i, i) = c*c*Si(i) - 2.0*s*c*Si(j) + s*s*Sj(j);
	    S(j, j) = s*s*Si(i) + 2.0*s*c*Si(j) + c*c*Sj(j);

      if(vectors) {
	        NumericVector Hi = H(_, i);
	        H(_, i) = c*Hi - s*H(_, j);
	        H(_, j) = s*Hi + c*H(_, j);
      }
    }
    if(vectors) {
      return List::create(_["values"] = diag(S),
                          _["vectors"] = H);
    } else {
      return List::create(_["values"] = diag(S),
                          _["vectors"] = R_NilValue);
    }
}

// [[Rcpp::export]]
List JacobiSCpp(NumericMatrix x, bool only_values = false, double eps = 0.0)
{
  NumericMatrix S(clone(x));
  int nr = S.nrow();
  bool vectors = !only_values;
  NumericMatrix H;
  
  if(vectors) {
    H = Ident(nr);
  }
  
  double eps0 =  as<double>((as<List>(Environment::base_env()[".Machine"]))["double.eps"]);
  double tol = eps > eps0 ? eps : eps0;  // i.e. no lower than .Machine$double.eps
  if(only_values & (eps == 0.0)) tol = sqrt(tol); // a lower accuracy is adequate here.

  double maxS = 0.0;
  for(int row = 1; row < nr; row++) {  // find value of maximum |off-diagonal|
    for(int col = 0; col < row; col++) {
      double val = fabs(S(row, col));
      maxS = maxS < val ? val : maxS;
    }
  }
  
  while(true) {
    if(maxS <= tol) break;
    double maxS0 = 0.0;
    for(int i = 1; i < nr; i++) {
      for(int j = 0; j < i; j++) {
        double val = fabs(S(i, j));
        maxS0 = maxS0 < val ? val : maxS0;
        if (val > maxS/2.0) {
          NumericVector Si = S(_, i), Sj = S(_, j);
          
          double theta = 0.5*atan2(2.0*Si(j), Sj(j) - Si(i));
          double s = sin(theta), c = cos(theta);
          
          S(i, _) = S(_, i) = c*Si - s*Sj;
          S(j, _) = S(_, j) = s*Si + c*Sj;
          S(i, j) = S(j, i) = 0.0;
          S(i, i) = c*c*Si(i) - 2.0*s*c*Si(j) + s*s*Sj(j);
          S(j, j) = s*s*Si(i) + 2.0*s*c*Si(j) + c*c*Sj(j);
          
          if(vectors) {
            NumericVector Hi = H(_, i);
            H(_, i) = c*Hi - s*H(_, j);
            H(_, j) = s*Hi + c*H(_, j);
          }
        }
      }
    }
    maxS = maxS0;
  }
  if(vectors) {
    return List::create(_["values"] = diag(S),
                        _["vectors"] = H);
  } else {
    return List::create(_["values"] = diag(S),
                        _["vectors"] = R_NilValue);
  }
}
