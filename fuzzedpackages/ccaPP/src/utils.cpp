/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#include <R.h>
#include "utils.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// C++ isnan() is not portable and gives error on Windows systems
// use R macro ISNAN() instead


// -------------
// MCD-estimator
// -------------

// MCD-estimator
mat covMCD(const mat& x) {
  // call R function from package robustbase
  Environment robustbase("package:robustbase");
	Function covMcd = robustbase["covMcd"];
  NumericMatrix Rcpp_x = wrap(x);             // does this reuse memory?
	List Rcpp_MCD = covMcd(Rcpp_x);             // call R function
  NumericMatrix Rcpp_Sigma = Rcpp_MCD["cov"]; // extract covariance matrix
	mat Sigma(Rcpp_Sigma.begin(), Rcpp_Sigma.nrow(), Rcpp_Sigma.ncol(), false); // reuse memory
	return Sigma;
}


// --------------
// median and MAD
// --------------

// L1 median
vec l1Median(const mat& x) {
	// call R function from package pcaPP
	// TODO: call underlying C++ code directly
	Environment pcaPP("package:pcaPP");
	Function l1median = pcaPP["l1median"];
	NumericMatrix Rcpp_x = wrap(x);					// does this reuse memory?
	NumericVector Rcpp_center = l1median(Rcpp_x);	// call R function
	vec center(Rcpp_center.begin(), Rcpp_center.size(), false);	// reuse memory
	return center;
}

// R interface to l1Median() (for testing)
SEXP R_l1Median(SEXP R_x) {
	NumericMatrix Rcpp_x(R_x);	// convert data to Rcpp type
	mat x(Rcpp_x.begin(), Rcpp_x.nrow(), Rcpp_x.ncol(), false);	// reuse memory
	vec center = l1Median(x);	// call arma version
	return wrap(center.memptr(), center.memptr() + center.n_elem);
}

// median using std::vector
// order of observations is messed up
double median(vector<double>& x) {
	int n = x.size();
	// find median
	int half = (n + 1) / 2;	// divide by integer to get truncated result
	half--;					// indices start with 0
	double med;
	if((n % 2) == 1) {
		// odd number of observations, take the middle one
		nth_element(x.begin(), x.begin()+half, x.end());
		med = x[half];
	} else {
		// even number of observations, take the mean of the two middle ones
		nth_element(x.begin(), x.begin()+half, x.end());
//		double tmp = x[half];
//		nth_element(x.begin(), x.begin()+half+1, x.end());
//  	med = 0.5 * (tmp + x[half+1]);
    med = 0.5 * (x[half] + *min_element(x.begin()+half+1, x.end()));
	}
	return med;
}

// median using arma::vec
double median(const vec& x) {
	uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i))) return(NA_REAL);
	}
	// copy data to std::vector to not mess up the order of observations
	vector<double> xx(n);
	for(uword i = 0; i < n; i++) {
		xx[i] = x(i);
	}
	return median(xx);
}

// R interface to median()
SEXP R_fastMedian(SEXP R_x) {
	NumericVector Rcpp_x(R_x);						// convert data to Rcpp type
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	return wrap(median(x));	// call arma version and wrap result
}

// MAD
// median is returned through the corresponding parameter
double mad(const vec& x, const double& constant, double& center) {
	uword n = x.n_elem;
	// make sure it returns NA in the presence of NA's
	for(uword i = 0; i < n; i++) {
		if(ISNAN(x(i))) return(NA_REAL);
	}
	// copy data to std::vector to not mess up the order of observations
	vector<double> xx(n);
	for(uword i = 0; i < n; i++) {
		xx[i] = x(i);
	}
	// find median
	center = median(xx);
	// compute MAD
	for(uword i = 0; i < n; i++) {
		xx[i] = abs(xx[i] - center);
	}
	return constant * median(xx);
}
double mad(const vec& x, double& center) {
	return mad(x, 1.4826, center);
}

// R interfaces to mad()
SEXP R_fastMAD(SEXP R_x, SEXP R_constant) {
	NumericVector Rcpp_x(R_x);						// convert data to Rcpp type
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	double constant = as<double>(R_constant), center;
	double MAD = mad(x, constant, center);			// call arma version
	return List::create(
			Named("center") = wrap(center),
			Named("MAD") = wrap(MAD)
			);
}


// ---------------
// order and ranks
// ---------------

// class definition
class SortData {
public:
	uword index;
	double value;

	// constructors
	SortData();
	SortData(uword&, const double&);
//	// overloaded < (is less) operator for sorting and ordering
//	bool operator< (const SortData&);
//	// overloaded > (is greater) operator for sorting and ordering
//	bool operator> (const SortData&);
};

// constructors
inline SortData::SortData() {}
inline SortData::SortData(uword& first, const double& second) {
	index = first;
	value = second;
}

//// overloaded < (is less) operator for sorting and ordering
//bool SortData::operator< (const SortData& other) {
//      return (this->value < other.value);
//}

//// overloaded > (is greater) operator for sorting and ordering
//bool SortData::operator> (const SortData& other) {
//      return (this->value > other.value);
//}

// compare two objects with < (is less) operator for sorting and ordering
bool sortDataIsLess(const SortData& left, const SortData& right) {
	return left.value < right.value;
}

// compare two objects with > (is greater) operator for sorting and ordering
bool sortDataIsGreater(const SortData& left, const SortData& right) {
	return left.value > right.value;
}

// compute order of observations
// stable sorting is not necessary since ties are broken by averaging
uvec order(const vec& x, const bool& decreasing) {
	// initialize data structure for sorting
	const uword n = x.n_elem;
	vector<SortData> foo(n);
	for(uword i = 0; i < n; i++) {
		foo[i] = SortData(i, x(i));
	}
	// call STL's sort()
	if(decreasing) {
		sort(foo.begin(), foo.end(), sortDataIsGreater);
	} else {
		sort(foo.begin(), foo.end(), sortDataIsLess);
	}
	// construct and return vector of indices
	uvec indices(n);
	for(uword i = 0; i < n; i++) {
		SortData bar = foo[i];
		indices(i) = bar.index;
	}
	return indices;
}

// compute increasing order of observations
uvec order(const vec& x) {
	return order(x, false);
}

// compute ranks of observations in a vector
// function name 'rank' causes error with clang++ on OS X Mavericks
vec rank_ccaPP(const vec& x) {
	const uword n = x.n_elem;
	uword i, j, k;
	// compute order of observations
	uvec ord = order(x);
	// compute ranks (break ties by taking average)
	vec ranks(n);
	for(i = 0; i < n; i = j+1) {
		j = i;
		// check if there is a series of equal values
		while((j < n-1) && (x(ord(j)) == x(ord(j+1)))) {
			j++;
		}
		// break ties by average rank, otherwise this gives just the rank
		for(k = i; k <= j; k++) {
			ranks(ord(k)) = (i + j + 2) / 2.0;
		}
	}
	// return ranks
	return ranks;
}

// R interface to rank_ccaPP() (for testing)
SEXP R_rank(SEXP R_x) {
  NumericVector Rcpp_x(R_x);                    // convert data to Rcpp type
  vec x(Rcpp_x.begin(), Rcpp_x.size(), false);  // convert data to arma type
  vec ranks = rank_ccaPP(x);                    // call arma version
  return wrap(ranks.memptr(), ranks.memptr() + ranks.n_elem);
}
