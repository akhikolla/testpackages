/*
 * Author: Andreas Alfons
 *         Erasmus University Rotterdam
 */

#include <R.h>
#include "cca.h"

using namespace Rcpp;
using namespace arma;
using namespace std;


// ******************************************
// control class definitions for correlations
// ******************************************

// The control classes that handle how the correlations are computed.  They
// store the values for the additional control parameters and have a cor()
// method to call the corresponding correlation function.

// TODO: use Rcpp function containsElementNamed() in constructors
//       to check whether the list contains an element with a certain name
//       (this is not available yet in the CRAN release Rcpp_0.9.10)

// -------------------------------------
// control class for Pearson correlation
// -------------------------------------

class CorPearsonControl {
public:
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// method to compute correlation
double CorPearsonControl::cor(const vec& x, const vec& y) {
	return corPearson(x, y);
}


// ---------------------------------------------
// control classes for nonparametric correlation
// ---------------------------------------------

// abstract class to share code for Spearman, Kendall and quadrant correlation

class CorNPControl {
public:
	bool consistent;
};


// --------------------------------------
// control class for Spearman correlation
// --------------------------------------

class CorSpearmanControl: public CorNPControl {
public:
	// constructors
	CorSpearmanControl();
	CorSpearmanControl(List&);
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// constructors
inline CorSpearmanControl::CorSpearmanControl() {
	consistent = false;
}
inline CorSpearmanControl::CorSpearmanControl(List& control) {
	consistent = as<bool>(control["consistent"]);
}

// method to compute correlation
double CorSpearmanControl::cor(const vec& x, const vec& y) {
	return corSpearman(x, y, consistent);
}


// -------------------------------------
// control class for Kendall correlation
// -------------------------------------

class CorKendallControl: public CorNPControl {
public:
	// constructors
	CorKendallControl();
	CorKendallControl(List&);
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// constructors
inline CorKendallControl::CorKendallControl() {
	consistent = false;
}
inline CorKendallControl::CorKendallControl(List& control) {
	consistent = as<bool>(control["consistent"]);
}

// method to compute correlation
double CorKendallControl::cor(const vec& x, const vec& y) {
	return corKendall(x, y, consistent);
}


// --------------------------------------
// control class for quadrant correlation
// --------------------------------------

class CorQuadrantControl: public CorNPControl {
public:
	// constructors
	CorQuadrantControl();
	CorQuadrantControl(List&);
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// constructors
inline CorQuadrantControl::CorQuadrantControl() {
	consistent = false;
}
inline CorQuadrantControl::CorQuadrantControl(List& control) {
	consistent = as<bool>(control["consistent"]);
}

// method to compute correlation
double CorQuadrantControl::cor(const vec& x, const vec& y) {
	return corQuadrant(x, y, consistent);
}


// -----------------------------
// control class for M-estimator
// -----------------------------

class CorMControl {
public:
	double prob;
	string initial;
	double tol;
	// constructors
	CorMControl();
	CorMControl(List&);
	// method to compute correlation
	double cor(const vec&, const vec&);
};

// constructors
inline CorMControl::CorMControl() {
	prob = 0.9;
	initial = "quadrant";
	tol = 0.000001;
}
inline CorMControl::CorMControl(List& control) {
	prob = as<double>(control["prob"]);
	initial = as<string>(control["initial"]);
	tol = as<double>(control["tol"]);
}

// method to compute correlation
double CorMControl::cor(const vec& x, const vec& y) {
	return corM(x, y, prob, initial, tol);
}


// -------------------------------
// correlation with control object
// -------------------------------

// for testing purposes

template <class CorControl>
double fastCor(const vec& x, const vec& y, CorControl control) {
	return control.cor(x, y);
}

// R interface
SEXP R_fastCor(SEXP R_x, SEXP R_y, SEXP R_method, SEXP R_control) {
	// convert data
	NumericVector Rcpp_x(R_x), Rcpp_y(R_y);
	vec x(Rcpp_x.begin(), Rcpp_x.size(), false);	// reuse memory
	vec y(Rcpp_y.begin(), Rcpp_y.size(), false);	// reuse memory
	string method = as<string>(R_method);	// convert character string
	List Rcpp_control(R_control);
	// define control object for the correlations and call arma version
	double r;
	if(method == "spearman") {
		CorSpearmanControl corControl(Rcpp_control);
		r = fastCor(x, y, corControl);
	} else if(method == "kendall") {
		CorKendallControl corControl(Rcpp_control);
		r = fastCor(x, y, corControl);
	} else if(method == "quadrant") {
		CorQuadrantControl corControl(Rcpp_control);
		r = fastCor(x, y, corControl);
	} else if(method == "M") {
		CorMControl corControl(Rcpp_control);
		r = fastCor(x, y, corControl);
	} else if(method == "pearson") {
		CorPearsonControl corControl;
		r = fastCor(x, y, corControl);
	} else {
		error("method not available");
	}
	// wrap and return result
	return wrap(r);
}


// ***********************************************************
// control class definitions for projection pursuit algorithms
// ***********************************************************

// The control classes that handle how the maximum correlations are computed at
// each step of CCA.  They store the values for the additional control
// parameters for the respective algorithm and have a maxCor() method to
// execute the algorithm.

// The data are assumed to be standardized for computing the maximum
// correlation.  Standardization and backtransformation of the canonical
// vectors is done in the workhorse function for canonical correlation
// analysis to share code between the projection pursuit algorithms.

// -----------------------------------------
// control class for alternate grid searches
// -----------------------------------------

class GridControl {
public:
	uword nIterations;	// number of iterations of alternate grid searches
	uword nAlternate;	// number of alternate grid searches in each iteration
	uword nGrid;		// number of grid points to be used for grid searches
	uvec selectX;		// x-variables to be used for finding order of y-variables
	uvec selectY;		// y-variables to be used for finding order of x-variables
	double tol;			// numeric tolerance for convergence
	// constructors
	GridControl();
	GridControl(List&);
	// get the order in which to update the elements of a weighting vector
	template <class CorControl>
	void findOrder(const mat&, const vec&, CorControl, uvec&, double&, vec&);
	// get the order in which to update the elements of the weighting vectors
	template <class CorControl>
	void findOrder(const mat&, const mat&, CorControl, uvec&, uvec&, double&,
			vec&, vec&, bool&);
	// get equispaced grid of angles in a plane
	vec getGrid(const uword&);
	// grid search to update one weighting vector
	template <class CorControl>
	void gridSearch(const mat&, const uvec&, const vec&, CorControl, vec&,
			double&, vec&);
	// set counter of consecutive iterations without improvement
	void setCounter(uword&, const double&, const double&);
	// maximum correlation between multivariate data sets x and y
	template <class CorControl>
	double maxCor(const mat&, const mat&, CorControl, vec&, vec&);
};

// constructors
inline GridControl::GridControl() {
	nIterations = 10;
	nAlternate = 10;
	nGrid = 25;
	tol = 0.000001;
}
inline GridControl::GridControl(List& control) {
	nIterations = as<uword>(control["nIterations"]);
	nAlternate = as<uword>(control["nAlternate"]);
	nGrid = as<uword>(control["nGrid"]);
	IntegerVector Rcpp_selectX = control["selectX"];
	int p = Rcpp_selectX.size();
	selectX.set_size(p);
	for(int j = 0; j < p; j++) {
		selectX(j) = Rcpp_selectX[j];
	}
	IntegerVector Rcpp_selectY = control["selectY"];
	int q = Rcpp_selectY.size();
	selectY.set_size(q);
	for(int j = 0; j < q; j++) {
		selectY(j) = Rcpp_selectY[j];
	}
	tol = as<double>(control["tol"]);
}

// get the order in which to update the elements of a weighting vector
// x ............ data matrix for which to get the order of updating the
//                weighting vector
// y ............ linear combination of the other data matrix according to the
//                other weighting vector, which is kept fixed
// corControl ... control object to compute correlation
// orderX ....... order of updating the weighting vector to be computed
// maxCor ....... maximum correlation to get initial value
// a ............ weighting vector to get initial value
template <class CorControl>
void GridControl::findOrder(const mat& x, const vec& y, CorControl corControl,
		uvec& orderX, double& maxCor, vec& a) {
	// compute columnwise absolute correlations of x with y
	const uword p = x.n_cols;
	vec corY(p);
	for(uword j = 0; j < p; j++) {
		corY(j) = abs(corControl.cor(x.unsafe_col(j), y));
	}
	// order columns of x according to absolute correlations with y
	orderX = order(corY, true);
	// store maximum correlation and set weight of corresponding variable to 1
	uword first = orderX(0);
	maxCor = corY(first);
	a(first) = 1;
}

// get the order in which to update the elements of the weighting vectors
// x ............ first data matrix for which to get the order of updating the
//                weighting vector
// y ............ second data matrix for which to get the order of updating the
//                weighting vector
// corControl ... control object to compute correlation
// orderX ....... order of updating the first weighting vector to be computed
// orderY ....... order of updating the first weighting vector to be computed
// maxCor ....... maximum correlation to get initial value
// a ............ first weighting vector to get initial value
// b ............ second weighting vector to get initial value
// startWithX ... logical to be computed that indicates whether to start with
//                the first data set in the alternate grid searches
template <class CorControl>
void GridControl::findOrder(const mat& x, const mat& y, CorControl corControl,
		uvec& orderX, uvec& orderY, double& maxCor, vec& a, vec& b,
		bool& startWithX) {
	const uword p = x.n_cols, q = y.n_cols;
	const uword pp = selectX.n_elem, qq = selectY.n_elem;
	mat corMat;
	vec avgCorX, avgCorY;
	if((pp > 0) && (qq > 0)) {
		// selected x- and y-variables for faster computation
		// sort indices
		selectX = sort(selectX);
		selectY = sort(selectY);
		// compute the absolute correlations of all x-variables with the
		// selected y-variables
		mat corMatX(p, qq);
		for(uword j = 0; j < qq; j++) {
			vec yy = y.unsafe_col(selectY(j));
			for(uword i = 0; i < p; i++) {
				corMatX(i, j) = abs(corControl.cor(x.unsafe_col(i), yy));
			}
		}
		// compute the absolute correlations of all y-variables with the
		// selected x-variables with
		// avoid recomputing already computed absolute correlations
		uword indexY = 0, nextY = selectY(0);	// already computed
		mat corMatY(q, pp);
		for(uword i = 0; i < q; i++) {
			if(i == nextY) {
				// use already computed absolute correlations
				for(uword j = 0; j < pp; j++) {
					corMatY(i, j) = corMatX(selectX(j), indexY);
				}
				indexY++; nextY = selectY(indexY);
			} else {
				// compute absolute correlations with selected x-variables
				vec yy = y.unsafe_col(i);
				for(uword j = 0; j < pp; j++) {
					vec xx = x.unsafe_col(selectX(j));
					corMatY(i, j) = abs(corControl.cor(yy, xx));
				}
			}
		}
		// compute average absolute correlations
		avgCorX = mean(corMatX, 1); avgCorY = mean(corMatY, 1);
	} else {
		// compute complete matrix of columnwise absolute correlations
		corMat.set_size(p, q);
		for(uword i = 0; i < p; i++) {
			vec xx = x.unsafe_col(i);
			for(uword j = 0; j < q; j++) {
				corMat(i, j) = abs(corControl.cor(xx, y.unsafe_col(j)));
			}
		}
		// compute average absolute correlations
		avgCorX = mean(corMat, 1); avgCorY = mean(corMat, 0).t();
	}
	// determine the order of the variables and the starting values for the
	// weighting vectors from average correlations
	orderX = order(avgCorX, true), orderY = order(avgCorY, true);
	// set weights of corresponding variables to 1 and determine with which
	// data set to start in alternate grid searches
	uword firstX = orderX(0),  firstY = orderY(0);
	a(firstX) = 1;
	b(firstY) = 1;
	startWithX = (avgCorX(firstX) >= avgCorY(firstY));
	// store maximum correlation
	if((pp > 0) && (qq > 0)) {
		// compute absolute correlation between variables with highest averages
		// (not all absolute correlations are precomuted)
		vec xx = x.unsafe_col(firstX), yy = y.unsafe_col(firstY);
		maxCor = abs(corControl.cor(xx, yy));
	} else {
		// take maximum correlation from full absolute correlation matrix
		maxCor = corMat(firstX, firstY);
	}
}

// get equispaced grid of angles in a plane
// in each iteration of alternate grid searches, the interval over which to
// search is cut in half
// i ....... iteration of grid search
vec GridControl::getGrid(const uword& i) {
	const int j = int(i);	// prevents error on OS X
	vec grid(nGrid);
	grid(0) = - M_PI / pow(2.0, j);	// lower end point of grid
	double step = M_PI / (nGrid * pow(2.0, j-1));	// step size
	for(uword k = 1; k < nGrid; k++) {
		grid(k) = grid(k-1) + step;	// iteratively add step size to define grid
	}
	// return equispaced grid of angles
	return grid;
}

// grid search to update one weighting vector
// x ............ data matrix for which to update the weighting vector
// y ............ linear combination of the other data matrix according to the
//                other weighting vector, which is kept fixed
// corControl ... control object to compute correlation
// grid ........  grid points to be used for grid search
// maxCor ....... maximum correlation to be updated
// a ............ weighting vector to be updated
template <class CorControl>
void GridControl::gridSearch(const mat& x, const uvec& orderX, const vec& y,
		CorControl corControl, vec& grid, double& maxCor, vec& a) {
	// initializations
	const uword p = x.n_cols, nGrid = grid.n_elem;
	// perform grid searches for each canonical basis vector
	for(uword j = 0; j < p; j++) {
		// define current canonical basis vector according to order of columns
		vec ej = zeros<vec>(p);
		ej(orderX(j)) = 1;
		// perform grid search for the current canonical basis vector
		vec corY(nGrid);
		for(uword k = 0; k < nGrid; k++) {
			double angle = grid(k);
			vec currentA = cos(angle) * a + sin(angle) * ej;
			corY(k) = abs(corControl.cor(x * currentA, y));
		}
		// find grid point that maximizes correlation functional and keep
		// maximum correlation of current grid search
		uword whichMax;
		double currentMaxCor = corY.max(whichMax);
		// update maximum correlation and weighting vector
		// if 0 degree angle is not part of the grid, the maximum correlation
		// of the current grid search may be smaller than the previous maximum
		if(currentMaxCor > maxCor) {
			maxCor = currentMaxCor;
			double optAngle = grid(whichMax);
			a = cos(optAngle) * a + sin(optAngle) * ej;
		}
	}
}

// set counter of consecutive iterations without improvement
// counter .... counter to be updated
// current .... value from current iteration
// previous ... value from previous iteration
void GridControl::setCounter(uword& counter, const double& current,
		const double& previous) {
	if((current - previous) > tol) {
		counter = 0;
	} else {
		counter++;
	}
}

// maximum correlation between multivariate data sets x and y based on
// alternate grid searches
// x ............. first data matrix
// y ............. second data matrix
// corControl .... control object to compute correlation
// a ............. first weighting vector to be updated
// b ............. second weighting vector to be updated
template <class CorControl>
double GridControl::maxCor(const mat& x, const mat& y, CorControl corControl,
		vec& a, vec& b) {
	// initializations
	uword p = x.n_cols, q = y.n_cols;
	// perform alternate grid searches if both data sets are multivariate
	// if one data set is univariate, alternate grid searches are not necessary
	double maxCor;	// initialize maximum correlation
	if((p == 1) && (q == 1)) {
		// both data sets are univariate
		a.ones(p); b.ones(q);
		vec xx = x.unsafe_col(0), yy = y.unsafe_col(0);	// reuse memory
		maxCor = abs(corControl.cor(xx, yy));			// compute correlation
	} else {
		double previousMaxCor = R_NegInf;
		if((p > 1) && (q == 1)) {
			// x is multivariate, y is univariate
			vec yy = y.unsafe_col(0);							// reuse memory
			uvec orderX(p);
			a.zeros(p); b.ones(q);
			findOrder(x, yy, corControl, orderX, maxCor, a);	// column order
			// stop if there are two consecutive iterations without improvement
			uword i = 0, convCounter = 0;
			while((i < nIterations) && (convCounter < 2)) {
				previousMaxCor = maxCor;
				vec grid = getGrid(i+1);		// define vector of grid points
				gridSearch(x, orderX, yy, corControl, grid, maxCor, a);	// grid search
				i++;
				setCounter(convCounter, maxCor, previousMaxCor);
			}
		} else if((p == 1) && (q > 1)) {
			// x is univariate, y is multivariate
			vec xx = x.unsafe_col(0);							// reuse memory
			uvec orderY(q);
			a.ones(p); b.zeros(q);
			findOrder(y, xx, corControl, orderY, maxCor, b);	// column order
			// stop if there are two consecutive iterations without improvement
			uword i = 0, convCounter = 0;
			while((i < nIterations) && (convCounter < 2)) {
				previousMaxCor = maxCor;
				vec grid = getGrid(i+1);		// define vector of grid points
				gridSearch(y, orderY, xx, corControl, grid, maxCor, b);	// grid search
				i++;
				setCounter(convCounter, maxCor, previousMaxCor);
			}
		} else if((p > 1) && (q > 1)) {
			// both data sets are multivariate
			// compute correlations between variables in x and y
			uvec orderX(p), orderY(q);
			a.zeros(p); b.zeros(q);
			bool startWithX;
			findOrder(x, y, corControl, orderX, orderY,
					maxCor, a, b, startWithX);	// column orders
			// perform alternate grid searches
			if(startWithX) {
				// start with grid search for x
				// stop if there are two consecutive iterations without improvement
				uword i = 0, convCounter = 0;
				while((i < nIterations) && (convCounter < 2)) {
					previousMaxCor = maxCor;
					vec grid = getGrid(i+1);	// define vector of grid points
					uword j = 0;
					double altMaxCor = R_NegInf;
					while((j < nAlternate) && ((maxCor - altMaxCor) > tol)) {
						altMaxCor = maxCor;
						// maximize correlation functional over a keeping b fixed
						vec yb = y * b;		// linear combination of columns of y
						gridSearch(x, orderX, yb, corControl, grid, maxCor, a);	// grid search
						// maximize correlation functional over b keeping a fixed
						vec xa = x * a;		// linear combination of columns of x
						gridSearch(y, orderY, xa, corControl, grid, maxCor, b);	// grid search
						j++;
					}
					i++;
					setCounter(convCounter, maxCor, previousMaxCor);
				}
			} else {
				// start with grid search for y
				// stop if there are two consecutive iterations without improvement
				uword i = 0, convCounter = 0;
				while((i < nIterations) && (convCounter < 2)) {
					previousMaxCor = maxCor;
					vec grid = getGrid(i+1);	// define vector of grid points
					uword j = 0;
					double altMaxCor = R_NegInf;
					while((j < nAlternate) && ((maxCor - altMaxCor) > tol)) {
						altMaxCor = maxCor;
						// maximize correlation functional over b keeping a fixed
						vec xa = x * a;		// linear combination of columns of x
						gridSearch(y, orderY, xa, corControl, grid, maxCor, b);	// grid search
						// maximize correlation functional over a keeping b fixed
						vec yb = y * b;		// linear combination of columns of y
						gridSearch(x, orderX, yb, corControl, grid, maxCor, a);	// grid search
						j++;
					}
					i++;
					setCounter(convCounter, maxCor, previousMaxCor);
				}
			}
		} else {
			return NA_REAL;	// should never happen
		}
	}
	// ensure that norm of weighting vectors is 1
	a = a / norm(a, 2);
	b = b / norm(b, 2);
	// check direction
	double r = corControl.cor(x * a, y * b);
	if(r < 0) {
		if((p > 1) && (q == 1)) {
			a = -a;
		} else {
			b = -b;
		}
	}
	// return maximum correlation
	return maxCor;
}


// -------------------------------------------------
// control class for projections through data points
// -------------------------------------------------

class ProjControl {
public:
	bool useL1Median;
	// constructors
	ProjControl();
	ProjControl(List&);
	// get matrix of directions through data points
	mat getDirections(const mat&);
	// maximum correlation between multivariate data sets x and y
	template <class CorControl>
	double maxCor(const mat&, const mat&, CorControl, vec&, vec&);
};

// constructors
inline ProjControl::ProjControl() {
	useL1Median = true;
}
inline ProjControl::ProjControl(List& control) {
	useL1Median = as<bool>(control["useL1Median"]);
}

// get matrix of directions through data points
mat ProjControl::getDirections(const mat& x) {
	// initializations
	uword n = x.n_rows, p = x.n_cols;
	mat A(p, n);
	if(useL1Median) {
		// fill columns of A with centered and normalized rows of x
		vec center = l1Median(x);
		for(uword i = 0; i < n; i++) {
			vec xi = x.row(i).t() - center;
			A.col(i) = xi / norm(xi, 2);
		}
	} else {
		// fill columns of A with normalized rows of x
		// (the data are assumed to be standardized)
		for(uword i = 0; i < n; i++) {
			vec xi = x.row(i).t();
			A.col(i) = xi / norm(xi, 2);
		}
	}
	return A;
}

// maximum correlation between multivariate data sets x and y based on
// projections through the data points
// x ............ first data matrix
// y ............ second data matrix
// corControl ... control object to compute correlation
// a ............ first weighting vector to be updated
// b ............ second weighting vector to be updated
template <class CorControl>
double ProjControl::maxCor(const mat& x, const mat& y, CorControl corControl,
		vec& a, vec& b) {
	// initializations
	uword n = x.n_rows, p = x.n_cols, q = y.n_cols;
	// explore projections through the data points
	double maxCor = R_NegInf;	// initialize maximum correlation
	mat A, B;
	if(p > 1) {
		A = getDirections(x);	// directions through data points of x
	} else {
		a.ones(p);
	}
	if(q > 1) {
		B = getDirections(y);	// directions through data points of y
	} else {
		b.ones(q);
	}
	// if one data set is univariate, projections through those data points
	// don't need to be explored
	if((p == 1) && (q == 1)) {
		// both data sets are univariate
		vec xx = x.unsafe_col(0), yy = y.unsafe_col(0);	// reuse memory
		maxCor = abs(corControl.cor(xx, yy));			// compute correlation
	} else if((p > 1) && (q == 1)) {
		// x is multivariate, y is univariate
		vec yy = y.unsafe_col(0);				// reuse memory
		double corY;
		uword whichMax = 0;
		for(uword i = 0; i < n; i++) {
			vec xa = x * A.unsafe_col(i);		// projection in current direction
			corY = abs(corControl.cor(xa, yy));	// absolute correlation with y
			if(corY > maxCor) {
				// update maximum correlation
				maxCor = corY;
				whichMax = i;
			}
		}
		// update directions corresponding to maximum correlation
		a = A.col(whichMax);
	} else if((p == 1) && (q > 1)) {
		// x is univariate, y is multivariate
		vec xx = x.unsafe_col(0);				// reuse memory
		double corX;
		uword whichMax = 0;
		for(uword i = 0; i < n; i++) {
			vec yb = y * B.unsafe_col(i);		// projection in current direction
			corX = abs(corControl.cor(xx, yb));	// absolute correlation with x
			if(corX > maxCor) {
				// update maximum correlation
				maxCor = corX;
				whichMax = i;
			}
		}
		// update directions corresponding to maximum correlation
		b = B.col(whichMax);
	} else if((p > 1) && (q > 1)) {
		// both data sets are multivariate
		// scan all n^2 possible combinations of directions
		uword whichMaxX = 0, whichMaxY = 0;
		double corXY;
		for(uword i = 0; i < n; i++) {
			vec xa = x * A.unsafe_col(i);				// projection in x space
			for(uword j = 0; j < n; j++) {
				vec yb = y * B.unsafe_col(j);			// projection in y space
				corXY = abs(corControl.cor(xa, yb));	// absolute correlation
				if(corXY > maxCor) {
					// update maximum correlation
					maxCor = corXY;
					whichMaxX = i;
					whichMaxY = j;
				}
			}
		}
		// update directions corresponding to maximum correlation
		a = A.col(whichMaxX);
		b = B.col(whichMaxY);
	} else {
		return NA_REAL;	// should never happen
	}
	// check direction
	double r = corControl.cor(x * a, y * b);
	if(r < 0) {
		b = -b;
	}
	// return maximum correlation
	return maxCor;
}


// *****************************************************
// canonical correlation analysis via projection pursuit
// *****************************************************

bool isDummy(const vec& x) {
	// initializations
	const uword n = x.n_elem;
	uword i = 0;
	bool dummy = true;
	// loop over vector elements to check whether they are either 0 or 1
	while(dummy && i < n) {
		dummy = dummy && ((x[i] == 0.0) || (x[i] == 1.0));
		i++;
	}
	// return whether the variable is a dummy variable
	return dummy;
}

// standardize data using median/MAD or mean/SD
// only scale is needed for backtransformation of canonical vectors
// x .......... data matrix
// robust ..... should the data be robustly standardized?
// fallback ... should the fallback mode for robust standardization be used?
// scale ...... scale estimates of the variables to be computed
mat standardize(const mat& x, const bool& robust,
		const bool& fallback, vec& center, vec& scale) {
	const uword n = x.n_rows, p = x.n_cols;
	mat xs(n, p);
  center.set_size(p);
  scale.set_size(p);
	if(robust) {
		// median and MAD
		for(uword j = 0; j < p; j++) {
			vec xj = x.unsafe_col(j);
      double med;
      scale(j) = mad(xj, med);  // compute median and MAD
      center(j) = med;
			if((scale(j) == 0.0) && (fallback || isDummy(xj))) {
				// compute mean and standard deviation
				center(j) = mean(xj);
				scale(j) = norm(xj - center(j), 2) / sqrt((double)(n-1));
			}
			if(scale(j) == 0.0) error("zero scale");
			xs.col(j) = (xj - center(j)) / scale(j);	// standardize variable
		}
	} else {
		// mean and standard deviation
		for(uword j = 0; j < p; j++) {
			// with unsafe_col(), the original data would be changed when
			// sweeping out the mean
			vec xj = x.col(j);
			center(j) = mean(xj);						// compute mean
			xj -= center(j);									// sweep out mean
			scale(j) = norm(xj, 2) / sqrt((double)(n-1));	// compute SD
			if(scale(j) == 0.0) error("zero scale");
			xs.col(j) = xj / scale(j);						// sweep out SD
		}
	}
	return xs;
}

// transform canonical vectors back to the original scale
// a ....... canonical vector
// scale ... scale estimates of the corresponding original variables
void backtransform(vec& a, const vec& scale) {
	a /= scale;			// divide by scale of corresponding variable
	a /= norm(a, 2);	// divide by norm
}

// compute rotation matrix for Householder transformation
// a ... canonical vector
mat householder(const vec& a) {
	const uword p = a.n_elem;
	vec e1 = zeros<vec>(p); e1(0) = 1;		// first basis vector
	vec n1 = e1 - a; n1 = n1 / norm(n1, 2);	// unit normal vector
	mat P = eye<mat>(p, p) - 2 * n1 * n1.t();
	return P;
}

// canonical correlation analysis via projection pursuit
// x ............ first data matrix
// y ............ second data matrix
// k ............ number of canonical variables to compute
// corControl ... control object to compute correlation
// ppControl .... control object for algorithm
// standard ..... should the data be standardized?
// robust ....... should robust standardization be used?
// fallback ..... should the fallback mode for robust standardization be used?
// A ............ matrix of canonical vectors for first matrix to be updated
// B ............ matrix of canonical vectors for second matrix to be updated
template <class CorControl, class PPControl>
vec ccaPP(const mat& x, const mat& y, const uword& k, CorControl corControl,
		PPControl ppControl, const bool& standard, const bool& robust, 
    const bool& fallback, mat& A, mat& B, vec& centerX, vec& centerY, 
    vec& scaleX, vec& scaleY) {
	// initializations
	uword p = x.n_cols, q = y.n_cols;
	A.set_size(p, k); B.set_size(q, k);
	vec r(k), a, b;
	// standardize the data if requested
  mat xs, ys;
  if(standard) {
    xs = standardize(x, robust, fallback, centerX, scaleX);
    ys = standardize(y, robust, fallback, centerY, scaleY);
    // compute first canonical correlation variables with standardized data
    r(0) = ppControl.maxCor(xs, ys, corControl, a, b);
  } else {
    centerX = zeros<vec>(p); centerY = zeros<vec>(q);
    scaleX = ones<vec>(p); scaleY = ones<vec>(q);
    // compute first canonical correlation variables with original data
    r(0) = ppControl.maxCor(x, y, corControl, a, b);
  }
	A.col(0) = a; B.col(0) = b;
  // compute higher order canonical correlations
	if(k > 1) {
    // data to be reduced in each step
    mat xl, yl;
    if(standard) {
      xl = xs;
      yl = ys;
    } else {
      xl = x;
      yl = y;
    }
    // compute covariance matrices
    mat SigmaX, SigmaY;
    if(robust) {
      SigmaX = covMCD(xl);
      SigmaY = covMCD(yl);
    } else {
      SigmaX = cov(xl);
      SigmaY = cov(yl);
    }
    // compute spectral decompositions
    vec eigValX, eigValY;
    mat eigVecX, eigVecY;
    eig_sym(eigValX, eigVecX, SigmaX); eig_sym(eigValY, eigVecY, SigmaY);
    eigValX = sqrt(eigValX); eigValY = sqrt(eigValY);
    // orthogonalize the data
    xl = xl * eigVecX; 
    for(uword j = 0; j < p; j++) xl.col(j) /= eigValX(j);
    yl = yl * eigVecY;
    for(uword j = 0; j < q; j++) yl.col(j) /= eigValY(j);
    // transform first canonical vectors accordingly and divide by norm
    a = eigValX % (eigVecX.t() * a); a /= norm(a, 2);
    b = eigValY % (eigVecY.t() * b); b /= norm(b, 2);
    // transform data to orthogonal subspaces and compute higher order 
    // canonical correlations
    mat P, Q;   // for backtransformation
    for(uword l = 1; l < k; l++) {
      // perform Householder transformation
      mat Pl = householder(a), Ql = householder(b);
      xl = xl * Pl; xl.shed_col(0);	// reduced x data
      yl = yl * Ql; yl.shed_col(0);	// reduced y data
      // compute canonical correlation and canonical vectors for reduced data
      r(l) = ppControl.maxCor(xl, yl, corControl, a, b);
      // transform canonical vectors back to original space
      if(l == 1) {
        P = Pl; Q = Ql;
      } else {
        // expand current Householder matrix for x and premultiply with
        // product of previous ones
        Pl.insert_rows(0, zeros<mat>(l-1, p-l+1));
        Pl.insert_cols(0, eye<mat>(p, l-1));
        P = P * Pl;
        // expand current Householder matrix for y and premultiply with
        // product of previous ones
        Ql.insert_rows(0, zeros<mat>(l-1, q-l+1));
        Ql.insert_cols(0, eye<mat>(q, l-1));
        Q = Q * Ql;
      }
      // expand canonical vectors and premultiply with product of
      // corresponding Householder matrices
      vec al = eigVecX * (P * join_cols(zeros<vec>(l), a) / eigValX);
      A.col(l) = al / norm(al, 2);
      vec bl = eigVecY * (Q * join_cols(zeros<vec>(l), b) / eigValY); 
      B.col(l) = bl / norm(bl, 2);
    }
  }
  // back-transform canonical vectors in case of standardization
  if(standard) {
    for(uword l = 0; l < k; l++) {
      vec al = A.unsafe_col(l), bl = B.unsafe_col(l);
      backtransform(al, scaleX); backtransform(bl, scaleY);
    }
  }
  // return canonical correlations
  return r;
}

// R interface
SEXP R_ccaPP(SEXP R_x, SEXP R_y, SEXP R_k, SEXP R_method, SEXP R_corControl,
		SEXP R_algorithm, SEXP R_ppControl, SEXP R_standardize, SEXP R_fallback) {
	// initializations
	NumericMatrix Rcpp_x(R_x), Rcpp_y(R_y);	// convert data to Rcpp types
	mat x(Rcpp_x.begin(), Rcpp_x.nrow(), Rcpp_x.ncol(), false);	// convert data
	mat y(Rcpp_y.begin(), Rcpp_y.nrow(), Rcpp_y.ncol(), false);	// to arma types
	uword k = as<uword>(R_k);
	string method = as<string>(R_method);       // convert character string
	List Rcpp_corControl(R_corControl);         // list of control parameters
	string algorithm = as<string>(R_algorithm); // convert character string
	List Rcpp_ppControl(R_ppControl);           // list of control parameters
  bool standard = as<bool>(R_standardize);    // convert to boolean
  bool fallback = as<bool>(R_fallback);       // convert to boolean
	// initialize results
	vec r, centerX, centerY, scaleX, scaleY;
	mat A, B;
	if(algorithm == "grid") {
		// define control object for alternate grid searches
		GridControl ppControl(Rcpp_ppControl);
		// define control object for the correlations and call the arma version
		if(method == "spearman") {
			CorSpearmanControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "kendall") {
			CorKendallControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "quadrant") {
			CorQuadrantControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "M") {
			CorMControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "pearson") {
			CorPearsonControl corControl;
      r = ccaPP(x, y, k, corControl, ppControl, standard, false, false, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else {
			error("method not available");
		}
	} else if(algorithm == "proj") {
		// define control object for projections through data points
		ProjControl ppControl(Rcpp_ppControl);
		// define control object for the correlations and call the arma version
		if(method == "spearman") {
			CorSpearmanControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "kendall") {
			CorKendallControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "quadrant") {
			CorQuadrantControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "M") {
			CorMControl corControl(Rcpp_corControl);
      r = ccaPP(x, y, k, corControl, ppControl, standard, true, fallback, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else if(method == "pearson") {
			CorPearsonControl corControl;
      r = ccaPP(x, y, k, corControl, ppControl, standard, false, false, 
          A, B, centerX, centerY, scaleX, scaleY);
		} else {
			error("method not available");
		}
	} else {
		error("algorithm not available");
	}
  // wrap and return result
  return List::create(
    Named("cor") = r,
    Named("A") = A,
    Named("B") = B,
    Named("centerX") = centerX,
    Named("centerY") = centerY,
    Named("scaleX") = scaleX,
    Named("scaleY") = scaleY
    );
}
