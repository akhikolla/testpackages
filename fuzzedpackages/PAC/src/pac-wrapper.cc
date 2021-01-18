#include <Rcpp.h>
#include <string>
#include "bsp.h"

// [[Rcpp::plugins(cpp11)]]
void printInfo(const int numpts, const int dim, std::string &method, const int maxlevel) {
	std::string method_name = "";
	if (method == "dsp") {
		method_name = "Discrepancy based partition";
	} else if (method == "ll") {
		method_name = "BSP with limited-lookahead";
	} else {
		method_name = "unknown method";
	}
	Rcpp::Rcout << "Input Data: " <<  numpts << " by " << dim << std::endl;
	Rcpp::Rcout << "Partition method: " <<  method_name << std::endl;
	Rcpp::Rcout << "Maximum level: " <<  maxlevel << std::endl;
}

void transformInput(Rcpp::NumericMatrix data, matrix &x, vector<double> &mmax, vector<double> &mmin) {
	int n = data.nrow();
	int dim = data.ncol();
	x.resize(n);
	mmax.resize(dim);
	mmin.resize(dim);
	for(int j = 0; j < dim; j++) {
		mmax[j] = data(0,j);
		mmin[j] = data(0,j);
	}
	for(int i = 0; i < n; i++) {
		x[i].resize(dim);
		for(int j = 0; j < dim; j++) {
			x[i][j] = data(i,j);
			if(data(i,j) > mmax[j])	 mmax[j] = data(i,j);
			if(data(i,j) < mmin[j])	 mmin[j] = data(i,j);
		}
	}
}

// [[Rcpp::export]]
Rcpp::NumericMatrix BSPLeaveCenterCpp(Rcpp::NumericMatrix data, SEXP numleaves, SEXP methodr){
	/* fill x with input data from R */
	matrix x;
	vector<double> mmax;
	vector<double> mmin;
	transformInput(data, x, mmax, mmin);
	std::string method = Rcpp::as<std::string>(methodr);
	int maxlevel = Rcpp::as<int>(numleaves);
    printInfo(x.size(), x[0].size(), method, maxlevel);

	bspTree T(x, mmax, mmin);
	if (method == "dsp") {
		T.dsp(NUM_CUT, maxlevel, THETA);
	} else if (method == "ll") {
		T.lla(maxlevel, MIN_PTS);
	} else {
		Rcpp::NumericMatrix dummy(0);
		Rcpp::Rcout << "unknown method" << std::endl;
		return dummy;
	}

	Rcpp::Rcout << "partition completed" << std::endl;
 	T.CalculateLeafCenter();

 	/* copy results back to R */
 	Rcpp::NumericMatrix out(T.leafctr.size(), T.leafctr[0].size());
 	
 	for(uint i = 0; i < T.leafctr.size(); i++) {
 		for(int j = 0; j < T.leafctr[0].size(); j++) {
 			out(i,j) = T.leafctr[i][j];
 		}
 	}
 	return out;
}
