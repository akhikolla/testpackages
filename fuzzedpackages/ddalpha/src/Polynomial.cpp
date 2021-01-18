/*
  File:             Polynomial.cpp
  Created by:       Oleksii Pokotylo
  First published:  07.05.2014
  Last revised:     07.05.2014

  Contains the polynomial classifier the DD-plot classification.

  For a description of the algorithm, see:
    Li, J., Cuesta-Albertos, J. A. and Liu, R. Y. (2012). DD-classifier: Nonparametric classification procedure based on
DD-plot, Journal of the American Statistical Association 107(498): 737 - 753.
*/

#include "stdafx.h"

#include <limits>
#include <boost/math/special_functions/binomial.hpp>

#include "asa047.h"

#ifndef _MSC_VER
#include <Rcpp.h>
using namespace Rcpp;
#endif

/**
Calculates the empirical risk for two classes on the basis of given depths

@param polynomial : Polynomial as a vector of coefficients starting with the first degree(a0 = 0 always)
@param points : nx2 matrix of depths, where each column contains the depths against the corresponding class
@param numClass1:  Number of points belonging to the first class
@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)

@throws empirical risk
*/
double GetEmpiricalRisk(TPoint& polynomial, TDMatrix points, unsigned numClass1, unsigned numClass2){
	unsigned degree = polynomial.size();
	unsigned n = numClass1 + numClass2;

	double risk = 0;
	int sign = 1;
	for (unsigned i = 0; i < n; i++){
		if (i >= numClass1)
			sign = -1;

		double val = points[i][0];
		double res = 0;
		for (unsigned j = 0; j<degree; j++){
			res += polynomial[j] * std::pow(val, j+1);
		}
		if ((points[i][1] - res)*sign > 0){    // for class1 depths[i,2] > res, for class 2 <
			risk++;
		}
	}

	return risk / n;
}

/**
Calculates the coefficients of the polynomial of a given degree
going through given points and the origin

@param degree : degree of the polynomial, should be equal the number of points
@param points : degreex2 points for the polynomial to path through

@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)

@throws runtime_error in case of singularity
*/
bool GetPolynomial(unsigned degree, TDMatrix points, TPoint& polynomial) {

	bMatrix A(degree, degree);
	for (unsigned i = 0; i < degree; i++){
		for (unsigned j = 0; j < degree; j++){
			A(i, j) = (std::pow(points[i][0], j + 1));
		}
	}

	bVector b(degree);
	for (unsigned i = 0; i < degree; i++){
		b[i] = points[i][1];
	}

	bPM pm(A.size1());
	bPM::size_type singular = boost::numeric::ublas::lu_factorize(A, pm);
	if (singular != 0) return false;
	boost::numeric::ublas::lu_substitute(A, pm, b);

	for (unsigned i = 0; i < degree; i++){
		if (!(b[i] < std::numeric_limits<double>::max()) || b[i] < -std::numeric_limits<double>::max()){  return false; }
		polynomial[i] = b[i];
	}
  
	return true;
}

/**
Chooses the best in classification sense polynomial among
"cardinality" randomly chosen polynomials, passing through
"degree" arbitrary points

@param points:       nx2 matrix of points where first column is an absciss,	n = numClass1 + numClass2
@param numClass1:    Number of points belonging to the first class
@param degree:       Degree of the polynomial
@param n_polynomials:  Number of randomly chosen polynomials

@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)
*/
TPoint GetRandomMinPolynomial(TDMatrix points, int numClass1, int numClass2, int degree, int n_polynomials){
  int n = numClass1 + numClass2;
	vector<int> usedIndexesX(n);
	vector<int> usedIndexesY(n);
	int nx = 0, ny = 0;

	for (int i = 0; i<n; i++){
		if (points[i][0] != 0){
			usedIndexesX[nx++] = i;
			if (points[i][1] != 0)
				usedIndexesY[ny++] = i;
		}
	}

	// 0.3 of all combination
	int numOfCombinations = boost::math::binomial_coefficient<double>(nx - 1, degree - 1) * ny * 0.3; 
	int numCandidates = (numOfCombinations > n_polynomials
		? n_polynomials
		: numOfCombinations);

	TPoint minPolynomial(degree);
	double minEmpRisk = 1;
	TDMatrix sample = new double*[degree];
	for (int i = 0; i < numCandidates; i++){
		// generate sample
		set<int> smp;
		smp.insert(usedIndexesY[random(ny)]);
		while ((int)smp.size() < degree){
			smp.insert(usedIndexesX[random(nx)]);
		}

		set <int>::const_iterator s = smp.begin();
		for (int j = 0; j < degree; j++, s++) {
			sample[j] = points[*s];
		}

		try{
			TPoint pol(degree);
			if(!GetPolynomial(degree, sample, pol))
				continue;
			double rsk = GetEmpiricalRisk(pol, points, numClass1, numClass2);
			if (rsk < minEmpRisk) {
  				minPolynomial = pol;
  				minEmpRisk = rsk;        
			}
		}
		catch (runtime_error &e){ /* singular matrix*/ }
		catch (...){ /* NA or inf */ }
	}
	delete[] sample;
	return minPolynomial;
}

static int _degree;
static TDMatrix _points;
static int _numClass1;
static int _numClass2;

/**
Calculates the empirical risk for two classes on the basis of given depths
and approximates it to get continuous derivative

@param polynomial : Polynomial as a vector of coefficients starting with the first degree(a0 = 0 always)
@param _points : nx2 matrix of depths, where each column contains the depths against the corresponding class
@param _numClass1:  Number of points belonging to the first class
@param _numClass2:  Number of points belonging to the first class
@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)
*/
double GetEmpiricalRiskSmoothed(double polynomial[]){
	const float smoothingConstant = 100;
	
	double risk = 0;
	int sign = 1;
	for (int i = 0; i < _numClass1 + _numClass2; i++){
		if (i >= _numClass1)
			sign = -1;

		double val = (_points)[i][0];
		double res = 0;
		for (int j = 0; j < _degree; j++){
			res += polynomial[j] * std::pow(val, j+1);
		}
		risk += 1 / (1 + exp(-smoothingConstant*((_points)[i][1] - res)*sign));
	}

	return risk / _numClass1 + _numClass2;
}


TPoint nlm_optimize(TDMatrix points, TPoint& minCandidate, int numClass1, int numClass2){
	/* static variables for GetEmpiricalRiskSmoothed */
	_points = points;
	_numClass1 = numClass1;
	_numClass2 = numClass2;
	_degree = minCandidate.size();

	double* start = new double[_degree];
	std::copy(minCandidate.begin(), minCandidate.end(), start);

	int icount;
	int ifault;
	int kcount;
	int konvge;
	int n = _degree;
	int numres;
	double reqmin;
	double *step;
	double *xmin;
	double ynewlo;

	step = new double[n];
	xmin = new double[n];

	reqmin = 1.0E-06;

	for (int i = 0; i < n; i++)
	{
		// determines the size and shape of the initial simplex.
		// The relative magnitudes of its elements should reflect the units of the variables.  
		step[i] = 1.0;
	}

	konvge = 10;
	kcount = 500;
	/*
	cout << "\n";
	cout << "  Starting point X:\n";
	cout << "\n";
	for (i = 0; i < n; i++)
	{
		cout << "  " << start[i] << "";
	}
	cout << "\n";
	
	ynewlo = GetEmpiricalRiskSmoothed(start);

	cout << "\n";
	cout << "  F(X) = " << ynewlo << "\n";
	*/
	nelmin(GetEmpiricalRiskSmoothed, n, start, xmin, &ynewlo, reqmin, step,
		konvge, kcount, &icount, &numres, &ifault);
	/*
	cout << "\n";
	cout << "  Return code IFAULT = " << ifault << "\n";
	cout << "\n";
	cout << "  Estimate of minimizing value X*:\n";
	cout << "\n";
	for (i = 0; i < n; i++)
	{
		cout << "  " << setw(14) << xmin[i] << "\n";
	}

	cout << "\n";
	cout << "  F(X*) = " << ynewlo << "\n";

	cout << "\n";
	cout << "  Number of iterations = " << icount << "\n";
	cout << "  Number of restarts =   " << numres << "\n";
*/
	TPoint minpol = TPoint(xmin, xmin + _degree);

	delete[] start;
	delete[] step;
	delete[] xmin;

	return minpol;
}

/**
Chooses the best in classification sense polynomial

@param points:       nx2 matrix of points where first column is an abscissa, n = numClass1 + numClass2
@param numClass1:    Number of points belonging to the first class
@param degree:       Degree of the polynomial
@param presize:		 if true - run evaluation 5 times

@return polynomial as a vector of coefficients starting with the first degree (a0 = 0 always)
*/
TPoint GetOptPolynomial(TDMatrix points, unsigned numClass1, unsigned numClass2, unsigned degree, bool presize /*default = FALSE*/){

	double minError = 100.1;
	TPoint minPol;
	
	for (int i = 0; i < (presize ? 3 : 1); i++){
		TPoint minCandidate = GetRandomMinPolynomial(points, numClass1, numClass2, degree, 10 ^ degree);
		double err = GetEmpiricalRisk(minCandidate, points, numClass1, numClass2);
		if (err < minError){
			minPol = minCandidate;
			minError = err;
		}

//#define DEBUG
#ifdef DEBUG
Rcpp::Rcout << "candminPol: ";
for (int i = 0; i< degree; i++){
  Rcpp::Rcout << minCandidate[i] << " ";
}
Rcpp::Rcout << " ; error = "<< err << " \n";
#endif

		TPoint optPolynomial = nlm_optimize(points, minCandidate, numClass1, numClass2);
		err = GetEmpiricalRisk(optPolynomial, points, numClass1, numClass2);
		if (err <= minError){
			minPol = optPolynomial;
			minError = err;
		}
  
#ifdef DEBUG
Rcpp::Rcout << "minPol: ";
for (int i = 0; i< minPol.size(); i++){
  Rcpp::Rcout << minPol[i] << " ";
}
Rcpp::Rcout << " ; error = "<< minError << " \n";
#endif		
	}  

	return(minPol);
}


/**
Calculates classification error of "degree" - degree polynomial using cross - validation approach

@param points:       nx2 matrix of points where first column is an absciss,	n = numClass1 + numClass2
@param numClass1:    Number of points belonging to the first class
@param degree:       Degree of the polynomial
@param chunkNumber:  Number of chunks in which data set should be splitted, chunkNumber should be a divider of n(n%%chunkNumber = 0)

@return Number of errors
*/
double GetCvError(TDMatrix points, int numClass1, int numClass2, int degree, int chunkNumber){

	int n = numClass1 + numClass2;
	int chunkSize = ceil((double)n / chunkNumber);

	TDMatrix learnpoints = new double*[n - chunkSize+1]; 
	TDMatrix checkpoints = new double*[chunkSize];

	int chunk = 0;
	int n1 = 0; // number of Class1 points in checkpoints
	for (int j = 0, l = 0, c = 0; j < n; j++){
		if (j%chunkNumber)
			learnpoints[l++] = points[j];
		else
			checkpoints[c++] = points[j];
		if (j < numClass1 && (j%chunkNumber == 0)) n1++;
	}

	double err = 0;
	bool bigch = true;
	for (; chunk < chunkNumber; chunk++){

		if (chunk > 0){
			if (bigch && (chunkNumber)*(chunkSize - 1) + chunk == n){
				bigch = false;
				chunkSize--;
				//checkpoints.resize(chunkSize);
				//learnpoints.resize(n - chunkSize);
				learnpoints[n - chunkSize - 1] = points[n - 1];
			}

			for (int i = 0; i < chunkSize; i++){
				checkpoints[i] = learnpoints[(chunkNumber - 1)*i + chunk - 1];
				learnpoints[(chunkNumber - 1)*i + chunk - 1] = points[chunkNumber*i + chunk - 1];
				if (chunkNumber*i + chunk == numClass1)
					n1--;
			}
		}

		TPoint minPolynomial = GetOptPolynomial(learnpoints, numClass1 - n1, numClass2 - chunkSize + n1, degree, false);
		double curErr = GetEmpiricalRisk(minPolynomial, checkpoints, n1, chunkSize - n1);
		err += curErr;//  chunkSize;
	}

	delete[] learnpoints;
	delete[] checkpoints;
	return err/n;
}

TPoint PolynomialLearnCV(TDMatrix input, int numClass1, int numClass2, int maxDegree, int chunkNumber, int *degree, int *axis){
  int numPoints = numClass1 + numClass2;

  int polOptDegree = 0;
	double polOptError = numPoints;
	int polOptAxis = 0;

	TDMatrix input2 = newM(numPoints, 2); // copy
	for (int i = 0, tmp; i < numPoints; i++){ input2[i][0] = input[i][1]; input2[i][1] = input[i][0]; } // swap columns

	for (int degree = 1; degree <= maxDegree; degree++){
		double polError = GetCvError(input, numClass1, numClass2, degree, chunkNumber);
		//cout << degree << " " << polError << "\n";
		if (polError < polOptError){
			polOptAxis = 0;
			polOptDegree = degree;
			polOptError = polError;
		}
		polError = GetCvError(input2, numClass1, numClass2, degree, chunkNumber);
		//cout << degree << " " << polError << "\n";
		if (polError < polOptError){
			polOptAxis = 1;
			polOptDegree = degree;
			polOptError = polError;
		}
	}

	//cout << polOptDegree << " " << polOptError << "\n";
	TPoint polynomial = polOptAxis == 0
		? GetOptPolynomial(input, numClass1, numClass2, polOptDegree, true)
		: GetOptPolynomial(input2, numClass1, numClass2, polOptDegree, true);

	deleteM(input2);

	*axis = polOptAxis;
	*degree = polOptDegree;
	return polynomial;
}
