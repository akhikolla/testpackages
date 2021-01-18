#ifndef IO_H
#define IO_H
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include "string.h"
#include <sstream>
#include "Eigen/SVD"
#include <RcppEigen.h>
#include "conversion.h"

using namespace std;
using namespace Rcpp;
using namespace Eigen;

class IO{
 public:
  int instantiated;
  int seed;
  string datafile;
  string prefix;
  int nItEM;
  int nItMC;
  int nBurn;
  int dp;
  
  bool sparse;
  int n;
  int p;
  int g;
  int nsample;
  int maxit;
  double tol;

  string analysis;
  string family;
  string algorithm;
  
  MatrixXd  x;
  MatrixXd  U;
  MatrixXd  V;
  VectorXd su;

  VectorXd  y;
  VectorXd  v;
  
  // Tools
  double   sy;
  VectorXd xTy;
  VectorXd xTx;
  VectorXd  sx;
 
  VectorXi Z0;
  double lambda;
  
  // Model parameters
  bool IsModelInitialized;
  VectorXd b;        // Centers
  VectorXd pi;       // Proportion
  double sigma2;     // Residual Variance
  double gamma2;     // Intra-Component Variance
  double intercept;  // Intercept
  double likelihood, entropy;
  MatrixXd P;        // Assignment probabilities
  MatrixXd theta;
  MatrixXi Zw;
  MatrixXd Bw;
  
  // Many different constructors
  // will be considered
  IO(Rcpp::S4 & obj);
  ~IO(){};
  void updateObj(Rcpp::S4 & obj);
};
#endif
