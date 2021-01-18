#pragma once

#include <stdio.h>
#include <iostream>
#include <string>
#include <vector>
#include <math.h> 
#include <iostream>
#include <RcppArmadillo.h>


//bool gomp = false;

long double mu(double t, arma::mat y1, arma::mat gamma1, arma::mat fH, arma::mat f1H, double mu0H, double thetaH, arma::mat QH);
void func1(arma::mat *res, double t, arma::mat *y, arma::mat fH, arma::mat f1H, arma::mat aH, arma::mat bH, arma::mat QH, double theta);
arma::mat Q(double t, arma::mat QH, double theta);
void func1(double t, arma::mat *y, arma::mat fH, arma::mat f1H, arma::mat aH, arma::mat bH, arma::mat QH, double theta);
