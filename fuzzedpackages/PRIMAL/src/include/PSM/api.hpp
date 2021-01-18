#ifndef API_HPP
#define API_HPP

#include <RcppEigen.h>

#include "PSM/PSM.hpp"

#include <iostream>


extern "C" void QuantileRegression_api(
	int *_n,
	int *_d,
	double *_x,
	double *_y,
 	double *_tau,
 	int *max_it,
 	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
);

extern "C" void SparseSVM_api(
 	int *_n,
	int *_d,
	double *_x,
	double *_y,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list,
 	double *x0_list
);

extern "C" void Dantzig_api(
	int *_n,
	int *_d,
	double *_x,
	double *_y,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
);

extern "C" void CompressedSensing_api(
 	int *_n,
	int *_d,
	double *_x,
	double *_y,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list
);

extern "C" void ParametricSimplexMethod_api(
	int *_m,
	int *_n,
	double *_a,
	double *_b,
	double *_b_bar,
	double *_c,
	double *_c_bar,
	int *max_it,
	double *lambda_threshold,
	int *T,
	double *lambda_list,
	double *x_list,
	double *y_list,
	int *B_init
);

#endif
