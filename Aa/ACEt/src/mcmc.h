#include "stdlib.h"
#include "math.h"

#include <fstream>
#include <iostream>
#include <boost/random.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/gamma_distribution.hpp>
#include <boost/random/bernoulli_distribution.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/uniform_real.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/generator_iterator.hpp>

/*
void ci_mh(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, 
	double *ph_m, double *ph_d, double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *G_a, double *G_c,
	double *var, double *var_b_a, double *var_b_c, int *D_a, int *D_c, int *iter_n, int *burn, double *sd_mcmc);

void ci_mh_atet(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_e, 
	double *ph_m, double *ph_d, double *B_des_a_m, double *B_des_a_d, double *B_des_e_m, double *B_des_e_d, double *G_a, double *G_e,
	double *var_b_a, double *var_b_e, int *D_a, int *D_e, int *iter_n, int *burn, double *sd_mcmc);
*/
void ci_mh_atctet(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, int * num_col_e, 
	double *ph_m, double *ph_d, double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d, double *G_a, double *G_c, double *G_e, double *ei_a, double *ei_c, double *ei_e,
	double *var_b_a, double *var_b_c, double *var_b_e, double *beta_a, double *beta_c, double *beta_e, int *D_a, int *D_c, int *D_e, int *iter_n, int *burn, double *sd_mcmc);

void ci_mh_atctet_2(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, int * num_col_e, 
  double *ph_m, double *ph_d, double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d,
	double *var_b_a, double *var_b_c, double *var_b_e, double *D_a, double *D_c, double *D_e, int *iter_n, int *burn, double *sd_mcmc);
