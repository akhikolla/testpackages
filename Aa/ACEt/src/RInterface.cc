#include "mcmc.h"

#include <R.h> // R functions
// #include "Rmath.h" // Rmath

void CppGibbs_mcmc_atctet(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, int * num_col_e, 
	double *ph_m, double *ph_d, double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d, double *G_a, double *G_c, double *G_e, double *ei_a, double *ei_c, double *ei_e,
	double *var_b_a, double *var_b_c, double *var_b_e, double *beta_a, double *beta_c, double *beta_e, int *D_a, int *D_c, int *D_e, int *iter_n, int *burn, double *sd_mcmc)
{

ci_mh_atctet(result,num_p_mz, num_p_dz, num_col_a, num_col_c, num_col_e, ph_m, ph_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, G_a, G_c, G_e, ei_a, ei_c, ei_e, var_b_a, var_b_c, var_b_e, beta_a, beta_c, beta_e, D_a, D_c, D_e, iter_n, burn, sd_mcmc);

}

void CppGibbs_mcmc_atctet_2(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, int * num_col_e, 
  double *ph_m, double *ph_d, double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d,
	double *var_b_a, double *var_b_c, double *var_b_e, double *D_a, double *D_c, double *D_e, int *iter_n, int *burn, double *sd_mcmc)
{

ci_mh_atctet_2(result,num_p_mz, num_p_dz, num_col_a, num_col_c, num_col_e, ph_m, ph_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e, iter_n, burn, sd_mcmc);

}


extern "C" {
	
void CWrapper_mcmc_atctet(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, int * num_col_e, double *ph_m, double *ph_d, 
	double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d, double *G_a, double *G_c, double *G_e, double *ei_a, double *ei_c, double *ei_e,
	double *var_b_a, double *var_b_c, double *var_b_e, double *beta_a, double *beta_c, double *beta_e, int *D_a, int *D_c, int *D_e, int *iter_n, int *burn, double *sd_mcmc)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_mcmc_atctet(result,num_p_mz, num_p_dz, num_col_a, num_col_c, num_col_e, ph_m, ph_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, G_a, G_c, G_e, ei_a, ei_c, ei_e, var_b_a, var_b_c, var_b_e, beta_a, beta_c, beta_e, D_a, D_c, D_e, iter_n, burn, sd_mcmc);
}

void CWrapper_mcmc_atctet_2(double *result,int * num_p_mz, int * num_p_dz, int * num_col_a, int * num_col_c, int * num_col_e, double *ph_m, double *ph_d, 
  double *B_des_a_m, double *B_des_a_d, double *B_des_c_m, double *B_des_c_d, double *B_des_e_m, double *B_des_e_d, 
	double *var_b_a, double *var_b_c, double *var_b_e, double *D_a, double *D_c, double *D_e, int *iter_n, int *burn, double *sd_mcmc)
{
//- Invoke second function which internally can do C++ things.
//
CppGibbs_mcmc_atctet_2(result,num_p_mz, num_p_dz, num_col_a, num_col_c, num_col_e, ph_m, ph_d, B_des_a_m, B_des_a_d, B_des_c_m, B_des_c_d, B_des_e_m, B_des_e_d, var_b_a, var_b_c, var_b_e, D_a, D_c, D_e, iter_n, burn, sd_mcmc);
}

}
