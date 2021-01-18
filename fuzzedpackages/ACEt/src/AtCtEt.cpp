#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>

//#include <Rcpp.h>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace std;

/*============GLOBAL VARIABLES===========*/
/*=====END OF GLOBAL VARIABLES===========*/

/*============FUNCTION DEFINITIONS*/

RcppExport SEXP loglik_AtCtEt_esp_c(SEXP b_a, SEXP b_c, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d, SEXP B_e_m, SEXP B_e_d) {
    
    arma::vec ba = as<arma::vec>(b_a); // Number of columns
    arma::vec bc = as<arma::vec>(b_c);
	arma::vec be = as<arma::vec>(b_e);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
    arma::mat b_e_d = as<arma::mat>(B_e_d);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = ba.n_elem;
	int num_c = bc.n_elem;
	int num_e = be.n_elem;
	arma::mat k_m(2,2);
	arma::mat k_d(2,2);
	arma::mat k_d_a(2,2);
	// k_d_a << 1 << 0.5 << arma::endr << 0.5 << 1 << arma::endr;
	k_d_a.fill(0.5);
	k_d_a(0,0) = 1;
	k_d_a(1,1) = 1;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_d.fill(1);
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;

	arma::mat V_m(2,2);
	arma::mat inv_V_m(2,2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double temp_com_a = exp(as_scalar(b_a_m.row(start)*ba));
		if(num_a==1)
		{
			temp_com_a = exp(ba(0));
		}
		double temp_com_c = exp(as_scalar(b_c_m.row(start)*bc));
		if(num_c==1)
		{
			temp_com_c = exp(bc(0));
		}
		double temp_com_e = exp(as_scalar(b_e_m.row(start)*be));
		if(num_e==1)
		{
			temp_com_e = exp(be(0));
		}
		V_m = temp_com_e*diag+temp_com_c*k_m+temp_com_a*k_m;
		double det_m = V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1);
		//inv_V_m = inv(V_m);
		inv_V_m = V_m/det_m;
		inv_V_m(0,1) = (-1)*inv_V_m(0,1);
		inv_V_m(1,0) = inv_V_m(0,1);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat temp = p_m_r*inv_V_m*p_m_v;
		YSY_m = YSY_m + as_scalar(temp);
		D_m = D_m + log(det_m);
	}

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double temp_com_a = exp(as_scalar(b_a_d.row(start)*ba));
		if(num_a==1)
		{
			temp_com_a = exp(ba(0));
		}
		double temp_com_c = exp(as_scalar(b_c_d.row(start)*bc));
		if(num_c==1)
		{
			temp_com_c = exp(bc(0));
		}
		double temp_com_e = exp(as_scalar(b_e_d.row(start)*be));
		if(num_e==1)
		{
			temp_com_e = exp(be(0));
		}
		arma::mat V_d = temp_com_e*diag+temp_com_c*k_d+temp_com_a*k_d_a;
		//arma::mat inv_V_d = inv(V_d);
		arma::mat inv_V_d(2,2);
		double det_d = V_d(0,0)*V_d(0,0)-V_d(0,1)*V_d(0,1);
		inv_V_d = V_d/det_d;
		inv_V_d(0,1) = (-1)*inv_V_d(0,1);
		inv_V_d(1,0) = inv_V_d(0,1);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat temp = p_d_r*inv_V_d*p_d_v;
		YSY_d = YSY_d + as_scalar(temp);
		D_d = D_d + log(det_d);
	}
	
	double res = (D_m + YSY_m + D_d + YSY_d)/2;
   
    return(Rcpp::wrap(res));
}


RcppExport SEXP gr_AtCtEt_esp_c(SEXP b_a, SEXP b_c, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d, SEXP B_e_m, SEXP B_e_d) {
    
    arma::vec ba = as<arma::vec>(b_a); // Number of columns
    arma::vec bc = as<arma::vec>(b_c);
	arma::vec be = as<arma::vec>(b_e);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
    arma::mat b_c_d = as<arma::mat>(B_c_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
    arma::mat b_e_d = as<arma::mat>(B_e_d);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = ba.n_elem;
	int num_c = bc.n_elem;
	int num_e = be.n_elem;
	arma::mat k_m(2,2);
	arma::mat k_d(2,2);
	arma::mat k_d_a(2,2);
	// k_d_a << 1 << 0.5 << arma::endr << 0.5 << 1 << arma::endr;
	k_d_a.fill(0.5);
	k_d_a(0,0) = 1;
	k_d_a(1,1) = 1;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_d.fill(1);
	
	arma::vec var_bc_m = arma::zeros<arma::vec>(num_c);
	arma::vec var_ba_m = arma::zeros<arma::vec>(num_a);
	arma::vec var_be_m = arma::zeros<arma::vec>(num_e);

	arma::mat V_m(2,2);
	arma::mat inv_V_m(2,2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double temp_com_a = exp(as_scalar(b_a_m.row(start)*ba));
		if(num_a==1)
		{
			temp_com_a = exp(ba(0));
		}
		double temp_com_c = exp(as_scalar(b_c_m.row(start)*bc));
		if(num_c==1)
		{
			temp_com_c = exp(bc(0));
		}
		double temp_com_e = exp(as_scalar(b_e_m.row(start)*be));
		if(num_e==1)
		{
			temp_com_e = exp(be(0));
		}
		V_m = temp_com_e*diag+temp_com_c*k_m+temp_com_a*k_m;
		inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat inv_ph_m = inv_V_m*p_m_v;
		arma::mat ph_inv_m = p_m_r*inv_V_m;
	
		arma::mat temp_prod = inv_V_m*k_m;
		//double temp_c = exp(as_scalar(b_c_m.row(start)*bc))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
		double temp_c = arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m);
		//if(num_c>1)
		//{
			temp_c = temp_com_c*temp_c;
		//}
		for(int j = 0; j < num_c; j++)
		{
			var_bc_m(j) = var_bc_m(j) + b_c_m(start,j)*temp_c;
		}
		double temp_a = arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m);
		//if(num_a>1)
		//{
			temp_a = temp_com_a*temp_a;
		//}
		for(int j = 0; j < num_a; j++)
		{
			var_ba_m(j) = var_ba_m(j) + b_a_m(start,j)*temp_a;
		}
		double temp_e = arma::sum(inv_V_m.diag())-as_scalar(ph_inv_m*inv_ph_m);
		//if(num_e>1)
		//{
			temp_e = temp_com_e*temp_e;
		//}
		for(int j = 0; j < num_e; j++)
		{
			var_be_m(j) = var_be_m(j) + b_e_m(start,j)*temp_e;
		}
		
	}

	arma::vec var_bc_d = arma::zeros<arma::vec>(num_c);
	arma::vec var_ba_d = arma::zeros<arma::vec>(num_a);
	arma::vec var_be_d = arma::zeros<arma::vec>(num_e);

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double temp_com_a = exp(as_scalar(b_a_d.row(start)*ba));
		if(num_a==1)
		{
			temp_com_a = exp(ba(0));
		}
		double temp_com_c = exp(as_scalar(b_c_d.row(start)*bc));
		if(num_c==1)
		{
			temp_com_c = exp(bc(0));
		}
		double temp_com_e = exp(as_scalar(b_e_d.row(start)*be));
		if(num_e==1)
		{
			temp_com_e = exp(be(0));
		}
		arma::mat V_d = temp_com_e*diag+temp_com_c*k_d+temp_com_a*k_d_a;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat inv_ph_d = inv_V_d*p_d_v;
		arma::mat ph_inv_d = p_d_r*inv_V_d;
		
		arma::mat temp_prod = inv_V_d*k_d;
		//double temp_c = exp(as_scalar(b_c_d.row(start)*bc))*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d*inv_ph_d));
		double temp_c = arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d*inv_ph_d);
		//if(num_c>1)
		//{
			temp_c = temp_com_c*temp_c;
		//}
		for(int j = 0; j < num_c; j++)
		{
			var_bc_d(j) = var_bc_d(j) + b_c_d(start,j)*temp_c;
		}
		temp_prod = inv_V_d*k_d_a;
		double temp_a = arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*k_d_a*inv_ph_d);
		//if(num_a>1)
		//{
			temp_a = temp_com_a*temp_a;
		//}
		for(int j = 0; j < num_a; j++)
		{
			var_ba_d(j) = var_ba_d(j) + b_a_d(start,j)*temp_a;
		}
		double temp_e = arma::sum(inv_V_d.diag())-as_scalar(ph_inv_d*inv_ph_d);
		//if(num_e>1)
		//{
			temp_e = temp_com_e*temp_e;
		//}
		for(int j = 0; j < num_e; j++)
		{
			var_be_d(j) = var_be_d(j) + b_e_d(start,j)*temp_e;
		}
	}

	arma::vec d_var_c = (var_bc_m + var_bc_d)/2;
	arma::vec d_var_a = (var_ba_m + var_ba_d)/2;
	arma::vec d_var_e = (var_be_m + var_be_d)/2;
	
	arma::vec res(num_e+num_c+num_a);
	for(int i = 0; i < num_a; i++)
	{
		res(i) = d_var_a(i);
	}
	for(int i = 0; i < num_c; i++)
	{
		res(num_a+i) = d_var_c(i);
	}
	for(int i = 0; i < num_e; i++)
	{
		res(num_a+num_c+i) = d_var_e(i);
	}
   
    return(Rcpp::wrap(res));
}


RcppExport SEXP hessian_AtCtEt_esp_c(SEXP b_a, SEXP b_c, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_c_m, SEXP B_c_d, SEXP B_e_m, SEXP B_e_d) 
{
    
  arma::vec p_m = as<arma::vec>(pheno_m);
  arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
  arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_c_m = as<arma::mat>(B_c_m);
  arma::mat b_c_d = as<arma::mat>(B_c_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
  arma::mat b_e_d = as<arma::mat>(B_e_d);
	
	arma::vec ba = as<arma::vec>(b_a);
	arma::vec bc = as<arma::vec>(b_c);
	arma::vec be = as<arma::vec>(b_e);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
	int num_c = b_c_m.n_cols;
	int num_e = b_e_m.n_cols;
	arma::mat k_m(2,2);
	arma::mat k_1d(2,2);
	arma::mat k_2d(2,2);
	arma::mat k_3d(2,2);
	// k_1d << 0.5 << 0 << arma::endr << 0 << 0.5 << arma::endr;
	k_1d.zeros();
	k_1d(0,0) = 0.5;
	k_1d(1,1) = 0.5;
	arma::vec  v = arma::ones<arma::vec>(2);
	arma::mat diag = arma::diagmat(v);
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	/*
	 double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;
	 */
	arma::rowvec  gsd_a2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2m2c = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_dd = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_2dd2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3dd3 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d2 = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_3d3 = arma::zeros<arma::rowvec>(num_d/2);

	arma::rowvec  gsd_e2mm2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_c2mm2e = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_emme = arma::zeros<arma::rowvec>(num_m);
	arma::rowvec  gsd_adda = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_adde = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_edde = arma::zeros<arma::rowvec>(num_d);
	//arma::rowvec  gsd_a2dd2a = arma::zeros<arma::rowvec>(num_d/2); = gsd_2dd2
	arma::rowvec  gsd_a2dd2e = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_c3dd3e = arma::zeros<arma::rowvec>(num_d/2);
	//arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(num_d/2); = gsd_2d2

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double bc_m = exp(as_scalar(b_c_m.row(start)*bc));
		double be_m = exp(as_scalar(b_e_m.row(start)*be));
		arma::mat V_m = be_m*diag + ba_m*k_m + bc_m*k_m;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat temp = p_m_r*inv_V_m*p_m_v;
		//YSY_m = YSY_m + as_scalar(temp);
		//D_m = D_m + log(V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1));
		double accu_temp = arma::accu(inv_V_m*k_m*inv_V_m);
		gsd_a2m2a(i) = ba_m*ba_m*accu_temp;
		gsd_c2m2a(i) = ba_m*accu_temp*bc_m;
		gsd_c2m2c(i) = bc_m*bc_m*accu_temp;
		arma::mat temp_mm = inv_V_m*inv_V_m;
		gsd_emme(start) = be_m*be_m*temp_mm(0,0);
		gsd_emme(end) = be_m*be_m*temp_mm(1,1);
		accu_temp = arma::accu(temp_mm);
		gsd_e2mm2a(i) = ba_m*be_m*accu_temp;
		gsd_c2mm2e(i) = bc_m*be_m*accu_temp;
	}

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double bc_d = exp(as_scalar(b_c_d.row(start)*bc));
		double be_d = exp(as_scalar(b_e_d.row(start)*be));
		double ba2_d = ba_d*ba_d;
		double bac_d = ba_d*bc_d;
		double bc2_d = bc_d*bc_d;
		double bae_d = ba_d*be_d;
		double be2_d = be_d*be_d;
		double bce_d = be_d*bc_d;
		arma::mat V_d = be_d*diag + ba_d*(k_1d+k_2d) + bc_d*k_3d;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat temp = p_d_r*inv_V_d*p_d_v;
		//YSY_d = YSY_d + as_scalar(temp);
		//D_d = D_d + log(V_d(0,0)*V_d(0,0)-V_d(0,1)*V_d(0,1));
		arma::mat temp_dd = inv_V_d*inv_V_d;
		gsd_dd(start) = ba2_d*temp_dd(0,0);
		gsd_dd(end) = ba2_d*temp_dd(1,1);
		gsd_adde(start) = bae_d*temp_dd(0,0);
		gsd_adde(end) = bae_d*temp_dd(1,1);
		gsd_edde(start) = be2_d*temp_dd(0,0);
		gsd_edde(end) = be2_d*temp_dd(1,1);
		double accu_temp = arma::accu(temp_dd);
		gsd_2dd2(i) = ba2_d*0.5*accu_temp;
		gsd_3dd3(i) = bac_d*accu_temp;
		gsd_a2dd2e(i) = bae_d*0.5*accu_temp;
		gsd_c3dd3e(i) = bce_d*accu_temp;
		double temp_si2i = arma::accu(inv_V_d*k_2d*inv_V_d);
		gsd_2d2(i) = ba2_d*0.5*temp_si2i;
		gsd_3d2(i) = bac_d*temp_si2i;
		gsd_3d3(i) = bc2_d*2*temp_si2i;
	}

	arma::mat gsd_max(num_a+num_c+num_e,num_a+num_c+num_e);
	arma::mat b_a_m_h(num_m/2,num_a);
	arma::mat b_c_m_h(num_m/2,num_c);
	arma::mat b_e_m_h(num_m/2,num_e);
	arma::mat b_a_d_h(num_d/2,num_a);
	arma::mat b_c_d_h(num_d/2,num_c);
	arma::mat b_e_d_h(num_d/2,num_e);
	for(int i = 0; i < (0.5*num_m); i++)
	{
		b_a_m_h.row(i) = b_a_m.row(2*i);
		b_c_m_h.row(i) = b_c_m.row(2*i);
		b_e_m_h.row(i) = b_e_m.row(2*i);
	}
	for(int i = 0; i < (0.5*num_d); i++)
	{
		b_a_d_h.row(i) = b_a_d.row(2*i);
		b_c_d_h.row(i) = b_c_d.row(2*i);
		b_e_d_h.row(i) = b_e_d.row(2*i);
	}
	//arma::mat b_a_m_t = b_a_m.t();
	//arma::mat b_c_m_t = b_c_m.t();
	arma::mat b_e_m_t = b_e_m.t();
	arma::mat b_a_d_t = b_a_d.t();
	arma::mat b_e_d_t = b_e_d.t();
	//arma::mat b_c_d_t = b_c_d.t();
	arma::mat b_a_m_ht = b_a_m_h.t();
	arma::mat b_c_m_ht = b_c_m_h.t();
	arma::mat b_e_m_ht = b_e_m_h.t();
	arma::mat b_a_d_ht = b_a_d_h.t();
	arma::mat b_c_d_ht = b_c_d_h.t();
	arma::mat b_e_d_ht = b_e_d_h.t();
	arma::mat temp1 = b_a_m_ht;
	arma::mat temp4 = b_c_m_ht;
	arma::mat temp6 = b_c_m_ht;
	arma::mat temp8 = b_e_m_ht;
	arma::mat temp13 = b_c_m_ht;
	for(int i = 0; i < (0.5*num_m); i++)
	{
		temp1.col(i) = gsd_a2m2a(i)*b_a_m_ht.col(i);
		temp4.col(i) = gsd_c2m2a(i)*b_c_m_ht.col(i);
		temp6.col(i) = gsd_c2m2c(i)*b_c_m_ht.col(i);
		temp8.col(i) = gsd_e2mm2a(i)*b_e_m_ht.col(i);
		temp13.col(i) = gsd_c2mm2e(i)*b_c_m_ht.col(i);
	}
	arma::mat temp3 = b_a_d_ht;
	arma::rowvec tempv1 = gsd_2dd2+gsd_2d2;
	arma::mat temp5 = b_c_d_ht;
	arma::rowvec tempv2 = 0.5*gsd_3dd3+gsd_3d2;
	arma::mat temp7 = b_c_d_ht;
	arma::mat temp10 = b_e_d_ht;
	arma::mat temp14 = b_c_d_ht;
	for(int i = 0; i < (0.5*num_d); i++)
	{
		temp3.col(i) = tempv1(i)*b_a_d_ht.col(i);
		temp5.col(i) = tempv2(i)*b_c_d_ht.col(i);
		temp7.col(i) = gsd_3d3(i)*b_c_d_ht.col(i);
		temp10.col(i) = gsd_a2dd2e(i)*b_e_d_ht.col(i);
		temp14.col(i) = gsd_c3dd3e(i)*b_c_d_ht.col(i);
	}
	arma::mat temp2 = b_a_d_t;
	arma::mat temp9 = b_e_d_t;
	arma::mat temp12 = b_e_d_t;
	for(int i = 0; i < num_d; i++)
	{
		temp2.col(i) = gsd_dd(i)*b_a_d_t.col(i);
		temp9.col(i) = gsd_adde(i)*b_e_d_t.col(i);
		temp12.col(i) = gsd_edde(i)*b_e_d_t.col(i);
	}
	arma::mat temp11 = b_e_m_t;
	for(int i = 0; i < num_m; i++)
	{
		temp11.col(i) = gsd_emme(i)*b_e_m_t.col(i);
	}
	gsd_max.submat(0,0,num_a-1,num_a-1) = temp1*b_a_m_h + 0.25*temp2*b_a_d + temp3*b_a_d_h;
	gsd_max.submat(num_a,0,num_a+num_c-1,num_a-1) = temp4*b_a_m_h + temp5*b_a_d_h;
	gsd_max.submat(0,num_a,num_a-1,num_a+num_c-1) = trans(gsd_max.submat(num_a,0,num_a+num_c-1,num_a-1));
	gsd_max.submat(num_a,num_a,num_a+num_c-1,num_a+num_c-1) = temp6*b_c_m_h + temp7*b_c_d_h;
	gsd_max.submat(num_a+num_c,0,num_a+num_c+num_e-1,num_a-1) = temp8*b_a_m_h + 0.5*temp9*b_a_d + temp10*b_a_d_h;
	gsd_max.submat(0,num_a+num_c,num_a-1,num_a+num_c+num_e-1) = trans(gsd_max.submat(num_a+num_c,0,num_a+num_c+num_e-1,num_a-1));
	gsd_max.submat(num_a+num_c,num_a+num_c,num_a+num_c+num_e-1,num_a+num_c+num_e-1) = temp11*b_e_m + temp12*b_e_d;
	gsd_max.submat(num_a,num_a+num_c,num_a+num_c-1,num_a+num_c+num_e-1) = temp13*b_e_m_h + temp14*b_e_d_h;
	gsd_max.submat(num_a+num_c, num_a, num_a+num_c+num_e-1,num_a+num_c-1) = trans(gsd_max.submat(num_a,num_a+num_c,num_a+num_c-1,num_a+num_c+num_e-1));

	gsd_max = 0.5*gsd_max;
  int num_t = num_a+num_c+num_e;
	arma::vec res((num_t+1)*num_t/2);
	int k = 0;
	for(int i=0; i<num_t; i++)
	  for(int j=i; j<num_t; j++)
	  {
	    res[k] = gsd_max(i,j);
	    k++;
	  }

	return(Rcpp::wrap(res));
}
