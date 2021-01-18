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

RcppExport SEXP loglik_AtEt_epsp_c(SEXP v_b_a, SEXP v_b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP b_a, SEXP D_a, SEXP B_e_m, SEXP B_e_d, SEXP b_e, SEXP D_e) {
    
    double va = as<double>(v_b_a); // Number of columns
    double ve = as<double>(v_b_e);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
    arma::mat b_e_d = as<arma::mat>(B_e_d);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_e = as<arma::mat>(D_e);
	arma::vec ba = as<arma::vec>(b_a);
	arma::vec be = as<arma::vec>(b_e);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
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
	int penal_a = 2;
	if(d_a(0,1)==-1)
	{penal_a = 1;}
	int penal_e = 2;
	if(d_e(0,1)==-1)
	{penal_e = 1;}
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;
	arma::rowvec  gsd_a2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_e2mm2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_emme = arma::zeros<arma::rowvec>(num_m);
	arma::rowvec  gsd_adda = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_adde = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_edde = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_a2dd2a = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_a2dd2e = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(num_d/2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double be_m = exp(as_scalar(b_e_m.row(start)*be));
		arma::mat V_m = be_m*diag + ba_m*k_m;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat temp = p_m_r*inv_V_m*p_m_v;
		YSY_m = YSY_m + as_scalar(temp);
		D_m = D_m + log(V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1));
		double accu_temp = arma::accu(inv_V_m*k_m*inv_V_m);
		gsd_a2m2a(i) = ba_m*ba_m*accu_temp;
		arma::mat temp_mm = inv_V_m*inv_V_m;
		gsd_emme(start) = be_m*be_m*temp_mm(0,0);
		gsd_emme(end) = be_m*be_m*temp_mm(1,1);
		accu_temp = arma::accu(temp_mm);
		gsd_e2mm2a(i) = ba_m*be_m*accu_temp;
	}

	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double be_d = exp(as_scalar(b_e_d.row(start)*be));
		double ba2_d = ba_d*ba_d;
		double bae_d = ba_d*be_d;
		double be2_d = be_d*be_d;
		arma::mat V_d = be_d*diag + ba_d*(k_1d+k_2d);
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat temp = p_d_r*inv_V_d*p_d_v;
		YSY_d = YSY_d + as_scalar(temp);
		D_d = D_d + log(V_d(0,0)*V_d(0,0)-V_d(0,1)*V_d(0,1));
		arma::mat temp_dd = inv_V_d*inv_V_d;
		gsd_adda(start) = ba2_d*temp_dd(0,0);
		gsd_adda(end) = ba2_d*temp_dd(1,1);
		gsd_adde(start) = bae_d*temp_dd(0,0);
		gsd_adde(end) = bae_d*temp_dd(1,1);
		gsd_edde(start) = be2_d*temp_dd(0,0);
		gsd_edde(end) = be2_d*temp_dd(1,1);
		double accu_temp = arma::accu(temp_dd);
		gsd_a2dd2a(i) = ba2_d*0.5*accu_temp;
		gsd_a2dd2e(i) = bae_d*0.5*accu_temp;
		double temp_si2i = arma::accu(inv_V_d*k_2d*inv_V_d);
		gsd_2d2(i) = ba2_d*0.5*temp_si2i;
	}

	arma::mat gsd_max(num_a+num_e,num_a+num_e);
	arma::mat b_a_m_h(num_m/2,num_a);
	arma::mat b_e_m_h(num_m/2,num_e);
	arma::mat b_a_d_h(num_d/2,num_a);
	arma::mat b_e_d_h(num_d/2,num_e);
	for(int i = 0; i < (0.5*num_m); i++)
	{
		b_a_m_h.row(i) = b_a_m.row(2*i);
		b_e_m_h.row(i) = b_e_m.row(2*i);
	}
	for(int i = 0; i < (0.5*num_d); i++)
	{
		b_a_d_h.row(i) = b_a_d.row(2*i);
		b_e_d_h.row(i) = b_e_d.row(2*i);
	}
	//arma::mat b_a_m_t = b_a_m.t();
	arma::mat b_e_m_t = b_e_m.t();
	arma::mat b_a_d_t = b_a_d.t();
	arma::mat b_e_d_t = b_e_d.t();
	arma::mat b_a_m_ht = b_a_m_h.t();
	arma::mat b_e_m_ht = b_e_m_h.t();
	arma::mat b_a_d_ht = b_a_d_h.t();
	arma::mat b_e_d_ht = b_e_d_h.t();
	arma::mat temp1 = b_a_m_ht;
	arma::mat temp4 = b_e_m_ht;
	arma::mat temp7 = b_e_m_t;
	for(int i = 0; i < (0.5*num_m); i++)
	{
		temp1.col(i) = gsd_a2m2a(i)*b_a_m_ht.col(i);
		temp4.col(i) = gsd_e2mm2a(i)*b_e_m_ht.col(i);
		//temp6.col(i) = gsd_c2m2c(i)*b_c_m_ht.col(i);
	}
	arma::mat temp3 = b_a_d_ht;
	arma::rowvec tempv1 = gsd_a2dd2a+gsd_2d2;
	arma::mat temp6 = b_e_d_ht;
	//arma::rowvec tempv2 = 0.5*gsd_3dd3+gsd_3d2;
	arma::mat temp8 = b_e_d_t;
	for(int i = 0; i < (0.5*num_d); i++)
	{
		temp3.col(i) = tempv1(i)*b_a_d_ht.col(i);
		temp6.col(i) = gsd_a2dd2e(i)*b_e_d_ht.col(i);
		// temp7.col(i) = gsd_3d3(i)*b_c_d_ht.col(i);
	}
	arma::mat temp2 = b_a_d_t;
	arma::mat temp5 = b_e_d_t;
	for(int i = 0; i < num_d; i++)
	{
		temp2.col(i) = gsd_adda(i)*b_a_d_t.col(i);
		temp5.col(i) = gsd_adde(i)*b_e_d_t.col(i);
		temp8.col(i) = gsd_edde(i)*b_e_d_t.col(i);
	}
	for(int i = 0; i < num_m; i++)
	{
		temp7.col(i) = gsd_emme(i)*b_e_m_t.col(i);
	}
	gsd_max.submat(0,0,num_a-1,num_a-1) = arma::diagmat(2*d_a.diag())/va + temp1*b_a_m_h + 0.25*temp2*b_a_d + temp3*b_a_d_h;
	gsd_max.submat(num_a,0,num_a+num_e-1,num_a-1) = temp4*b_a_m_h + 0.5*temp5*b_a_d + temp6*b_a_d_h;
	gsd_max.submat(0,num_a,num_a-1,num_a+num_e-1) = trans(gsd_max.submat(num_a,0,num_a+num_e-1,num_a-1));
	gsd_max.submat(num_a,num_a,num_a+num_e-1,num_a+num_e-1) = arma::diagmat(2*d_e.diag())/ve + temp7*b_e_m + temp8*b_e_d;

	double gsd = log(fabs(det(0.5*gsd_max)));

	double res = D_m + YSY_m + D_d + YSY_d + (num_a-penal_a)*log(va) + (num_e-penal_e)*log(ve) + as_scalar(trans(ba)*d_a*ba)/va + as_scalar(trans(be)*d_e*be)/ve + gsd;

	return(Rcpp::wrap(res));
}

RcppExport SEXP loglik_AtEt_epsp_g_c(SEXP b_a, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_e_m, SEXP B_e_d, SEXP v_b_a, SEXP v_b_e, SEXP D_a, SEXP D_e) {

    arma::vec ba = as<arma::vec>(b_a); // Number of columns
    arma::vec be = as<arma::vec>(b_e);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
    arma::mat b_e_d = as<arma::mat>(B_e_d);

	double va = as<double>(v_b_a);
	double ve = as<double>(v_b_e);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_e = as<arma::mat>(D_e);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
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
	
	double YSY_m = 0;
	double YSY_d = 0;
	double D_m = 0;
	double D_d = 0;

	for(int i = 0; i < (0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double be_m = exp(as_scalar(b_e_m.row(start)*be));
		arma::mat V_m = be_m*diag+exp(as_scalar(b_a_m.row(start)*ba))*k_m;
		double det_m = V_m(0,0)*V_m(0,0)-V_m(0,1)*V_m(0,1);
		//inv_V_m = inv(V_m);
		arma::mat inv_V_m = V_m/det_m;
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

	for(int i = 0; i < (0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double be_d = exp(as_scalar(b_e_d.row(start)*be));
		arma::mat V_d = be_d*diag+exp(as_scalar(b_a_d.row(start)*ba))*(k_1d+k_2d);
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

	double res = D_m + YSY_m + D_d + YSY_d;
	if(va>0)
		res += as_scalar(trans(ba)*d_a*ba)/va;
	if(ve>0)
		res += as_scalar(trans(be)*d_e*be)/ve;

	return(Rcpp::wrap(res));
}

RcppExport SEXP gr_AtEt_epsp_g_c(SEXP b_a, SEXP b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP B_e_m, SEXP B_e_d, SEXP v_b_a, SEXP v_b_e, SEXP D_a, SEXP D_e)
{
	arma::vec ba = as<arma::vec>(b_a); // Number of columns
    arma::vec be = as<arma::vec>(b_e);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
    arma::mat b_e_d = as<arma::mat>(B_e_d);

	double va = as<double>(v_b_a);
	double ve = as<double>(v_b_e);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_e = as<arma::mat>(D_e);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
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
	
	arma::vec v_b_a_m = arma::ones<arma::vec>(num_a);
	arma::vec v_b_e_m = arma::ones<arma::vec>(num_e);
	
	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double be_m = exp(as_scalar(b_e_m.row(start)*be));
		arma::mat V_m = ba_m*k_m + be_m*diag;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat inv_ph_m = inv_V_m*p_m_v;
		arma::mat ph_inv_m = p_m_r*inv_V_m;
		arma::mat temp_prod = inv_V_m*k_m;
		double temp_a = ba_m*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_m*k_m*inv_ph_m));
		for(int j = 0; j < num_a; j++)
		{
			v_b_a_m(j) = v_b_a_m(j) + b_a_m(start,j)*temp_a;
		}
		double temp_e = be_m*(arma::sum(inv_V_m.diag())-as_scalar(ph_inv_m*inv_ph_m));
		for(int j = 0; j < num_e; j++)
		{
			v_b_e_m(j) = v_b_e_m(j) + b_e_m(start,j)*temp_e;
		}
	}

	arma::vec v_b_a_d = arma::ones<arma::vec>(num_a);
	arma::vec v_b_e_d = arma::ones<arma::vec>(num_e);
	
	for(int i = 0; i<(0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double be_d = exp(as_scalar(b_e_d.row(start)*be));
		arma::mat V_d = ba_d*(k_1d+k_2d) + be_d*diag;
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat inv_ph_d = inv_V_d*p_d_v;
		arma::mat ph_inv_d = p_d_r*inv_V_d;
		arma::mat temp_prod = inv_V_d*(k_1d+k_2d);
		double temp_a = ba_d*(arma::sum(temp_prod.diag())-as_scalar(ph_inv_d*(k_1d+k_2d)*inv_ph_d));
		for(int j = 0; j < num_a; j++)
		{
			v_b_a_d(j) = v_b_a_d(j) + b_a_d(start,j)*temp_a;
		}
		double temp_e = be_d*(arma::sum(inv_V_d.diag())-as_scalar(ph_inv_d*inv_ph_d));
		for(int j = 0; j < num_e; j++)
		{
			v_b_e_d(j) = v_b_e_d(j) + b_e_d(start,j)*temp_e;
		}
	}

	arma::vec res(num_a+num_e);
	arma::vec res_a = v_b_a_m + v_b_a_d;
	if(va>0)
		res_a += 2*d_a*ba/va;
	arma::vec res_e = v_b_e_m + v_b_e_d;
	if(ve>0)
		res_e += 2*d_e*be/ve;
	for(int i = 0; i < num_a; i++)
	{
		res(i) = res_a(i);
	}
	for(int i = 0; i < num_e; i++)
	{
		res(i+num_a) = res_e(i);
	}
	return(Rcpp::wrap(res));
}

RcppExport SEXP gr_AtEt_epsp_c(SEXP v_b_a, SEXP v_b_e, SEXP pheno_m, SEXP pheno_d, SEXP B_a_m, SEXP B_a_d, SEXP b_a, SEXP D_a, SEXP B_e_m, SEXP B_e_d, SEXP b_e, SEXP D_e)
{
    double va = as<double>(v_b_a); // Number of columns
    double ve = as<double>(v_b_e);
	
    arma::vec p_m = as<arma::vec>(pheno_m);
    arma::vec p_d = as<arma::vec>(pheno_d);
	arma::mat b_a_m = as<arma::mat>(B_a_m);
    arma::mat b_a_d = as<arma::mat>(B_a_d);
	arma::mat b_e_m = as<arma::mat>(B_e_m);
    arma::mat b_e_d = as<arma::mat>(B_e_d);
	arma::mat d_a = as<arma::mat>(D_a);
	arma::mat d_e = as<arma::mat>(D_e);
	arma::vec ba = as<arma::vec>(b_a);
	arma::vec be = as<arma::vec>(b_e);
	
	int num_m = p_m.n_elem;
	int num_d = p_d.n_elem;
	int num_a = b_a_m.n_cols;
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
	int penal_a = 2;
	if(d_a(0,1)==-1)
	{penal_a = 1;}
	int penal_e = 2;
	if(d_e(0,1)==-1)
	{penal_e = 1;}
	
	k_m.fill(1);
	k_2d.fill(0.5);
	k_3d.fill(1);
	
	arma::rowvec  gsd_a2m2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_e2mm2a = arma::zeros<arma::rowvec>(num_m/2);
	arma::rowvec  gsd_emme = arma::zeros<arma::rowvec>(num_m);
	
	arma::rowvec  gsd_adda = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_adde = arma::zeros<arma::rowvec>(num_d);
	arma::rowvec  gsd_edde = arma::zeros<arma::rowvec>(num_d);
	
	arma::rowvec  gsd_a2dd2a = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_a2dd2e = arma::zeros<arma::rowvec>(num_d/2);
	arma::rowvec  gsd_2d2 = arma::zeros<arma::rowvec>(num_d/2);

	for(int i = 0; i<(0.5*num_m); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_m = exp(as_scalar(b_a_m.row(start)*ba));
		double be_m = exp(as_scalar(b_e_m.row(start)*be));
		double ba2_m = ba_m*ba_m;
		double be2_m = be_m*be_m;
		double bae_m = ba_m*be_m;
		arma::mat V_m = be_m*diag + ba_m*k_m;
		arma::mat inv_V_m = inv(V_m);
		arma::vec p_m_v(2);
		arma::rowvec p_m_r(2);
		p_m_v << p_m[start] << p_m[end];
		p_m_r << p_m[start] << p_m[end];
		arma::mat inv_ph_m = inv_V_m*p_m_v;
		arma::mat ph_inv_m = p_m_r*inv_V_m;
		
		arma::mat temp_prod = inv_V_m*k_m*inv_V_m;
		gsd_a2m2a(i) = ba2_m*arma::accu(temp_prod);
		temp_prod = inv_V_m*inv_V_m;
		gsd_emme(start) = be2_m*temp_prod(0,0);
		gsd_emme(end) = be2_m*temp_prod(1,1);
		gsd_e2mm2a(i) = bae_m*arma::accu(temp_prod);
	}

	for(int i = 0; i < (0.5*num_d); i++)
	{
		int start = 2*i;
		int end = start + 1;
		double ba_d = exp(as_scalar(b_a_d.row(start)*ba));
		double be_d = exp(as_scalar(b_e_d.row(start)*be));
		double ba2_d = ba_d*ba_d;
		double be2_d = be_d*be_d;
		double bae_d = ba_d*be_d;
		arma::mat V_d = be_d*diag + ba_d*(k_1d+k_2d);
		arma::mat inv_V_d = inv(V_d);
		arma::vec p_d_v(2);
		arma::rowvec p_d_r(2);
		p_d_v << p_d[start] << p_d[end];
		p_d_r << p_d[start] << p_d[end];
		arma::mat inv_ph_d = inv_V_d*p_d_v;
		arma::mat ph_inv_d = p_d_r*inv_V_d;
		
		arma::mat temp_dd = inv_V_d*inv_V_d;
		gsd_adda(start) = ba2_d*temp_dd(0,0);
		gsd_adda(end) = ba2_d*temp_dd(1,1);
		gsd_adde(start) = bae_d*temp_dd(0,0);
		gsd_adde(end) = bae_d*temp_dd(1,1);
		gsd_edde(start) = be2_d*temp_dd(0,0);
		gsd_edde(end) = be2_d*temp_dd(1,1);
		gsd_a2dd2a(i) = ba2_d*0.5*arma::accu(temp_dd);
		gsd_a2dd2e(i) = bae_d*0.5*arma::accu(temp_dd);
		arma::mat temp_i2i = inv_V_d*k_2d*inv_V_d;
		gsd_2d2(i) = ba2_d*0.5*arma::accu(temp_i2i);
	}

	arma::mat gsd_max(num_a+num_e,num_a+num_e);
	arma::mat gsd_b_a(num_a+num_e,num_a+num_e);
	arma::mat gsd_b_e(num_a+num_e,num_a+num_e);
	
	gsd_b_a.fill(0);
	gsd_b_e.fill(0);

	arma::mat b_a_m_h(num_m/2,num_a);
	arma::mat b_e_m_h(num_m/2,num_e);
	arma::mat b_a_d_h(num_d/2,num_a);
	arma::mat b_e_d_h(num_d/2,num_e);
	for(int i = 0; i < (0.5*num_m); i++)
	{
		b_a_m_h.row(i) = b_a_m.row(2*i);
		b_e_m_h.row(i) = b_e_m.row(2*i);
	}
	for(int i = 0; i < (0.5*num_d); i++)
	{
		b_a_d_h.row(i) = b_a_d.row(2*i);
		b_e_d_h.row(i) = b_e_d.row(2*i);
	}
	//arma::mat b_a_m_t = b_a_m.t();
	arma::mat b_e_m_t = b_e_m.t();
	arma::mat b_a_d_t = b_a_d.t();
	arma::mat b_e_d_t = b_e_d.t();
	arma::mat b_a_m_ht = b_a_m_h.t();
	arma::mat b_e_m_ht = b_e_m_h.t();
	arma::mat b_a_d_ht = b_a_d_h.t();
	arma::mat b_e_d_ht = b_e_d_h.t();
	arma::mat temp1 = b_a_m_ht;
	arma::mat temp4 = b_e_m_ht;
	arma::mat temp7 = b_e_m_t;
	for(int i = 0; i < (0.5*num_m); i++)
	{
		temp1.col(i) = gsd_a2m2a(i)*b_a_m_ht.col(i);
		temp4.col(i) = gsd_e2mm2a(i)*b_e_m_ht.col(i);
		//temp6.col(i) = gsd_c2m2c(i)*b_c_m_ht.col(i);
	}
	arma::mat temp3 = b_a_d_ht;
	arma::rowvec tempv1 = gsd_a2dd2a+gsd_2d2;
	arma::mat temp6 = b_e_d_ht;
	//arma::rowvec tempv2 = 0.5*gsd_3dd3+gsd_3d2;
	arma::mat temp8 = b_e_d_t;
	for(int i = 0; i < (0.5*num_d); i++)
	{
		temp3.col(i) = tempv1(i)*b_a_d_ht.col(i);
		temp6.col(i) = gsd_a2dd2e(i)*b_e_d_ht.col(i);
		// temp7.col(i) = gsd_3d3(i)*b_c_d_ht.col(i);
	}
	arma::mat temp2 = b_a_d_t;
	arma::mat temp5 = b_e_d_t;
	for(int i = 0; i < num_d; i++)
	{
		temp2.col(i) = gsd_adda(i)*b_a_d_t.col(i);
		temp5.col(i) = gsd_adde(i)*b_e_d_t.col(i);
		temp8.col(i) = gsd_edde(i)*b_e_d_t.col(i);
	}
	for(int i = 0; i < num_m; i++)
	{
		temp7.col(i) = gsd_emme(i)*b_e_m_t.col(i);
	}
	gsd_max.submat(0,0,num_a-1,num_a-1) = arma::diagmat(2*d_a.diag())/va + temp1*b_a_m_h + 0.25*temp2*b_a_d + temp3*b_a_d_h;
	gsd_max.submat(num_a,0,num_a+num_e-1,num_a-1) = temp4*b_a_m_h + 0.5*temp5*b_a_d + temp6*b_a_d_h;
	gsd_max.submat(0,num_a,num_a-1,num_a+num_e-1) = trans(gsd_max.submat(num_a,0,num_a+num_e-1,num_a-1));
	gsd_max.submat(num_a,num_a,num_a+num_e-1,num_a+num_e-1) = arma::diagmat(2*d_e.diag())/ve + temp7*b_e_m + temp8*b_e_d;

	arma::mat inv_gsd_max = arma::inv(gsd_max);

	gsd_b_a.submat(0,0,num_a-1,num_a-1) = arma::diagmat(2*d_a.diag());
	gsd_b_e.submat(num_a,num_a,num_a+num_e-1,num_a+num_e-1) = arma::diagmat(2*d_e.diag());

	arma::mat temp_m = inv_gsd_max*gsd_b_a;
	double d_v_b_a = ((-1)*as_scalar(trans(ba)*d_a*ba)-0.5*arma::sum(temp_m.diag()))/(va*va) + (num_a-penal_a)/va;
	temp_m = inv_gsd_max*gsd_b_e;
	double d_v_b_e = ((-1)*as_scalar(trans(be)*d_e*be)-0.5*arma::sum(temp_m.diag()))/(ve*ve) + (num_e-penal_e)/ve;

	arma::vec res(2);
	res(0) = d_v_b_a;
	res(1) = d_v_b_e;
   
    return(Rcpp::wrap(res));

}
