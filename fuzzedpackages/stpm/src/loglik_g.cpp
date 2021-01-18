#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

double f_j1_g(double a0, double a2, double b0, double b2, double q0, double q2, double r0, double t, double t0, double geno_a, double geno_b, double geno_q)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double rs = (a + sqrt(pow(a, 2) + 2 * q*pow(b, 2))) / (2 * q);
	double u = 2 * sqrt(pow(a, 2) + 2 * q*pow(b, 2));

	double f = (-0.5)*log(rs + u / ((u / (r0 - rs) + 2 * q)*exp(u*(t - t0)) - 2 * q));

	return(f);
}

double f_j2_g(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f1, double m0, double r0, double t, double y, double t0, double geno_a, double geno_b, double geno_q, double geno_f)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;

	double rs = (a + sqrt(pow(a, 2) + 2 * q*pow(b, 2))) / (2 * q);
	double u = 2 * sqrt(pow(a, 2) + 2 * q*pow(b, 2));
	double c1 = (u / (r0 - rs) + 2 * q)*exp(-1 * u*t0);
	double w_u = (-1)*q*f + q*a*(f - f1) / sqrt(pow(a, 2) + 2 * q*pow(b, 2));
	double c2 = (4 * pow(a, 2)*(f - f1) / u - u*(f - m0))*exp(-0.5*u*t0) / rs;

	// c1*exp(u*tij)
	double temp1 = (u / (r0 - rs) + 2 * q)*exp(u*(t - t0));

	double vf = (-0.5)*pow(y - ((-2)*w_u - (2 * a*(f - f1) / u + f)*temp1 + c2*exp(u*t / 2)) / (2 * q - temp1), 2) / (rs + u / (temp1 - 2 * q));

	return(vf);
}

double f_i1_g(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f1, double mu0, double theta, double m0, double r0, double tau, double t0, double geno_a, double geno_b, double geno_q, double geno_f)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;

	double rs = (a + sqrt(pow(a, 2) + 2 * q*pow(b, 2))) / (2 * q);
	double u = 2 * sqrt(pow(a, 2) + 2 * q*pow(b, 2));
	double c1 = (u / (r0 - rs) + 2 * q)*exp((-1) * u*t0);
	double w_u = (-1)*q*f + q*a*(f - f1) / sqrt(pow(a, 2) + 2 * q*pow(b, 2));
	double c2 = (4 * pow(a, 2)*(f - f1) / u - u*(f - m0))*exp((-0.5)*u*t0) / rs;

	// 2*q-c1*exp(u*tau)
	double temp1 = 2 * q - (u / (r0 - rs) + 2 * q)*exp(u*(tau - t0));

	double vf = log(mu0*exp(theta*tau) + q*pow(2 * a*(f - f1)*(1 - 4 * q / temp1) / u + c2*exp(0.5*u*tau) / temp1, 2) + 0.5*(a + sqrt(pow(a, 2) + 2 * q*pow(b, 2))) + u*q / ((-1)*temp1));

	return(vf);
}

double mu_int_g(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f1, double mu0, double theta, double m0, double r0, double tij, double t0, double geno_a, double geno_b, double geno_q, double geno_f)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;

	double rs = (a + sqrt(pow(a, 2) + 2 * q*pow(b, 2))) / (2 * q);
	double u = 2 * sqrt(pow(a, 2) + 2 * q*pow(b, 2));
	double c1 = (u / (r0 - rs) + 2 * q)*exp(-1 * u*t0);
	double w_u = (-1)*q*f + q*a*(f - f1) / sqrt(pow(a, 2) + 2 * q*pow(b, 2));
	double c2 = (4 * pow(a, 2)*(f - f1) / u - u*(f - m0))*exp((-0.5)*u*t0) / rs;

	// 2*q - c1*exp(u*t0)
	double temp1 = (-1)*(u) / (r0 - rs);
	// 2*q - c1*exp(u*tij)
	double temp2 = 2 * q - (u / (r0 - rs) + 2 * q)*exp(u*(tij - t0));

	double f_1 = (a + sqrt(pow(a, 2) + 2 * q*pow(b, 2)))*(tij - t0) / 2;
	double f_2 = 0.5*(log(temp2 / temp1) - u*(tij - t0));
	double f_3 = 4 * (tij - t0)*pow(a*(f - f1) / u, 2);
	double f_4 = (32 * q*pow(a*(f - f1) / u, 2))*((tij - t0) / (2 * q) - 0.5*log(temp2 / temp1) / (u*q) + 1 / (u*temp2) - 1 / (u*temp1));
	double f_5 = pow(c2, 2) / (c1*u)*(1 / temp2 - 1 / temp1);
	double f_6 = (-8)*a*c2*(f - f1)*sqrt(-0.5*c1 / q) / (c1*pow(u, 2))*(atan(sqrt((-0.5)*c1 / q)*exp(0.5*u*tij)) - atan(sqrt(-0.5*c1 / q)*exp(0.5*u*t0)));
	double f_7 = (-16)*pow(a, 2)*pow(f - f1, 2) / pow(u, 3)*(log(temp2 / temp1) - u*(tij - t0));
	double f_8 = 16 * a*q*c2*(f - f1) / (c1*pow(u, 2))*(exp((-0.5)*u*tij) / (temp2)-exp((-0.5)*u*t0) / (temp1)-exp((-0.5)*u*tij) / (2 * q) + exp((-0.5)*u*t0) / (2 * q) - sqrt((-0.5)*c1 / q) / (2 * q)*(atan(sqrt((-0.5)*c1 / q)*exp(0.5*u*tij)) - atan(sqrt((-0.5)*c1 / q)*exp(0.5*u*t0))));
	double f_9 = 0;
	if (theta != 0)
	{
		f_9 = mu0*(exp(theta*tij) - exp(theta*t0)) / theta;
	}
	else
	{
		f_9 = mu0*(tij - t0);
	}

	double vf = f_1 + f_2 + f_9 + q*(f_3 + f_4 + f_5 + f_6 - f_7 - f_8);

	return(vf);
}


RcppExport SEXP mloglik_g(SEXP param, SEXP m0, SEXP r0, SEXP tau, SEXP yij, SEXP delta, SEXP tij, SEXP n_j, SEXP geno_a, SEXP geno_b, SEXP geno_q, SEXP geno_f, SEXP t0)
{

	arma::vec param_c = as<arma::vec>(param); //
	double lik = 0;
	double a0 = param_c[0];
	double a2 = param_c[1];
	double b0 = param_c[2];
	double b2 = param_c[3];
	double q0 = param_c[4];
	double q2 = param_c[5];
	double f0 = param_c[6];
	double f2 = param_c[7];
	double f1 = param_c[8];
	double mu0 = param_c[9];
	double theta = param_c[10];

	arma::vec m0_v = as<arma::vec>(m0);
	arma::vec r0_v = as<arma::vec>(r0);
	arma::vec tau_v = as<arma::vec>(tau);
	arma::vec y_v = as<arma::vec>(yij);
	arma::vec delta_v = as<arma::vec>(delta);
	arma::vec t_v = as<arma::vec>(tij);
	arma::vec n_v = as<arma::vec>(n_j);
	arma::vec t0_v = as<arma::vec>(t0);
	arma::vec geno_a_v = as<arma::vec>(geno_a);
  arma::vec geno_b_v = as<arma::vec>(geno_b);
  arma::vec geno_q_v = as<arma::vec>(geno_q);
  arma::vec geno_f_v = as<arma::vec>(geno_f);

	int i_n = n_v.n_elem;
	int t_n = m0_v.n_elem;
	int index = 0;

	for (int i = 0; i < i_n; i++)
	{
		int start = index;
		int end = index + n_v[i];
		double temp1 = 0;
		double geno_a_i = geno_a_v[i];
    double geno_b_i = geno_b_v[i];
    double geno_q_i = geno_q_v[i];
    double geno_f_i = geno_f_v[i];

		for (int j = start; j < end; j++)
		{
			temp1 = temp1 + f_j1_g(a0, a2, b0, b2, q0, q2, r0_v[j], t_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i) + f_j2_g(a0, a2, b0, b2, q0, q2, f0, f2, f1, m0_v[j], r0_v[j], t_v[j], y_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, geno_f_i) - mu_int_g(a0, a2, b0, b2, q0, q2, f0, f2, f1, mu0, theta, m0_v[j], r0_v[j], t_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, geno_f_i);
			// temp1 = temp1 + f_j1(a, b, q, r0_v[j], t_v[j], t0_v[j]) + f_j2(a, b, q, f, f1, m0_v[j], r0_v[j], t_v[j], y_v[j], t0_v[j]);
		}

		lik += temp1;

		int pre_end = end - 1;
		lik += delta_v[i] * f_i1_g(a0, a2, b0, b2, q0, q2, f0, f2, f1, mu0, theta, m0_v[pre_end], r0_v[pre_end], tau_v[i], t0_v[pre_end], geno_a_i, geno_b_i, geno_q_i, geno_f_i);

		index = end;

	}

	lik *= -1;

	return(Rcpp::wrap(lik));

}
