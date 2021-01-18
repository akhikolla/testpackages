#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

void d_f_j1_g_2(double a0, double a2, double b0, double b2, double q0, double q2, double r0, double tij, double t0, double geno_a, double geno_b, double geno_q, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;

	double e1 = pow(b, 2);
	double e4 = sqrt(2 * (e1 * q) + pow(a, 2));
	double e5 = 2 * q;
	double e6 = a + e4;
	double e7 = e6 / e5;
	double e8 = r0 - e7;
	double e9 = tij - t0;
	double e12 = 2 * (e4 / e8) + e5;
	double e13 = exp(2 * (e4 * e9));
	double e15 = e12 * e13 - e5;
	double e16 = a / e4;
	double e17 = e7 + 2 * (e4 / e15);
	double e20 = 1 + e16;
	double e21 = 2 * (e6 / pow(e5, 2));
	double e22 = 2 * e16;
	double e23 = 2 * (e1 / e4);
	double e25 = 4 * (q / e4);
	double e26 = e1 / (2 * (q * e4));

	// double df[3];
	double d_a = (-1)*(0.5 * ((e20 / e5 + (e22 - 2 * (((e20 * e4 / (q * e8) + e22) / e8 + 2 * (a * e12 * e9 / e4)) * e13 *e4 / e15)) / e15) / e17));
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = (-1)*(0.5 * (b * ((e25 - 2 * (((2 / e8 + e25) / e8 + 4 * (q * e12 * e9 / e4)) * e13 *e4 / e15)) / e15 + 1 / e4) / e17));
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = (-1)*(0.5 * (((e23 - 2 * ((((2 * ((e26 - e21) * e4 / e8) + e23) / e8 + 2 + 2 * (e1 * e12 * e9 / e4)) * e13 - 2) * e4 / e15)) / e15 + e26 - e21) / e17));
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;

	// return(df);
}

void d_f_j2_g_2(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f10, double f12, double m0, double r0, double tij, double yij, double t0, double geno_a, double geno_b, double geno_q, double geno_f, double geno_f1, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;
	double f1 = (f12 - f10) / 2 * geno_f1 + f10;

	double e1 = pow(a, 2);
	double e2 = pow(b, 2);
	double e3 = e2 * q;
	double e5 = 2 * e3 + e1;
	double e6 = sqrt(e5);
	double e7 = 2 * q;
	double e8 = a + e6;
	double e9 = e8 / e7;
	double e10 = r0 - e9;
	double e11 = tij - t0;
	double e13 = 2 * (e6 / e10) + e7;
	double e14 = f - f1;
	double e16 = exp(2 * (e6 * e11));
	double e17 = e13 * e16;
	double e18 = a * e14;
	double e19 = e18 / e6;
	double e20 = e1 * e14;
	double e24 = 2 * (e20 / e6) - 2 * ((f - m0) * e6);
	double e26 = e7 - e17;
	double e27 = exp((-1)*t0 * e6);
	double e28 = exp(tij * e6);
	double e29 = e17 - e7;
	double e30 = e19 + f;
	double e39 = q * (2 * (e24 * e27 * e28 / e8) - 2 * (e19 - f)) - e13 * e30 * e16;
	double e40 = e39 / e26;
	double e41 = a / e6;
	double e42 = e9 + 2 * (e6 / e29);
	double e43 = yij - e40;
	double e46 = 1 + e41;
	double e47 = 2 * e8 / pow(e7, 2);
	double e49 = 2 * e6;
	double e50 = e2 / (2 * (q * e6));
	double e51 = pow(e49, 2);
	double e52 = 2 * (e2 / e6);
	double e53 = e50 - e47;
	double e56 = (2 * (e53 * e6 / e10) + e52) / e10 + 2 + 2 * (e2 * e13 * e11 / e6);
	double e57 = 0.5 * (t0 * e24);
	double e58 = 2 * e41;
	double e59 = 4 * (q / e6);
	double e60 = e1 / e5;
	double e62 = (e46 * e6 / (q * e10) + e58) / e10 + 2 * (a * e13 * e11 / e6);
	double e63 = e56 * e16;
	double e64 = e42 * e26;
	double e65 = e51 * e6;
	double e67 = (2 / e10 + e59) / e10 + 4 * (q * e13 * e11 / e6);
	double e68 = e40 - e30;
	double e71 = e57 + 4 * (e20 / e51) + f - m0;

	//double df[5];
	double d_a = ((e62 * e68 * e16 + q * (2 * (a * ((8 - 4 * e60) * e14 / 2 - 2 * (e57 + f - m0)) / e6 - e46 * e24 / e8) + 2 * (a * tij * e24 / e6)) * e27 * e28 / e8 - ((2 - 2 * e60) * e13 * e16 / 2 + 2 * (q * (1 - e60))) * e14 / e6) / e26 + 0.5 * ((e46 / e7 + (e58 - 2 * (e62 * e16 * e6 / e29)) / e29) * e43 / e42)) * e43 / e42;
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = b * (((e67 * e68 + 8 * (a * q * e13 * e14 / e65)) * e16 + pow(q, 2) * ((4 * (tij * e24) - 2 * (2 * (e24 / e8) + 4 * e71)) * e27 * e28 / e8 + 4 * (e18 / e5)) / e6) / e26 + 0.5 * (((e59 - 2 * (e67 * e16 * e6 / e29)) / e29 + 1 / e6) * e43 / e42)) * e43 / e42;
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = (((4 * (a * e2 * e13 * e14 / e65) - e56 * e30) * e16 + q * (2 * (e2 * tij * e24 / e6) - 2 * (2 * (e2 * e71 / e6) + 2 * (q * e24 * e53 / e8))) * e27 * e28 / e8 - ((2 - e63) * e39 / e26 + 2 * (a * (1 - e3 / e5) * e14 / e6 - f))) / e26 + 0.5 * (((e52 - 2 * ((e63 - 2) * e6 / e29)) / e29 + e50 - e47) * e43 / e42)) * e43 / e42;
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;
	double d_f = (q * (2 * ((2 * (e1 / e6) - e49) * e27 * e28 / e8) - 2 * (e41 - 1)) - e46 * e13 * e16) * e43 / e64;
	df[6] = d_f*(1 - 0.5*geno_f);
	df[7] = d_f*0.5*geno_f;
	double d_f1 = a * (e17 + q * (2 - 4 * (a * e27 * e28 / e8))) * e43 / (e64 * e6);
	df[8] = d_f1*(1 - 0.5*geno_f1);
	df[9] = d_f1*0.5*geno_f1;

	//return(df);

}

void d_f_i1_g_2(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f10, double f12, double mu00, double mu02, double theta, double m0, double r0, double tau, double t0, double geno_a, double geno_b, double geno_q, double geno_f, double geno_f1, double geno_mu, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;
	double f1 = (f12 - f10) / 2 * geno_f1 + f10;
	double mu0 = (mu02 - mu00) / 2 * geno_mu + mu00;

	double e1 = pow(a, 2);
	double e2 = pow(b, 2);
	double e3 = e2 * q;
	double e5 = 2 * e3 + e1;
	double e6 = sqrt(e5);
	double e7 = 2 * q;
	double e8 = a + e6;
	double e10 = r0 - e8 / e7;
	double e11 = tau - t0;
	double e14 = 2 * (e6 / e10) + e7;
	double e15 = exp(2 * (e6 *e11));
	double e17 = e7 - e14 *e15;
	double e18 = f - f1;
	double e19 = e1 *e18;
	double e23 = 2 * (e19 / e6) - 2 * ((f - m0) *e6);
	double e25 = exp((-1)*t0 *e6);
	double e26 = exp(tau *e6);
	double e27 = e17 *e8;
	double e28 = 1 - 4 * (q / e17);
	double e29 = a *e28;
	double e31 = e29 *e18 / e6;
	double e32 = q *e23;
	double e34 = 2 * (e32 *e25 *e26 / e27) + e31;
	double e35 = exp(tau * theta);
	double e36 = a / e6;
	double e41 = 0.5 *e8 + mu0 *e35 + q * (pow(e34, 2) - 2 * (e6 / e17));
	double e42 = 2 * e6;
	double e43 = 1 + e36;
	double e49 = e2 / (2 * (q *e6)) - 2 * (e8 / pow(e7, 2));
	double e50 = e2 / e6;
	double e52 = pow((-1)*e17, 2);
	double e54 = (e43 *e6 / (q *e10) + 2 * e36) / e10 + 2 * (a *e14 *e11 / e6);
	double e56 = ((2 / e10 + 4 * (q / e6)) / e10 + 4 * (q *e14 *e11 / e6)) *e15;
	double e57 = 0.5 * t0 *e23;
	double e58 = 2 - ((2 * (e49 *e6 / e10) + 2 * e50) / e10 + 2 + 2 * (e2 *e14 *e11 / e6)) *e15;
	double e59 = e54 *e15;
	double e60 = e17 *e6;
	double e61 = pow(e17, 2);
	double e65 = e57 + 4 * (e19 / pow(e42, 2)) + f - m0;
	double e66 = a * q;
	double e67 = e1 / e5;
	double e68 = q *e58;

	// double df[11];
	double d_a = (0.5 *e43 + q * (2 * (((e28 * (2 - 2 * e67) - 8 * (e66 *e54 *e15 / e61)) *e18 / e42 + q * (2 * ((e59 / e17 + a * tau / e6) *e23) + 2 * (a * ((8 - 4 * e67) *e18 / 2 - 2 * (e57 + f - m0)) / e6 - e43 *e23 / e8)) *e25 *e26 / e27) *e34) - (2 * (e59 *e6 / e52) + 2 * (a / e60)))) / e41;
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = b * q * (1 / e6 + q * (2 * (((2 * ((e56 / e17 + 2 * (q * tau / e6)) *e23) - 2 * (q * (2 * (e23 / e8) + 4 * e65) / e6)) *e25 *e26 / e27 - a * (4 * (e28 / e5) + 8 * (e56 / e61)) *e18 / e42) *e34) - 4 / e60) - 2 * (e56 *e6 / e52)) / e41;
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = (e34 * (e31 + q * (2 * (e23 *e25 *e26 / e27) + 2 * (q * (2 * (e23 * (e2 * tau / e6 - e58 / e17)) - 2 * (2 * (e2 *e65 / e6) + 2 * (e32 *e49 / e8))) *e25 *e26 / e27 - a * (2 * ((4 - 4 * (e68 / e17)) / e17) + 2 * (e2 *e28 / e5)) *e18 / e42))) + 0.5 *e50 + 2 * (e68 *e6 / e52) - (2 * (e3 / e6) + e42) / e17) / e41;
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;
	double d_f = 2 * (q *e34 * (2 * (q * (2 * (e1 / e6) - e42) *e25 *e26 / e27) + e29 / e6) / e41);
	df[6] = d_f*(1 - 0.5*geno_f);
	df[7] = d_f*0.5*geno_f;
	double d_f1 = -(2 * (e66 *(1 + q * (4 * (a *e25 *e26 / e8) - 4) / e17) *e34 / (e41 *e6)));
	df[8] = d_f1*(1 - 0.5*geno_f1);
	df[9] = d_f1*0.5*geno_f1; 
	//df[9] = e35 / e41;
	double d_mu = e35 / e41;
	df[10] = d_mu*(1 - 0.5*geno_mu);
	df[11] = d_mu*0.5*geno_mu;
	df[12] = mu0 * tau *e35 / e41;

	// return(df);
}

void dev_mu_int_g_2(double a0, double a2, double b0, double b2, double q0, double q2, double f0, double f2, double f10, double f12, double mu00, double mu02, double theta, double tij, double r0, double m0, double t0, double geno_a, double geno_b, double geno_q, double geno_f, double geno_f1, double geno_mu, double * df)
{
	double a = (a2 - a0) / 2 * geno_a + a0;
	double b = (b2 - b0) / 2 * geno_b + b0;
	double q = (q2 - q0) / 2 * geno_q + q0;
	double f = (f2 - f0) / 2 * geno_f + f0;
	double f1 = (f12 - f10) / 2 * geno_f1 + f10;
	double mu0 = (mu02 - mu00) / 2 * geno_mu + mu00;

	double e1 = pow(a,2);
    double e2 = pow(b,2);
    double e4 = 2 * (e2 * q) + e1;
    double e5 = sqrt(e4);
    double e6 = 2 * q;
    double e7 = a + e5;
    double e8 = e7/e6;
    double e9 = r0 - e8;
    double e10 = 2 * (e5/e9);
    double e11 = e10 + e6;
    double e12 = f - f1;
    double e13 = 2 * e5;
    double e14 = tij - t0;
    double e15 = t0 * e5;
    double e16 = e1 * e12;
    double e18 = exp(2 * (e5 * e14));
    double e19 = e16/e5;
    double e23 = 2 * e19 - 2 * ((f - m0) * e5);
    double e25 = e6 - e11 * e18;
    double e27 = exp(-(2 * e15));
    double e29 = exp(-e15);
    double e30 = pow(e13,2);
    double e31 = a/e5;
    double e32 = q * e9;
    double e33 = e11 * e30;
    double e36 = 1 + e31;
    double e37 = 2 * (e7/pow(e6,2));
    double e38 = e33 * e27;
    double e41 = e2/(2 * (q * e5)) - e37;
    double e43 = 0.5 + e32/e13;
    double e44 = e2/e5;
    double e45 = exp(-(tij * e5));
    double e46 = q * e23;
    double e48 = q/e25 - 0.5;
    double e51 = e36 * e5/e32 + 2 * e31;
    double e55 = 2 * (e41 * e5/e9) + 2 * e44;
    double e57 = 2/e9 + 4 * (q/e5);
    double e58 = a * e12;
    double e60 = e43 * e29 + e45 * e48;
    double e61 = pow(e58/e13,2);
    double e62 = pow(e12,2);
    double e64 = e9/e13 + 1/e25;
    double e65 = 0.5 * (t0 * e23);
    double e66 = e51/e9;
    double e68 = e55/e9 + 2;
    double e70 = e11 * pow(e7,2) * e27;
    double e71 = e57/e9;
    double e72 = e16/e30;
    double e73 = pow(e29,2);
    double e76 = pow(2 * (e46 * e29/e7),2)/(e11 * e27) + 32 * (q * e61);
    double e77 = pow(e13,3);
    double e80 = e65 + 4 * e72 + f - m0;
    double e81 = a * e11;
    double e82 = e1/e4;
    double e83 = e2 * e11;
    double e84 = q * e11;
    double e87 = pow(-(4 * (e4/e9)),2) * e9;
    double e88 = e66 + 2 * (e81 * e14/e5);
    double e89 = (e68 + 2 * (e83 * e14/e5)) * e18;
    double e90 = pow(e38,2);
    double e91 = (e71 + 4 * (e84 * e14/e5)) * e18;
    double e92 = pow(2 * (e25 * e5),2);
    double e93 = 32 * e23;
    double e94 = a * t0;
    double e95 = e2 * t0;
    double e96 = exp(t0 * theta);
    double e97 = exp(theta * tij);
    double e100 = pow(-e10,2) * e9;
    double e101 = e88 * e18;
    double e102 = e66 - 2 * (e94 * e11/e5);
    double e105 = e68 - 2 * (e95 * e11/e5);
    double e108 = e33 * e7 * e27;
    double e109 = pow(e25,2);
    double e110 = e77 * e5;
    double e111 = e71 - 4 * (q * t0 * e11/e5);
    double e115 = (e97 - e96)/theta;
    double e117 = 2 - e89;
    double e121 = 2 * (e1/e5) - e13;
    double e123 = 2 * (e2 * e80/e5) + 2 * (e46 * e41/e7);
    double e126 = a * ((8 - 4 * e82) * e12/2 - 2 * (e65 + f - m0))/e5 - e36 * e23/e7;
    double e129 = a * (1 - e82) * e62/e30;
    double e130 = e1 * e62;
    double e131 = q * (2 * (e23/e7) + 4 * e80);
	
	double d_a = (q * (e36/e6 + 8 * e129) - e31) * e14 + q * (q * (e64 * (2 * (e46 * (4 * e126 - 2 * (e102 * e23/e11)) * e73/e70) + 64 * e129)/e13 - (e60 * ((e93 + 32 * (a * e126))/e38 - 64 * (a * (2 * (e102 * e5) + 4 * (e81/e5)) * e23 * e27 * e5/e90)) - 32 * (a * ((e94 * e43/e5 + q * e51/e100) * e29 + (a * tij * e48/e5 - q * e88 * e18/e109) * e45) * e23/e38)) * e29 * e12/e7) - ((2 * (e51 * e5) + 4 * a)/e87 + (2 * (a * e25/e5) - 2 * (e101 * e5))/e92) * e76) - 0.5 * ((e101 + e51 * e25/e13)/e25);
	df[0] = d_a*(1 - 0.5*geno_a);
	df[1] = d_a*0.5*geno_a;
	double d_b = b * (q * (q * (a * (e60 * (32 * (e131/(e38 * e5)) + 64 * ((2 * (e111 * e5) + 8 * (e84/e5)) * e23 * e27 * e5/e90)) + 32 * (q * ((e57/e100 + 2 * (t0 * e43/e5)) * e29 + (2 * (tij * e48/e5) - e91/e109) * e45) * e23/e38)) * e29 * e12/e7 - q * e64 * (2 * ((2 * (e111 * e23/e11) + 4 * (e131/e5)) * e23 * e73/e70) + 256 * (e130/e110))/e13) - (((2 * (e57 * e5) + 8 * q)/e87 + (4 * (q * e25/e5) - 2 * (e91 * e5))/e92) * e76 + (1 + 32 * (e1 * q * e62/e77)) * e14/e5)) - 0.5 * ((e91 + e25 * e57/e13)/e25));
	df[2] = d_b*(1 - 0.5*geno_b);
	df[3] = d_b*0.5*geno_b;
	double d_q = e76 * e64/e13 + (e8 + 4 * e61 + q * (e2 * (1/e6 - 16 * (e130/e77))/e5 - e37) - e44) * e14 + 0.5 * ((2 - (e89 + e55 * e25/e13))/e25) + q * (e64 * (32 * e61 - q * (128 * (e1 * e2 * e62/e110) + 2 * (q * (2 * (e105 * e23/e11) + 4 * e123) * e23 * e73/e70)))/e13 - (((2 * (e117 * e5) + 2 * (e2 * e25/e5))/e92 + (2 * (e55 * e5) + 4 * e2)/e87) * e76 + a * (32 * (e60 * e23/e38) + q * (32 * ((((1 - q * e117/e25)/e25 - e2 * tij * e48/e5) * e45 - (e95 * e43 - (1 - q * e55/e13) * e9/2) * e29/e5) * e23/e38) - e60 * (32 * (e123/e38) + 64 * ((2 * (e105 * e5) + 4 * (e83/e5)) * e23 * e27 * e5/e90)))) * e29 * e12/e7));
	df[4] = d_q*(1 - 0.5*geno_q);
	df[5] = d_q*0.5*geno_q;
	double d_f = q * (8 * (e16 * e14/e30) + q * (e64 * (64 * e72 + 8 * (e46 * e121 * e73/e70))/e13 - a * e60 * (32 * (e121 * e12) + e93) * e29/e108));
	df[6] = d_f*(1 - 0.5*geno_f);
	df[7] = d_f*0.5*geno_f;
	double d_f1 = a * q * (q * (e60 * (e93 + 64 * e19) * e29/e108 - a * e64 * (16 * (e46 * e73/(e70 * e5)) + 64 * (e12/e30))/e13) - 8 * (e58 * e14/e30));
	df[8] = d_f1*(1 - 0.5*geno_f1);
	df[9] = d_f1*0.5*geno_f1; 
	//df[9] = e115;
	double d_mu = e115;
	df[10] = d_mu*(1 - 0.5*geno_mu);
	df[11] = d_mu*0.5*geno_mu;
	df[12] = 0;
	if (theta != 0)
	{
		df[12] = mu0 * (tij * e97 - (e115 + t0 * e96))/theta;
	}

}


RcppExport SEXP devlik_g_2(SEXP param, SEXP m0, SEXP r0, SEXP tau, SEXP yij, SEXP delta, SEXP tij, SEXP n_j, SEXP geno_a, SEXP geno_b, SEXP geno_q, SEXP geno_f, SEXP geno_f1, SEXP geno_mu, SEXP t0)
{

	arma::vec param_c = as<arma::vec>(param);
	// double lik = 0;

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
	arma::vec geno_f1_v = as<arma::vec>(geno_f1);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu);

	double a0 = param_c(0);
	double a2 = param_c(1);
	double b0 = param_c(2);
	double b2 = param_c(3);
	double q0 = param_c(4);
	double q2 = param_c(5);
	double f0 = param_c(6);
	double f2 = param_c(7);
	double f10 = param_c(8);
	double f12 = param_c(9);
	double mu00 = param_c(9);
	double mu02 = param_c(10);
	double theta = param_c(11);

	int i_n = n_v.n_elem;
	int t_n = m0_v.n_elem;
	int index = 0;

	double dev_a0 = 0;
	double dev_a2 = 0;
	double dev_b0 = 0;
	double dev_b2 = 0;
	double dev_q0 = 0;
	double dev_q2 = 0;
	double dev_f0 = 0;
	double dev_f2 = 0;
	double dev_f10 = 0;
	double dev_f12 = 0;
	double dev_mu00 = 0;
	double dev_mu02 = 0;
	double dev_theta = 0;

	for (int i = 0; i < i_n; i++)
	{
		int start = index;
		int end = index + n_v(i);
		double geno_a_i = geno_a_v(i);
		double geno_b_i = geno_b_v(i);
		double geno_q_i = geno_q_v(i);
		double geno_f_i = geno_f_v(i);
		double geno_f1_i = geno_f1_v(i);
		double geno_mu_i = geno_mu_v(i);

		for (int j = start; j < end; j++)
		{
			double re_1[6];
			for (int k = 0; k < 6; k++)
			{
				re_1[k] = 0;
			}
			d_f_j1_g_2(a0, a2, b0, b2, q0, q2, r0_v[j], t_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, re_1);
			double re_2[10];
			for (int k = 0; k < 10; k++)
			{
				re_2[k] = 0;
			}
			d_f_j2_g_2(a0, a2, b0, b2, q0, q2, f0, f2, f10, f12, m0_v[j], r0_v[j], t_v[j], y_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, geno_f_i, geno_f1_i, re_2);
			double re_4[13];
			for (int k = 0; k < 13; k++)
			{
				re_4[k] = 0;
			}

			dev_mu_int_g_2(a0, a2, b0, b2, q0, q2, f0, f2, f10, f12, mu00, mu02, theta, t_v[j], r0_v[j], m0_v[j], t0_v[j], geno_a_i, geno_b_i, geno_q_i, geno_f_i, geno_f1_i, geno_mu_i, re_4);

			dev_a0 = dev_a0 + re_1[0] + re_2[0] - re_4[0];
			dev_a2 = dev_a2 + re_1[1] + re_2[1] - re_4[1];
			dev_b0 = dev_b0 + re_1[2] + re_2[2] - re_4[2];
			dev_b2 = dev_b2 + re_1[3] + re_2[3] - re_4[3];
			dev_q0 = dev_q0 + re_1[4] + re_2[4] - re_4[4];
			dev_q2 = dev_q2 + re_1[5] + re_2[5] - re_4[5];
			dev_f0 = dev_f0 + re_2[6] - re_4[6];
			dev_f2 = dev_f2 + re_2[7] - re_4[7];
			dev_f10 = dev_f10 + re_2[8] - re_4[8];
			dev_f12 = dev_f12 + re_2[9] - re_4[9];
			dev_mu00 -= re_4[10];
			dev_mu02 -= re_4[11];
			dev_theta -= re_4[12];
		}

		int pre_end = end - 1;
		double re_3[13];

		if (delta_v(i)>0)
		{
			d_f_i1_g_2(a0, a2, b0, b2, q0, q2, f0, f2, f10, f12, mu00, mu02, theta, m0_v(pre_end), r0_v(pre_end), tau_v(i), t0_v(pre_end), geno_a_i, geno_b_i, geno_q_i, geno_f_i, geno_f1_i, geno_mu_i, re_3);

			dev_a0 += re_3[0];
			dev_a2 += re_3[1];
			dev_b0 += re_3[2];
			dev_b2 += re_3[3];
			dev_q0 += re_3[4];
			dev_q2 += re_3[5];
			dev_f0  += re_3[6];
			dev_f2 += re_3[7];
			dev_f10 += re_3[8];
			dev_f12 += re_3[9];
			dev_mu00 += re_3[10];
			dev_mu02 += re_3[11];
			dev_theta += re_3[12];

		}

		index = end;
	}

	arma::vec res = arma::zeros<arma::vec>(13);
	res(0) = dev_a0*(-1);
	res(1) = dev_a2*(-1);
	res(2) = dev_b0*(-1);
	res(3) = dev_b2*(-1);
	res(4) = dev_q0*(-1);
	res(5) = dev_q2*(-1);
	res(6) = dev_f0*(-1);
	res(7) = dev_f2*(-1);
	res(8) = dev_f10*(-1);
	res(9) = dev_f12*(-1);
	res(10) = dev_mu00*(-1);
	res(11) = dev_mu02*(-1);
	res(12) = dev_theta*(-1);

	return(Rcpp::wrap(res));

}
