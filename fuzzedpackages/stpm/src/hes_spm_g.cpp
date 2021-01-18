#include <stdio.h>
#include <string>
#include <vector>
#include <math.h>
#include <iostream>

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;

RcppExport SEXP d_f_i1_a0_abqff1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
	double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);
	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e6 = b0 + geno_b * (b2 - b0)/2;
	double    e7 = pow(e6,2);
	double    e8 = e5 + 2 * (e7 * e3);
	double    e9 = sqrt(e8);
	double    e10 = 2 * e3;
	double    e11 = e4 + e9;
	double    e12 = e11/e10;
	double    e13 = r0 - e12;
	double    e14 = tau - t0;
	double    e15 = e10 + 2 * (e9/e13);
	double    e18 = geno_f * (f2 - f0)/2;
	double    e20 = exp(2 * (e9 * e14));
	double    e21 = f0 + e18;
	double    e23 = e10 - e15 * e20;
	double    e24 = e4/e9;
	double    e25 = e21 - f1;
	double    e26 = e21 - m0;
	double    e27 = e24 + 1;
	double    e28 = e5 * e25;
	double    e29 = e3 * e13;
	double    e33 = 2 * (e28/e9) - 2 * (e26 * e9);
	double    e34 = e27 * e9;
	double    e36 = e34/e29 + 2 * e24;
	double    e37 = e36/e13;
	double    e40 = e15 * e4 * e14/e9;
	double    e41 = 0.5 * geno_a;
	double    e42 = e37 + 2 * e40;
	double    e43 = e23 * e11;
	double    e44 = e5/e8;
	double    e46 = 1 - 4 * (e3/e23);
	double    e48 = exp(-(t0 * e9));
	double    e49 = exp(tau * e9);
	double    e50 = 1 - e41;
	double    e51 = e7/e9;
	double    e52 = e46 * e4;
	double    e53 = e42 * e20;
	double    e54 = 2 * e29;
	double    e57 = e52 * e25/e9 + 2 * (e33 * e48 * e49 * e3/e43);
	double    e58 = tau * e4;
	double    e59 = t0 * e33;
	double    e60 = pow(e23,2);
	double    e61 = 4 * e44;
	double    e62 = e58/e9;
	double    e64 = e53/e23 + e62;
	double    e65 = 8 - e61;
	double    e68 = 0.5 * e59 + f0 + e18;
	double    e70 = 2 * (e68 - m0);
	double    e71 = e50 * e4;
	double    e72 = 2 * e9;
	double    e75 = e65 * e25/2 - e70;
	double    e76 = e71/e9;
	double    e77 = 0.5 * e24;
	double    e78 = 2 - 2 * e44;
	double    e79 = 0.5 * e51;
	double    e84 = 2 * (e64 * e33);
	double    e85 = 2 * (e75 * e4/e9 - e27 * e33/e11);
	double    e86 = e42 * e4;
	double    e87 = e84 + e85;
	double    e92 = e46 * e78 - 8 * (e86 * e20 * e3/e60);
	double    e93 = e3/e9;
	double    e94 = e15 * e7;
	double    e95 = e15 * e3;
	double    e97 = pow(-e23,2);
	double    e98 = 2 * e26;
	double    e99 = 0.5 * geno_q;
	double    e100 = exp(tau * theta);
	double    e103 = e92 * e25/e72 + e87 * e48 * e49 * e3/e43;
	double    e104 = pow(e57,2);
	double    e105 = 2 * (e9/e23);
	double    e106 = e104 - e105;
	double    e108 = e76 + 1 - e41;
	double    e110 = e94 * e14/e9;
	double    e112 = e95 * e14/e9;
	double    e113 = e28/e8;
	double    e114 = 0.5 + e77;
	double    e117 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e120 = e106 * e3 + 0.5 * e11 + e100 * e117;
	double    e121 = e23 * e9;
	double    e122 = e4/e8;
	double    e123 = 2 * e113;
	double    e125 = e51 - e11/e3;
	double    e126 = e79 - e12;
	double    e127 = e27/e9;
	double    e128 = pow(e29,2);
	double    e129 = 1 - e99;
	double    e130 = e123 + e98;
	double    e131 = 0.5 * e44;
	double    e132 = 2 * e93;
	double    e140 = e108 * e9/e54 + e76;
	double    e144 = e114 * e9/e54 + e77;
	double    e146 = 2 * (e53 * e9/e97);
	double    e147 = 2 * (e103 * e57);
	double    e149 = 2 * ((e126 * e9/e54 + e79)/e13);
	double    e150 = 2 * ((e93 + 1/(2 * e13))/e13);
	double    e151 = 2 * ((1/e13 + e132)/e13);
	double    e152 = 2 + 2 * ((e125 * e9/e54 + e51)/e13);
	double    e153 = e8 * e9;
	double    e156 = e147 - (e146 + 2 * (e4/e121));
	double    e161 = e156 * e3 + 0.5 * e27;
	double    e162 = 0.5 - e131;
	double    e163 = 1 - (e50 * e5/e8 + e41);
	double    e171 = e40 + 2 * (e144/e13);
	double    e173 = 1 - (e110 + 1 + e149) * e20;
	double    e174 = 2 - (e152 + 2 * e110) * e20;
	double    e176 = 2 * (e140/e13) + 2 * (e50 * e15 * e4 * e14/e9);
	double    e177 = e150 + 2 * e112;
	double    e178 = e151 + 4 * e112;
	double    e181 = 2 * (e5/e9) - e72;
	double    e182 = e37 + e27 * e3 * e9/e128;
	double    e183 = e127 - e122;
	double    e185 = e23 * e7/e9;
	double    e187 = e23 * e3/e9;
	double    e188 = 0.5 * geno_f;
	double    e190 = e4 * e7/e153;
	double    e191 = 1 - e188;
	double    e192 = e42 * e7;
	double    e193 = pow(e121,2);
	double    e194 = e94/e8;
	double    e195 = e95/e8;
	double    e196 = pow(e72,2);
	double    e198 = 0.5 * e127 - 0.5 * e122;
	double    e200 = tau * e3/e9;
	double    e203 = e177 * e20;
	double    e204 = e178 * e20;
	double    e206 = 2 * ((1 - e131) * e25) + m0;
	double    e207 = 2 * ((2 - e44) * e25);
	double    e208 = e64 * e130;
	double    e209 = e42 * e3;
	double    e210 = pow(e8,2);
	double    e211 = e162/e9;
	double    e212 = e163/e9;
	double    e214 = e173 * e3/e23;
	double    e217 = e129 * e174 * e3/e23;
	double    e218 = e52/e9;
	double    e219 = e46/e8;
	double    e221 = e33 * e4/e8;
	double    e222 = e33/e11;
	double    e226 = e181 * e48 * e49 * e3/e43;
	double    e227 = 0.5 * e185;
	double    e228 = 1 - e44;
	double    e229 = 2 * e187;
	double    e230 = t0 * e87;
	double    e232 = t0 * e130/e9;
	double    e233 = ((e182 * e108/2 + (e183 * e50 * e4 + 1 - e41)/e13)/e3 + 2 * e212)/e13;
	double    e234 = ((e182 * e114/2 + (e198 * e4 + 0.5)/e13)/e3 + 2 * e211)/e13;
	double    e235 = ((e182/2 - 2 * (e4 * e3/e8))/e9 + e183/e13)/e13;
	double    e237 = ((e36 * e125/2 + e183 * e7)/e29 - (e27 * (r0 - (e125/2 + e12)) * e9/e128 + 2 * e190))/e13 + (2 * e192 + 2 * ((e152 - e194) * e4)) * e14/e9;
	double    e239 = ((e36 * e126/2 + e198 * e7)/e29 - (e27 * (0.5 * e13 - e126/2) * e9/e128 + e190))/e13 + (e192 + 2 * ((1 + e149 - 0.5 * e194) * e4)) * e14/e9;
	double    e240 = ((e37 + (e34/e128 - 4 * e122) * e3)/e9 + (2 * e127 - 2 * e122)/e13)/e13;
	double    e241 = pow(e120,2);
	double    e242 = (e108 * e23 - e176 * e11 * e20)/e43;
	double    e245 = (e114 * e23 - e171 * e11 * e20)/e43;
	double    e249 = e171 * e20;
	double    e251 = e162 * e15 + 2 * (e144 * e4/e13);
	double    e253 = e163 * e15 + 2 * (e140 * e4/e13);
	double    e255 = e173 * e11 + e227;
	double    e256 = e129 * e46;
	double    e259 = e174 * e11 + e185;
	double    e260 = e176 * e20;
	double    e261 = (e150 - e195) * e4;
	double    e262 = (e151 - 2 * e195) * e4;
	double    e266 = 0.5 - e214;
	double    e267 = 0.5 * e130;
	double    e270 = 0.5 * tau - 0.5 * t0;
	double    e271 = 1 - (e217 + e99);
	double    e272 = e206 - e21;
	double    e273 = e207 - e98;
	double    e274 = 2 * e130;
	double    e276 = 2 * e200 - (e229 - e178 * e11 * e20)/e43;
	double    e277 = e200 - (e187 - e177 * e11 * e20)/e43;
	double    e278 = e234 + (e86 + 2 * e251) * e14/e9;
	double    e280 = e208 * e7/e9;
	double    e282 = e208 * e3/e9;
	double    e295 = (e162 * e46 - 4 * (e171 * e4 * e20 * e3/e60)) * e25/e9 + 2 * (((0.5 * e62 - e245) * e33 + (e206 - e68) * e4/e9) * e48 * e49 * e3/e43);
	double    e298 = (e163 * e46 - 4 * (e176 * e4 * e20 * e3/e60)) * e25/e9 + 2 * ((e50 * (e207 - (e98 + e59)) * e4/e9 + e33 * (tau * e50 * e4/e9 - e242)) * e48 * e49 * e3/e43);
	double    e302 = e92/e196;
	double    e303 = e57 * e161;
	double    e313 = e249 * e9;
	double    e318 = e173 * e9;
	double    e320 = e228 * e25/e8;
	double    e322 = e71 * e14/e9;
	double    e323 = e52/e210;
	double    e324 = e218 + 2 * e226;
	double    e325 = e174 * e9;
	double    e327 = (2 * (e42 * e50 * e4) + 2 * e253) * e14/e9;
	double    e329 = (2 * e209 + 2 * e261) * e14/e9;
	double    e330 = e260 * e9;
	double    e331 = e203 * e9;
	double    e335 = e204 * e9;
	double    e337 = (2 * e262 + 4 * e209) * e14/e9;
	double    e338 = e226 + 0.5 * e218;
	double    e341 = (4 * (e4 * e48 * e49/e11) - 4) * e3/e23 + 1;
	double    e343 = e3 * e14/e9;
	double    e345 = 2 * ((((e65/2 - e61) * e25 - e70)/e8 - e232) * e4 - ((e222 + e123 + e98) * e27 + e221)/e11) + e230;
	double    e347 = 2 * (((e270 * e33 - e267) * e7 * e3/e9 + (0.5 - e255 * e3/e43) * e33) * e48 * e49/e43) - (0.5 * (e46 * e7/e8) + 4 * (e266/e23)) * e4 * e25/e9;
	double    e349 = 2 * (((e33 * e14 - e130) * e129 * e7 * e3/e9 + (1 - (e259 * e129 * e3/e43 + e99)) * e33) * e48 * e49/e43) - (e256 * e7/e8 + 4 * (e271/e23)) * e4 * e25/e9;
	double    e351 = 2 * ((e33 * e276 - (e274 + 2 * e59) * e3/e9) * e48 * e49/e43) - (2 * e219 + 4 * (e204/e60)) * e4 * e25/e9;
	double    e353 = 2 * ((e33 * e277 - (e130 + e59) * e3/e9) * e48 * e49/e43) - (e219 + 4 * (e203/e60)) * e4 * e25/e9;
	double    e355 = e58 * e7/e153;
	double    e357 = e58 * e3/e153;
	double    e359 = tau * e7/e9;
	
	arma::vec df = arma::zeros<arma::vec>(12);
	
	df(0) = ((2 * ((((e322 - e242) * e87 + 2 * (((e233 + e42 * 
        e176 * e20/e23 + e327) * e20/e23 + tau * e163/e9) * 
        e33 + e64 * e50 * e273 * e4/e9) + 2 * ((e75 * 
        e163 - e50 * (4 * e320 + t0 * e273/e9) * e5)/e9 - 
        ((e50 * e273 * e4/e9 - e108 * e33/e11) * e27 + 
            e163 * e33/e9)/e11)) * e48 * e49 * e3/e43 - 
        (((4 * (e78 * e176) + 8 * ((e233 + 2 * (e253 * e14/e9)) * 
            e4 + e42 * ((2 * e322 + 2 * (e260/e23)) * e4 + 
            1 - e41))) * e20 * e3/e60 + 4 * (e228 * e50 * 
            e46 * e4/e8))/2 + 2 * (e92 * e50 * e4/e196)) * 
            e25/e9) * e57 + e298 * e103) - (2 * (((e233 + 
        e327) * e9 + e42 * (e76 + 2 * (e176 * e23 * e20 * 
        e9/e97))) * e20/e97) + 2 * (e50/e121 - (e50 * 
        e23 * e4/e9 - e330) * e4/e193))) * e3 + 0.5 * 
        e212 - ((2 * (e298 * e57) - 2 * ((e76 + e330/e23)/e23)) * 
        e3 + 0.5 * e108) * e161/e120) * e50/e120;
	
	df(1) = geno_a * 
        ((2 * ((((e270 * e4/e9 - e245) * e87 + 2 * (((e278 + 
            e42 * e171 * e20/e23) * e20/e23 + tau * e162/e9) * 
            e33 + e64 * e272 * e4/e9) + 2 * ((e75 * e162 - 
            (2 * e320 + t0 * e272/e9) * e5)/e9 - ((e272 * 
            e4/e9 - e114 * e33/e11) * e27 + e162 * e33/e9)/e11)) * 
            e48 * e49 * e3/e43 - (e92 * e4/e196 + ((4 * 
            (e171 * e78) + 8 * ((e234 + 2 * (e251 * e14/e9)) * 
            e4 + e42 * ((e4 * e14/e9 + 2 * (e249/e23)) * 
            e4 + 0.5))) * e20 * e3/e60 + 2 * (e228 * e46 * 
            e4/e8))/2) * e25/e9) * e57 + e295 * e103) - 
            (2 * ((e278 * e9 + e42 * (e77 + 2 * (e171 * 
                e23 * e20 * e9/e97))) * e20/e97) + 2 * 
                (0.5/e121 - (0.5 * (e23 * e4/e9) - e313) * 
                  e4/e193))) * e3 + 0.5 * e211 - ((2 * (e295 * 
            e57) - 2 * ((e313/e23 + e77)/e23)) * e3 + 0.5 * 
            e114) * e161/e120) * e50/e120;
	
	df(2) = ((2 * ((e229 - 
        e335)/e193) - 1/e153) * e4 + 2 * ((e103 * e351 + 
        ((e87 * e276 + 2 * (((e240 + e42 * e178 * e20/e23 + 
            e337) * e20/e23 - 2 * e357) * e33 - 2 * e282) - 
            (2 * (((2 * e75 - 8 * e113)/e8 - 2 * e232) * 
                e4 - (e27 * (2 * e222 + e274) + 2 * e221)/e11) + 
                2 * e230) * e3/e9) * e48 * e49/e43 - (((8 * 
            ((e240 + e42 * (2 * (e204/e23) + 4 * e343) + 
                2 * (e262 * e14/e9)) * e20/e60) - 8 * e323) * 
            e4 + 4 * (e78 * e178 * e20/e60))/2 + 4 * e302) * 
            e25/e9) * e57) * e3) - (e161 * (1/e9 + 2 * 
        (e57 * e351 * e3) - 2 * ((e335/e23 + e132)/e23))/e120 + 
        2 * (((e240 + e337) * e9 + e42 * (2 * (e178 * e23 * 
            e20 * e9/e97) + e132)) * e20/e97))) * e50 * 
        (1 - 0.5 * geno_b) * e6 * e3/e120;
	
	df(3) = geno_b * ((2 * 
        ((e187 - e331)/e193) - 0.5/e153) * e4 + 2 * ((e103 * 
        e353 + ((e87 * e277 + 2 * (((e235 + e42 * e177 * 
        e20/e23 + e329) * e20/e23 - e357) * e33 - e282) - 
        e345 * e3/e9) * e48 * e49/e43 - (((8 * ((e235 + 
        e42 * (2 * (e203/e23) + 2 * e343) + 2 * (e261 * 
        e14/e9)) * e20/e60) - 4 * e323) * e4 + 4 * (e78 * 
        e177 * e20/e60))/2 + 2 * e302) * e25/e9) * e57) * 
        e3) - (e161 * (0.5/e9 + 2 * (e57 * e353 * e3) - 
        2 * ((e331/e23 + e93)/e23))/e120 + 2 * (((e235 + 
        e329) * e9 + e42 * (e93 + 2 * (e177 * e23 * e20 * 
        e9/e97))) * e20/e97))) * e50 * e6 * e3/e120;
		
	df(4) = ((e129 * (2 * ((e325 + e185) * e4/e193) - 2 * 
            ((e237 * e9 + e42 * (e51 - 2 * (e174 * e23 * 
                e9/e97))) * e20/e97)) + 2 * ((((e87 * (e359 - 
            e259/e43) + 2 * (((e237 - e42 * e174/e23) * 
            e20/e23 - e355) * e33 - e280) - e345 * e7/e9) * 
            e3 + e84 + e85) * e129 * e48 * e49/e43 - (((8 * 
            ((e237 * e129 * e3 + e42 * (1 - (e99 + 2 * e217))) * 
                e20/e60) - 4 * (e256 * e4 * e7/e210)) * 
            e4 + 4 * (e271 * e78/e23))/2 + 2 * (e92 * e129 * 
            e7/e196)) * e25/e9) * e57 + e103 * e349)) * 
            e3 + e129 * (e147 - ((0.5 * (e7/e8) + 2/e23) * 
            e4/e9 + e146)) - ((e104 + e79 - e105) * e129 + 
            (2 * (e57 * e349) - 2 * ((e51 - e325/e23) * 
                e129/e23)) * e3) * e161/e120) * e50/e120;
				
	df(5) = geno_q * ((2 * (((((0.5 * e359 - e255/e43) * 
            e87 + 2 * (((e239 - e42 * e173/e23) * e20/e23 - 
            0.5 * e355) * e33 - 0.5 * e280) - (0.5 * e230 + 
            2 * (((0.5 * e75 - e123)/e8 - 0.5 * e232) * e4 - 
                (e27 * (0.5 * e222 + e267) + 0.5 * e221)/e11)) * 
            e7/e9) * e3 + 0.5 * e87) * e48 * e49/e43 - 
            (e92 * e7/e196 + ((8 * ((e239 * e3 + e42 * 
                (0.5 - 2 * e214)) * e20/e60) - 2 * (e52 * 
                e7/e210)) * e4 + 4 * (e266 * e78/e23))/2) * 
                e25/e9) * e57 + e103 * e347) + 2 * ((e318 + 
            e227) * e4/e193) - 2 * ((e239 * e9 + e42 * 
            (e79 - 2 * (e173 * e23 * e9/e97))) * e20/e97)) * 
            e3 + 0.5 * e156 - (e161 * ((2 * (e57 * e347) - 
            2 * ((e79 - e318/e23)/e23)) * e3 + 0.25 * e51 + 
            0.5 * e106)/e120 + 0.25 * e190)) * e50/e120;
			
	df(6) = e50 * (2 * ((e92 * e191/e72 + (2 * (e64 * e191 * 
            e181) + 2 * ((0.5 * (e191 * e65) - 2 * (0.5 * 
            (t0 * e191 * e181) + 1 - e188)) * e4/e9 - e27 * 
            e191 * e181/e11)) * e48 * e49 * e3/e43) * 
            e57 + e103 * e324 * e191) - 2 * (e57 * e324 * 
            e161 * e191/e120)) * e3/e120;
			
	df(7) = geno_f * 
            e50 * (2 * (((e64 * e181 + 2 * ((0.25 * e65 - 
            2 * (0.25 * (t0 * e181) + 0.5)) * e4/e9 - 0.5 * 
            (e27 * e181/e11))) * e48 * e49 * e3/e43 + 
            e92/(4 * e9)) * e57 + e103 * e338) - 2 * (e303 * 
            e338/e120)) * e3/e120;
			
	df(8) = e50 * (2 * (e303 * 
            e341 * e4/e120) - 2 * (e103 * e341 * e4 + (e92/2 + 
            (2 * (0.5 * e65 - (2 * (e27/e11) + 2 * (t0 * e4/e9)) * 
                e4) + 4 * (e64 * e4)) * e4 * e48 * e49 * 
                e3/e43) * e57)) * e3/(e120 * e9);
				
	df(9) = -(e161 * e50 * (1 - 0.5 * geno_mu) * e100/e241);
	df(10) = -(0.5 * (geno_mu * e161 * e50 * e100/e241));
	df(11) = -(tau * e161 * e50 * e100 * e117/e241);
	
	return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_a2_abqff1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
double    e3 = geno_q * (q2 - q0)/2 + q0;
double    e4 = a0 + geno_a * (a2 - a0)/2;
double    e5 = pow(e4,2);
double    e6 = b0 + geno_b * (b2 - b0)/2;
double    e7 = pow(e6,2);
double    e8 = e5 + 2 * (e7 * e3);
double    e9 = sqrt(e8);
double    e10 = 2 * e3;
double    e11 = e4 + e9;
double    e12 = e11/e10;
double    e13 = r0 - e12;
double    e14 = tau - t0;
double    e15 = e10 + 2 * (e9/e13);
double    e18 = geno_f * (f2 - f0)/2;
double    e20 = exp(2 * (e9 * e14));
double    e21 = f0 + e18;
double    e23 = e10 - e15 * e20;
double    e24 = e4/e9;
double    e25 = e21 - f1;
double    e26 = e24 + 1;
double    e27 = e5 * e25;
double    e28 = e21 - m0;
double    e29 = e3 * e13;
double    e33 = 2 * (e27/e9) - 2 * (e28 * e9);
double    e34 = 2 * e29;
double    e35 = e26 * e9;
double    e37 = e35/e34 + e24;
double    e40 = e15 * e4 * e14/e9;
double    e41 = e37/e13;
double    e42 = e23 * e11;
double    e43 = e41 + e40;
double    e45 = 1 - 4 * (e3/e23);
double    e47 = exp(-(t0 * e9));
double    e48 = exp(tau * e9);
double    e49 = e5/e8;
double    e50 = e7/e9;
double    e51 = e45 * e4;
double    e52 = e43 * e20;
double    e55 = e51 * e25/e9 + 2 * (e33 * e47 * e48 * e3/e42);
double    e56 = tau * e4;
double    e57 = pow(e23,2);
double    e58 = 2 * e11;
double    e59 = t0 * e33;
double    e60 = 0.5 * (e56/e9);
double    e61 = 2 * e49;
double    e63 = e52/e23 + e60;
double    e64 = 4 - e61;
double    e67 = 0.5 * e59 + f0 + e18;
double    e68 = 0.5 * e24;
double    e69 = 0.5 * e50;
double    e70 = 2 * e9;
double    e71 = e26 * e33;
double    e75 = e64 * e25/2 + m0 - e67;
double    e76 = 1 - e49;
double    e81 = e3/e9;
double    e82 = 2 * (e63 * e33);
double    e83 = 2 * (e75 * e4/e9 - e71/e58);
double    e84 = e82 + e83;
double    e90 = e76 * e45 - 8 * (e43 * e4 * e20 * e3/e57);
double    e91 = e15 * e7;
double    e92 = e15 * e3;
double    e93 = 0.5 * geno_q;
double    e95 = pow(-e23,2);
double    e96 = 2 * e13;
double    e98 = e91 * e14/e9;
double    e100 = e92 * e14/e9;
double    e101 = e27/e8;
double    e102 = 0.5 + e68;
double    e103 = exp(tau * theta);
double    e106 = e90 * e25/e70 + e84 * e47 * e48 * e3/e42;
double    e107 = pow(e55,2);
double    e108 = 2 * (e9/e23);
double    e109 = e107 - e108;
double    e111 = e50 - e11/e3;
double    e112 = e69 - e12;
double    e114 = 2 * e101 + 2 * e28;
double    e117 = geno_mu * (mu02 - mu00)/2 + mu00;
double    e120 = e109 * e3 + 0.5 * e11 + e103 * e117;
double    e121 = e23 * e9;
double    e122 = e4/e8;
double    e123 = 1 - e93;
double    e124 = 0.5 * e49;
double    e125 = 2 * e81;
double    e126 = e26/e9;
double    e127 = pow(e34,2);
double    e133 = e102 * e9/e34 + e68;
double    e136 = 2 * ((e112 * e9/e34 + e69)/e13);
double    e137 = 2 * ((e81 + 1/e96)/e13);
double    e138 = 2 * ((1/e13 + e125)/e13);
double    e139 = 2 + 2 * ((e111 * e9/e34 + e50)/e13);
double    e143 = e8 * e9;
double    e144 = 2 * (e52 * e9/e95);
double    e145 = 2 * (e106 * e55);
double    e148 = 0.5 - e124;
double    e149 = e145 - (e4/e121 + e144);
double    e153 = e149 * e3 + 0.25 * e26;
double    e154 = e40 + 2 * (e133/e13);
double    e156 = 1 - (e98 + 1 + e136) * e20;
double    e157 = 2 - (e139 + 2 * e98) * e20;
double    e158 = e137 + 2 * e100;
double    e159 = e138 + 4 * e100;
double    e162 = 2 * (e5/e9) - e70;
double    e164 = e23 * e7/e9;
double    e166 = e23 * e3/e9;
double    e168 = e4 * e7/e143;
double    e169 = 0.5 * geno_f;
double    e170 = e43 * e7;
double    e171 = e41 + 2 * (e26 * e3 * e9/e127);
double    e172 = e126 - e122;
double    e173 = e91/e8;
double    e174 = e92/e8;
double    e176 = 0.5 * e126 - 0.5 * e122;
double    e177 = 1 - e169;
double    e180 = tau * e3/e9;
double    e181 = pow(e121,2);
double    e183 = e158 * e20;
double    e184 = e159 * e20;
double    e185 = pow(e58,2);
double    e186 = pow(e70,2);
double    e188 = 2 * ((1 - e124) * e25) + m0;
double    e189 = e63 * e114;
double    e190 = e43 * e3;
double    e191 = e154 * e20;
double    e192 = e26 * e114;
double    e193 = e71/e185;
double    e194 = pow(e8,2);
double    e195 = e148/e9;
double    e197 = e156 * e3/e23;
double    e200 = e123 * e157 * e3/e23;
double    e201 = e51/e9;
double    e202 = e45/e8;
double    e204 = e33 * e4/e8;
double    e208 = e162 * e47 * e48 * e3/e42;
double    e209 = 0.5 * e164;
double    e210 = 2 * e166;
double    e211 = t0 * e84;
double    e213 = t0 * e114/e9;
double    e214 = ((e171 * e102/2 + (e176 * e4 + 0.5)/e96)/e3 + e195)/e13;
double    e215 = ((e171/2 - e4 * e3/e8)/e9 + e172/e96)/e13;
double    e217 = ((e37 * e111/2 + e172 * e7/2)/e29 - (e168 + 2 * (e26 * (r0 - (e111/2 + e12)) * e9/e127)))/e13 + ((e139 - e173) * e4 + 2 * e170) * e14/e9;
double    e219 = ((e37 * e112/2 + e176 * e7/2)/e29 - (0.5 * e168 + 2 * (e26 * (0.5 * e13 - e112/2) * e9/e127)))/e13 + (e170 + (1 + e136 - 0.5 * e173) * e4) * e14/e9;
double    e220 = ((e41 + (2 * (e35/e127) - 2 * e122) * e3)/e9 + (2 * e126 - 2 * e122)/e96)/e13;
double    e221 = pow(e120,2);
double    e222 = (e102 * e23 - e154 * e11 * e20)/e42;
double    e224 = e148 * e15;
double    e226 = e156 * e11 + e209;
double    e227 = e123 * e45;
double    e230 = e157 * e11 + e164;
double    e231 = (e137 - e174) * e4;
double    e232 = (e138 - 2 * e174) * e4;
double    e236 = 0.5 - e197;
double    e237 = 0.5 * e64;
double    e240 = 0.5 * tau - 0.5 * t0;
double    e241 = 1 - (e200 + e93);
double    e242 = e188 - e21;
double    e244 = 2 * e180 - (e210 - e159 * e11 * e20)/e42;
double    e245 = e180 - (e166 - e158 * e11 * e20)/e42;
double    e246 = e214 + (((e37 + 2 * e133)/e13 + e40) * e4 + e224) * e14/e9;
double    e248 = e189 * e7/e9;
double    e250 = e189 * e3/e9;
double    e263 = (e148 * e45 - 4 * (e154 * e4 * e20 * e3/e57)) * e25/e9 + 2 * (((e60 - e222) * e33 + (e188 - e67) * e4/e9) * e47 * e48 * e3/e42);
double    e264 = e90/e186;
double    e268 = e55 * e153;
double    e273 = (e231 + 2 * e190) * e14/e9;
double    e275 = (e232 + 4 * e190) * e14/e9;
double    e282 = e191 * e9;
double    e287 = e156 * e9;
double    e288 = e51/e194;
double    e289 = e201 + 2 * e208;
double    e290 = e157 * e9;
double    e291 = e183 * e9;
double    e295 = e184 * e9;
double    e296 = e208 + 0.5 * e201;
double    e299 = (4 * (e4 * e47 * e48/e11) - 4) * e3/e23 + 1;
double    e301 = e3 * e14/e9;
double    e303 = 2 * (((e240 * e33 - 0.5 * e114) * e7 * e3/e9 + (0.5 - e226 * e3/e42) * e33) * e47 * e48/e42) - (0.5 * (e45 * e7/e8) + 4 * (e236/e23)) * e4 * e25/e9;
double    e305 = 2 * (((e33 * e14 - e114) * e123 * e7 * e3/e9 + (1 - (e230 * e123 * e3/e42 + e93)) * e33) * e47 * e48/e42) - (e227 * e7/e8 + 4 * (e241/e23)) * e4 *  e25/e9;
double    e307 = 2 * ((e192 + e204)/e58 + 2 * e193 - (((e64/2 - e61) * e25 + m0 - e67)/e8 - 0.5 * e213) * e4) -  e211;
double    e309 = 2 * ((e33 * e244 - (2 * e114 + 2 * e59) * e3/e9) * e47 * e48/e42) - (2 * e202 + 4 * (e184/e57)) * e4 * e25/e9;
double    e311 = 2 * ((e33 * e245 - (e114 + e59) * e3/e9) * e47 * e48/e42) - (e202 + 4 * (e183/e57)) * e4 * e25/e9;
double    e312 = t0 * e162;
double    e314 = e56 * e7/e143;
double    e316 = e56 * e3/e143;
double    e318 = tau * e7/e9;

arma::vec df = arma::zeros<arma::vec>(11);
	
	df(0) = pow(geno_a,2) * (((0.5 * (e23 * e4/e9) - e282) * e4/e181 + 
        2 * ((((e240 * e4/e9 - e222) * e84 + 2 * (((e246 + 
            e43 * e154 * e20/e23) * e20/e23 + 0.5 * (tau * 
            e148/e9)) * e33 + e63 * e242 * e4/e9) + 2 * 
            ((e75 * e148 - ((e26 * e242 * e4 + e148 * e33)/e58 + 
                (e76 * e25/e8 + 0.5 * (t0 * e242/e9)) * 
                  e5))/e9 + 2 * (e26 * e102 * e33/e185))) * 
            e47 * e48 * e3/e42 - (((e51/e8 + 4 * (e191 * 
            e3/e57)) * e76 + 8 * (((e214 + (e224 + 2 * (e133 * 
            e4/e13)) * e14/e9) * e4 + e43 * ((e4 * e14/e9 + 
            2 * (e191/e23)) * e4 + 0.5)) * e20 * e3/e57))/2 + 
            e90 * e4/e186) * e25/e9) * e55 + e263 * e106) - 
        (0.5/e121 + 2 * ((e246 * e9 + e43 * (e68 + 2 * (e154 * 
            e23 * e20 * e9/e95))) * e20/e95))) * e3 + 
        0.25 * e195 - ((2 * (e263 * e55) - 2 * ((e282/e23 + 
        e68)/e23)) * e3 + 0.5 * e102) * e153/e120)/e120;
		
	df(1) = geno_a * (((e210 - e295)/e181 - 0.5/e143) * 
            e4 + 2 * ((e106 * e309 + ((e84 * e244 + (2 * 
            ((2 * e192 + 2 * e204)/e58 + 4 * e193 - ((2 * 
                e75 - 4 * e101)/e8 - e213) * e4) - 2 * e211) * 
            e3/e9 + 2 * (((e220 + e43 * e159 * e20/e23 + 
            e275) * e20/e23 - e316) * e33 - 2 * e250)) * 
            e47 * e48/e42 - (((8 * ((e220 + e43 * (2 * (e184/e23) + 
            4 * e301) + e232 * e14/e9) * e20/e57) - 4 * 
            e288) * e4 + 4 * (e76 * e159 * e20/e57))/2 + 
            4 * e264) * e25/e9) * e55) * e3) - (e153 * 
            (1/e9 + 2 * (e55 * e309 * e3) - 2 * ((e295/e23 + 
                e125)/e23))/e120 + 2 * (((e220 + e275) * 
            e9 + e43 * (2 * (e159 * e23 * e20 * e9/e95) + 
            e125)) * e20/e95))) * (1 - 0.5 * geno_b) * e6 * 
            e3/e120;
			
	df(2) = geno_a * geno_b * (((e166 - e291)/e181 - 
            0.25/e143) * e4 + 2 * ((e106 * e311 + ((e84 * 
            e245 + e307 * e3/e9 + 2 * (((e215 + e43 * e158 * 
            e20/e23 + e273) * e20/e23 - 0.5 * e316) * e33 - 
            e250)) * e47 * e48/e42 - (((8 * ((e215 + e43 * 
            (2 * (e183/e23) + 2 * e301) + e231 * e14/e9) * 
            e20/e57) - 2 * e288) * e4 + 4 * (e76 * e158 * 
            e20/e57))/2 + 2 * e264) * e25/e9) * e55) * 
            e3) - (e153 * (0.5/e9 + 2 * (e55 * e311 * e3) - 
            2 * ((e291/e23 + e81)/e23))/e120 + 2 * (((e215 + 
            e273) * e9 + e43 * (e81 + 2 * (e158 * e23 * 
            e20 * e9/e95))) * e20/e95))) * e6 * e3/e120;
			
	df(3) = geno_a * ((((e290 + e164) * e4/e181 - 2 * ((e217 * 
            e9 + e43 * (e50 - 2 * (e157 * e23 * e9/e95))) * 
            e20/e95)) * e123 + 2 * ((((e84 * (e318 - e230/e42) + 
            e307 * e7/e9 + 2 * (((e217 - e43 * e157/e23) * 
            e20/e23 - 0.5 * e314) * e33 - e248)) * e3 + 
            e82 + e83) * e123 * e47 * e48/e42 - (((8 * 
            ((e217 * e123 * e3 + e43 * (1 - (e93 + 2 * e200))) * 
                e20/e57) - 2 * (e227 * e4 * e7/e194)) * 
            e4 + 4 * (e241 * e76/e23))/2 + 2 * (e90 * e123 * 
            e7/e186)) * e25/e9) * e55 + e106 * e305)) * 
            e3 + e123 * (e145 - ((0.25 * (e7/e8) + 1/e23) * 
            e4/e9 + e144)) - ((e107 + e69 - e108) * e123 + 
            (2 * (e55 * e305) - 2 * ((e50 - e290/e23) * 
                e123/e23)) * e3) * e153/e120)/e120;
				
	df(4) = geno_a * 
            geno_q * (((e287 + e209) * e4/e181 + 2 * (((((0.5 * 
            e318 - e226/e42) * e84 + (2 * (e193 + (0.5 * 
            e192 + 0.5 * e204)/e58 - ((0.5 * e75 - e101)/e8 - 
            0.25 * e213) * e4) - 0.5 * e211) * e7/e9 + 2 * 
            (((e219 - e43 * e156/e23) * e20/e23 - 0.25 * 
                e314) * e33 - 0.5 * e248)) * e3 + 0.5 * e84) * 
            e47 * e48/e42 - (e90 * e7/e186 + ((8 * ((e219 * 
            e3 + e43 * (0.5 - 2 * e197)) * e20/e57) - e51 * 
            e7/e194) * e4 + 4 * (e236 * e76/e23))/2) * 
            e25/e9) * e55 + e106 * e303) - 2 * ((e219 * 
            e9 + e43 * (e69 - 2 * (e156 * e23 * e9/e95))) * 
            e20/e95)) * e3 + 0.5 * e149 - (e153 * ((2 * 
            (e55 * e303) - 2 * ((e69 - e287/e23)/e23)) * 
            e3 + 0.25 * e50 + 0.5 * e109)/e120 + 0.125 * 
            e168))/e120;
			
	df(5) = geno_a * (2 * ((e90 * e177/e70 + 
            (2 * (e63 * e177 * e162) + 2 * (((e237 - 0.5 * 
                e312) * e177 + e169 - 1) * e4/e9 - e26 * 
                e177 * e162/e58)) * e47 * e48 * e3/e42) * 
            e55 + e106 * e289 * e177) - 2 * (e55 * e289 * 
            e153 * e177/e120)) * e3/e120;
			
	df(6) = geno_a * 
            geno_f * (2 * (((e63 * e162 + 2 * ((0.5 * (e237 - 
            1) - 0.25 * e312) * e4/e9 - e26 * e162/(4 * 
            e11))) * e47 * e48 * e3/e42 + e90/(4 * e9)) * 
            e55 + e106 * e296) - 2 * (e268 * e296/e120)) * 
            e3/e120;
	
	df(7) = geno_a * (2 * (e268 * e299 * e4/e120) - 
            2 * (e106 * e299 * e4 + (e90/2 + (2 * (e237 - 
                (e26/e11 + t0 * e4/e9) * e4) + 4 * (e63 * 
                e4)) * e4 * e47 * e48 * e3/e42) * e55)) * 
            e3/(e120 * e9);
			
	df(8) = -(geno_a * e153 * (1 - 0.5 * geno_mu) * e103/e221);
			
	df(9) = -(0.5 * (geno_a * geno_mu * e153 * e103/e221));
			
	df(10) = -(geno_a * tau * e153 * e103 * e117/e221);
    
	return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_b0_bqff1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
	double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
    double e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e6 = b0 + geno_b * (b2 - b0)/2;
	double    e7 = pow(e6,2);
	double    e8 = e7 * e3;
	double    e9 = e5 + 2 * e8;
	double    e10 = sqrt(e9);
	double    e11 = 2 * e3;
	double    e12 = e4 + e10;
	double    e13 = e12/e11;
	double    e14 = r0 - e13;
	double    e15 = tau - t0;
	double    e18 = geno_f * (f2 - f0)/2;
	double    e19 = e11 + 2 * (e10/e14);
	double    e21 = exp(2 * (e10 * e15));
	double    e22 = f0 + e18;
	double    e24 = e11 - e19 * e21;
	double    e25 = e22 - f1;
	double    e26 = e5 * e25;
	double    e27 = e22 - m0;
	double    e31 = 2 * (e26/e10) - 2 * (e27 * e10);
	double    e32 = e3/e10;
	double    e33 = e19 * e3;
	double    e35 = e33 * e15/e10;
	double    e36 = 4 * e32;
	double    e37 = 4 * e35;
	double    e39 = (2/e14 + e36)/e14 + e37;
	double    e40 = e24 * e12;
	double    e42 = 1 - 4 * (e3/e24);
	double    e44 = exp(-(t0 * e10));
	double    e45 = exp(tau * e10);
	double    e46 = e39 * e21;
	double    e47 = 2 * e10;
	double    e48 = e7/e10;
	double    e49 = e42 * e4;
	double    e50 = 0.5 * geno_q;
	double    e53 = e49 * e25/e10 + 2 * (e31 * e44 * e45 * e3/e40);
	double    e55 = tau * e3/e10;
	double    e56 = e31/e12;
	double    e57 = 2 * e55;
	double    e58 = t0 * e31;
	double    e59 = pow(e24,2);
	double    e60 = pow(e47,2);
	double    e61 = 1 - e50;
	double    e63 = e46/e24 + e57;
	double    e64 = 2 * e56;
	double    e71 = e64 + 4 * (0.5 * e58 + 4 * (e26/e60) + f0 + e18 - m0);
	double    e72 = e42/e9;
	double    e74 = 0.5 * e48;
	double    e75 = 2 * (e3 * e14);
	double    e80 = 2 * (e63 * e31) - 2 * (e71 * e3/e10);
	double    e81 = 4 * e72;
	double    e83 = 8 * (e46/e59);
	double    e84 = e81 + e83;
	double    e87 = e19 * e7 * e15/e10;
	double    e88 = e84 * e4;
	double    e92 = e80 * e44 * e45/e40 - e88 * e25/e47;
	double    e94 = pow(-e24,2);
	double    e96 = 2 * (e26/e9);
	double    e97 = 2 * e27;
	double    e98 = exp(tau * theta);
	double    e99 = pow(e53,2);
	double    e100 = 2 * (e10/e24);
	double    e101 = e99 - e100;
	double    e102 = e96 + e97;
	double    e105 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e108 = e101 * e3 + 0.5 * e12 + e98 * e105;
	double    e109 = e24 * e10;
	double    e111 = e48 - e12/e3;
	double    e112 = e74 - e13;
	double    e113 = 2 * e14;
	double    e114 = 2 * e32;
	double    e121 = 2 * ((e112 * e10/e75 + e74)/e14);
	double    e122 = 2 * (e53 * e92);
	double    e123 = 2 * ((e32 + 1/e113)/e14);
	double    e124 = 2 * ((1/e14 + e114)/e14);
	double    e125 = 2 + 2 * ((e111 * e10/e75 + e48)/e14);
	double    e126 = 0.5 * geno_b;
	double    e127 = 1/e10;
	double    e128 = e122 - 4/e109;
	double    e133 = e128 * e3 + e127 - 2 * (e46 * e10/e94);
	double    e134 = 1 - e126;
	double    e135 = e36 + 4/e14;
	double    e138 = pow(e3,2);
	double    e139 = 1 - (e87 + 1 + e121) * e21;
	double    e140 = 2 - (e125 + 2 * e87) * e21;
	double    e141 = e123 + 2 * e35;
	double    e142 = e124 + e37;
	double    e143 = e5/e10;
	double    e145 = e24 * e7/e10;
	double    e147 = e24 * e3/e10;
	double    e153 = 0.5 - 0.5 * (e8/e9);
	double    e154 = 0.5 * geno_f;
	double    e155 = 1 - (e61 * e7 * e3/e9 + e50);
	double    e156 = e33/e9;
	double    e157 = e138/e9;
	double    e158 = 1 - e154;
	double    e159 = e39 * e61;
	double    e160 = e9 * e10;
	double    e161 = e141 * e21;
	double    e162 = e142 * e21;
	double    e163 = e143 - e10;
	double    e165 = 2 * e143 - e47;
	double    e166 = pow(e109,2);
	double    e167 = e49/e10;
	double    e168 = e71/e9;
	double    e170 = e26/pow(e47,3);
	double    e171 = 0.5 * e145;
	double    e172 = 2 * e147;
	double    e173 = t0 * e80;
	double    e174 = t0 * e102;
	double    e175 = pow(e108,2);
	double	e178 = (e39 * e7 + 4 * (e153 * e19 + (1 + e121) * e3)) * e15/e10 + (e112 * e135/e75 + 4 * (e153/e10))/e14;
	double    e179 = e63 * e102;
	double    e181 = (e111 * e61 * e135/e75 + 4 * (e155/e10))/e14 + (2 * (e159 * e7) + 4 * (e155 * e19 + e61 * e125 * e3)) * e15/e10;
	double    e186 = e139 * e12 + e171;
	double    e189 = e61 * e42 * e7/e9;
	double    e191 = e42 * e7/e9;
	double    e193 = e140 * e12 + e145;
	double    e200 = (2 * e39 + 4 * (e123 - e156)) * e3 * e15;
	double    e202 = (4 * e39 + 4 * (e124 - 2 * e156)) * e3 * e15;
	double    e205 = (4 * e157 - e135/e113)/e14;
	double    e206 = (8 * e157 - e135/e14)/e14;
	double    e207 = 0.5 * e102;
	double    e208 = 2 * e72;
	double    e209 = 2 * e102;
	double    e210 = e57 - (e172 - e142 * e12 * e21)/e40;
	double    e211 = 4 * ((0.5 - e139 * e3/e24)/e24);
	double    e212 = 4 * ((1 - (e61 * e140 * e3/e24 + e50))/e24);
	double    e213 = 4 * (e161/e59);
	double    e214 = 4 * (e162/e59);
	double    e215 = e55 - (e147 - e141 * e12 * e21)/e40;
	double    e223 = e179 * e3/e10;
	double    e227 = e53 * e133;
	double    e231 = (e200 - e205)/e10;
	double    e240 = e39 * e139/e24;
	double    e242 = e159 * e140/e24;
	double    e245 = e39 * e141 * e21/e24;
	double    e248 = e39 * e142 * e21/e24;
	double    e249 = (e202 - e206)/e10;
	double    e257 = e139 * e10;
	double    e258 = e167 + 2 * (e165 * e44 * e45 * e3/e40);
	double    e259 = e140 * e10;
	double    e261 = (2 * ((e56 + e96 + e97)/e12) + 4 * (0.5 * e174 + 16 * e170))/e10 + e168;
	double    e262 = e161 * e10;
	double    e266 = e162 * e10;
	double    e269 = (4 * (e4 * e44 * e45/e12) - 4) * e3/e24 + 1;
	double    e270 = e5/e60;
	double    e272 = 0.5 * e167 + 2 * (e163 * e44 * e45 * e3/e40);
	double    e274 = 2 * ((((0.5 * tau - 0.5 * t0) * e31 - e207) * e7 * e3/e10 + (0.5 - e186 * e3/e40) * e31) * e44 * e45/e40) - (0.5 * e191 + e211) * e4 * e25/e10;
	double    e276 = 2 * (((e31 * e15 - e102) * e61 * e7 * e3/e10 + (1 - (e193 * e61 * e3/e40 + e50)) * e31) * e44 * e45/e40) - (e189 + e212) * e4 * e25/e10;
	double    e278 = 2 * ((e31 * e210 - (e209 + 2 * e58) * e3/e10) * e44 * e45/e40) - (e208 + e214) * e4 * e25/e10;
	double    e280 = 2 * ((e31 * e215 - (e102 + e58) * e3/e10) * e44 * e45/e40) - (e72 + e213) * e4 * e25/e10;
	double    e281 = e173 * e7;
	double    e283 = tau * e7/e10;
	double    e285 = tau * e138/e160;
	
	arma::vec df = arma::zeros<arma::vec>(10);

	df(0) = (((2 * (((e80 * e210 + (2 * (((2 * ((e64 + e209)/e12) + 
        4 * (32 * e170 + e174))/e10 + 2 * e168) * e3) - 
        2 * e173) * e3/e10 + 2 * (((e248 + e249) * e21/e24 - 
        4 * e285) * e31 - 2 * e223)) * e44 * e45/e40 - 
        (8 * ((e249 + 2 * e248) * e21/e59) - (2 * e84 + 
            4 * (e81 + e214)) * e3/e9) * e4 * e25/e47) * 
        e53 + e92 * e278 * e3) + 4 * ((e172 - e266)/e166) - 
        2/e160) * e3 - 2 * ((e39 * (2 * (e142 * e24 * e21 * 
        e10/e94) + e114) + e202 - e206) * e21/e94)) * 
        e134 * e7 + e133 * (1 - (e134 * (e127 + 2 * (e53 * 
        e278 * e3) - 2 * ((e266/e24 + e114)/e24)) * e7 * 
        e3/e108 + e126))) * e134 * e3/e108;
		
	df(1) = geno_b * 
        (((2 * (((e80 * e215 + (2 * (e261 * e3) - e173) * 
            e3/e10 + 2 * (((e231 + e245) * e21/e24 - 2 * 
            e285) * e31 - e223)) * e44 * e45/e40 - (8 * 
            ((e231 + 2 * e245) * e21/e59) - (e81 + 4 * (e208 + 
            e213) + e83) * e3/e9) * e4 * e25/e47) * e53 + 
            e92 * e280 * e3) + 4 * ((e147 - e262)/e166) - 
            1/e160) * e3 - 2 * ((e39 * (e32 + 2 * (e141 * 
            e24 * e21 * e10/e94)) + e200 - e205) * e21/e94)) * 
            e7 + e133 * (0.5 - (0.5/e10 + 2 * (e53 * e280 * 
            e3) - 2 * ((e262/e24 + e32)/e24)) * e7 * e3/e108)) * 
        e134 * e3/e108;
		
	df(2) = ((e61 * (e122 - (e7/e9 + 
        4/e24)/e10) + (2 * ((((e80 * (e283 - e193/e40) - 
        (2 * (e71 - e261 * e7 * e3) + e281)/e10) * e61 + 
        2 * (((e181 - e242) * e21/e24 + 2 * (tau * e155/e10)) * 
            e31 - e63 * e61 * e102 * e7/e10)) * e44 * 
        e45/e40 - (8 * ((e181 - 2 * e242) * e21/e59) - 
        (e61 * e84 * e7 + 4 * (2 * e189 + e212))/e9) * 
        e4 * e25/e47) * e53 + e92 * e276) + 4 * ((e259 + 
        e145) * e61/e166)) * e3 - (((e99 + e74 - e100) * 
        e61 + (2 * (e53 * e276) - 2 * ((e48 - e259/e24) * 
        e61/e24)) * e3) * e133/e108 + 2 * ((e181 * e10 + 
        e39 * (e48 - 2 * (e140 * e24 * e10/e94)) * e61) * 
        e21/e94))) * e3 + e133 * e61) * e134 * e6/e108;
		
	df(3) = geno_q * (((2 * ((((0.5 * e283 - e186/e40) * 
            e80 + 2 * (((e178 - e240) * e21/e24 + 2 * (tau * 
            e153/e10)) * e31 - 0.5 * (e179 * e7/e10)) - 
            (0.5 * e281 + 2 * (0.5 * e71 - ((2 * ((0.5 * e56 + 
                e207)/e12) + 4 * (0.25 * e174 + 8 * e170))/e10 + 
                0.5 * e168) * e7 * e3))/e10) * e44 * e45/e40 - 
            (8 * ((e178 - 2 * e240) * e21/e59) - (e84 * 
                e7/2 + 4 * (e191 + e211))/e9) * e4 * e25/e47) * 
            e53 + e92 * e274) + 4 * ((e257 + e171)/e166)) * 
            e3 + 0.5 * e128 - (e133 * ((2 * (e53 * e274) - 
            2 * ((e74 - e257/e24)/e24)) * e3 + 0.25 * e48 + 
            0.5 * e101)/e108 + 0.5 * (e7/e160) + 2 * ((e178 * 
            e10 + e39 * (e74 - 2 * (e139 * e24 * e10/e94))) * 
            e21/e94))) * e3 + 0.5 * e133) * e134 * e6/e108;
			
	df(4) = e134 * (2 * (e53 * ((2 * (e63 * e158 * e165) - 
            2 * ((2 * (e158 * e165/e12) + 4 * ((0.5 * (t0 * 
                e165) + 4 * e270) * e158 + 1 - e154)) * e3/e10)) * 
            e44 * e45/e40 - e158 * e84 * e4/e47) + e258 * 
            e92 * e158) - 2 * (e53 * e258 * e133 * e158/e108)) * 
            e6 * e138/e108;
			
	df(5) = geno_f * e134 * (2 * (e53 * 
            ((2 * (e63 * e163) - 2 * ((2 * (e163/e12) + 4 * 
                (0.5 + 0.5 * (t0 * e163) + 2 * e270)) * e3/e10)) * 
                e44 * e45/e40 - e88/(4 * e10)) + e92 * 
            e272) - 2 * (e227 * e272/e108)) * e6 * e138/e108;
			
	df(6) = e134 * (2 * (e227 * e269/e108) + 2 * (e53 * 
            ((2 * ((4 * (4/e60 + t0/e10) + 4/(e12 * e10)) * 
                e3) - 4 * e63) * e4 * e44 * e45/e40 + e84/2) - 
            e92 * e269)) * e4 * e6 * e138/(e108 * e10);
			
	df(7) = -(e133 * e134 * (1 - 0.5 * geno_mu) * e6 * e98 * e3/e175);
			
	df(8) = -(0.5 * (geno_mu * e133 * e134 * e6 * e98 * e3/e175));
			
	df(9) = -(tau * e133 * e134 * e6 * e98 * e105 * e3/e175);
	
	return(Rcpp::wrap(df));
	
}

RcppExport SEXP d_f_i1_b2_bqff1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
	double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e6 = b0 + geno_b * (b2 - b0)/2;
	double    e7 = pow(e6,2);
	double    e8 = e7 * e3;
	double    e9 = e5 + 2 * e8;
	double    e10 = sqrt(e9);
	double    e11 = 2 * e3;
	double    e12 = e4 + e10;
	double    e13 = e12/e11;
	double    e14 = r0 - e13;
	double    e15 = tau - t0;
	double    e16 = e11 + 2 * (e10/e14);
	double    e19 = geno_f * (f2 - f0)/2;
	double    e21 = exp(2 * (e10 * e15));
	double    e22 = f0 + e19;
	double    e24 = e11 - e16 * e21;
	double    e25 = e22 - f1;
	double    e26 = e5 * e25;
	double    e27 = e22 - m0;
	double    e31 = 2 * (e26/e10) - 2 * (e27 * e10);
	double    e32 = e3/e10;
	double    e33 = e16 * e3;
	double    e36 = 2 * (e33 * e15/e10);
	double    e37 = 2 * e32;
	double    e38 = e24 * e12;
	double    e40 = (1/e14 + e37)/e14 + e36;
	double    e42 = 1 - 4 * (e3/e24);
	double    e44 = exp(-(t0 * e10));
	double    e45 = exp(tau * e10);
	double    e46 = e7/e10;
	double    e47 = e40 * e21;
	double    e48 = 2 * e10;
	double    e49 = 0.5 * geno_q;
	double    e50 = e42 * e4;
	double    e53 = e50 * e25/e10 + 2 * (e31 * e44 * e45 * e3/e38);
	double    e54 = 1 - e49;
	double    e56 = 0.5 * e46;
	double    e57 = 2 * (e3 * e14);
	double    e59 = tau * e3/e10;
	double    e61 = e47/e24 + e59;
	double    e62 = e31/e12;
	double    e63 = pow(e24,2);
	double    e64 = pow(e48,2);
	double    e65 = t0 * e31;
	double    e66 = e42/e9;
	double    e67 = e62 + 2 * (0.5 * e65 + 4 * (e26/e64) + f0 + e19 - m0);
	double    e68 = 2 * e66;
	double    e74 = 2 * (e61 * e31) - 2 * (e67 * e3/e10);
	double    e75 = 8 * (e47/e63);
	double    e76 = e68 + e75;
	double    e79 = e16 * e7 * e15/e10;
	double    e80 = e76 * e4;
	double    e84 = e74 * e44 * e45/e38 - e80 * e25/e48;
	double    e85 = exp(tau * theta);
	double    e87 = pow(-e24,2);
	double    e88 = pow(e53,2);
	double    e90 = e46 - e12/e3;
	double    e91 = e56 - e13;
	double    e92 = 2 * e14;
	double    e93 = 2 * (e10/e24);
	double    e94 = e88 - e93;
	double    e97 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e100 = e94 * e3 + 0.5 * e12 + e85 * e97;
	double    e101 = e24 * e10;
	double    e103 = 2 * (e26/e9);
	double    e104 = 2 * e27;
	double    e109 = 2 * ((e91 * e10/e57 + e56)/e14);
	double    e110 = 2 * ((e32 + 1/e92)/e14);
	double    e111 = 2 + 2 * ((e90 * e10/e57 + e46)/e14);
	double    e113 = 2 * (e53 * e84);
	double    e114 = e103 + e104;
	double    e115 = 0.5/e10;
	double    e116 = e113 - 2/e101;
	double    e122 = e116 * e3 + e115 - 2 * (e47 * e10/e87);
	double    e124 = 1 - (e79 + 1 + e109) * e21;
	double    e125 = 2 - (e111 + 2 * e79) * e21;
	double    e126 = e110 + e36;
	double    e127 = e5/e10;
	double    e128 = e37 + 2/e14;
	double    e134 = e24 * e7/e10;
	double    e136 = 0.5 - 0.5 * (e8/e9);
	double    e137 = 1 - (e54 * e7 * e3/e9 + e49);
	double    e138 = pow(e3,2);
	double    e139 = 0.5 * geno_f;
	double    e140 = 1 - e139;
	double    e141 = e40 * e54;
	double    e142 = e126 * e21;
	double    e143 = e127 - e10;
	double    e145 = 2 * e127 - e48;
	double    e146 = e50/e10;
	double    e148 = e24 * e3/e10;
	double    e149 = 0.5 * e134;
	double    e150 = pow(e100,2);
	double    e153 = (e40 * e7 + 2 * (e136 * e16 + (1 + e109) * e3)) * e15/e10 + (e91 * e128/e57 + 2 * (e136/e10))/e14;
	double    e155 = (e90 * e54 * e128/e57 + 2 * (e137/e10))/e14 + (2 * (e141 * e7) + 2 * (e137 * e16 + e54 * e111 * e3)) * e15/e10;
	double    e157 = pow(e101,2);
	double    e158 = e9 * e10;
	double    e162 = e124 * e12 + e149;
	double    e165 = e54 * e42 * e7/e9;
	double    e167 = e42 * e7/e9;
	double    e169 = e125 * e12 + e134;
	double    e172 = (2 * e40 + 2 * (e110 - e33/e9)) * e3 * e15;
	double    e174 = (2 * (e138/e9) - e128/e92)/e14;
	double    e176 = e26/pow(e48,3);
	double    e177 = 0.5 * e114;
	double    e178 = 4 * ((0.5 - e124 * e3/e24)/e24);
	double    e179 = 4 * ((1 - (e54 * e125 * e3/e24 + e49))/e24);
	double    e180 = 4 * (e142/e63);
	double    e181 = t0 * e114;
	double    e182 = e59 - (e148 - e126 * e12 * e21)/e38;
	double    e186 = e61 * e114;
	double    e193 = e53 * e122;
	double    e198 = e40 * e124/e24;
	double    e200 = e141 * e125/e24;
	double    e203 = e40 * e126 * e21/e24;
	double    e204 = (e172 - e174)/e10;
	double    e209 = (e62 + e103 + e104)/e12 + 2 * (0.5 * e181 + 16 * e176);
	double    e217 = e124 * e10;
	double    e218 = e146 + 2 * (e145 * e44 * e45 * e3/e38);
	double    e219 = e125 * e10;
	double    e220 = e142 * e10;
	double    e223 = (4 * (e4 * e44 * e45/e12) - 4) * e3/e24 + 1;
	double    e224 = e5/e64;
	double    e226 = 0.5 * e146 + 2 * (e143 * e44 * e45 * e3/e38);
	double    e228 = 2 * ((((0.5 * tau - 0.5 * t0) * e31 - e177) * e7 * e3/e10 + (0.5 - e162 * e3/e38) * e31) * e44 * e45/e38) - (0.5 * e167 + e178) * e4 * e25/e10;
	double    e230 = 2 * (((e31 * e15 - e114) * e54 * e7 * e3/e10 + (1 - (e169 * e54 * e3/e38 + e49)) * e31) * e44 * e45/e38) - (e165 + e179) * e4 * e25/e10;
	double    e232 = 2 * ((e31 * e182 - (e114 + e65) * e3/e10) * e44 * e45/e38) - (e66 + e180) * e4 * e25/e10;
	double    e233 = t0 * e74;
	double    e235 = tau * e7/e10;
	
	arma::vec df = arma::zeros<arma::vec>(9);

	df(0) = pow(geno_b,2) * (((2 * ((((2 * ((e209/e10 + e67/e9) * 
        e3) - e233) * e3/e10 + e74 * e182 + 2 * (((e203 + 
        e204) * e21/e24 - tau * e138/e158) * e31 - e186 * 
        e3/e10)) * e44 * e45/e38 - (8 * ((e204 + 2 * e203) * 
        e21/e63) - (e68 + 2 * (e68 + e180) + e75) * e3/e9) * 
        e4 * e25/e48) * e53 + e84 * e232 * e3) + 2 * ((e148 - 
        e220)/e157) - 0.5/e158) * e3 - 2 * ((e40 * (e32 + 
        2 * (e126 * e24 * e21 * e10/e87)) + e172 - e174) * 
        e21/e87)) * e7 + e122 * (0.5 - (e115 + 2 * (e53 * 
        e232 * e3) - 2 * ((e220/e24 + e32)/e24)) * e7 * 
        e3/e100)) * e3/e100;
	
	df(1) = geno_b * ((e54 * (e113 - 
        (0.5 * (e7/e9) + 2/e24)/e10) + (2 * (((e54 * e74 * 
        (e235 - e169/e38) + 2 * (((e155 - e200) * e21/e24 + 
        tau * e137/e10) * e31 - e61 * e54 * e114 * e7/e10) - 
        (2 * (e67 * e137 - e209 * e54 * e7 * e3/e10) + 
            t0 * e54 * e74 * e7)/e10) * e44 * e45/e38 - 
        (8 * ((e155 - 2 * e200) * e21/e63) - (e54 * e76 * 
            e7 + 2 * (2 * e165 + e179))/e9) * e4 * e25/e48) * 
        e53 + e84 * e230) + 2 * ((e219 + e134) * e54/e157)) * 
        e3 - (((e88 + e56 - e93) * e54 + (2 * (e53 * e230) - 
        2 * ((e46 - e219/e24) * e54/e24)) * e3) * e122/e100 + 
        2 * ((e155 * e10 + e40 * (e46 - 2 * (e125 * e24 * 
            e10/e87)) * e54) * e21/e87))) * e3 + e122 * 
        e54) * e6/e100;
		
	df(2) = geno_b * geno_q * (((2 * ((((0.5 * 
        e235 - e162/e38) * e74 + 2 * (((e153 - e198) * 
        e21/e24 + tau * e136/e10) * e31 - 0.5 * (e186 * 
        e7/e10)) - (0.5 * (e233 * e7) + 2 * (e67 * e136 - 
        ((0.5 * e62 + e177)/e12 + 2 * (0.25 * e181 + 8 * 
            e176)) * e7 * e3/e10))/e10) * e44 * e45/e38 - 
        (8 * ((e153 - 2 * e198) * e21/e63) - (e76 * e7/2 + 
            2 * (e167 + e178))/e9) * e4 * e25/e48) * e53 + 
        e84 * e228) + 2 * ((e217 + e149)/e157)) * e3 + 
        0.5 * e116 - (e122 * ((2 * (e53 * e228) - 2 * ((e56 - 
        e217/e24)/e24)) * e3 + 0.25 * e46 + 0.5 * e94)/e100 + 
        0.25 * (e7/e158) + 2 * ((e153 * e10 + e40 * (e56 - 
        2 * (e124 * e24 * e10/e87))) * e21/e87))) * e3 + 
        0.5 * e122) * e6/e100;
		
	df(3) = geno_b * (2 * (e53 * 
        ((2 * (e61 * e140 * e145) - 2 * ((e140 * e145/e12 + 
            2 * ((0.5 * (t0 * e145) + 4 * e224) * e140 + 1 - 
                e139)) * e3/e10)) * e44 * e45/e38 - e140 * 
            e76 * e4/e48) + e218 * e84 * e140) - 2 * (e53 * 
        e218 * e122 * e140/e100)) * e6 * e138/e100;
		
	df(4) = geno_b * 
        geno_f * (2 * (e53 * ((2 * (e61 * e143) - 2 * ((e143/e12 + 
        2 * (0.5 + 0.5 * (t0 * e143) + 2 * e224)) * e3/e10)) * 
        e44 * e45/e38 - e80/(4 * e10)) + e84 * e226) - 
        2 * (e193 * e226/e100)) * e6 * e138/e100;
		
	df(5) = geno_b * 
        (2 * (e193 * e223/e100) + 2 * (e53 * (e76/2 + (2 * 
            ((2 * (4/e64 + t0/e10) + 2/(e12 * e10)) * e3) - 
            4 * e61) * e4 * e44 * e45/e38) - e84 * e223)) * 
        e4 * e6 * e138/(e100 * e10);
		
	df(6) = -(geno_b * e122 * (1 - 0.5 * geno_mu) * e6 * e85 * e3/e150);
		
	df(7) = -(0.5 * (geno_b * geno_mu * e122 * e6 * e85 * e3/e150));
		
	df(8) = -(geno_b * tau * e122 * e6 * e85 * e97 * e3/e150);
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_q0_qff1mtb_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
	double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = b0 + geno_b * (b2 - b0)/2;
	double    e6 = pow(e5,2);
	double    e7 = pow(e4,2);
	double    e8 = e6 * e3;
	double    e9 = e7 + 2 * e8;
	double    e10 = sqrt(e9);
	double    e11 = 2 * e3;
	double    e12 = e4 + e10;
	double    e13 = e12/e11;
	double    e14 = r0 - e13;
	double    e15 = tau - t0;
	double    e16 = e11 + 2 * (e10/e14);
	double    e18 = exp(2 * (e10 * e15));
	double    e21 = geno_f * (f2 - f0)/2;
	double    e22 = f0 + e21;
	double    e24 = e11 - e16 * e18;
	double    e25 = e22 - f1;
	double    e26 = e7 * e25;
	double    e27 = e6/e10;
	double    e28 = e22 - m0;
	double    e31 = 2 * (e26/e10);
	double    e32 = 2 * (e28 * e10);
	double    e33 = e31 - e32;
	double    e34 = pow(e11,2);
	double    e36 = 2 * (e3 * e10);
	double    e39 = e6/e36 - 2 * (e12/e34);
	double    e40 = e16 * e6;
	double    e42 = e40 * e15/e10;
	double    e43 = e24 * e12;
	double    e44 = 2 * e42;
	double    e48 = 2 * (e39 * e10/e14) + 2 * e27;
	double    e51 = 1 - 4 * (e3/e24);
	double    e52 = exp(-(t0 * e10));
	double    e53 = exp(tau * e10);
	double    e54 = e48/e14;
	double    e56 = e54 + 2 + e44;
	double    e58 = 2 - e56 * e18;
	double    e59 = 2 * e10;
	double    e60 = e51 * e4;
	double    e62 = e33 * e52 * e53;
	double    e64 = e60 * e25/e10;
	double    e65 = 0.5 * e27;
	double    e66 = e3/e10;
	double    e67 = e64 + 2 * (e62 * e3/e43);
	double    e68 = e58 * e3;
	double    e70 = 2 * (e3 * e14);
	double    e71 = 0.5 * geno_q;
	double    e73 = pow(e59,2);
	double    e76 = tau * e6/e10 - e58/e24;
	double    e77 = e12/e3;
	double    e83 = 0.5 * (t0 * e33) + 4 * (e26/e73) + f0 + e21 - m0;
	double    e84 = 1 - e71;
	double    e86 = e51 * e6/e9;
	double    e88 = 4 - 4 * (e68/e24);
	double    e91 = e16 * e3 * e15/e10;
	double    e101 = 2 * (e33 * e76) - 2 * (2 * (e39 * e33 * e3/e12) + 2 * (e83 * e6/e10));
	double    e104 = 2 * e86 + 2 * (e88/e24);
	double    e106 = e8/e9;
	double    e108 = 2 * (e26/e9) + 2 * e28;
	double    e109 = e27 - e77;
	double    e110 = e65 - e13;
	double    e111 = e8/e10;
	double    e112 = 2 * e66;
	double    e115 = e109 * e10/e70 + e27;
	double    e118 = e110 * e10/e70 + e65;
	double    e119 = e104 * e4;
	double    e120 = e66 + 1/(2 * e14);
	double    e122 = 1/e14 + e112;
	double    e123 = exp(tau * theta);
	double    e124 = pow(e67,2);
	double    e134 = 2 * (e101 * e52 * e53 * e3/e43 - e119 * e25/e59) + 2 * (e62/e43);
	double    e135 = 2 * (e10/e24);
	double    e137 = pow(-e24,2);
	double    e138 = e124 - e135;
	double    e141 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e144 = e138 * e3 + 0.5 * e12 + e123 * e141;
	double    e148 = 2 * (e118/e14);
	double    e149 = 2 + 2 * (e115/e14);
	double    e150 = pow(e36,2);
	double    e153 = e64 + e134 * e3;
	double    e156 = 1 - (e42 + 1 + e148) * e18;
	double    e157 = 2 - (e149 + e44) * e18;
	double    e159 = 2 * (e120/e14) + 2 * e91;
	double    e161 = 2 * (e122/e14) + 4 * e91;
	double    e163 = 2 * e111 + e59;
	double    e164 = pow(e24,2);
	double    e165 = pow(e5,4);
	double    e166 = 1 - e106;
	double    e167 = 2 * e106;
	double    e171 = e153 * e67 + e65 + 2 * (e68 * e10/e137) - e163/e24;
	double    e172 = e9 * e10;
	double    e173 = e7/e10;
	double    e174 = 2 - e167;
	double    e176 = e24 * e6/e10;
	double    e178 = e24 * e3/e10;
	double    e180 = e165/e172;
	double    e181 = 0.5 * geno_f;
	double    e182 = 2 * (e8/e150);
	double    e183 = 2/e34;
	double    e184 = 1 - e181;
	double    e185 = e182 + e183;
	double    e186 = e56 * e3;
	double    e187 = e51/e9;
	double    e189 = e33 * e15 - e108;
	double    e190 = e185 * e3;
	double    e191 = e40/e9;
	double    e194 = e173 - e10;
	double    e195 = 0.5 - e156 * e3/e24;
	double    e197 = 2 * e173 - e59;
	double    e200 = e39 * e108;
	double    e208 = e156 * e12 + 0.5 * e176;
	double    e209 = e166/e10;
	double    e212 = e84 * e51 * e6/e9;
	double    e214 = e157 * e12 + e176;
	double    e215 = e174/e10;
	double    e218 = e159 * e18;
	double    e221 = e161 * e18;
	double    e222 = e178 - e159 * e12 * e18;
	double    e224 = e26/pow(e59,3);
	double    e227 = 0.5 * tau - 0.5 * t0;
	double    e228 = 1/e11;
	double    e229 = 1/e3;
	double    e231 = 2 * ((e111 + e10) * e6/e150) + 2 * ((e27 - 2 * e77)/e34);
	double    e233 = 2 * ((0.5 * e111 + 0.5 * e10) * e6/e150) + 2 * ((e65 - e77)/e34);
	double    e235 = 2 * e178 - e161 * e12 * e18;
	double    e236 = 4 * (e195/e24);
	double    e237 = 4 * ((1 - (e84 * e157 * e3/e24 + e71))/e24);
	double    e238 = t0 * e108;
	double    e239 = (((e110 * e48/e11 + 2 * (e118 * e39 - e233 * e10))/e14 - e180)/e14 + (e56 + 2 * (1 + e148 - 0.5 * e191)) * e6 * e15/e10) * e18;
	double    e240 = pow(e144,2);
	double    e242 = ((e109 * e48/e11 + 2 * (e115 * e39 - e231 * e10))/e14 - 2 * e180)/e14 + (2 * e56 + 2 * (e149 - e191)) * e6 * e15/e10;
	double    e246 = ((e48/e59 + 2 * (e39 * e120 + e228 - e190))/e14 + 2 * e209)/e14 + (2 * e186 + 2 * (e166 * e16 + 2 * (e120 * e6/e14))) * e15/e10;
	double    e248 = ((e48/e10 + 2 * (e39 * e122 + e229 - 2 * e190))/e14 + 2 * e215)/e14 + (2 * (e174 * e16 + 2 * (e122 * e6/e14)) + 4 * e186) * e15/e10;
	double    e257 = (e227 * e33 - 0.5 * e108) * e6/e10 - e208 * e33/e43;
	double    e260 = (e212 + e237) * e4 * e25/e10;
	double    e263 = (e187 + 4 * (e218/e164)) * e4 * e25/e10;
	double    e268 = e189 * e6/e10 - e214 * e33/e43;
	double    e271 = (0.5 * e86 + e236) * e4 * e25/e10;
	double    e272 = e64 + (e134 - 2 * (e171 * e67/e144)) * e3;
	double    e273 = e60/e10;
	double    e276 = (2 * e187 + 4 * (e221/e164)) * e4 * e25/e10;
	double    e277 = 2 * (((e33 * (2 * tau - 2 * t0) - 2 * e108) * e3/e10 - e235 * e33/e43) * e52 * e53/e43);
	double    e278 = 2 * ((e189 * e3/e10 - e222 * e33/e43) * e52 * e53/e43);
	double    e279 = (e242 * e3 + e54 + 2 + e44) * e18;
	double    e280 = e239 * e3;
	double    e289 = e194 * e52 * e53;
	double    e291 = e200 * e6/e10;
	double    e292 = e200 * e3;
	double    e293 = e39/e12;
	double    e294 = e83/e9;
	double    e295 = e184 * e197;
	double    e296 = e58 * e157;
	double    e299 = e58 * e159/e24 - e246;
	double    e302 = e58 * e161/e24 - e248;
	double    e303 = e104 * e3;
	double    e304 = e104/2;
	double    e307 = e108 * e6 * e76/e10;
	double    e310 = e108 * e3 * e76/e10;
	double    e311 = e12 * e10;
	double    e315 = e7/e73;
	double    e316 = 0.5 * e273;
	double    e318 = 0.5 * e238 + 16 * e224;
	double    e320 = 2 * ((e257 * e3 + 0.5 * e33) * e52 * e53/e43) - e271;
	double    e322 = 2 * ((e268 * e3 + e31 - e32) * e84 * e52 * e53/e43) - e260;
	double    e323 = e277 - e276;
	double    e324 = e278 - e263;
	double    e326 = 2 + 2 * e166;
	double    e327 = 4 * (e4 * e52 * e53/e12);
	double    e329 = t0 * e101 * e3;
	double    e331 = tau * e165/e172;
	double    e333 = tau * e3/e10;
	
	arma::vec df = arma::zeros<arma::vec>(10);

	df(0) = ((e84 * e134 + (2 * (((((e6 * e15/e10 - e214/e43) * 
        e101 - 2 * (e307 + e33 * (e331 - (e242 * e18 + 
        e296/e24)/e24))) * e84 - 2 * (2 * ((e39 * (1 - (e84 * 
        e6 * e3/e311 + e71)) * e33 - (e291 + e231 * e33) * 
        e84 * e3)/e12) - 2 * ((e318/e10 + e294) * e84 * 
        e165/e10))) * e3 + e84 * e101) * e52 * e53/e43 + 
        ((e104 * e6/e9 + 2 * ((e157 * e88 + 4 * (2 - (e279 + 
            e296 * e3/e24)))/e164)) * e84 + 2 * ((2 * e212 + 
            e237) * e6/e9)) * e4 * e25/e59) + 2 * (e268 * 
        e84 * e52 * e53/e43)) * e3 - e260) * e67 + e153 * 
        e322 + e84 * (2 * (((e27 - 2 * (e157 * e24 * e10/e137)) * 
        e58 * e3 + (2 - e279) * e10)/e137) - ((e326 * e6/e10 - 
        e157 * e163/e24)/e24 + 0.5 * e180)) - ((e124 + 
        e65 - e135) * e84 + (2 * (e67 * e322) - 2 * ((e27 - 
        e157 * e10/e24) * e84/e24)) * e3) * e171/e144) * 
        e84/e144;
		
	df(1) = geno_q * (e153 * e320 + e67 * ((2 * 
        ((((e227 * e6/e10 - e208/e43) * e101 - (2 * ((0.5 * 
            e331 - (e239 + e156 * e58/e24)/e24) * e33 + 
            0.5 * e307) + 2 * (2 * ((e39 * (0.5 - 0.5 * (e8/e311)) * 
            e33 - (e233 * e33 + 0.5 * e291) * e3)/e12) - 
            2 * (((0.25 * e238 + 8 * e224)/e10 + 0.5 * e294) * 
                e165/e10)))) * e3 + 0.5 * e101) * e52 * 
            e53/e43 + ((e304 + 2 * (e86 + e236)) * e6/e9 + 
            2 * ((e156 * e88 + 4 * (e195 * e58 - e280))/e164)) * 
            e4 * e25/e59) + 2 * (e257 * e52 * e53/e43)) * 
        e3 + 0.5 * e134 - e271) + 2 * (((e65 - 2 * (e156 * 
        e24 * e10/e137)) * e58 * e3 + (0.5 * e58 - e280) * 
        e10)/e137) - (e171 * ((2 * (e67 * e320) - 2 * ((e65 - 
        e156 * e10/e24)/e24)) * e3 + 0.25 * e27 + 0.5 * 
        e138)/e144 + ((1 + 2 * (0.5 - 0.5 * e106)) * e6/e10 - 
        e156 * e163/e24)/e24 + 0.25 * e180)) * e84/e144;
		
	df(2) = ((e184 * e51 * e4/e10 + (2 * (e295 * e52 * 
            e53/e43) + 2 * ((2 * (e295 * e76) - 2 * (2 * 
            (((0.5 * (t0 * e197) + 4 * e315) * e184 + 1 - 
                e181) * e6/e10) + 2 * (e39 * e184 * e197 * 
            e3/e12))) * e52 * e53 * e3/e43 - e184 * e104 * 
            e4/e59)) * e3) * e67 + e272 * (e273 + 2 * (e197 * 
            e52 * e53 * e3/e43)) * e184) * e84/e144;
			
	df(3) = geno_f * 
            (e272 * (e316 + 2 * (e289 * e3/e43)) + e67 * 
                ((2 * (e289/e43) + 2 * ((2 * (e194 * e76) - 
                  2 * (2 * (e194 * e39 * e3/e12) + 2 * ((0.5 + 
                    0.5 * (t0 * e194) + 2 * e315) * e6/e10))) * 
                  e52 * e53 * e3/e43 - e119/(4 * e10))) * 
                  e3 + e316)) * e84/e144;
				  
	df(4) = ((((4 - e327)/e24 + 
            2 * (e304 + (2 * (2 * ((4/e73 + t0/e10) * e6) + 
                4 * (e39 * e3/e12)) - 4 * e76) * e4 * e52 * 
                e53 * e3/e43)) * e3 - 1) * e67 - e272 * 
            ((e327 - 4) * e3/e24 + 1)) * e84 * e4/(e144 * 
            e10);
			
	df(5) = -(e171 * (1 - 0.5 * geno_mu) * e84 * e123/e240);
			
	df(6) = -(0.5 * (geno_mu * e171 * e84 * e123/e240));
			
	df(7) = -(tau * e171 * e84 * e123 * e141/e240);
			
	df(8) = ((e153 * e323 + e67 * (e277 + 
            2 * ((e101 * (2 * e333 - e235/e43) + 2 * (e33 * 
                (tau * e174/e10 - e302 * e18/e24) - 2 * 
                e310) - (2 * (2 * (((e229 - (2 * e293 + 2 * 
                e185) * e3) * e33 - 2 * e292) * e3/e12) + 
                2 * (e83 * e174 - (32 * e224 + e238) * e6 * 
                  e3/e10)) + 2 * e329)/e10) * e52 * e53 * 
                e3/e43 - ((2 * (e51 * (2 - 4 * e106) - 4 * 
                (e161 * e6 * e18 * e3/e164)) - 2 * e303)/e9 - 
                2 * ((4 * (e302 * e3) - e161 * e88) * e18/e164)) * 
                e4 * e25/e59) - e276) + 2 * ((e58 * (2 * 
            (e161 * e24 * e18 * e10/e137) + e112) - e248 * 
            e18 * e10)/e137) - e171 * (1/e10 + 2 * (e67 * 
            e323 * e3) - 2 * ((e221 * e10/e24 + e112)/e24))/e144) * 
            e3 + 0.5 * e215 - (e161 * e163 * e18/e24 + 
            (2 * e174 + 4) * e3/e10)/e24) * (1 - 0.5 * geno_b) * 
            e84 * e5/e144;
			
	df(9) = geno_b * ((e153 * e324 + 
            e67 * (e278 + 2 * ((e101 * (e333 - e222/e43) + 
                2 * (e33 * (tau * e166/e10 - e299 * e18/e24) - 
                  e310) - (2 * (2 * (((e228 - (e293 + e182 + 
                e183) * e3) * e33 - e292) * e3/e12) + 2 * 
                (e83 * e166 - e318 * e6 * e3/e10)) + e329)/e10) * 
                e52 * e53 * e3/e43 - ((2 * ((1 - e167) * 
                e51 - 4 * (e159 * e6 * e18 * e3/e164)) - 
                e303)/e9 - 2 * ((4 * (e299 * e3) - e159 * 
                e88) * e18/e164)) * e4 * e25/e59) - e263) + 
            2 * (((e66 + 2 * (e159 * e24 * e18 * e10/e137)) * 
                e58 - e246 * e18 * e10)/e137) - e171 * 
            (0.5/e10 + 2 * (e67 * e324 * e3) - 2 * ((e218 * 
                e10/e24 + e66)/e24))/e144) * e3 + 0.5 * 
            e209 - (e159 * e163 * e18/e24 + e326 * e3/e10)/e24) * 
            e84 * e5/e144;
	
    return(Rcpp::wrap(df));
}


RcppExport SEXP d_f_i1_q2_qff1mtb_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = b0 + geno_b * (b2 - b0)/2;
	double    e6 = pow(e5,2);
	double    e7 = pow(e4,2);
	double    e8 = e6 * e3;
	double    e9 = e7 + 2 * e8;
	double    e10 = sqrt(e9);
	double    e11 = 2 * e3;
	double    e12 = e4 + e10;
	double    e13 = e12/e11;
	double    e14 = r0 - e13;
	double    e15 = tau - t0;
	double    e16 = e11 + 2 * (e10/e14);
	double    e18 = exp(2 * (e10 * e15));
	double    e21 = geno_f * (f2 - f0)/2;
	double    e22 = f0 + e21;
	double    e24 = e11 - e16 * e18;
	double    e25 = e22 - f1;
	double    e26 = e7 * e25;
	double    e27 = e22 - m0;
	double    e28 = e6/e10;
	double    e32 = 2 * (e26/e10) - 2 * (e27 * e10);
	double    e33 = pow(e11,2);
	double    e35 = 2 * (e3 * e10);
	double    e38 = e6/e35 - 2 * (e12/e33);
	double    e39 = e24 * e12;
	double    e42 = e16 * e6 * e15/e10;
	double    e45 = e38 * e10/e14 + e28;
	double    e47 = 1 - 4 * (e3/e24);
	double    e49 = exp(-(t0 * e10));
	double    e50 = exp(tau * e10);
	double    e53 = e45/e14 + e42 + 1;
	double    e55 = 1 - e53 * e18;
	double    e56 = 2 * e10;
	double    e57 = e47 * e4;
	double    e60 = e57 * e25/e10 + 2 * (e32 * e49 * e50 * e3/e39);
	double    e61 = 0.5 * e28;
	double    e62 = e3/e10;
	double    e63 = e55 * e3;
	double    e64 = t0 * e32;
	double    e68 = e16 * e3 * e15/e10;
	double    e69 = pow(e56,2);
	double    e71 = 0.5 * (tau * e6/e10) - e55/e24;
	double    e74 = e47 * e6/e9;
	double    e80 = 0.5 * e64 + 4 * (e26/e69) + f0 + e21 - m0;
	double    e81 = 2 - 4 * (e63/e24);
	double    e82 = e8/e9;
	double    e83 = e61 - e13;
	double    e92 = 2 * (e71 * e32) - 2 * (e38 * e32 * e3/e12 + e80 * e6/e10);
	double    e93 = 2 * (e81/e24);
	double    e94 = 2 * e62;
	double    e95 = e74 + e93;
	double    e98 = e83 * e10/(2 * (e3 * e14)) + e61;
	double    e99 = e62 + 1/(2 * e14);
	double    e101 = 1/e14 + e94;
	double    e104 = 2 * (e26/e9) + 2 * e27;
	double    e105 = exp(tau * theta);
	double    e106 = e95 * e4;
	double    e108 = pow(e60,2) - 2 * (e10/e24);
	double    e115 = e92 * e49 * e50 * e3/e39 - e106 * e25/e56;
	double    e116 = e8/e10;
	double    e119 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e121 = pow(-e24,2);
	double    e124 = e108 * e3 + 0.5 * e12 + e105 * e119;
	double    e125 = 0.5 * e60;
	double    e130 = 1 - (e42 + 1 + 2 * (e98/e14)) * e18;
	double    e132 = 2 * (e99/e14) + 2 * e68;
	double    e134 = 2 * (e101/e14) + 4 * e68;
	double    e136 = pow(e35,2);
	double    e137 = e116 + e10;
	double    e138 = e125 + 2 * (e115 * e3);
	double    e139 = pow(e24,2);
	double    e140 = 2 * e82;
	double    e141 = e7/e10;
	double    e142 = 0.25 * e28;
	double    e146 = e60 * e138 + e142 + 2 * (e63 * e10/e121) - e137/e24;
	double    e147 = 1 - e82;
	double    e148 = 2 - e140;
	double    e150 = e24 * e3/e10;
	double    e152 = 2 * (e8/e136);
	double    e153 = 2/e33;
	double    e155 = tau * e3/e10;
	double    e156 = e47/e9;
	double    e157 = 0.5 * geno_f;
	double    e158 = e152 + e153;
	double    e159 = e53 * e3;
	double    e160 = e158 * e3;
	double    e161 = pow(e5,4);
	double    e162 = 1 - e157;
	double    e163 = e9 * e10;
	double    e166 = e132 * e18;
	double    e167 = e141 - e10;
	double    e168 = 0.5 - e130 * e3/e24;
	double    e170 = 2 * e141 - e56;
	double    e177 = e130 * e12 + 0.5 * (e24 * e6/e10);
	double    e178 = e147/e10;
	double    e179 = e57/e10;
	double    e180 = e156 + 4 * (e166/e139);
	double    e181 = e148/e10;
	double    e184 = e134 * e18;
	double    e186 = e161/e163;
	double    e189 = 0.5 * tau - 0.5 * t0;
	double    e190 = 1/e11;
	double    e191 = 1/e3;
	double    e193 = 2 * ((0.5 * e116 + 0.5 * e10) * e6/e136) + 2 * ((e61 - e12/e3)/e33);
	double    e195 = 2 * e155 - (2 * e150 - e134 * e12 * e18)/e39;
	double    e196 = 4 * (e168/e24);
	double    e197 = e155 - (e150 - e132 * e12 * e18)/e39;
	double    e198 = (((e45 * e83/e11 + e98 * e38 - e193 * e10)/e14 - 0.5 * e186)/e14 + ((e45 + 2 * e98)/e14 + (e15/e10 - 0.5/e9) * e16 * e6 + 2) * e6 * e15/e10) * e18;
	double    e200 = ((e45/e56 + e38 * e99 + e190 - e160)/e14 + e178)/e14 + (e147 * e16 + 2 * e159 + 2 * (e99 * e6/e14)) * e15/e10;
	double    e202 = ((e45/e10 + e38 * e101 + e191 - 2 * e160)/e14 + e181)/e14 + (e148 * e16 + 2 * (e101 * e6/e14) + 4 * e159) * e15/e10;
	double    e203 = pow(e124,2);
	double    e218 = e38 * e104;
	double    e222 = e71 * e104;
	double    e227 = (2 * e115 - 2 * (e146 * e60/e124)) * e3 + e125;
	double    e229 = e26/pow(e56,3);
	double    e231 = 2 * (((e189 * e32 - 0.5 * e104) * e6 * e3/e10 + (0.5 - e177 * e3/e39) * e32) * e49 * e50/e39) - (0.5 * e74 + e196) * e4 * e25/e10;
	double    e233 = 2 * ((e32 * e195 - (2 * e104 + 2 * e64) * e3/e10) * e49 * e50/e39) - (2 * e156 + 4 * (e184/e139)) * e4 * e25/e10;
	double    e235 = 2 * ((e32 * e197 - (e104 + e64) * e3/e10) * e49 * e50/e39) - e180 * e4 * e25/e10;
	double    e236 = t0 * e104;
	double    e237 = e198 * e3;
	double    e238 = e95/2;
	double    e243 = e218 * e3;
	double    e244 = e38/e12;
	double    e246 = e222 * e3/e10;
	double    e249 = e55 * e132/e24 - e200;
	double    e252 = e55 * e134/e24 - e202;
	double    e253 = e179 + 2 * (e170 * e49 * e50 * e3/e39);
	double    e256 = (4 * (e4 * e49 * e50/e12) - 4) * e3/e24 + 1;
	double    e257 = e7/e69;
	double    e259 = 0.5 * e179 + 2 * (e167 * e49 * e50 * e3/e39);
	double    e261 = t0 * e92 * e3;
	
	arma::vec df = arma::zeros<arma::vec>(9);

	df(0) = pow(geno_q,2) * (e60 * (0.5 * e231 + 2 * (((((e189 * 
        e6/e10 - e177/e39) * e92 - (2 * ((e38 * (0.5 - 
        0.5 * (e8/(e12 * e10))) * e32 - (e193 * e32 + 0.5 * 
        (e218 * e6/e10)) * e3)/e12 - ((0.25 * e236 + 8 * 
        e229)/e10 + 0.5 * (e80/e9)) * e161/e10) + 2 * ((0.25 * 
        (tau * e161/e163) - (e198 + e55 * e130/e24)/e24) * 
        e32 + 0.5 * (e222 * e6/e10)))) * e3 + 0.5 * e92) * 
        e49 * e50/e39 + ((e238 + e74 + e196) * e6/e9 + 
        2 * ((e130 * e81 + 4 * (e168 * e55 - e237))/e139)) * 
        e4 * e25/e56) * e3 + 0.5 * e115)) + e138 * e231 + 
        2 * (((e61 - 2 * (e130 * e24 * e10/e121)) * e55 * 
            e3 + (0.5 * e55 - e237) * e10)/e121) - (e146 * 
        ((2 * (e60 * e231) - 2 * ((e61 - e130 * e10/e24)/e24)) * 
            e3 + e142 + 0.5 * e108)/e124 + ((1 - 0.5 * e82) * 
        e6/e10 - e137 * e130/e24)/e24 + 0.125 * e186))/e124;
		
	df(1) = geno_q * (e60 * (0.5 * (e253 * e162) + 2 * (((2 * 
            (e71 * e162 * e170) - 2 * (((0.5 * (t0 * e170) + 
            4 * e257) * e162 + 1 - e157) * e6/e10 + e38 * 
            e162 * e170 * e3/e12)) * e49 * e50 * e3/e39 - 
            e95 * e162 * e4/e56) * e3)) + e253 * e227 * 
            e162)/e124;
			
	df(2) = geno_f * geno_q * (e60 * (0.5 * 
            e259 + 2 * (((2 * (e167 * e71) - 2 * (e167 * 
            e38 * e3/e12 + (0.5 + 0.5 * (t0 * e167) + 2 * 
            e257) * e6/e10)) * e49 * e50 * e3/e39 - e106/(4 * 
            e10)) * e3)) + e227 * e259)/e124;
	
	df(3) = geno_q * 
            (e60 * (2 * ((e238 + (2 * ((4/e69 + t0/e10) * 
                e6 + 2 * (e38 * e3/e12)) - 4 * e71) * e4 * 
                e49 * e50 * e3/e39) * e3) - 0.5 * e256) - 
                e227 * e256) * e4/(e124 * e10);
				
	df(4) = -(geno_q * e146 * (1 - 0.5 * geno_mu) * e105/e203);
	
	df(5) = -(0.5 * (geno_mu * geno_q * e146 * e105/e203));
	
	df(6) = -(geno_q * tau * e146 * e105 * e119/e203);
	
	df(7) = geno_q * 
            ((e60 * (0.5 * e233 + 2 * ((e92 * e195 + 2 * 
                ((0.5 * (tau * e148/e10) - e252 * e18/e24) * 
                  e32 - 2 * e246) - (2 * ((((e191 - (2 * e244 + 
                2 * e158) * e3) * e32 - 2 * e243)/e12 - 
                (32 * e229 + e236) * e6/e10) * e3 + e80 * 
                e148) + 2 * e261)/e10) * e49 * e50 * e3/e39 - 
                ((e47 * (2 - 4 * e82) - (2 * e95 + 4 * (e134 * 
                  e6 * e18/e139)) * e3)/e9 - 2 * ((4 * (e252 * 
                  e3) - e81 * e134) * e18/e139)) * e4 * 
                  e25/e56)) + e138 * e233 + 2 * ((e55 * 
                (2 * (e134 * e24 * e18 * e10/e121) + e94) - 
                e202 * e18 * e10)/e121) - e146 * (1/e10 + 
                2 * (e60 * e233 * e3) - 2 * ((e184 * e10/e24 + 
                e94)/e24))/e124) * e3 + 0.25 * e181 - (e137 * 
                e134 * e18/e24 + (4 - e140) * e3/e10)/e24) * 
            (1 - 0.5 * geno_b) * e5/e124;
	
	df(8) = geno_b * geno_q * 
            ((e60 * (0.5 * e235 + 2 * ((e92 * e197 + 2 * 
                ((0.5 * (tau * e147/e10) - e249 * e18/e24) * 
                  e32 - e246) - (2 * ((((e190 - (e244 + e152 + 
                e153) * e3) * e32 - e243)/e12 - (0.5 * e236 + 
                16 * e229) * e6/e10) * e3 + e80 * e147) + 
                e261)/e10) * e49 * e50 * e3/e39 - (((1 - 
                e140) * e47 - (e180 * e6 + e93) * e3)/e9 - 
                2 * ((4 * (e249 * e3) - e81 * e132) * e18/e139)) * 
                e4 * e25/e56)) + e138 * e235 + 2 * (((e62 + 
                2 * (e132 * e24 * e18 * e10/e121)) * e55 - 
                e200 * e18 * e10)/e121) - e146 * (0.5/e10 + 
                2 * (e60 * e235 * e3) - 2 * ((e166 * e10/e24 + 
                e62)/e24))/e124) * e3 + 0.25 * e178 - (e137 * 
                e132 * e18/e24 + (2 - e82) * e3/e10)/e24) * 
            e5/e124;
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_f0_ff1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e7 = sqrt(e5 + 2 * (pow(b0 + geno_b * (b2 - b0)/2,2) * e3));
	double    e8 = 2 * e3;
	double    e9 = e4 + e7;
	double    e11 = e8 - (e8 + 2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * (tau - t0)));
	double    e12 = f0 + geno_f * (f2 - f0)/2;
	double    e13 = e12 - f1;
	double    e15 = exp(-(t0 * e7));
	double    e16 = exp(tau * e7);
	double    e17 = (1 - 4 * (e3/e11)) * e4;
	double    e18 = e11 * e9;
	double    e21 = e17 * e13/e7 + 2 * ((2 * (e5 * e13/e7) - 2 * ((e12 - m0) * e7)) * e15 * e16 * e3/e18);
	double    e22 = pow(e21,2);
	double    e23 = exp(tau * theta);
	double    e26 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e29 = (e22 - 2 * (e7/e11)) * e3 + 0.5 * e9 + e23 * e26;
	double    e31 = e17/e7 + 2 * ((2 * (e5/e7) - 2 * e7) * e15 * e16 * e3/e18);
	double    e33 = 1 - 0.5 * geno_f;
	double    e34 = pow(e29,2);
	double    e37 = 1 - 2 * (e22 * e3/e29);
	double    e38 = pow(e31,2);
	
	arma::vec df = arma::zeros<arma::vec>(6);

	df(0) = 2 * (e38 * pow(e33,2) * e37 * e3/e29);
		
	df(1) = geno_f * e38 * e33 * e37 * e3/e29;
			
	df(2) = -(2 * (e31 * ((4 * (e4 * e15 * e16/e9) - 4) * e3/e11 + 1) * e33 * e37 * e4 * e3/(e29 * e7)));
	
	df(3) = -(2 * (e21 * e31 * e33 * (1 - 0.5 * geno_mu) * e23 * e3/e34));
				
	df(4) = -(geno_mu * e21 * e31 * e33 * e23 * e3/e34);
	
	df(5) = -(2 * (tau * e21 * e31 * e33 * e23 * e26 * e3/e34));
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_f2_ff1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e7 = sqrt(e5 + 2 * (pow(b0 + geno_b * (b2 - b0)/2,2) * e3));
	double    e8 = 2 * e3;
	double    e9 = e4 + e7;
	double    e11 = e8 - (e8 + 2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * (tau - t0)));
	double    e12 = f0 + geno_f * (f2 - f0)/2;
	double    e13 = e12 - f1;
	double    e15 = exp(-(t0 * e7));
	double    e16 = exp(tau * e7);
	double    e17 = (1 - 4 * (e3/e11)) * e4;
	double    e18 = e11 * e9;
	double    e21 = e17 * e13/e7 + 2 * ((2 * (e5 * e13/e7) - 2 * ((e12 - m0) * e7)) * e15 * e16 * e3/e18);
	double    e22 = 2 * e7;
	double    e23 = exp(tau * theta);
	double    e24 = pow(e21,2);
	double    e27 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e30 = (e24 - 2 * (e7/e11)) * e3 + 0.5 * e9 + e23 * e27;
	double    e34 = (2 * (e5/e7) - e22) * e15 * e16 * e3/e18;
	double    e36 = e17/e22 + e34;
	double    e37 = pow(e30,2);
	double    e40 = 1 - 2 * (e24 * e3/e30);
	
	arma::vec df = arma::zeros<arma::vec>(5);

	df(0) = pow(geno_f,2) * e36 * (e17/e7 + 2 * e34) * e40 * e3/e30;
		
	df(1) = -(2 * (geno_f * e36 * ((4 * (e4 * e15 * e16/e9) - 4) * e3/e11 + 1) * e40 * e4 * e3/(e30 * e7)));
			
	df(2) = -(2 * (geno_f * e21 * e36 * (1 - 0.5 * geno_mu) * e23 * e3/e37));
	
	df(3) = -(geno_f * geno_mu * e21 * e36 * e23 * e3/e37);
				
	df(4) = -(2 * (geno_f * tau * e21 * e36 * e23 * e27 * e3/e37));
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_f1_f1mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e7 = sqrt(e5 + 2 * (pow(b0 + geno_b * (b2 - b0)/2,2) * e3));
	double    e8 = 2 * e3;
	double    e9 = e4 + e7;
	double    e11 = e8 - (e8 + 2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * (tau - t0)));
	double    e12 = f0 + geno_f * (f2 - f0)/2;
	double    e13 = e12 - f1;
	double    e15 = exp(-(t0 * e7));
	double    e16 = exp(tau * e7);
	double    e20 = (1 - 4 * (e3/e11)) * e4 * e13/e7 + 2 * ((2 * (e5 * e13/e7) - 2 * ((e12 - m0) * e7)) * e15 * e16 * e3/(e11 * e9));
	double    e21 = exp(tau * theta);
	double    e22 = pow(e20,2);
	double    e25 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e26 = ((e22 - 2 * (e7/e11)) * e3 + 0.5 * e9 + e21 * e25) * e7;
	double    e27 = pow(e26,2);
	double    e30 = (4 * (e4 * e15 * e16/e9) - 4) * e3/e11 + 1;
		
	arma::vec df = arma::zeros<arma::vec>(4);

	df(0) = 2 * (pow(e30,2) * (1/e26 - 2 * (e22 * e3 * e7/e27)) * e5 * e3/e7);
		
	df(1) = 2 * (e20 * e30 * (1 - 0.5 * geno_mu) * e4 * e21 * e3 * e7/e27);
			
	df(2) = geno_mu * e20 * e30 * e4 * e21 * e3 * e7/e27;
	
	df(3) = 2 * (tau * e20 * e30 * e4 * e21 * e25 * e3 * e7/e27);
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_m0_mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e7 = sqrt(e5 + 2 * (pow(b0 + geno_b * (b2 - b0)/2,2) * e3));
	double    e8 = 2 * e3;
	double    e9 = e4 + e7;
	double    e12 = f0 + geno_f * (f2 - f0)/2;
	double    e13 = exp(tau * theta);
	double    e15 = e13 * (geno_mu * (mu02 - mu00)/2 + mu00);
	double    e19 = (pow((1 - 4 * (e3/(e8 - (e8 + 2 * (e7/(r0 - e9/e8))) * 
        exp(2 * (e7 * (tau - t0)))))) * e4 * (e12 - f1)/e7 + 
        2 * ((2 * (e5 * (e12 - f1)/e7) - 2 * ((e12 - m0) * 
            e7)) * exp(-(t0 * e7)) * exp(tau * e7) * e3/((e8 - 
            (e8 + 2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * 
                (tau - t0)))) * e9)),2) - 2 * (e7/(e8 - (e8 + 
        2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * (tau - t0)))))) * 
        e3 + 0.5 * e9 + e15;
	double    e21 = 1 - 0.5 * geno_mu;
	double    e22 = pow(e19,2);
	double    e23 = pow(e13,2);
		
	arma::vec df = arma::zeros<arma::vec>(3);

	df(0) = -(pow(e21,2) * e23/e22);
		
	df(1) = -(0.5 * (geno_mu * e21 * e23/e22));
			
	df(2) = tau * e21 * (1 - e15/e19) * e13/e19;
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_m2_mt_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
		double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);

	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e7 = sqrt(e5 + 2 * (pow(b0 + geno_b * (b2 - b0)/2,2) * e3));
	double    e8 = 2 * e3;
	double    e9 = e4 + e7;
	double    e12 = f0 + geno_f * (f2 - f0)/2;
	double    e13 = exp(tau * theta);
	double    e15 = e13 * (geno_mu * (mu02 - mu00)/2 + mu00);
    double e19 = (pow((1 - 4 * (e3/(e8 - (e8 + 2 * (e7/(r0 - e9/e8))) * 
        exp(2 * (e7 * (tau - t0)))))) * e4 * (e12 - f1)/e7 + 
        2 * ((2 * (e5 * (e12 - f1)/e7) - 2 * ((e12 - m0) * 
            e7)) * exp(-(t0 * e7)) * exp(tau * e7) * e3/((e8 - 
            (e8 + 2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * 
                (tau - t0)))) * e9)),2) - 2 * (e7/(e8 - (e8 + 
        2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * (tau - t0)))))) * e3 + 0.5 * e9 + e15;
		
	arma::vec df = arma::zeros<arma::vec>(2);

	df(0) = -(0.25 * (pow(geno_mu,2) * pow(e13,2)/pow(e19,2)));
		
	df(1) = 0.5 * (geno_mu * tau * (1 - e15/e19) * e13/e19);
	
    return(Rcpp::wrap(df));
}

RcppExport SEXP d_f_i1_t_t_g_c(SEXP a0_p, SEXP a2_p, SEXP b0_p, SEXP b2_p, SEXP q0_p, SEXP q2_p, SEXP f0_p, SEXP f2_p, SEXP f1_p, SEXP mu00_p, SEXP mu02_p, SEXP theta_p,
SEXP m0_p, SEXP r0_p, SEXP tau_p, SEXP t0_p, SEXP geno_a_p, SEXP geno_b_p, SEXP geno_q_p, SEXP geno_f_p, SEXP geno_mu_p)
{
	arma::vec a0_v = as<arma::vec>(a0_p);
	arma::vec a2_v = as<arma::vec>(a2_p);
	arma::vec b0_v = as<arma::vec>(b0_p);
	arma::vec b2_v = as<arma::vec>(b2_p);
	arma::vec q0_v = as<arma::vec>(q0_p);
	arma::vec q2_v = as<arma::vec>(q2_p);
	arma::vec f0_v = as<arma::vec>(f0_p);
	arma::vec f2_v = as<arma::vec>(f2_p);
	arma::vec f1_v = as<arma::vec>(f1_p);
	arma::vec mu00_v = as<arma::vec>(mu00_p);
	arma::vec mu02_v = as<arma::vec>(mu02_p);
	arma::vec theta_v = as<arma::vec>(theta_p);
	arma::vec m0_v = as<arma::vec>(m0_p);
	arma::vec r0_v = as<arma::vec>(r0_p);
	arma::vec tau_v = as<arma::vec>(tau_p);
	arma::vec t0_v = as<arma::vec>(t0_p);
	arma::vec geno_a_v = as<arma::vec>(geno_a_p);
	arma::vec geno_b_v = as<arma::vec>(geno_b_p);
	arma::vec geno_q_v = as<arma::vec>(geno_q_p);
	arma::vec geno_f_v = as<arma::vec>(geno_f_p);
	arma::vec geno_mu_v = as<arma::vec>(geno_mu_p);
	
	double	a0 = a0_v(0);
	double	a2 = a2_v(0);
	double	b0 = b0_v(0);
	double	b2 = b2_v(0);
	double	q0 = q0_v(0);
	double	q2 = q2_v(0);
	double	f0 = f0_v(0);
	double	f2 = f2_v(0);
	double	f1 = f1_v(0);
	double	mu00 = mu00_v(0);
	double	mu02 = mu02_v(0);
	double	theta = theta_v(0);
	double	m0 = m0_v(0);
	double	r0 = r0_v(0);
	double	tau = tau_v(0);
	double	t0 = t0_v(0);
	double	geno_a = geno_a_v(0);
	double	geno_b = geno_b_v(0);
	double	geno_q = geno_q_v(0);
	double	geno_f = geno_f_v(0);
	double	geno_mu = geno_mu_v(0);
	
	double	e3 = geno_q * (q2 - q0)/2 + q0;
	double    e4 = a0 + geno_a * (a2 - a0)/2;
	double    e5 = pow(e4,2);
	double    e7 = sqrt(e5 + 2 * (pow(b0 + geno_b * (b2 - b0)/2,2) * e3));
	double    e8 = 2 * e3;
	double    e9 = e4 + e7;
	double    e12 = f0 + geno_f * (f2 - f0)/2;
	double    e13 = exp(tau * theta);
	double    e17 = geno_mu * (mu02 - mu00)/2 + mu00;
	double    e18 = e13 * e17;
	double    e22 = (pow((1 - 4 * (e3/(e8 - (e8 + 2 * (e7/(r0 - e9/e8))) * 
        exp(2 * (e7 * (tau - t0)))))) * e4 * (e12 - f1)/e7 + 
        2 * ((2 * (e5 * (e12 - f1)/e7) - 2 * ((e12 - m0) * 
            e7)) * exp(-(t0 * e7)) * exp(tau * e7) * e3/((e8 - 
            (e8 + 2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * 
                (tau - t0)))) * e9)),2) - 2 * (e7/(e8 - (e8 + 
        2 * (e7/(r0 - e9/e8))) * exp(2 * (e7 * (tau - t0)))))) * e3 + 0.5 * e9 + e18;
		
	arma::vec df = arma::zeros<arma::vec>(1);

	df(0) = pow(tau,2) * (1 - e18/e22) * e13 * e17/e22;
	
    return(Rcpp::wrap(df));
}

