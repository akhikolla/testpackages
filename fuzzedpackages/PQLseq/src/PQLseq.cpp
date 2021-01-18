#include <fstream>
#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

#include <R.h>
#include <Rmath.h>
#include <cmath>
#include <stdio.h>
#include <stdlib.h> 
#include <cstring>

using namespace std;
using namespace arma;
using namespace Rcpp;


#define ARMA_DONT_PRINT_ERRORS


void invTransformH( vec eigval, mat &H ){

	const size_t num = sum( eigval < 1e-8 );
	const uvec idx = find( eigval < 1e-8 );
	for(size_t i=0; i<num; i++){ eigval[idx[i]] = 1e-5; }

	eigval += 0.5*randu( H.n_rows );

	mat Q(H.n_rows,H.n_cols), R(H.n_rows,H.n_cols), O(H.n_rows, H.n_cols);
	qr( Q, R, H );
	O = Q * diagmat( R.diag()/abs(R.diag()) );

	H = O.t() * diagmat(1.0/eigval) * O;
	return;
}


//*************************************************************************************//
//                             AVERAGE INFORMATION METHOD                              //
//*************************************************************************************//
// [[Rcpp::export]]
SEXP AI(SEXP Yin, SEXP Xin, SEXP numKin, SEXP Phiin, SEXP Din, SEXP tauin, SEXP fixtauin, SEXP tolin)
{/*Average Information*/
	try {
		vec Y = as<vec>(Yin);
		mat X = as<mat>(Xin);
		int numK = Rcpp::as<int>(numKin);
		const Rcpp::List Phi(Phiin);
		vec D = as<vec>(Din);
		vec tau = as<vec>(tauin);
		const uvec fixtau = as<uvec>(fixtauin);
		int numK2 = sum( fixtau == 0 );
		const double tol = Rcpp::as<double>(tolin);
		uvec ZERO = (tau < tol);
		size_t numIDV = X.n_rows, numCVT = X.n_cols;
		mat Hinv(numCVT, numCVT), XtHinvX_inv(numCVT, numCVT), P(numIDV, numIDV);
		vec alpha(numCVT), eta(numIDV);
		cube PHI(numIDV, numIDV, numK);

		mat H = tau[0] * diagmat(1.0 / D);

		for(size_t i=1; i<=numK; ++i) {
		    stringstream kins;
			kins << "kins" << i;
			PHI.slice(i-1) = symmatl(as<mat>(Phi[kins.str()]));
			H = H + tau[i] * PHI.slice(i-1);
		}
		
		mat U;
		vec eigval;
		
		// double time_start=clock();
		eig_sym(eigval, U, H, "dc" );

		if(any(eigval < 1e-8)){
			invTransformH( eigval, H );
			Hinv = H;
		}else{
			Hinv = U * diagmat(1.0/eigval) * U.t();
		}

		mat HinvX = Hinv * X;
		mat XtHinvX = X.t() * HinvX;
		
		mat U2;
		vec eigval2;

		eig_sym( eigval2, U2, XtHinvX, "dc" );

		if(any(eigval2 < 1e-8)){
			invTransformH( eigval2, XtHinvX );
			XtHinvX_inv = XtHinvX;
		}else{
			XtHinvX_inv = U2 * diagmat( 1.0/eigval2 ) * U2.t();
		}
		

		// double time_mv =(clock()-time_start)/(double(CLOCKS_PER_SEC));


		P = Hinv - HinvX * XtHinvX_inv * HinvX.t();
		alpha = XtHinvX_inv * HinvX.t() * Y;
		eta = Y - tau[0] * (Hinv * (Y - X * alpha)) / D;

		if(numK2 > 0) {
		    const uvec idxtau = find(fixtau == 0);
		    mat AImat(numK2, numK2);//average information matrix
		    vec PY = P * Y;
			vec score(numK2), PAPY;
			for(size_t i=0; i<numK2; ++i) {
			    if(i == 0 && idxtau[0] == 0) {
				    PAPY = PY / D;
					score[0] = dot(PAPY, PY) - sum(diagvec(P) / D);
					AImat(0, 0) = dot(PAPY, P * PAPY);
				} else {
				    PAPY = P * PHI.slice(idxtau[i]-1) * PY;
					score[i] = dot(Y, PAPY) - accu(P % PHI.slice(idxtau[i]-1));
					for(size_t j=0; j<=i; ++j) {
					    if(j == 0 && idxtau[0] == 0) {
						    AImat(i,0) = dot(PY / D, PAPY);
							AImat(0,i) = AImat(i,0);
						} else {
						    AImat(i,j) = dot(PY, PHI.slice(idxtau[j]-1) * PAPY);
							if(j!=i) {AImat(j,i) = AImat(i,j);}
						}
					}//end for j
				}
			}// end for i

			vec Dtau = solve(AImat, score);
			vec tau0 = tau;

			tau.elem( idxtau ) = tau0.elem( idxtau ) + Dtau;
		
			tau.elem( find(ZERO % (tau < tol)) ).zeros();
			double step = 1.0;
			while(any(tau < 0.0)) {
			    step *= 0.5;
				tau.elem( idxtau ) = tau0.elem( idxtau ) + step * Dtau;
				tau.elem( find(ZERO % (tau < tol)) ).zeros();
			}
			tau.elem( find(tau < tol) ).zeros();
		}
		// return values
		return List::create(Named("tau") = tau, Named("P") = P, Named("cov") = XtHinvX_inv, 
							Named("alpha") = alpha, Named("eta") = eta, Named("H")=H);
	} catch( std::exception &ex ) {
		forward_exception_to_r( ex );
	} catch(...) {
		::Rf_error( "C++ exception (unknown reason)..." );
	}
	return R_NilValue;
}


/////////////////////////////////////////////////////////////////////////////////////////
//                             CODE END HERE                                           //
/////////////////////////////////////////////////////////////////////////////////////////
