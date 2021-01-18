#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;

// Begin function prototypes
// follwing used for downloaded functions for  solution of cubic and quartic equation
#include <math.h>

#define TwoPi  6.28318530717958648
const double ceps=1e-14;
// poly34.h : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com


// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]±i*x[1]
int   SolveP2(double *x, double a, double b);	 // solve equation x^2 + a*x + b = 0

// x - array of size 3
// return 3: 3 real roots x[0], x[1], x[2]
// return 1: 1 real root x[0] and pair of complex roots: x[1]±i*x[2]
int   SolveP3(double *x, double a, double b, double c);			// solve cubic equation x^3 + a*x^2 + b*x + c = 0

// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]±i*x[3], 
// return 0: two pair of complex roots: x[0]±i*x[1],  x[2]±i*x[3], 
int   SolveP4(double *x,double a,double b,double c,double d);	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d = 0 by Dekart-Euler method

// x - array of size 5
// return 5: 5 real roots x[0], x[1], x[2], x[3], x[4], possible multiple roots
// return 3: 3 real roots x[0], x[1], x[2] and complex x[3]±i*x[4], 
// return 1: 1 real root x[0] and two pair of complex roots: x[1]±i*x[2],  x[3]±i*x[4], 
int   SolveP5(double *x, double a, double b, double c, double d, double e);	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0

//-----------------------------------------------------------------------------
// And some additional functions for internal use.
// Your may remove this definitions from here
int   SolveP4Bi(double *x, double b, double d);				// solve equation x^4 + b*x^2 + d = 0
int   SolveP4De(double *x, double b, double c, double d);	// solve equation x^4 + b*x^2 + c*x + d = 0
void  CSqrt( double x, double y, double &a, double &b);		// returns as a+i*s,  sqrt(x+i*y)
double N4Step(double x, double a,double b,double c,double d);// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
double SolveP5_1(double a,double b,double c,double d,double e);	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0

// end of comment for downloaded functions

// More function prototypes
double logdet(mat lowerTri);
double sumsqr(mat x);
vec  real_polyroots2(vec coef);
double min_root(double quartic, double cubic, double constant, int index, 
		const List& Kerns, vec sigma, mat res, double lambda1, double lambda2,
		vec lambda_factor);

void check_interrupt_fn(void *dummy);
int pending_interrupt();
double signif_dig(double x, int digits);
double check_poly(double quartic, double cubic, double constant, double sigma);

vec sqrt_vec(vec x);
arma::mat gchol(arma::mat matrix);
arma::mat gchol_inv(arma::mat matrix);

// Main function

//[[Rcpp::export]]
List vcpen_Rcpp(arma::vec& y, arma::mat& X, const List& Kerns, arma::vec lambda_factor,  arma::vec lambda_grid,  
		double frac1, arma::vec sigma2_init, int maxiter, bool print_iter){

  // frac1 is fraction of lambda_grid for lambda1 (lambda1 = frac1*lambda_grid)
  //       so lambda2 = (1-frac) * lambda_grid

  // lambda1 is penalty multiplier for L1 penalty ( lambda1 * sum(sigma) )
  // lambda2 is penalty multiplier for L2 penalty ( lambda2 * sum(sigma^2) )


  // logl will be scaled by 1/n, so penalty lambdas do not scale with n.
  // this makes lambdas to have small values. 


  double pi = 3.141593;
  double tol = 0.0000001;
  double eps8 = 0.000000001;

  double lambda1 = 0.0;
  double lambda2 = 0.0;
  double logllassoold;

  int n_vc = Kerns.size();
  int n_subj = y.size();
  int n_beta = X.n_cols;
  int n_grid = lambda_grid.size();


  Col <int> iter_grid(n_grid, fill::zeros);
  vec logl_grid(n_grid, fill::zeros);
  vec logllasso_grid(n_grid, fill::zeros);
  mat beta_grid(n_beta, n_grid, fill::zeros);
  mat sigma2_grid(n_vc, n_grid, fill::zeros);

  double logllasso_final;
  double logl_final;
  int iter_final;
  mat Xnew;
  mat ynew;

  mat xx =  trans(X) * X;
  mat xy =  trans(X) * y;
  mat beta = solve(xx, xy);

  vec sigma2 = sigma2_init;

  for(int i = 0; i < sigma2.size(); i ++){
    if(sigma2(i) < 0.00001){
      sigma2(i) = 0.1;
    }
  }

  vec sigma  = sqrt_vec(sigma2_init);
 
  mat omega = zeros(n_subj, n_subj);
  mat omega_chol = omega;
  mat chol_inv = omega;

  // initialization

  for (int i = 0; i < n_vc; ++i) {
    omega = omega + sigma2(i) * as<mat>(Kerns[i]);
  }

  omega_chol = gchol(omega);
  chol_inv = gchol_inv(omega_chol);

  mat Vinv = trans(chol_inv) * chol_inv;

  // res = Linv (y- X*beta), where V = L L'

  mat res = chol_inv * (y - (X * beta));

  // res_b = Linv'Linv * res = Vinv * res
 
  mat res_b = trans(chol_inv) * res;
  
  double loglConst = - 0.5 * n_subj * log(2.0 * pi);
   
  double logl = loglConst  - 0.5 * logdet(omega_chol)  - 0.5 * sumsqr(res);
 
  // scale logl by 1/n
  logl = logl / double(n_subj);

  double quad = 0.0;

  // scaled trace and quad
  double quad_n = 0.0;
  double trace_n = 0.0;

  int n_vc_penalized = n_vc-1;
  double sigma_sum = 0.0;
  double sigma2_sum = 0.0;
  double quartic = 0.0;
  double cubic = 0.0;
  double logllasso = 0.0;
  mat q(1,1,fill::zeros);
  double constant = 0.0;
     
  // loop over lambda grid values
 
  for(int index_grid = 0; index_grid < n_grid; index_grid ++){

    // if initial sigma2 = 0, then the estimation will be stuck at 0,
    // because updates to sigma2 multiply current estimate by a factor, 
    // so if current estimate = 0, updates will also be 0.
    // For this reason, if initial sigma2 near 0, move it away for 
    // its inital value
    for(int i = 0; i < sigma2.size(); i ++){
      if(sigma2(i) < 0.00001){
	sigma2(i) = 0.1;
      }
    }

    if(print_iter){
      Rcout << "index_grid = " << index_grid;
      Rcout << ", sigma2 = ";
      for(int i = 0; i < sigma2.size(); i++){
	Rcout  << sigma2(i) << ", ";
      }
      Rcout << endl;
    }

  
    lambda1 = frac1 * lambda_grid(index_grid);
    lambda2 = (1.0 - frac1) * lambda_grid(index_grid);
   
    sigma_sum = 0.0;
    sigma2_sum = 0.0;
 
    for(int i = 0; i < n_vc_penalized; i++){
	sigma_sum  += fabs(lambda_factor(i) * sigma(i));
	sigma2_sum +=  lambda_factor(i) * sigma2(i);   
    }

    
    logllasso = -logl + lambda1*sigma_sum + lambda2*sigma2_sum;
   
    for(int iter = 0; iter < maxiter; iter++){

    
      if(print_iter){
	if(iter % 100 == 0) Rcout << "... iter = " << iter << endl;
      }

      // update variance components
      omega.fill(0.0);

      for(int index_vc = 0; index_vc < n_vc_penalized; index_vc++){

	trace_n = trace(Vinv *  as<mat>(Kerns[index_vc]) )/double(n_subj);
	q =  trans(res_b) * as<mat>(Kerns[index_vc]) * res_b;
	quad_n = q(0,0) / double(n_subj);
	
	if(lambda1 < tol){
	  // L2 penalty alone has solution to quadratic equation
	  sigma2(index_vc) = sigma2(index_vc) * 
	                 sqrt(quad_n) / sqrt(trace_n + 2.0*lambda2*lambda_factor(index_vc));
	  sigma(index_vc) = sqrt(sigma2(index_vc));
	} else{
	  // L1 penalty alone or L1 & L2 penalties jointly
	  // has solution to polynomial equation of 4th degree

	  quartic = trace_n + 2.0 * lambda2 * lambda_factor(index_vc);
	  cubic = lambda1 * lambda_factor(index_vc);
	  constant = -pow(sigma(index_vc), 4) * quad_n;

	  sigma(index_vc) = min_root(quartic, cubic, constant, index_vc, Kerns, sigma, res, 
				     lambda1, lambda2, lambda_factor);
	  sigma2(index_vc) = sigma(index_vc) * sigma(index_vc);
	}

	omega = omega + sigma2(index_vc) * as<mat>(Kerns[index_vc]);
      }

     
      // update last vc, for residual error vc
      sigma2(n_vc-1) =  sigma2(n_vc-1) * sqrt( sumsqr(res_b) /  trace(Vinv) ) ;
      sigma(n_vc-1)  = sqrt(sigma2(n_vc-1));
      // problem: when resid vc -> 0, varmat becomes non-singular,
      // so below bounds to a small non-zero value

      if(sigma2(n_vc-1) < eps8){sigma2(n_vc-1) =  eps8; };

      // add last vc to omega
      omega =  omega + sigma2(n_vc-1)  * as<mat>(Kerns[n_vc-1]);
        
      // now with updated omega, update cholesky, inv,  beta, res, res_b

      omega_chol = gchol(omega);
      chol_inv = gchol_inv(omega_chol);


      Vinv = trans(chol_inv) * chol_inv;
      Xnew  =  chol_inv *  X; 
      ynew =   chol_inv * y;
      xx =  trans(Xnew) * Xnew;
      xy =  trans(Xnew) * ynew;
      beta = solve(xx, xy);
      res = chol_inv * (y - (X * beta));
      res_b = trans(chol_inv) * res;

 
      //  check convergence
      logllassoold =  logllasso;

      logl = loglConst  - 0.5 * logdet(omega_chol)  - 0.5 * sumsqr(res);
      // scaled logl by 1/n
      logl = logl / double(n_subj);

      sigma_sum = 0.0;
      sigma2_sum = 0.0;
      
      for(int i = 0; i < n_vc_penalized; i++){
	sigma_sum  += fabs(lambda_factor(i) * sigma(i));
	sigma2_sum += lambda_factor(i) * sigma2(i);
      }


      logllasso = -logl + lambda1 * sigma_sum + lambda2 * sigma2_sum;
       

      // save values in case of convergenece
      logllasso_final = logllasso;
      logl_final = logl;
      iter_final = iter;


      if(print_iter){
	Rcout << "     iter = " << iter;
	Rcout << ", logl = " << logl;
	Rcout << ", logllasso = " << logllasso;
	Rcout << ", sigma2 = ";
	for(int i = 0; i < sigma2.size(); i++){
	  Rcout  << sigma2(i) << ", ";
	}
	Rcout << endl;
      }

    
      if (fabs(logllasso - logllassoold) < tol * (fabs(logllassoold) + 1.0) ) {
	break; // this breaks out of iter loop since converged
      }

    } // end of iter loop

    // save values for a grid point
   
    iter_grid(index_grid) = iter_final;
    logl_grid(index_grid) = logl_final * double(n_subj);
    logllasso_grid(index_grid) = logllasso_final * double(n_subj);

    for(int i = 0; i < n_beta; i++){
      beta_grid(i,index_grid) = beta(i);
    }
    for(int i = 0; i < n_vc; i++){
      sigma2_grid(i, index_grid) = sigma2(i);
    }
  
  } // end of k loop over grid values
 

  return Rcpp::List::create(Rcpp::Named("n_vc") = n_vc,
			    Rcpp::Named("n_subj") = n_subj,
			    Rcpp::Named("n_beta")  = n_beta,
			    Rcpp::Named("n_grid") = n_grid,
			    Rcpp::Named("frac1") = frac1,
			    Rcpp::Named("lambda_grid") = lambda_grid,
			    Rcpp::Named("beta_grid") = beta_grid,
			    Rcpp::Named("iter_grid") = iter_grid,
			    Rcpp::Named("logl_grid") = logl_grid,
			    Rcpp::Named("logllasso_grid") = logllasso_grid,
 			    Rcpp::Named("vc_grid") = sigma2_grid);

}


// Utility functions
double signif_dig(double x, int digits){
  // round to specified significant digits
  return round(x*pow(10.0, digits))/ pow(10.0, digits);
}


double logdet(mat lowerTri){
  // assume input is square matrix L from LL' = V,
  // L is lower triangle matrix from chol(V)
  int nrow = lowerTri.n_rows;
  double lndet = 0.0;
  for(int i = 0; i < nrow; i++){
    if(lowerTri(i,i) > 0.0){
      lndet += log(lowerTri(i,i));
    }
  }
  lndet = 2.0 * lndet;
  return lndet;
}

vec sqrt_vec(vec x2){
  vec x(x2.size(), fill::zeros);
  for(int i =0; i < x2.size(); i++){
    x(i) = sqrt(x2(i));
  }
  return x;
}

  
double sumsqr(mat x){
  double sum = 0.0;
  int ncol = x.n_cols;
  int nrow = x.n_rows;
  for(int i = 0; i < nrow; i++){
    for(int j = 0; j < ncol; j++){
      sum  += x(i,j)*x(i,j);
    }
  }
  return sum;
}

double min_root(double quartic, double cubic, double constant, int index, 
		const List& Kerns, vec sigma, mat res, double lambda1, double lambda2,
		vec lambda_factor){

  // this function returns the sigma that minimizes the objective function: 
  // negative loglikelihood plus penalty term
  // find all real roots, evaluate pl at each root, 
  // and choose root that minimizes pl 


  // % is Schur product: element-wise multiplication of two objects
  vec sigma2 = sigma % sigma;

  double sigma_sum = 0.0;
  double sigma2_sum = 0.0;
  double d_subj = res.size();

  vec coef(5, fill::zeros);
  coef(0) = constant;
  coef(3) = cubic;
  coef(4) = quartic;

  vec roots = real_polyroots2(coef);
  for(int i = 0; i < roots.size();i++){
    roots(i) = fabs(roots(i));
  }

  int n_vc = sigma.size();
  int n_vc_penalized = n_vc - 1;
  int n_subj =  as<mat>(Kerns[0]).n_rows;
  vec obj(roots.size(), fill::zeros);
  int n_roots = roots.size();
  mat omega(n_subj, n_subj, fill::zeros);
  mat omega_chol(n_subj, n_subj, fill::zeros);
  mat chol_inv(n_subj, n_subj, fill::zeros);
  mat resnew(n_subj,1, fill::zeros);

  for(int i = 0; i < n_roots; i++){
    sigma(index) = roots(i);
    sigma2(index) =  sigma(index) * sigma(index);
    omega.fill(0.0);
    for(int j = 0; j< n_vc; j++){
      omega = omega + sigma2(j) * as<mat>(Kerns[j]);
    }

    omega_chol = gchol(omega);
    chol_inv = gchol_inv(omega_chol);
  
    resnew = chol_inv * res;

    // obj is  penalized lnlike (eq 4)
    // note that cubic here is lambda

    sigma_sum = 0.0;
    sigma2_sum = 0.0;
    for(int j = 0; j < n_vc_penalized; j++){
      sigma_sum += fabs(lambda_factor(j) * sigma(j));
      sigma2_sum += lambda_factor(j) * sigma2(j);
    }

    
    obj(i) =  (0.5 * logdet(omega_chol) + 0.5 * sumsqr(resnew)) / d_subj 
                 + lambda1*sigma_sum + lambda2*sigma2_sum;
  }

  // index of min element in array.
  uword index_min = obj.index_min();
 
  return roots(index_min);
}

vec  real_polyroots2(vec coef){
  // input vector of coefficients for polynomial and 
  // return vector of unique real roots. If none, return 0.
  // uses  Descartes-Euler method

  vec roots(4, fill::zeros);
  double *x = (double *)malloc( (size_t) (4 * sizeof(double)) );

  double a = coef(3) / coef(4);
  double b = coef(2) / coef(4);
  double c = coef(1) / coef(4);
  double d = coef(0) / coef(4);
			    
  int solType = SolveP4(x, a, b, c, d);
  
  if(solType == 2){
    roots.resize(2);
    roots(0) = x[0];
    roots(1) = x[1];
  }
  if(solType == 4){
    for(int i = 0; i < roots.size(); i++){
      roots(i) = x[i];
    }
  }
  vec uroots = unique(roots);

  free(x);
  return(uroots);

}

 

double check_poly(double quartic, double cubic, double constant, double sigma){
  double poly;
  poly = pow(sigma, 4.0)*quartic + cubic*pow(sigma, 3.0) + constant;
  return poly ;
}

/* Check for interrupt without long jumping */
void check_interrupt_fn(void *dummy) {
  R_CheckUserInterrupt();
}

int pending_interrupt() {
  return !(R_ToplevelExec(check_interrupt_fn, NULL));
} 

// poly34.cpp : solution of cubic and quartic equation
// (c) Khashin S.I. http://math.ivanovo.ac.ru/dalgebra/Khashin/index.html
// khash2 (at) gmail.com
// Thanks to Alexandr Rakhmanin <rakhmanin (at) gmail.com>
// public domain
// DJS downloaded from http://math.ivanovo.ac.ru/dalgebra/Khashin/poly/index.html
// on Apr 5, 2019


//=============================================================================
// _root3, root3 from http://prografix.narod.ru
//=============================================================================
static double _root3 ( double x )
{
    double s = 1.;
    while ( x < 1. )
    {
        x *= 8.;
        s *= 0.5;
    }
    while ( x > 8. )
    {
        x *= 0.125;
        s *= 2.;
    }
    double r = 1.5;
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    r -= 1./3. * ( r - x / ( r * r ) );
    return r * s;
}

double root3 ( double x )
{
    if ( x > 0 ) return _root3 ( x ); else
    if ( x < 0 ) return-_root3 (-x ); else
    return 0.;
}


// x - array of size 2
// return 2: 2 real roots x[0], x[1]
// return 0: pair of complex roots: x[0]±i*x[1]
int   SolveP2(double *x, double a, double b) {			// solve equation x^2 + a*x + b = 0
	double D = 0.25*a*a - b;
	if (D >= 0) {
		D = sqrt(D);
		x[0] = -0.5*a + D;
		x[1] = -0.5*a - D;
		return 2;
	}
	x[0] = -0.5*a;
	x[1] = sqrt(-D);
	return 0;
}
//---------------------------------------------------------------------------
// x - array of size 3
// In case 3 real roots: => x[0], x[1], x[2], return 3
//         2 real roots: x[0], x[1],          return 2
//         1 real root : x[0], x[1] ± i*x[2], return 1

int SolveP3(double *x, double a, double b, double c) { // solve cubic equation x^3 + a*x^2 + b*x + c = 0
    double a2 = a*a;
    double q = (a2 - 3 * b) / 9;
    double r = (a * (2 * a2 - 9 * b) + 27 * c) / 54;
    // equation x^3 + q*x + r = 0
    double r2 = r*r;
    double q3 = q * q*q;
    double A, B;
    if (r2 <= (q3 + ceps)) {//<<-- FIXED!
        double t = r / sqrt(q3);
        if (t<-1) t = -1;
        if (t > 1) t = 1;
        t = acos(t);
        a /= 3;
        q = -2 * sqrt(q);
        x[0] = q * cos(t / 3) - a;
        x[1] = q * cos((t + TwoPi) / 3) - a;
        x[2] = q * cos((t - TwoPi) / 3) - a;
        return (3);
    } else {
        //A =-pow(fabs(r)+sqrt(r2-q3),1./3); 
        A = -root3(fabs(r) + sqrt(r2 - q3));
        if (r < 0) A = -A;
        
        if(A == 0){
            B = 0;
        } else {
            B = q / A;
        }
        // original code that is confusing and caused warning
        // B = (A == 0 ? 0 : B = q / A);
        // so replaced with above if/else
        
        
        a /= 3;
        x[0] = (A + B) - a;
        x[1] = -0.5 * (A + B) - a;
        x[2] = 0.5 * sqrt(3.)*(A - B);
        if (fabs(x[2]) < ceps) {
            x[2] = x[1];
            return (2);
        }
        return (1);
    }
}// SolveP3(double *x,double a,double b,double c) {	
//---------------------------------------------------------------------------
// a>=0!
void  CSqrt( double x, double y, double &a, double &b) // returns:  a+i*s = sqrt(x+i*y)
{
	double r  = sqrt(x*x+y*y);
	if( y==0 ) { 
		r = sqrt(r);
		if(x>=0) { a=r; b=0; } else { a=0; b=r; }
	} else {		// y != 0
		a = sqrt(0.5*(x+r));
		b = 0.5*y/a;
	}
}
//---------------------------------------------------------------------------
int   SolveP4Bi(double *x, double b, double d)	// solve equation x^4 + b*x^2 + d = 0
{
	double D = b*b-4*d;
	if( D>=0 ) 
	{
		double sD = sqrt(D);
		double x1 = (-b+sD)/2;
		double x2 = (-b-sD)/2;	// x2 <= x1
		if( x2>=0 )				// 0 <= x2 <= x1, 4 real roots
		{
			double sx1 = sqrt(x1);
			double sx2 = sqrt(x2);
			x[0] = -sx1;
			x[1] =  sx1;
			x[2] = -sx2;
			x[3] =  sx2;
			return 4;
		}
		if( x1 < 0 )				// x2 <= x1 < 0, two pair of imaginary roots
		{
			double sx1 = sqrt(-x1);
			double sx2 = sqrt(-x2);
			x[0] =    0;
			x[1] =  sx1;
			x[2] =    0;
			x[3] =  sx2;
			return 0;
		}
		// now x2 < 0 <= x1 , two real roots and one pair of imginary root
			double sx1 = sqrt( x1);
			double sx2 = sqrt(-x2);
			x[0] = -sx1;
			x[1] =  sx1;
			x[2] =    0;
			x[3] =  sx2;
			return 2;
	} else { // if( D < 0 ), two pair of compex roots
		double sD2 = 0.5*sqrt(-D);
		CSqrt(-0.5*b, sD2, x[0],x[1]);
		CSqrt(-0.5*b,-sD2, x[2],x[3]);
		return 0;
	} // if( D>=0 ) 
} // SolveP4Bi(double *x, double b, double d)	// solve equation x^4 + b*x^2 d
//---------------------------------------------------------------------------
#define SWAP(a,b) { t=b; b=a; a=t; }
static void  dblSort3( double &a, double &b, double &c) // make: a <= b <= c
{
	double t;
	if( a>b ) SWAP(a,b);	// now a<=b
	if( c<b ) {
		SWAP(b,c);			// now a<=b, b<=c
		if( a>b ) SWAP(a,b);// now a<=b
	}
}
//---------------------------------------------------------------------------
int   SolveP4De(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
{
	//if( c==0 ) return SolveP4Bi(x,b,d); // After that, c!=0
	if( fabs(c)<1e-14*(fabs(b)+fabs(d)) ) return SolveP4Bi(x,b,d); // After that, c!=0

	int res3 = SolveP3( x, 2*b, b*b-4*d, -c*c);	// solve resolvent
	// by Viet theorem:  x1*x2*x3=-c*c not equals to 0, so x1!=0, x2!=0, x3!=0
	if( res3>1 )	// 3 real roots, 
	{				
		dblSort3(x[0], x[1], x[2]);	// sort roots to x[0] <= x[1] <= x[2]
		// Note: x[0]*x[1]*x[2]= c*c > 0
		if( x[0] > 0) // all roots are positive
		{
			double sz1 = sqrt(x[0]);
			double sz2 = sqrt(x[1]);
			double sz3 = sqrt(x[2]);
			// Note: sz1*sz2*sz3= -c (and not equal to 0)
			if( c>0 )
			{
				x[0] = (-sz1 -sz2 -sz3)/2;
				x[1] = (-sz1 +sz2 +sz3)/2;
				x[2] = (+sz1 -sz2 +sz3)/2;
				x[3] = (+sz1 +sz2 -sz3)/2;
				return 4;
			}
			// now: c<0
			x[0] = (-sz1 -sz2 +sz3)/2;
			x[1] = (-sz1 +sz2 -sz3)/2;
			x[2] = (+sz1 -sz2 -sz3)/2;
			x[3] = (+sz1 +sz2 +sz3)/2;
			return 4;
		} // if( x[0] > 0) // all roots are positive
		// now x[0] <= x[1] < 0, x[2] > 0
		// two pair of comlex roots
		double sz1 = sqrt(-x[0]);
		double sz2 = sqrt(-x[1]);
		double sz3 = sqrt( x[2]);

		if( c>0 )	// sign = -1
		{
			x[0] = -sz3/2;					
			x[1] = ( sz1 -sz2)/2;		// x[0]±i*x[1]
			x[2] =  sz3/2;
			x[3] = (-sz1 -sz2)/2;		// x[2]±i*x[3]
			return 0;
		}
		// now: c<0 , sign = +1
		x[0] =   sz3/2;
		x[1] = (-sz1 +sz2)/2;
		x[2] =  -sz3/2;
		x[3] = ( sz1 +sz2)/2;
		return 0;
	} // if( res3>1 )	// 3 real roots, 
	// now resoventa have 1 real and pair of compex roots
	// x[0] - real root, and x[0]>0, 
	// x[1]±i*x[2] - complex roots, 
	// x[0] must be >=0. But one times x[0]=~ 1e-17, so:
	if (x[0] < 0) x[0] = 0;
	double sz1 = sqrt(x[0]);
	double szr, szi;
	CSqrt(x[1], x[2], szr, szi);  // (szr+i*szi)^2 = x[1]+i*x[2]
	if( c>0 )	// sign = -1
	{
		x[0] = -sz1/2-szr;			// 1st real root
		x[1] = -sz1/2+szr;			// 2nd real root
		x[2] = sz1/2; 
		x[3] = szi;
		return 2;
	}
	// now: c<0 , sign = +1
	x[0] = sz1/2-szr;			// 1st real root
	x[1] = sz1/2+szr;			// 2nd real root
	x[2] = -sz1/2;
	x[3] = szi;
	return 2;
} // SolveP4De(double *x, double b, double c, double d)	// solve equation x^4 + b*x^2 + c*x + d
//-----------------------------------------------------------------------------
double N4Step(double x, double a,double b,double c,double d)	// one Newton step for x^4 + a*x^3 + b*x^2 + c*x + d
{
	double fxs= ((4*x+3*a)*x+2*b)*x+c;	// f'(x)
	if (fxs == 0) return x;	//return 1e99; <<-- FIXED!
	double fx = (((x+a)*x+b)*x+c)*x+d;	// f(x)
	return x - fx/fxs;
} 
//-----------------------------------------------------------------------------
// x - array of size 4
// return 4: 4 real roots x[0], x[1], x[2], x[3], possible multiple roots
// return 2: 2 real roots x[0], x[1] and complex x[2]±i*x[3], 
// return 0: two pair of complex roots: x[0]±i*x[1],  x[2]±i*x[3], 
int   SolveP4(double *x,double a,double b,double c,double d) {	// solve equation x^4 + a*x^3 + b*x^2 + c*x + d by Dekart-Euler method
	// move to a=0:
	double d1 = d + 0.25*a*( 0.25*b*a - 3./64*a*a*a - c);
	double c1 = c + 0.5*a*(0.25*a*a - b);
	double b1 = b - 0.375*a*a;
	int res = SolveP4De( x, b1, c1, d1);
	if( res==4) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; x[3]-= a/4; }
	else if (res==2) { x[0]-= a/4; x[1]-= a/4; x[2]-= a/4; }
	else             { x[0]-= a/4; x[2]-= a/4; }
	// one Newton step for each real root:
	if( res>0 )
	{
		x[0] = N4Step(x[0], a,b,c,d);
		x[1] = N4Step(x[1], a,b,c,d);
	}
	if( res>2 )
	{
		x[2] = N4Step(x[2], a,b,c,d);
		x[3] = N4Step(x[3], a,b,c,d);
	}
	return res;
}
//-----------------------------------------------------------------------------
#define F5(t) (((((t+a)*t+b)*t+c)*t+d)*t+e)
//-----------------------------------------------------------------------------
double SolveP5_1(double a,double b,double c,double d,double e)	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
	int cnt;
	if( fabs(e)< ceps ) return 0;

	double brd =  fabs(a);			// brd - border of real roots
	if( fabs(b)>brd ) brd = fabs(b);
	if( fabs(c)>brd ) brd = fabs(c);
	if( fabs(d)>brd ) brd = fabs(d);
	if( fabs(e)>brd ) brd = fabs(e);
	brd++;							// brd - border of real roots

	double x0, f0;					// less than root
	double x1, f1;					// greater than root
	double x2, f2, f2s;				// next values, f(x2), f'(x2)
	double dx;

	if( e<0 ) { x0 =   0; x1 = brd; f0=e; f1=F5(x1); x2 = 0.01*brd; }	// positive root
	else	  { x0 =-brd; x1 =   0; f0=F5(x0); f1=e; x2 =-0.01*brd; }	// negative root

	if( fabs(f0)< ceps ) return x0;
	if( fabs(f1)< ceps ) return x1;

	// now x0<x1, f(x0)<0, f(x1)>0
	// Firstly 10 bisections
	for( cnt=0; cnt<10; cnt++)
	{
		x2 = (x0 + x1) / 2;					// next point
		//x2 = x0 - f0*(x1 - x0) / (f1 - f0);		// next point
		f2 = F5(x2);				// f(x2)
		if( fabs(f2)< ceps ) return x2;
		if( f2>0 ) { x1=x2; f1=f2; }
		else       { x0=x2; f0=f2; }
	}

	// At each step:
	// x0<x1, f(x0)<0, f(x1)>0.
	// x2 - next value
	// we hope that x0 < x2 < x1, but not necessarily
	do {
		if(cnt++>50) break;
		if( x2<=x0 || x2>= x1 ) x2 = (x0 + x1)/2;	// now  x0 < x2 < x1
		f2 = F5(x2);								// f(x2)
		if( fabs(f2)< ceps ) return x2;
		if( f2>0 ) { x1=x2; f1=f2; }
		else       { x0=x2; f0=f2; }
		f2s= (((5*x2+4*a)*x2+3*b)*x2+2*c)*x2+d;		// f'(x2)
		if( fabs(f2s)< ceps ) { x2=1e99; continue; }
		dx = f2/f2s;
		x2 -= dx;
	} while(fabs(dx)> ceps);
	return x2;
} // SolveP5_1(double a,double b,double c,double d,double e)	// return real root of x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------
int   SolveP5(double *x,double a,double b,double c,double d,double e)	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
{
	double r = x[0] = SolveP5_1(a,b,c,d,e);
	double a1 = a+r, b1=b+r*a1, c1=c+r*b1, d1=d+r*c1;
	return 1+SolveP4(x+1, a1,b1,c1,d1);
} // SolveP5(double *x,double a,double b,double c,double d,double e)	// solve equation x^5 + a*x^4 + b*x^3 + c*x^2 + d*x + e = 0
//-----------------------------------------------------------------------------


arma::mat gchol(arma::mat matrix) {

/*

A symmetric matrix A can be decomposed as LDL', where L is a lower triangular matrix 
with 1's on the diagonal, L' is the transpose of L, and D is diagonal. The inverse of 
L is also lower-triangular, with 1's on the diagonal. If all elements of D are positive, 
then A must be symmetric positive definite (SPD), and the solution can be reduced 
the usual Cholesky decomposition U'U where U is upper triangular and U = sqrt(D) L'.

The main advantage of the generalized form is that it admits of matrices that are not of 
full rank: D will contain zeros marking the redundant columns, and the rank of A is the 
number of non-zero columns. If all elements of D are zero or positive, then A is a 
non-negative definite (NND) matrix. The generalized form also has the (quite minor) 
numerical advantage of not requiring square roots during its calculation. 
subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**

subroutine to do Cholesky decompostion on a matrix: C = FDF'
**   where F is lower triangular with 1's on the diagonal, and D is diagonal
**
** arguments are:
**     n         the size of the matrix to be factored
**     **matrix  a ragged array containing an n by n submatrix to be factored
**
**  The factorization is returned in the lower triangle, D occupies the
**    diagonal and the upper triangle is left undisturbed.
**    The lower triangle need not be filled in at the start.
**
**
**  If a column is deemed to be redundant, then that diagonal is set to zero.
**
**   Terry Therneau

DJS: return  L * diag(sqrt(D))

*/

  double EPSILON = 0.000000001;     /* <= EPS is considered a zero */
    
  double temp;   /* took out register definition from original code from Terry */
  int  i,j,k;
  double eps, pivot;
  int rank;
  int n =  matrix.n_rows;
  
  eps =0;
  for (i=0; i<n; i++) {
    if (matrix(i,i) > eps)  eps = matrix(i,i);
    for (j=(i+1); j<n; j++)  matrix(j,i) = matrix(i,j);
  }
  eps *= EPSILON;
  
  rank =0;
  for (i=0; i<n; i++) {
    pivot = matrix(i,i);
    if (pivot < eps) matrix(i,i) =0;
    else  {
      rank++;
      for (j=(i+1); j<n; j++) {
	temp = matrix(j,i)/pivot;
	matrix(j,i) = temp;
	matrix(j,j) -= temp*temp*pivot;
	for (k=(j+1); k<n; k++) matrix(k,j) -= temp*matrix(k,i);
      }
    }
  }
  
  vec diag(n);
  for(i = 0; i < n;i++){
    diag(i) = sqrt(matrix(i,i));
    matrix(i,i) = 1.0;
  }
  
  for(i=0; i<(n-1); i++){
    for(j=i+1;j<n; j++){
      matrix(i,j) = 0.0;
    }
  }

  matrix = matrix * diagmat(diag);

  return  matrix;
}


arma::mat gchol_inv(arma::mat matrix) {
  // find inverse of lower triangular matrix that
  // is a generalized cholesky  L * diag(sqrt(D))

  int n = matrix.n_rows;
  arma::mat inv = matrix;

  int i,j,k;
  double diag;

  for (k=0; k<n; k++){
  
    if (inv(k,k) > 0.0) {

      diag = inv(k,k);
      for (i=0; i < k; i++){
	inv(k,i) = inv(k,i)/diag;
      }
      for(i=k; i< n; i++){
	inv(i,k) = inv(i,k)/diag;
      }
      
      for(i=0; i < n; i++){
	if(i == k) continue;
	for(j=0; j < n; j++){
	  if(j == k) continue;
	  inv(i,j) = inv(i,j) - inv(i,k) * inv(k,j) * diag;
	}
      }
      inv(k,k) = - 1.0/diag;
    }
  }
  for(i=0; i < n; i++){
    for(j=0; j <=i; j++){
      inv(i,j) = - inv(i,j);
    }
  }

  return(inv);
}
